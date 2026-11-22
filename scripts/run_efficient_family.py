"""VNTR effective dictionary/HMM/families; STR legacy frozen-dictionary graphs."""
import argparse,json,subprocess,time
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
from vntrgraph import Locus,IndexedFastaReference
from vntrpathgraph import assign_reference_paths
from vntrlocal import classify,assign_slots,path_view,topology,plan_edits
from vntrlocal_export import export_budget
from vntrdisplay_export import dump,table
from vntrpathstructure import annotate_paths
from vntreffmotifs import effective_model,require_original_model,select_motifs
from vntrfamilies import select_families,apply_families
from vntrpartition import repeat_type,str_policy_path

VERSION='efficient-motif-family-v2'


def annotate(records,model_path,locus):
    model,model_hash,features=annotate_paths(records,model_path,locus,decoder='hmm_single')
    tokens={p:classify(f['tokens'],[]) for p,f in features.items()}
    slots={p:assign_slots(ts,tokens['P0001'])[0] for p,ts in tokens.items()}
    return model,model_hash,features,tokens,slots


def member_identity(records, model, fasta_sequence=None):
    """Compare independent member evidence before constructing paths from the model."""
    if len({r['id'] for r in records}) != len(records):
        raise ValueError('Member identity: duplicate path ID')
    refs = [r for r in records if r.get('role') == 'reference']
    if len(refs) > 1 or (refs and refs[0]['id'] != 'P0001'):
        raise ValueError('Member identity: expected a unique reference named P0001')
    checked = set()
    expected = {k: model[k] for k in ('locus', 'assembly', 'chrom', 'start', 'end')}

    def compare(key, value, source):
        if value != expected[key]:
            raise ValueError(f'Member identity mismatch in {source}: {key}={value!r}, expected {expected[key]!r}')
        checked.add(key)

    for row in records:
        if row.get('sha256') != sha256(row['sequence'].encode('ascii')).hexdigest():
            raise ValueError(f"Member identity: sequence hash mismatch in {row['id']}")
        for item in (row, *row.get('observations', [])):
            for key in ('locus','assembly','chrom'):
                if item.get(key) is not None:
                    compare(key, item[key], row['id'])
            interval = item.get('reference_interval')
            if isinstance(interval, dict):
                if 'chrom' in interval: compare('chrom', interval['chrom'], row['id'])
                bounds = (interval.get('start'),interval.get('end'))
            elif isinstance(interval, (tuple,list)) and len(interval) in (2,3):
                if len(interval)==3: compare('chrom',interval[0],row['id'])
                bounds = interval[-2:]
            elif interval is None:
                continue
            else:
                raise ValueError(f"Member identity: malformed reference_interval in {row['id']}")
            if any(type(x) is not int for x in bounds) or bounds[0]>bounds[1]:
                raise ValueError(f"Member identity: invalid reference_interval in {row['id']}")
            if item is row and row.get('role')=='reference':
                compare('start', bounds[0], row['id'])
                compare('end', bounds[1], row['id'])
            elif bounds[1]<model['start'] or bounds[0]>model['end']:
                raise ValueError(f"Member identity: disjoint reference interval in {row['id']}: {bounds}")
    if refs:
        ref = refs[0]
        if len(ref['sequence']) != model['end']-model['start']:
            raise ValueError('Member identity: reference length differs from the model interval')
        if ref.get('molecules') or ref.get('observations'):
            raise ValueError('Member identity: reference has observed support')
        if fasta_sequence is not None:
            if ref['sequence'] != fasta_sequence:
                raise ValueError('Member identity: reference bases differ from supplied FASTA')
            checked.add('reference_bases')
    missing = sorted((set(expected) | {'reference_bases'}) - checked)
    return dict(status='partially_verified' if missing else 'verified',
                checked=sorted(checked), missing=missing, reference_in_members=bool(refs))


def run(input_path,model_path,out,q=None,fraction=.03,reference=None,render=True,supernode_policy=None):
    start=time.monotonic()
    records=json.loads(input_path.read_text());original_model=json.loads(model_path.read_text())
    require_original_model(original_model)
    if supernode_policy is not None and q is not None and not 0 <= q <= 1:
        raise ValueError('q must be in [0,1]')
    locus=Locus(original_model['locus'],original_model['chrom'],original_model['start'],original_model['end'],(),boundary_status=original_model['boundary_status'])
    reference_sequence=None
    if reference is not None:
        with IndexedFastaReference(reference) as fasta:
            reference_sequence=fasta.fetch(locus.chrom,locus.start,locus.end).upper()
    identity=member_identity(records,original_model,reference_sequence)
    for r in records:r.setdefault('source_path_id',r['id'])
    if not any(r.get('role')=='reference' for r in records):
        if reference is None:raise ValueError('An observation-only input needs its reference FASTA')
        records=assign_reference_paths(records,reference_sequence,[locus.chrom,locus.start,locus.end])
    records.sort(key=lambda r:r['id'])
    if records[0]['id']!='P0001' or records[0].get('role')!='reference':
        raise ValueError('Member identity: invalid reference path ordering')
    out.mkdir(parents=True,exist_ok=False)
    dump(out/'input_identity.json',identity)
    dump(out/'source_members.json',records);dump(out/'locus_model.original.json',original_model)
    table(out/'input_path_mapping.tsv',[dict(path=r['id'],source_path=r.get('source_path_id',''),role=r['role'],sha256=r['sha256']) for r in records])
    kind=repeat_type(original_model)
    if supernode_policy is not None:
        if kind!='VNTR': raise ValueError('Supernode policy requires a VNTR locus')
        from vntrsupernode_export import export_supernodes
        return export_supernodes(out,locus,records,model_path,supernode_policy,render,input_path,q=q)
    if q is None: q=.1
    processing=dict(repeat_type=kind,decoder='seed_dp' if kind=='STR' else 'hmm_single',
        hmm_skip=kind=='VNTR',efficient_dictionary=kind=='VNTR',family_display=kind=='VNTR',
        q=q if kind=='VNTR' else None,loss_fraction=fraction if kind=='VNTR' else 0)
    dump(out/'processing.json',processing)
    if kind=='STR':
        from vntrpathstructure_export import export_structure
        folder=out/'str_legacy';folder.mkdir()
        metadata=dict(backbone_path='P0001',processing=processing,input_identity=identity)
        summary=export_structure(folder,locus,records,[],metadata,model_path,
            str_policy_path(model_path,locus.name),0,0,'dot' if render else None)
        stats=summary['graph_stats']['structure']
        comparison=[dict(stage='str_legacy',motifs=len(original_model['motifs']),N_d=stats['nodes'],
            E_d=stats['edges'],C=stats['nodes']+stats['edges'],paths=summary['paths'],molecules=summary['molecules'],
            total_edits_bp=summary['loss_bp'],lossy_paths=summary['lossy_paths'],max_path_loss_fraction=0)]
        dump(out/'comparison.json',comparison);table(out/'comparison.tsv',comparison)
        dump(out/'provenance.json',dict(version=VERSION,locus=locus.name,**processing,
            input_sha256=sha256(input_path.read_bytes()).hexdigest(),full_model_sha256=sha256(model_path.read_bytes()).hexdigest(),
            source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[1]).decode().strip(),
            source_input=str(input_path),source_model=str(model_path),reference=str(reference) if reference else None,
            seconds=round(time.monotonic()-start,3)))
        print('FINISHED',locus.name,json.dumps(comparison),flush=True)
        return comparison
    full,full_hash,features,tokens,slots=annotate(records,model_path,locus)
    dump(out/'full_dictionary_annotations.json',dict(features=features,tokens=tokens,slots=slots))
    choice=select_motifs(full,tokens,q);dump(out/'dictionary_selection.json',choice)
    effective=effective_model(full,choice)
    effective_path=out/'locus_model.efficient.json';dump(effective_path,effective)
    print('DICTIONARY',locus.name,json.dumps(choice,sort_keys=True),flush=True)
    def base_builder(ts,ss):
        return lambda p,current:path_view(current,ss[p],ts['P0001'],ts[p],parameterize=True)
    builder=base_builder(tokens,slots)
    full_rep,full_views,full_sequences,full_losses,_,full_decisions=plan_edits(records,tokens,slots,full,fraction,view_builder=builder)
    baseline=export_budget(out/'baseline_3pct',records,tokens,full_rep,slots,full_views,full_sequences,
        full_losses,fraction,full_decisions,render,full)
    model,model_hash,features,tokens,slots=annotate(records,effective_path,locus)
    dump(out/'selected_annotations.json',dict(tokens=tokens,slots=slots));dump(out/'original_annotations.json',features)
    builder=base_builder(tokens,slots);base_views={p:builder(p,ts) for p,ts in tokens.items()}
    family_zero,rules,zero_report=select_families(base_views,model)
    def family_builder(p,current):return apply_families(builder(p,current),rules)
    rep,_,sequences,losses,_,decisions=plan_edits(records,tokens,slots,model,fraction,view_builder=family_builder)
    expanded={p:builder(p,ts) for p,ts in rep.items()}
    final_views,_,family_report=select_families(expanded,model)
    frozen_views={p:family_builder(p,ts) for p,ts in rep.items()};frozen_metrics=topology(frozen_views)[3]
    if frozen_metrics['C']<family_report['metrics']['C']:
        final_views=frozen_views;family_report=deepcopy(zero_report)
        family_report.update(metrics=frozen_metrics,selection='pre-loss families retained because their final graph is smaller')
    decisions.update(effective_dictionary=choice,family_selection_before_loss=zero_report,
                     family_selection_after_loss=family_report,loss_objective='actual parameterized family N+E')
    expanded_summary=export_budget(out/'expanded_3pct',records,tokens,rep,slots,expanded,sequences,losses,fraction,decisions,render,model)
    final=export_budget(out/'family_3pct',records,tokens,rep,slots,final_views,sequences,losses,fraction,decisions,render,model)
    # The family layer must not change any occurrence, visit payload or loss ledger.
    for name in ['paths.original.fa','paths.represented.fa','losses.jsonl',locus.name+'.original.gfa',locus.name+'.represented.gfa',locus.name+'.path_sequences.txt']:
        assert (out/'expanded_3pct'/name).read_bytes()==(out/'family_3pct'/name).read_bytes(),name
    comparison=[]
    for name,r,size in [('full_dictionary_3pct',baseline,len(full['motifs'])),
                        ('efficient_dictionary_expanded_3pct',expanded_summary,len(model['motifs'])),
                        ('efficient_dictionary_family_3pct',final,len(model['motifs']))]:
        comparison.append(dict(stage=name,motifs=size,**{k:r[k] for k in ['N_d','E_d','C','paths','molecules','total_edits_bp','lossy_paths','max_path_loss_fraction']}))
    dump(out/'comparison.json',comparison);table(out/'comparison.tsv',comparison)
    dump(out/'family_selection.json',family_report)
    if render:
        from vntrfamily_viewer import write_viewer
        write_viewer(out,locus.name)
    dump(out/'provenance.json',dict(version=VERSION,locus=locus.name,**processing,
        input_sha256=sha256(input_path.read_bytes()).hexdigest(),full_model_sha256=full_hash,effective_model_sha256=model_hash,
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[1]).decode().strip(),
        source_input=str(input_path),source_model=str(model_path),reference=str(reference) if reference else None,
        seconds=round(time.monotonic()-start,3)))
    print('FINISHED',locus.name,json.dumps(comparison),flush=True)
    return comparison


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True);parser.add_argument('--model',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--reference',type=Path)
    parser.add_argument('--q',type=float);parser.add_argument('--loss-fraction',type=float)
    parser.add_argument('--supernode-policy',type=Path)
    parser.add_argument('--no-render',action='store_true');args=parser.parse_args()
    if args.supernode_policy and args.loss_fraction is not None:
        parser.error('Supernode loss budget comes from --supernode-policy')
    if args.supernode_policy and args.q is not None and not 0 <= args.q <= 1:
        parser.error('q must be in [0,1]')
    run(args.input,args.model,args.out,args.q,.03 if args.loss_fraction is None else args.loss_fraction,
        args.reference,not args.no_render,args.supernode_policy)


if __name__=='__main__':main()
