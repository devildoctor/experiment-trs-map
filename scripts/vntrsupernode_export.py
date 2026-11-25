"""Supernode bundle, on-disk validation, and annotation model provenance."""
from collections import Counter,defaultdict
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import csv
import shutil
import subprocess
import time
import warnings

from vntrdisplay_export import dump,table
from vntreffmotifs import effective_model,require_original_model,select_motifs
from vntrlocal import VERSION as LOCAL_VERSION,assign_slots
from vntrlocal_export import literal_graph,seq_id
from vntrpathgraph import validate_gfa,write_gfa
from vntrpathstructure import annotate_paths,restore_original,unit_key
from vntrsupernode import make_structure
from vntrsupernode_loss import plan_substitutions

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    checksum=sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):
            checksum.update(block)
    return checksum.hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def jsonlines(path,rows):
    with path.open('w',encoding='utf8') as f:
        for row in rows:
            f.write(json.dumps(row,ensure_ascii=False,separators=(',',':'))+'\n')


def load_policy(path, *, dictionary='full'):
    policy=read(path)
    if 'q' in policy or 'dictionary_q' in policy:
        raise ValueError('Supply q as an argument, not a policy field')
    rules = dict(policy_id='structure-supernode-v1',annotation_dictionary=dictionary,
                 loss_objective='expanded_literal_N_plus_E',allowed_operations=['S'],
                 allow_cross_motif_edits=False,require_assessed_structure_for_loss=True)
    if any(policy.get(k)!=v for k,v in rules.items()):
        raise ValueError('Unsupported structure-supernode policy')
    for k,low,high in [('max_module_width',1,12),('max_unit_substitutions',0,2)]:
        if type(policy.get(k))!=int or not low<=policy[k]<=high:
            raise ValueError('Invalid '+k)
    if not 0<=policy['loss_fraction']<=.03:
        raise ValueError('loss_fraction must be between 0 and 0.03')
    return policy


def unique_rows(rows,key,what):
    result={}
    for row in rows:
        identity=key(row)
        if identity in result: raise ValueError('Duplicate '+what+': '+str(identity))
        result[identity]=row
    return result


def path_budget(record,policy):
    return 0 if record['role']=='reference' else int(len(record['sequence'])*Fraction(str(policy['loss_fraction'])))


def verify_support(row,paths,records):
    ids=set(paths)
    molecules={m for p in ids if records[p]['role']!='reference' for m in records[p]['molecules']}
    if (sorted(row['path_ids'])!=sorted(ids) or row['path_count']!=len(ids)
            or row['visit_count']!=len(paths) or row['molecule_count']!=len(molecules)):
        raise ValueError('Incorrect path, visit or molecule support')


def validate_dictionary(out, prov):
    files=prov['files']
    original=read(out/'locus_model.original.json')
    require_original_model(original)
    original_hash=files['locus_model.original.json']
    if prov.get('full_model_sha256',original_hash)!=original_hash:
        raise ValueError('Full model hash mismatch')
    mode=prov['dictionary']
    q=prov['dictionary_q']
    policy=load_policy(out/'supernode_policy.json',dictionary=mode)
    if prov['policy']!=policy or prov['policy_sha256']!=files['supernode_policy.json']:
        raise ValueError('Resolved policy mismatch')
    if prov['decoder']!='hmm_assessed':
        raise ValueError('Expected assessed supernode annotation')
    if mode=='full':
        if q is not None or prov['model_sha256']!=original_hash or 'dictionary_selection.json' in files:
            raise ValueError('Invalid full dictionary provenance')
        return policy,original
    if mode!='effective' or type(q) not in (float,int) or not 0 <= q <= 1:
        raise ValueError('Invalid effective dictionary provenance')
    if not {'dictionary_selection.json','locus_model.efficient.json'} <= files.keys():
        raise ValueError('Incomplete effective dictionary manifest')
    choice=read(out/'dictionary_selection.json')
    model=read(out/'locus_model.efficient.json')
    if prov.get('bundle_version',1)>=3 or 'retained_interruptions' in choice:
        frozen=set(original.get('interruptions',[])).intersection(choice['original_ids'])
        if set(choice.get('retained_interruptions',[]))!=frozen or any(choice['mapping'].get(m)!=m for m in frozen):
            raise ValueError('Frozen interruptions missing or remapped in effective dictionary')
    if (choice['q']!=q or prov.get('full_model_sha256')!=original_hash
            or prov['model_sha256']!=files['locus_model.efficient.json']
            or prov.get('effective_model_sha256')!=prov['model_sha256']
            or choice['original_ids']!=sorted(m['id'] for m in original['motifs'])
            or sorted(m['id'] for m in model['motifs'])!=choice['selected_ids']
            or model!=effective_model(original,choice)):
        raise ValueError('Effective dictionary selection or hash mismatch')
    return policy,model


def validate_bundle(out,source_root=None):
    prov=read(out/'provenance.json')
    modes=('original','represented')
    required={'display_graph.json','sequence_variants.json','source_members.json','path_visits.jsonl',
              'supernode_policy.json','annotation_audit.json','losses.jsonl','paths.tsv','locus_model.original.json'}
    if prov.get('bundle_version',1)>=2:
        required.update({'input_manifest.json',prov['locus']+'.path_sequences.txt'})
        if prov.get('input_manifest_sha256')!=prov['files'].get('input_manifest.json'):
            raise ValueError('Input manifest digest mismatch')
        inp=read(out/'input_manifest.json')
        expected=inp.get('file_sha256') if prov.get('input_sha256_scope')=='input_file' else prov['input_manifest_sha256']
        if not expected or prov.get('input_sha256')!=expected:
            raise ValueError('Input identity digest mismatch')
        if 'scripts/supernode_projection.js' not in prov['source_files']:
            raise ValueError('Missing projection source hash')
    if source_root is not None:
        for name,expected in prov['source_files'].items():
            path=source_root/name
            if not path.is_file() or digest(path)!=expected:
                raise ValueError('Source file hash mismatch: '+name)
    for mode in modes:
        required.update((prov['locus']+'.'+mode+'.gfa','paths.'+mode+'.fa'))
    if not required.issubset(prov['files']): raise ValueError('Incomplete bundle manifest')
    for name,h in prov['files'].items():
        if digest(out/name)!=h: raise ValueError('Bundle hash mismatch: '+name)
    display=read(out/'display_graph.json')
    if display['locus']!=prov['locus']: raise ValueError('Locus mismatch')
    nodes=display['nodes']
    variants=read(out/'sequence_variants.json')
    records=unique_rows(read(out/'source_members.json'),lambda r:r['id'],'source path')
    if prov.get('bundle_version',1)>=2:
        with (out/(prov['locus']+'.path_sequences.txt')).open(encoding='utf-8') as handle:
            exported=unique_rows(csv.DictReader((line for line in handle if not line.startswith('#')),delimiter='\t'),
                                 lambda r:r['name'],'path in TXT')
        if set(exported)!=set(records):
            raise ValueError('TXT path names mismatch')
        for p,r in records.items():
            row=exported[p]
            if (row['role'],row['sequence'],int(row['length_bp']),int(row['molecules']))!=(
                    r['role'],r['sequence'],len(r['sequence']),len(set(r['molecules']))):
                raise ValueError('TXT path content mismatch: '+p)
    rows=unique_rows(map(json.loads,(out/'path_visits.jsonl').read_bytes().splitlines()),lambda r:r['path'],'path')
    paths={p:row['visits'] for p,row in rows.items()}
    unique_rows((v for vs in paths.values() for v in vs),lambda v:v['visit_id'],'visit')
    if set(paths)!=set(records): raise ValueError('Path identity mismatch')
    policy,annotation_model=validate_dictionary(out,prov)
    motif_ids={m['id'] for m in annotation_model['motifs']+annotation_model.get('terminal_templates',[])}
    annotations=read(out/'annotation_audit.json')
    if any(t['motif'] not in motif_ids for row in annotations.values()
           for t in row['tokens'] if t['status']=='unit'):
        raise ValueError('Annotated unit outside the recorded dictionary')
    events=defaultdict(list)
    ledger=unique_rows(map(json.loads,(out/'losses.jsonl').read_bytes().splitlines()),lambda e:e['id'],'event')
    for e in ledger.values():
        if e['path'] not in records: raise ValueError('Unknown event path')
        events[e['path']].append(e)
    walks,sequences,segments={},{},{}
    for mode in modes:
        gfa=out/(display['locus']+'.'+mode+'.gfa')
        walks[mode]={}
        segments[mode]={}
        for line in gfa.read_text(encoding='ascii').splitlines():
            fields=line.split('\t')
            if fields[0]=='S': segments[mode][int(fields[1])]=fields[2]
            if fields[0]=='P': walks[mode][fields[1]]=[int(x[:-1]) for x in fields[2].split(',')]
        lines=(out/('paths.'+mode+'.fa')).read_text(encoding='ascii').splitlines()
        pairs=unique_rows(zip(lines[::2],lines[1::2]),lambda pair:pair[0][1:],'FASTA path')
        sequences[mode]={p:pair[1] for p,pair in pairs.items()}
        if set(sequences[mode])!=set(paths): raise ValueError('FASTA path identities differ')
        validate_gfa(gfa,sequences[mode])
    occurrences={}
    node_visits=defaultdict(list)
    edge_paths=defaultdict(list)
    for p,vs in paths.items():
        es=events.get(p,[])
        by_occ=defaultdict(list)
        for e in es: by_occ[e['occurrence_id']].append(e)
        r=records[p]
        if (r['sequence']!=sequences['original'][p] or
                sha256(r['sequence'].encode()).hexdigest()!=r['sha256']):
            raise ValueError('Original sequence or source hash mismatch '+p)
        seen_positions=set()
        for e in es:
            a=e['ostart']
            if (e['op']!='S' or not 0<=a<len(r['sequence']) or a in seen_positions
                or (e['oend'],e['rstart'],e['rend'])!=(a+1,a,a+1)
                or e['original']!=r['sequence'][a] or e['represented']!=sequences['represented'][p][a:a+1]
                or e['original']==e['represented']):
                raise ValueError('Invalid single-base substitution '+e['id'])
            seen_positions.add(a)
        if restore_original(sequences['represented'][p],es)!=r['sequence']:
            raise ValueError('Residual mismatch '+p)
        budget=path_budget(r,policy)
        if len(es)>budget: raise ValueError('Path budget exceeded '+p)
        cursor=0
        for v in vs:
            node=nodes[v['supernode_id']]
            i,j=v['token_interval']
            ann=annotations[p]['tokens'][i:j]
            if not (v['path_id']==p and 0<=i<j<=len(annotations[p]['tokens'])):
                raise ValueError('Invalid visit coordinates or owner '+v['visit_id'])
            interval=[ann[0]['start'],ann[-1]['end']]
            for mode in modes:
                walk=v[mode+'_gfa_walk']
                ids=v[mode+'_variant_ids']
                if (v[mode+'_interval']!=interval or walk!=walks[mode][p][i:j] or len(ids)!=j-i
                    or any(segments[mode][n]!=variants[k] for n,k in zip(walk,ids))):
                    raise ValueError('Invalid local '+mode+' expansion '+v['visit_id'])
            keys=[list(unit_key(t)) for t in ann]
            if (cursor!=i or j<=i or node['ordered_pattern']*v['repeat_count']!=keys
                or v['ordered_occurrence_ids']!=[f'{p}:u{k}' for k in range(i,j)]
                or v['visit_id'] not in node['member_visit_ids']):
                raise ValueError('Invalid structural coverage '+v['visit_id'])
            cursor=j
            visit_loss=[]
            for k,oid in enumerate(v['ordered_occurrence_ids']):
                original,rep=[variants[v[mode+'_variant_ids'][k]] for mode in modes]
                t=ann[k]
                if t['exact_bases']!=original or len(original)!=len(rep):
                    raise ValueError('Frozen occurrence mismatch '+oid)
                changes=[a for a,(x,y) in enumerate(zip(original,rep)) if x!=y]
                unit_events=by_occ[oid]
                if (changes!=[e['ostart']-t['start'] for e in unit_events]
                        or len(changes)>policy['max_unit_substitutions']):
                    raise ValueError('Unit substitutions mismatch '+oid)
                if changes and (node['protected'] or v['structure_assessment']!='invariant_within_model'
                                or v['boundary_assessment']!='invariant_within_model' or v['position_assessment']!='aligned'):
                    raise ValueError('Edited protected unit '+oid)
                slot=annotations[p]['slot_audit']['selected_slots'][i+k]
                occurrences[oid]=(p,t,original,rep,node['id'],slot)
                visit_loss.extend(e['id'] for e in unit_events)
            if visit_loss!=v['loss_event_ids']: raise ValueError('Local loss mismatch '+v['visit_id'])
            node_visits[v['supernode_id']].append(v)
        if cursor!=len(annotations[p]['tokens']): raise ValueError('Incomplete path '+p)
        for a,b in zip(vs,vs[1:]): edge_paths[a['supernode_id'],b['supernode_id']].append(p)
        for mode in modes:
            spelled=''.join(variants[k] for v in vs for k in v[mode+'_variant_ids'])
            walked=[k for v in vs for k in v[mode+'_gfa_walk']]
            if walked!=walks[mode][p] or spelled!=sequences[mode][p]:
                raise ValueError('Expansion mismatch '+mode+' '+p)
    if set(nodes)!=set(node_visits): raise ValueError('Display node identities differ')
    for n,vs in node_visits.items():
        node=nodes[n]
        identity=[display['locus'],prov['model_sha256'],prov['policy_sha256'],
                  node['region_id'],node['ordered_pattern'],node['protected']]
        expected_id='SN_'+sha256(json.dumps(identity,separators=(',',':')).encode()).hexdigest()[:24]
        if n!=expected_id:
            raise ValueError('Node identity differs from the recorded model or policy')
        if (n!=node['id'] or sorted(node['member_visit_ids'])!=sorted(v['visit_id'] for v in vs)
            or sorted(node['observed_repeat_counts'])!=sorted({v['repeat_count'] for v in vs})
            or sorted(node.get('lossy_visit_ids',[]))!=sorted(v['visit_id'] for v in vs if v['loss_event_ids'])):
            raise ValueError('Incorrect node membership '+n)
        verify_support(node,[v['path_id'] for v in vs],records)
    edges=unique_rows(display['edges'],lambda e:(e['source'],e['target']),'edge')
    if set(edges)!=set(edge_paths): raise ValueError('Display adjacency mismatch')
    for edge,ps in edge_paths.items(): verify_support(edges[edge],ps,records)
    for e in ledger.values():
        if e['occurrence_id'] not in occurrences or e['representative_occurrence'] not in occurrences:
            raise ValueError('Unknown occurrence in event '+e['id'])
        rp,rt,donor,_,rnode,rslot=occurrences[e['representative_occurrence']]
        owner,t,_,rep,node,slot=occurrences[e['occurrence_id']]
        if (owner!=e['path'] or records[rp]['role']=='reference' or unit_key(t)!=unit_key(rt) or rnode!=node
            or rslot!=slot or donor!=rep or donor!=e['representative_sequence'] or e['supernode_id']!=node):
            raise ValueError('Invalid frozen representative '+e['id'])
    return dict(paths=len(paths),occurrences=len(occurrences),original_gfa=True,original_expansion=True,
                represented_gfa_and_expansion=True,dual_walks=True,loss_restore=True,
                structure_and_protected_intervals=True,local_loss_and_donors=True,local_walk_slices=True,
                node_and_edge_support=True,manifest_coverage=True,hashes=True)


def dictionary_q_for_size(size):
    if size >= 9:
        return 0.2
    if size >= 4:
        return 0.1
    return None


def export_supernodes(out,locus,records,model_path,policy_path,render=True,input_path=None,*,q=None,metadata=None):
    start=time.monotonic()
    policy=load_policy(policy_path)
    records=sorted(records,key=lambda r:r['id'])
    metadata=metadata or {}
    original=read(model_path)
    require_original_model(original)
    standard_q=dictionary_q_for_size(len(original['motifs']))
    if q is not None:
        if type(q) not in (int,float) or not 0 <= q <= 1:
            raise ValueError('q must be in [0,1]')
        if q != standard_q:
            warnings.warn(f'VNTR dictionary size determines q. Requested {q}, using {standard_q!r}.', stacklevel=2)
    q=standard_q
    resolved_policy=out/'supernode_policy.json'
    annotation_model=model_path
    if q is None:
        shutil.copyfile(policy_path,resolved_policy)
    else:
        _,_,full_features=annotate_paths(records,model_path,locus,decoder='hmm_single')
        # select_motifs counts only status and motif, before slot assignment or edits.
        choice=select_motifs(original,{p:f['tokens'] for p,f in full_features.items()},q)
        annotation_model=out/'locus_model.efficient.json'
        dump(annotation_model,effective_model(original,choice))
        dump(out/'dictionary_selection.json',choice)
        policy['annotation_dictionary']='effective'
        dump(resolved_policy,policy)
    model,mhash,features=annotate_paths(records,annotation_model,locus,decoder='hmm_assessed')
    tokens={p:f['tokens'] for p,f in features.items()}
    slots,audits={},{}
    for p,ts in tokens.items(): slots[p],audits[p]=assign_slots(ts,tokens['P0001'])
    dump(out/'source_members.json',records)
    shutil.copyfile(model_path,out/'locus_model.original.json')
    display,visits=make_structure(records,tokens,audits,model,mhash,policy,digest(resolved_policy))
    rep,events,decisions=plan_substitutions(records,tokens,slots,display,visits,policy)
    graphs={}; sequences={}; variants={}
    ledger=[dict(path=p,substitutions=len(es),loss_bp=len(es)) for p,es in events.items()]
    for mode,ts in [('original',tokens),('represented',rep)]:
        g=literal_graph(records,ts,slots)
        graphs[mode]=g
        seq={p:''.join(t['exact_bases'] for t in row) for p,row in ts.items()}
        sequences[mode]=seq
        write_gfa(out/(locus.name+'.'+mode+'.gfa'),g,records,(LOCAL_VERSION if mode=='original' else 'supernodes_v1')+'_'+mode,seq,ledger if mode=='represented' else [])
        (out/('paths.'+mode+'.fa')).write_text(''.join('>'+p+'\n'+s+'\n' for p,s in seq.items()),encoding='ascii')
        for p,vs in visits.items():
            for v in vs:
                i,j=v['token_interval']
                strings=[t['exact_bases'] for t in ts[p][i:j]]
                variants.update((seq_id(s),s) for s in strings)
                v[mode+'_variant_ids']=[seq_id(s) for s in strings]
                v[mode+'_gfa_walk']=g['paths'][p][i:j]
    observed={r['sequence'] for r in records if r['role']!='reference'}
    rows=[]
    for r in records:
        p=r['id']
        es=events.get(p,[])
        by_occ=defaultdict(list)
        for e in es: by_occ[e['occurrence_id']].append(e['id'])
        for v in visits[p]:
            v['loss_event_ids']=[e for oid in v['ordered_occurrence_ids'] for e in by_occ[oid]]
            if v['loss_event_ids']:
                display['nodes'][v['supernode_id']].setdefault('lossy_visit_ids',[]).append(v['visit_id'])
        rows.append(dict(path=p,role=r['role'],molecules=len(r['molecules']),length=len(r['sequence']),
                         edited_bp=len(es),fraction=len(es)/len(r['sequence']),
                         represented_is_observed_whole_sequence=sequences['represented'][p] in observed))
    table(out/'paths.tsv',rows)
    jsonlines(out/'losses.jsonl',(e for p,es in sorted(events.items()) for e in es))
    jsonlines(out/'path_visits.jsonl',(dict(path=p,visits=vs) for p,vs in sorted(visits.items())))
    dump(out/'sequence_variants.json',variants)
    dump(out/'display_graph.json',display)
    dump(out/'edit_decisions.json',decisions)
    (out/'annotation_audit.json').write_text(json.dumps({p:dict(tokens=ts,slot_audit=audits[p],
        annotation=features[p]['annotation']) for p,ts in tokens.items()},separators=(',',':')),encoding='utf8')
    metrics={mode:dict(nodes=len(g['nodes']),edges=len(g['edges'])) for mode,g in graphs.items()}
    if sum(metrics['original'].values())-sum(metrics['represented'].values())!=sum(d['delta_sequence_complexity'] for d in decisions):
        raise ValueError('Literal graph complexity accounting mismatch')
    metric=dict(nodes=len(display['nodes']),edges=len(display['edges']))
    nt=sum(map(len,tokens.values()))
    folded=sum(v['token_interval'][1]-v['token_interval'][0] for vs in visits.values() for v in vs if v['repeat_count']>1)
    summary=dict(locus=locus.name,paths=len(records),molecules=len({m for r in records if r['role']!='reference' for m in r['molecules']}),
        full_dictionary=len(original['motifs']),tokens=nt,folded_tokens=folded,folded_token_fraction=folded/nt,
        visits=sum(map(len,visits.values())),unknown_tokens=sum(t['status']=='unknown' for ts in tokens.values() for t in ts),
        structure_assessments=dict(Counter(t.get('structure_assessment','not_assessed') for ts in tokens.values() for t in ts)),
        budget=policy['loss_fraction'],edited_bp=sum(len(es) for es in events.values()),lossy_paths=len(events),
        max_path_fraction=max(row['fraction'] for row in rows),structure_graph=metric,literal_graphs=metrics,
        stages=[dict(stage='unfolded_original',display=metrics['original'],literal=metrics['original'],edited_bp=0),
                dict(stage='zero_loss_supernodes',display=metric,literal=metrics['original'],edited_bp=0),
                dict(stage='constrained_supernodes',display=metric,literal=metrics['represented'],edited_bp=sum(map(len,events.values())))])
    if q is not None:
        summary['effective_dictionary']=len(model['motifs'])
    from vntrpathstructure_export import write_path_sequences_txt
    write_path_sequences_txt(out/(locus.name+'.path_sequences.txt'),records,metadata)
    if metadata.get('core_length_evidence') is not None:
        summary['core_length_evidence']=metadata['core_length_evidence']
    if 'input_identity' in metadata:
        summary['input_identity']=metadata['input_identity']
    dump(out/'summary.json',summary)
    if render:
        from vntrsupernode_viewer import write_viewer
        write_viewer(out,display,visits,variants,records,events,policy,rows,display_model=model)
    manifest={k:metadata[k] for k in ('reference','reference_core_sha256','catalog_sha256','input_files',
              'extraction_version','boundary_model_sha256','input_identity','admission','bam_count','bam_audits',
              'reference_index_sha256','sources_sha256','previous_members_sha256') if k in metadata}
    if input_path:
        manifest.update(kind=metadata.get('input_kind','members'),file=str(input_path),file_sha256=digest(input_path))
    elif 'input_files' in manifest:
        manifest['kind']='bam'
    else:
        manifest.update(kind='in_memory',members_sha256=digest(out/'source_members.json'))
    dump(out/'input_manifest.json',manifest)
    files={f.name:digest(f) for f in sorted(out.iterdir()) if f.is_file() and f.name not in ('provenance.json','validation.json')}
    commit=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,capture_output=True,text=True)
    sources=[*ROOT.glob('scripts/vntr*.py'),ROOT/'scripts/run_efficient_family.py',ROOT/'scripts/supernodes.html',ROOT/'scripts/supernode_compare.js',
             ROOT/'scripts/supernode_projection.js',ROOT/'third_party/cytoscape-3.30.4/cytoscape.min.js']
    bundle_version=3 if q is not None and choice.get('retained_interruptions') else 2
    provenance=dict(bundle_version=bundle_version,locus=locus.name,decoder='hmm_assessed',dictionary=policy['annotation_dictionary'],dictionary_q=q,policy=policy,
        full_model_sha256=digest(model_path),model_sha256=mhash,policy_sha256=digest(resolved_policy),boundary_model_hashes=sorted({o['boundary_model_sha256']
            for r in records for o in r['observations'] if o.get('boundary_model_sha256')}),
        source_commit=commit.stdout.strip() if commit.returncode==0 else None,
        source_files={f.relative_to(ROOT).as_posix():digest(f) for f in sources},
        input_manifest_sha256=files['input_manifest.json'],parameters=metadata.get('parameters',{}),
        input_sha256_scope='input_file' if input_path else 'input_manifest',
        input_sha256=manifest.get('file_sha256',files['input_manifest.json']),files=files)
    if q is not None:
        provenance['effective_model_sha256']=mhash
    dump(out/'provenance.json',provenance)
    result=validate_bundle(out)
    result.update(seconds=round(time.monotonic()-start,3),hashed_payload_bytes=sum((out/name).stat().st_size for name in files))
    dump(out/'validation.json',result)
    return summary


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description='Validate a saved supernode bundle')
    parser.add_argument('bundle',type=Path)
    parser.add_argument('--source-root',type=Path,help='Additionally compare the recorded source file hashes with a checkout')
    args=parser.parse_args()
    print(json.dumps(validate_bundle(args.bundle,args.source_root)))

