"""Full-dictionary supernode bundle, on-disk validation, and provenance."""
from collections import Counter,defaultdict
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import shutil
import subprocess
import time

from vntrdisplay_export import dump,table
from vntrlocal import VERSION as LOCAL_VERSION,assign_slots
from vntrlocal_export import literal_graph,seq_id
from vntrpathgraph import validate_gfa,write_gfa
from vntrpathstructure import annotate_paths,restore_original,unit_key
from vntrsupernode import make_structure
from vntrsupernode_loss import plan_substitutions

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def jsonlines(path,rows):
    with path.open('w',encoding='utf8') as f:
        for row in rows:
            f.write(json.dumps(row,ensure_ascii=False,separators=(',',':'))+'\n')


def load_policy(path):
    policy=read(path)
    rules = dict(policy_id='structure-supernode-v1',annotation_dictionary='full',
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


def validate_bundle(out):
    prov=read(out/'provenance.json')
    for name,h in prov['files'].items():
        if digest(out/name)!=h: raise ValueError('Bundle hash mismatch: '+name)
    display=read(out/'display_graph.json')
    nodes=display['nodes']
    variants=read(out/'sequence_variants.json')
    records={r['id']:r for r in read(out/'source_members.json')}
    paths={row['path']:row['visits'] for row in map(json.loads,(out/'path_visits.jsonl').read_bytes().splitlines())}
    if set(paths)!=set(records): raise ValueError('Path identity mismatch')
    policy=load_policy(out/'supernode_policy.json')
    annotations=read(out/'annotation_audit.json')
    events=defaultdict(list)
    for e in map(json.loads,(out/'losses.jsonl').read_bytes().splitlines()):
        events[e['path']].append(e)
    modes=('original','represented')
    walks,sequences={},{}
    for mode in modes:
        gfa=out/(display['locus']+'.'+mode+'.gfa')
        walks[mode]={}
        for line in gfa.read_text(encoding='ascii').splitlines():
            fields=line.split('\t')
            if fields[0]=='P': walks[mode][fields[1]]=[int(x[:-1]) for x in fields[2].split(',')]
        lines=(out/('paths.'+mode+'.fa')).read_text(encoding='ascii').splitlines()
        sequences[mode]=dict(zip((s[1:] for s in lines[::2]),lines[1::2]))
        validate_gfa(gfa,sequences[mode])
    occurrences={}
    visit_edges=set()
    for p,vs in paths.items():
        es=events.get(p,[])
        by_occ=defaultdict(list)
        for e in es: by_occ[e['occurrence_id']].append(e)
        r=records[p]
        if (r['sequence']!=sequences['original'][p] or
                sha256(r['sequence'].encode()).hexdigest()!=r['sha256']):
            raise ValueError('Original sequence or source hash mismatch '+p)
        if restore_original(sequences['represented'][p],es)!=r['sequence']:
            raise ValueError('Residual mismatch '+p)
        budget=int(len(r['sequence'])*Fraction(str(policy['loss_fraction']))) if r['role']!='reference' else 0
        if len(es)>budget: raise ValueError('Path budget exceeded '+p)
        cursor=0
        for v in vs:
            node=nodes[v['supernode_id']]
            i,j=v['token_interval']
            ann=annotations[p]['tokens'][i:j]
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
                                or v['boundary_assessment']!='invariant_within_model'):
                    raise ValueError('Edited protected unit '+oid)
                occurrences[oid]=(p,t,original,node['id'])
                visit_loss.extend(e['id'] for e in unit_events)
            if visit_loss!=v['loss_event_ids']: raise ValueError('Local loss mismatch '+v['visit_id'])
        if cursor!=len(annotations[p]['tokens']): raise ValueError('Incomplete path '+p)
        visit_edges.update((a['supernode_id'],b['supernode_id']) for a,b in zip(vs,vs[1:]))
        for mode in modes:
            spelled=''.join(variants[k] for v in vs for k in v[mode+'_variant_ids'])
            walked=[k for v in vs for k in v[mode+'_gfa_walk']]
            if walked!=walks[mode][p] or spelled!=sequences[mode][p]:
                raise ValueError('Expansion mismatch '+mode+' '+p)
    if visit_edges!={(e['source'],e['target']) for e in display['edges']}:
        raise ValueError('Display adjacency mismatch')
    for p,es in events.items():
        for e in es:
            rp,rt,donor,rnode=occurrences[e['representative_occurrence']]
            _,t,source,node=occurrences[e['occurrence_id']]
            if (records[rp]['role']=='reference' or unit_key(t)!=unit_key(rt) or rnode!=node
                 or donor!=e['representative_sequence'] or e['op']!='S'
                 or e['supernode_id']!=node or donor[e['ostart']-t['start']]!=e['represented']):
                raise ValueError('Invalid frozen representative '+e['id'])
    return dict(paths=len(paths),occurrences=len(occurrences),original_gfa=True,original_expansion=True,
                represented_gfa_and_expansion=True,dual_walks=True,loss_restore=True,
                structure_and_protected_intervals=True,local_loss_and_donors=True,hashes=True)


def export_supernodes(out,locus,records,model_path,policy_path,render=True,input_path=None):
    start=time.monotonic()
    policy=load_policy(policy_path)
    records=sorted(records,key=lambda r:r['id'])
    model,mhash,features=annotate_paths(records,model_path,locus,decoder='hmm_assessed')
    tokens={p:f['tokens'] for p,f in features.items()}
    slots,audits={},{}
    for p,ts in tokens.items(): slots[p],audits[p]=assign_slots(ts,tokens['P0001'])
    dump(out/'source_members.json',records)
    shutil.copyfile(model_path,out/'locus_model.original.json')
    shutil.copyfile(policy_path,out/'supernode_policy.json')
    display,visits=make_structure(records,tokens,audits,model,mhash,policy,digest(policy_path))
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
        full_dictionary=len(model['motifs']),tokens=nt,folded_tokens=folded,folded_token_fraction=folded/nt,
        visits=sum(map(len,visits.values())),unknown_tokens=sum(t['status']=='unknown' for ts in tokens.values() for t in ts),
        structure_assessments=dict(Counter(t.get('structure_assessment','not_assessed') for ts in tokens.values() for t in ts)),
        budget=policy['loss_fraction'],edited_bp=sum(len(es) for es in events.values()),lossy_paths=len(events),
        max_path_fraction=max(row['fraction'] for row in rows),structure_graph=metric,literal_graphs=metrics,
        stages=[dict(stage='unfolded_original',display=metrics['original'],literal=metrics['original'],edited_bp=0),
                dict(stage='zero_loss_supernodes',display=metric,literal=metrics['original'],edited_bp=0),
                dict(stage='constrained_supernodes',display=metric,literal=metrics['represented'],edited_bp=sum(map(len,events.values())))])
    dump(out/'summary.json',summary)
    if render:
        from vntrsupernode_viewer import write_viewer
        write_viewer(out,display,visits,variants,records,events,policy)
    files={f.name:digest(f) for f in sorted(out.iterdir()) if f.is_file() and f.name not in ('provenance.json','validation.json')}
    commit=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,capture_output=True,text=True)
    sources=[*ROOT.glob('scripts/vntr*.py'),ROOT/'scripts/run_efficient_family.py',ROOT/'scripts/supernodes.html',
             ROOT/'third_party/cytoscape-3.30.4/cytoscape.min.js']
    provenance=dict(locus=locus.name,decoder='hmm_assessed',dictionary='full',dictionary_q=None,policy=policy,
        model_sha256=mhash,policy_sha256=digest(policy_path),boundary_model_hashes=sorted({o['boundary_model_sha256']
            for r in records for o in r['observations'] if o.get('boundary_model_sha256')}),
        source_commit=commit.stdout.strip() if commit.returncode==0 else None,
        source_files={str(f.relative_to(ROOT)):digest(f) for f in sources if f.exists()},
        input_sha256=digest(input_path) if input_path else None,files=files)
    dump(out/'provenance.json',provenance)
    result=validate_bundle(out)
    result.update(seconds=round(time.monotonic()-start,3),hashed_payload_bytes=sum((out/name).stat().st_size for name in files))
    dump(out/'validation.json',result)
    return summary


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description='Validate a saved supernode bundle')
    parser.add_argument('bundle',type=Path)
    print(json.dumps(validate_bundle(parser.parse_args().bundle)))

