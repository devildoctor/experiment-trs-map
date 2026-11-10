"""Auditable local graph artifacts; variable nodes are never literal GFA S records."""
from collections import Counter,defaultdict
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import csv,json,subprocess
from vntrpathgraph import validate_graph,write_gfa,validate_gfa
from vntrdisplay import add_support
from vntrdisplay_export import dump,table
from vntrlocal import topology,category,VERSION
from vntrcore import edit_distance


def literal_graph(records,tokens,slots):
    sequences={};walks={}
    for p,ts in sorted(tokens.items()):
        keys=[(slots[p][i],t['exact_bases']) for i,t in enumerate(ts)]
        walks[p]=keys
        for key in keys:sequences[key]=key[1]
    ids={k:i+1 for i,k in enumerate(sorted(sequences))}
    paths={p:[ids[k] for k in w] for p,w in walks.items()}
    result=dict(nodes={ids[k]:s for k,s in sequences.items()},paths=paths,
        edges={e for w in paths.values() for e in zip(w,w[1:])},backbone_path='P0001',
        coordinates={ids[k]:[list(k[0])] for k in ids})
    validate_graph(result,{p:''.join(t['exact_bases'] for t in ts) for p,ts in tokens.items()})
    return result


def seq_id(s):return 'V'+sha256(s.encode()).hexdigest()[:16]


def annotate_unknown_hints(display,model):
    """Keep nearest-template hints distinct from primary-optimal annotation choices."""
    hints={}
    for visits in display['paths'].values():
        for visit in visits:
            for part in visit['parts']:
                for variant,a in zip(part['variant_ids'],part['annotation']):
                    if a['status']!='unknown':continue
                    if variant not in hints:
                        sequence=display['sequence_variants'][variant]
                        distances={m['id']:edit_distance(sequence,m['sequence']) for m in model['motifs']}
                        best=min(distances.values())
                        hints[variant]=dict(minimum_edit_distance=best,nearest_types=sorted(k for k,v in distances.items() if v==best),
                            distances=distances,scope='fixed original fragment versus frozen forward templates; hint, not a selected annotation')
                    a['unknown_template_hint']=hints[variant]
    display['unknown_template_hints']=hints
    return hints


def materialize(records,views,graph,losses,original_tokens):
    chains,mapping,paths,metrics=topology(views);bykey=defaultdict(list);by_id={r['id']:r for r in records}
    variants={};variant_members=defaultdict(set);transition_members=defaultdict(set)
    for p,vs in views.items():
        previous=None
        for v in vs:
            bykey[v['key']].append((p,v))
            for s in v['variant_vector']:variants[seq_id(s)]=s
            vv=tuple(seq_id(s) for s in v['variant_vector'])
            variant_members[v['key'],vv].update(by_id[p]['molecules'])
            if previous is not None:transition_members[previous,(v['key'],vv)].update(by_id[p]['molecules'])
            previous=(v['key'],vv)
    atoms={}
    for key,members in bykey.items():
        v=members[0][1];ts=[t for p,w in members for t in w['tokens']]
        atoms[key]=dict(category=list(v['category']),position=v['position'],
            types=sorted({t['motif'] for t in ts}),counts=sorted({w['count'] for p,w in members}),
            variant_ids=sorted({seq_id(s) for p,w in members for s in w['variant_vector']}),
            lengths=sorted({sum(map(len,w['variant_vector'])) for p,w in members}),
            ambiguity_kinds=sorted({t['ambiguity_kind'] for t in ts}),
            candidate_types=sorted({x for t in ts for x in t['candidate_types']}),
            parameterized_weak=v.get('parameterized_weak',False),
            approximate=any(t['match_kind']=='approximate' for t in ts),
            weak_approximate=any(t['weak_approximate'] for t in ts))
    ids={i:'D%04d'%(i+1) for i in range(len(chains))};nodes={}
    for i,chain in enumerate(chains):
        parts=[atoms[k] for k in chain];flags=[]
        if any(p['category'][0]=='unknown' for p in parts):flags.append('unknown_content')
        if any(p['ambiguity_kinds']!=['none'] for p in parts):flags.append('annotation_ambiguity')
        if any(p['weak_approximate'] for p in parts):flags.append('weak_approximate')
        nodes[ids[i]]=dict(id=ids[i],kind='module' if len(chain)>1 else parts[0]['category'][0],
            parts=parts,uncertainty_flags=flags)
    path_visits={};occurrences=[];local_losses=defaultdict(set)
    for p,vs in views.items():
        out=[];position=0
        for v in vs:
            n=ids[mapping[v['key']]];ns=[graph['paths'][p][i] for i in v['occurrence_indices']]
            length=sum(map(len,v['variant_vector']))
            if any(original_tokens[p][i]['exact_bases']!=s for i,s in zip(v['occurrence_indices'],v['variant_vector'])):
                local_losses[n].add(p)
            part=dict(count=v['count'],variant_ids=[seq_id(s) for s in v['variant_vector']],
                types=[t['motif'] for t in v['tokens']],occurrence_indices=v['occurrence_indices'],
                original_interval=v['original_interval'],represented_interval=[position,position+length],
                gfa_nodes=ns,annotation=[dict(motif=t['motif'],status=t['status'],
                    match_kind=t['match_kind'],ambiguity_kind=t['ambiguity_kind'],candidate_types=t['candidate_types'],
                    ambiguity_regions=t['ambiguity_regions'],edit_distance=t.get('edit_distance',0),
                    weak_approximate=t['weak_approximate']) for t in v['tokens']])
            if v.get('skip_context') is not None:part['skip_context']=v['skip_context']
            if out and out[-1]['node']==n:
                out[-1]['parts'].append(part);out[-1]['gfa_nodes']+=ns
            else:out.append(dict(node=n,parts=[part],gfa_nodes=list(ns)))
            for i,t in zip(v['occurrence_indices'],v['tokens']):
                occurrences.append(dict(path=p,index=i,node=n,slot=list(v['slots'][i-v['occurrence_indices'][0]]),
                    original_start=t['start'],original_end=t['end'],represented_variant=seq_id(t['exact_bases']),
                    motif=t['motif'],status=t['status'],match_kind=t['match_kind'],ambiguity_kind=t['ambiguity_kind'],
                    candidate_types=t['candidate_types'],ambiguity_regions=t['ambiguity_regions'],
                    weak_approximate=t['weak_approximate']))
            position+=length
        path_visits[p]=out
    edges={(a['node'],b['node']) for vs in path_visits.values() for a,b in zip(vs,vs[1:])}
    display=dict(schema_version=3,version=VERSION,nodes=nodes,paths=path_visits,
        edges=[dict(source=a,target=b) for a,b in sorted(edges)],metrics=metrics,
        sequence_variants=variants,backbone_path='P0001',
        scope='All registered P paths with their own ordered visit parameters; topology alone does not enumerate observations.')
    groups={p:p for p in views};add_support(display,graph,records,losses,groups)
    # Compare the actual strings of each occurrence, including boundary deletions.
    # A loss elsewhere on the same P must not color this node.
    for n,node in display['nodes'].items():node['lossy_member_ids']=sorted(local_losses[n])
    display['loss_marker_scope']='node_local_string_change'
    display['weak_skip_parameterized']=any(v.get('parameterized_weak',False) for vs in views.values() for v in vs)
    variant_rows=[dict(node=ids[mapping[k]],family_position=repr(k),variant_vector=list(v),
        route_molecules=len(ms),molecule_ids=sorted(ms)) for (k,v),ms in variant_members.items()]
    transition_rows=[dict(source=ids[mapping[a[0]]],target=ids[mapping[b[0]]],source_variants=list(a[1]),
        target_variants=list(b[1]),route_molecules=len(ms),molecule_ids=sorted(ms)) for (a,b),ms in transition_members.items()]
    return display,occurrences,variant_rows,transition_rows


def validate_parameters(display,graph,represented):
    for p,vs in display['paths'].items():
        expanded=''.join(display['sequence_variants'][v] for visit in vs for part in visit['parts'] for v in part['variant_ids'])
        assert expanded==represented[p],p
        assert [n for visit in vs for n in visit['gfa_nodes']]==graph['paths'][p],p
        for visit in vs:
            specs=display['nodes'][visit['node']]['parts']
            assert len(specs)==len(visit['parts'])
            for spec,part in zip(specs,visit['parts']):
                assert len(part['variant_ids'])==part['count'] and part['count'] in spec['counts']
                assert all(x in spec['variant_ids'] for x in part['variant_ids'])
    return dict(paths=len(represented),parameter_sequence_reconstruction=True,gfa_walk_reconstruction=True)


def atom_label(part):
    cat=part['category'][0];counts=part['counts']
    values=lambda ns: str(ns[0]) if len(ns)==1 else '{'+','.join(map(str,ns))+'}'
    if cat=='unknown':label='X? '+values(part['lengths'])+'bp'
    elif cat=='partial':label='PART '+','.join(part['types'])+' '+values(part['lengths'])+'bp'
    else:label='|'.join(part['types'])+' x'+values(counts)
    if cat in ('type_ambiguous','partition_ambiguous') or part['ambiguity_kinds']!=['none']:label+=' ?'
    if part['weak_approximate'] and not part.get('parameterized_weak',False):label+=' weak'
    elif part['approximate']:label+=' ~'
    if len(part['variant_ids'])>1:label+=' v'+str(len(part['variant_ids']))
    return label


def render_graph(out,display,dot='dot'):
    lines=['digraph D1S80 {','graph [rankdir=LR,bgcolor="white",nodesep=0.10,ranksep=0.45,splines=polyline,outputorder=edgesfirst,nslimit=2,nslimit1=2,mclimit=0.2,searchsize=5];',
        'node [shape=box,style="rounded,filled",fontname="DejaVu Sans",fontsize=10];',
        'edge [color="#a5b4c7",arrowsize=0.35,penwidth=0.6,fontsize=7];']
    legend=('D1S80 | All 307 P paths (1 REF + 306 OBS) | Per-P loss cap '+str(round(display.get('loss_budget_fraction',0)*100))+'%'
        '\nNODE: unit/run or linear module; x{...}=observed count set; v=exact sequence variants; X?=unresolved fragment (actual bp retained)'
        '\n?=annotation ambiguity; ~=approximate match (per-visit weak/distance in parameters); orange LOSS*=edits inside this node; blue border=REF path'
        '\nEDGE: registered adjacency; P=path count; M=distinct molecule provenance keys. Only registered paths plus their visit parameters are observations.')
    lines.append('graph [labelloc="t",labeljust="l",fontname="DejaVu Sans",fontsize=14,label='+json.dumps(legend)+'];')
    ref={v['node'] for v in display['paths']['P0001']}
    for n,node in display['nodes'].items():
        label='\n'.join(atom_label(x) for x in node['parts'])+'\n'+n+' | P'+str(len(node['path_ids']))+' M'+str(node['route_molecules'])
        if node['lossy_member_ids']:label+=' LOSS*'
        flags=node['uncertainty_flags']
        fill='#e9e5f6' if 'annotation_ambiguity' in flags else '#eceff2' if 'unknown_content' in flags else '#e8f4f1'
        if node['lossy_member_ids']:fill='#fff0cc'
        color='#2563eb' if n in ref else '#64748b'
        tooltip=dict(node=n,paths=node['path_ids'],lengths=node['represented_lengths'],parts=node['parts'])
        lines.append('%s [label=%s,fillcolor="%s",color="%s",penwidth=%s,tooltip=%s];'%(n,json.dumps(label),fill,color,2 if n in ref else 1,json.dumps(json.dumps(tooltip))))
    for e in display['edges']:
        lines.append('%s -> %s [label="P%d M%d",tooltip=%s];'%(e['source'],e['target'],len(e['path_ids']),e['route_molecules'],json.dumps(','.join(e['path_ids']))))
    lines.append('}');f=out/'D1S80.display.dot';f.write_text('\n'.join(lines)+'\n',encoding='ascii')
    subprocess.run([dot,'-Tsvg',str(f),'-o',str(out/'D1S80.display.svg')],check=True,timeout=300)


def export_budget(out,records,original_tokens,rep_tokens,slots,views,represented,losses,fraction,decisions,render=True,model=None):
    out.mkdir(parents=True,exist_ok=True);originals={r['id']:r['sequence'] for r in records}
    exact=literal_graph(records,original_tokens,slots);graph=literal_graph(records,rep_tokens,slots)
    display,occurrences,variants,transitions=materialize(records,views,graph,losses,original_tokens)
    display['loss_budget_fraction']=fraction
    if model:annotate_unknown_hints(display,model)
    validation=validate_parameters(display,graph,represented)
    for name,g,seqs,loss in [('original',exact,originals,[]),('represented',graph,represented,losses)]:
        f=out/('D1S80.'+name+'.gfa');write_gfa(f,g,records,VERSION+'_'+name,seqs,loss)
        validation[name+'_gfa']=validate_gfa(f,seqs)
        (out/('paths.'+name+'.fa')).write_text(''.join('>'+p+'\n'+s+'\n' for p,s in sorted(seqs.items())),encoding='ascii')
    loss_by_id={v['path']:v for v in losses};observed={r['sequence'] for r in records if r.get('role')!='reference'};rows=[]
    for r in sorted(records,key=lambda x:x['id']):
        p=r['id'];loss=loss_by_id.get(p,{})
        rows.append(dict(path=p,role=r.get('role','observation'),molecules=len(r['molecules']),
            original_length=len(originals[p]),represented_length=len(represented[p]),loss_bp=loss.get('loss_bp',0),
            loss_fraction=loss.get('fraction',0),substitutions=loss.get('substitutions',0),
            insertions=loss.get('insertions',0),deletions=loss.get('deletions',0),
            represented_is_observed_whole_sequence=represented[p] in observed,
            original_sequence=originals[p],represented_sequence=represented[p]))
    table(out/'D1S80.path_sequences.txt',rows);table(out/'paths.tsv',rows)
    with (out/'D1S80.path_sequences.txt').open(encoding='utf-8') as h:readback=list(csv.DictReader(h,delimiter='\t'))
    assert {r['path']:r['original_sequence'] for r in readback}==originals
    assert {r['path']:r['represented_sequence'] for r in readback}==represented
    validation.update(txt_roundtrip=True,all_P_preserved=len(rows)==len(records),original_reconstruction=True,
        unknown_sequences_unchanged=True,unit_count_preserved=True,reference_unchanged=True)
    for p,ts in original_tokens.items():
        assert len(ts)==len(rep_tokens[p])
        for i,t in enumerate(ts):
            if t['status']!='unit' or t['ambiguity_kind']!='none':assert t['exact_bases']==rep_tokens[p][i]['exact_bases']
    dump(out/'display_graph.json',display);dump(out/'local_decisions.json',decisions)
    table(out/'unit_occurrences.tsv',occurrences);table(out/'variant_support.tsv',variants);table(out/'variant_transitions.tsv',transitions)
    dump(out/'sequence_variants.json',display['sequence_variants'])
    (out/'path_visits.jsonl').write_text(''.join(json.dumps(dict(path=p,visits=vs))+'\n' for p,vs in display['paths'].items()),encoding='utf-8')
    (out/'losses.jsonl').write_text(''.join(json.dumps(x)+'\n' for x in losses),encoding='utf-8')
    table(out/'nodes.tsv',[dict(id=n,expression=' / '.join(atom_label(x) for x in v['parts']),
        lengths=v['represented_lengths'],paths=v['path_ids'],molecules=v['route_molecules'],flags=v['uncertainty_flags'],
        lossy_paths=v['lossy_member_ids']) for n,v in display['nodes'].items()])
    table(out/'edges.tsv',display['edges']);dump(out/'validation.json',validation)
    summary=dict(loss_budget_fraction=fraction,**display['metrics'],paths=len(records),
        molecules=len({m for r in records for m in r['molecules']}),lossy_paths=len(losses),
        total_edits_bp=sum(x['loss_bp'] for x in losses),
        substitutions=sum(x['substitutions'] for x in losses),insertions=sum(x['insertions'] for x in losses),deletions=sum(x['deletions'] for x in losses),
        max_path_loss_fraction=max((x['fraction'] for x in losses),default=0),
        mean_path_loss_fraction=sum(x['fraction'] for x in losses)/sum(r.get('role')!='reference' for r in records),
        total_original_observed_bp=sum(len(r['sequence']) for r in records if r.get('role')!='reference'),
        literal_nodes=len(graph['nodes']),literal_edges=len(graph['edges']),
        sequence_variant_options=len(display['sequence_variants']),
        parameterized_parts=sum(len(n['parts']) for n in display['nodes'].values()),
        path_visits=sum(len(vs) for vs in display['paths'].values()),
        unknown_nodes=sum('unknown_content' in n['uncertainty_flags'] for n in display['nodes'].values()),
        ambiguous_nodes=sum('annotation_ambiguity' in n['uncertainty_flags'] for n in display['nodes'].values()),
        synthetic_represented_paths=sum(not r['represented_is_observed_whole_sequence'] and r['role']!='reference' for r in rows),
        validation=validation)
    dump(out/'summary.json',summary)
    if render:render_graph(out,display)
    return summary
