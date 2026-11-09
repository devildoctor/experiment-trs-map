"""Parameterized display only; literal GFA segments remain fixed length."""
from collections import Counter,defaultdict
from copy import deepcopy
from vntrdisplay_loss import alignment_key


def graph_metrics(nodes,edges):
    before=defaultdict(set);after=defaultdict(set);neighbors=defaultdict(set)
    for a,b in edges:
        after[a].add(b);before[b].add(a);neighbors[a].add(b);neighbors[b].add(a)
    pending=set(nodes);components=0
    while pending:
        components+=1;stack=[pending.pop()]
        while stack:
            for v in neighbors[stack.pop()] & pending:pending.remove(v);stack.append(v)
    return dict(N_d=len(nodes),E_d=len(edges),C=len(nodes)+len(edges),
        branch_nodes=sum(len(v)>1 for v in after.values()),
        branch_exits=sum(len(v) for v in after.values() if len(v)>1),
        merge_nodes=sum(len(v)>1 for v in before.values()),independent_blocks=components)


def build_display(graph,features,representatives,model,compress=True):
    ref=graph['backbone_path']
    # Reference coordinates are fixed across every budget, including when only
    # a subset of active representatives is evaluated for a candidate batch.
    rt=features[ref]['tokens']; ref_run={};run=0
    for i,t in enumerate(rt):
        if i==0 or t['exact_bases']!=rt[i-1]['exact_bases'] or alignment_key(t)!=alignment_key(rt[i-1]):run+=1
        ref_run[i]=run
    keyed={};paths={}
    for p,walk in sorted(graph['paths'].items()):
        tokens=features[representatives[p]]['tokens'];visits=[];i=0
        assert len(tokens)==len(walk)
        while i<len(tokens):
            t=tokens[i];j=i+1
            if t['status']=='unit':
                while j<len(tokens) and tokens[j]['status']=='unit' and alignment_key(tokens[j])==alignment_key(t) and tokens[j]['exact_bases']==t['exact_bases']:j+=1
            expansion=walk[i:j];seq=t['exact_bases'];count=j-i
            coords=[graph['coordinates'][n][0] for n in expansion]
            matching=sorted({ref_run[c[0]] for c in coords if c[1]==1 and c[0]<len(rt) and rt[c[0]]['exact_bases']==seq})
            position=('ref_run',tuple(matching)) if matching and t['status']=='unit' else ('slot',tuple(coords[0]))
            key=(position,seq,'unit' if t['status']=='unit' else t['status'])
            if key not in keyed:
                keyed[key]=dict(kind='repeat_run' if t['status']=='unit' else 'segment',
                    unit_sequence=seq,protected_types=set(),counts=set(),uncertainty_flags=set(),position=position)
            node=keyed[key];node['protected_types'].add(t['motif']);node['counts'].add(count)
            if any(x.get('annotation_ambiguity') for x in tokens[i:j]):node['uncertainty_flags'].add('annotation_ambiguity')
            if t['status']=='unknown':node['uncertainty_flags'].add('unknown_content')
            visits.append(dict(key=key,count=count,gfa_nodes=expansion,
                protected_types=[x['motif'] for x in tokens[i:j]],unit_range=[i,j]))
            i=j
        paths[p]=visits
    ids={k:f'D{i:04d}' for i,k in enumerate(sorted(keyed,key=repr),1)}
    nodes={ids[k]:dict(v,id=ids[k],protected_types=sorted(v['protected_types']),
        counts=sorted(v['counts']),uncertainty_flags=sorted(v['uncertainty_flags'])) for k,v in keyed.items()}
    for visits in paths.values():
        for v in visits:v['node']=ids[v.pop('key')]
    edges={(a['node'],b['node']) for visits in paths.values() for a,b in zip(visits,visits[1:])}
    if compress:
        before=defaultdict(set);after=defaultdict(set)
        for a,b in edges:after[a].add(b);before[b].add(a)
        starts={v[0]['node'] for v in paths.values()};ends={v[-1]['node'] for v in paths.values()}
        variable={n for n,v in nodes.items() if len(v['counts'])>1 or v['counts']!=[1]}
        def successor(n):
            if n in ends or n in variable or len(after[n])!=1:return None
            b=next(iter(after[n]))
            return b if b not in starts|variable and len(before[b])==1 and b!=n else None
        incoming={b:a for a in nodes if (b:=successor(a)) is not None}
        chains=[];seen=set()
        for head in sorted(nodes,key=lambda n:(n in incoming,n)):
            if head in seen:continue
            chain=[];n=head
            while n is not None and n not in seen:
                seen.add(n);chain.append(n);n=successor(n)
            chains.append(chain)
        mapping={n:f'D{i:04d}' for i,c in enumerate(chains,1) for n in c};new_nodes={}
        for chain in chains:
            name=mapping[chain[0]]
            if len(chain)==1:new_nodes[name]=dict(nodes[chain[0]],id=name)
            else:new_nodes[name]=dict(id=name,kind='module',parts=[nodes[n] for n in chain],
                protected_types=[t for n in chain for t in nodes[n]['protected_types']],
                uncertainty_flags=sorted({x for n in chain for x in nodes[n]['uncertainty_flags']}),counts=[1])
        new_paths={}
        for p,visits in paths.items():
            out=[]
            for v in visits:
                name=mapping[v['node']]
                if out and name==out[-1]['node'] and new_nodes[name]['kind']=='module':
                    out[-1]['gfa_nodes']+=v['gfa_nodes'];out[-1]['parts'].append(v)
                else:out.append(dict(node=name,count=v['count'],gfa_nodes=list(v['gfa_nodes']),parts=[v]) if new_nodes[name]['kind']=='module' else dict(v,node=name))
            new_paths[p]=out
        nodes,paths=new_nodes,new_paths
        edges={(a['node'],b['node']) for visits in paths.values() for a,b in zip(visits,visits[1:])}
    return dict(schema_version=2,nodes=nodes,edges=[dict(source=a,target=b) for a,b in sorted(edges)],
                paths=paths,metrics=graph_metrics(nodes,edges),backbone_path=ref,
                scope='Only registered full P walks and per-visit parameters are observations.')


def validate_display(display,graph,represented):
    def unit_string(node,count):
        if node['kind']=='module':return ''.join(unit_string(p,1) for p in node['parts'])
        if count not in node['counts']:raise ValueError('Unobserved repeat count')
        return node['unit_sequence']*count
    for p,visits in display['paths'].items():
        nodes=[n for v in visits for n in v['gfa_nodes']]
        if nodes!=graph['paths'][p]:raise ValueError('Display expansion changed GFA walk')
        literal=''.join(graph['nodes'][n] for n in nodes)
        parameterized=''.join(unit_string(display['nodes'][v['node']],v['count']) for v in visits)
        if literal!=represented[p] or parameterized!=literal:raise ValueError('Display parameter expansion changed bases')
    return dict(paths=len(represented),gfa_expansion=True,parameter_expansion=True)


def add_support(display,graph,records,losses,groups):
    """Local exact support is computed from the represented-coordinate edit map.

    A connection is exact when neither adjacent base, nor its interbase
    boundary, was changed. Edits elsewhere on the same P do not remove support.
    """
    members={r['id']:set(r['molecules']) for r in records};loss_by_id={x['path']:x for x in losses}
    node_ps=defaultdict(set);edge_ps=defaultdict(set);node_exact=defaultdict(set);edge_exact=defaultdict(set)
    lengths=defaultdict(set)
    def clean(p,a,b):
        for e in loss_by_id.get(p,{}).get('events',[]):
            if e['rstart']==e['rend']:
                if a<e['rstart']<b:return False
            elif e['rstart']<b and a<e['rend']:return False
        return True
    for p,visits in display['paths'].items():
        position=0;previous=None
        for v in visits:
            n=v['node'];length=sum(len(graph['nodes'][x]) for x in v['gfa_nodes'])
            lengths[n].add(length)
            node_ps[n].add(p)
            if clean(p,position,position+length):node_exact[n].add(p)
            if previous is not None:
                edge=(previous,n);edge_ps[edge].add(p)
                if clean(p,max(0,position-1),position+1):edge_exact[edge].add(p)
            v['represented_interval']=[position,position+length];position+=length;previous=n
    def stats(ps,raw):
        mols=set().union(*(members[p] for p in ps));exact=set().union(*(members[p] for p in raw)) if raw else set()
        return dict(path_ids=sorted(ps),structure_group_ids=sorted({groups[p] for p in ps}),
            route_molecules=len(mols),raw_exact_molecules=len(exact),
            lossy_member_ids=sorted(p for p in ps if p in loss_by_id))
    for n,v in display['nodes'].items():v.update(stats(node_ps[n],node_exact[n]),represented_lengths=sorted(lengths[n]))
    for e in display['edges']:e.update(stats(edge_ps[e['source'],e['target']],edge_exact[e['source'],e['target']]))
    return display
