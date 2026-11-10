"""Expandable, position-matched families with a direct one-base representative."""
from collections import Counter,defaultdict
from itertools import combinations
from vntrcore import edit_distance
from vntrlocal import topology


def apply_families(visits,rules):
    targets=[rules[v['key']]['key'] if v['key'] in rules else v['key'] for v in visits]
    repeated={k for k,n in Counter(targets).items() if n>1}
    out=[]
    for v,key in zip(visits,targets):
        if v['key'] in rules and key not in repeated:
            rule=rules[v['key']];v=dict(v,key=key,display_family=rule['family'])
            v['category']=('unit_family',rule['family']['id'])+v['category'][2:]
        out.append(v)
    return out


def select_families(views,model):
    """Search global partitions; merge only mutually exclusive single-part nodes.

    Sequence payloads, selected types and run boundaries are never rewritten.
    This searches the stated candidate space, not all possible graph designs.
    """
    sequences={m['id']:m['sequence'] for m in model['motifs']};ids=sorted(sequences)
    counts=Counter(t['motif'] for vs in views.values() for v in vs for t in v['tokens'])
    distances={(a,b):edit_distance(sequences[a],sequences[b]) for a in ids for b in ids}
    allowed={}
    for size in range(1,len(ids)+1):
        for block in combinations(ids,size):
            representatives=[a for a in block if all(distances[a,b]<=1 for b in block)]
            if representatives:allowed[frozenset(block)]=min(representatives,key=lambda a:(-counts[a],a))
    def partitions(rest):
        if not rest:yield [];return
        first=min(rest)
        for block in sorted(allowed,key=lambda b:tuple(sorted(b))):
            if first in block and block<=rest:
                for tail in partitions(rest-block):yield [tuple(sorted(block))]+tail
    chains,_,_,_=topology(views);single={chain[0] for chain in chains if len(chain)==1}
    atoms={v['key']:v for vs in views.values() for v in vs};members=defaultdict(set)
    for p,vs in views.items():
        for v in vs:members[v['key']].add(p)
    positions=defaultdict(list)
    for key in single:
        v=atoms[key];cat=v['category']
        if cat[0]=='unit' and cat[-1]=='internal' and cat[1] in sequences:
            positions[(v['position'],cat[2:])].append(key)
    best=None;trials=0
    for partition in partitions(frozenset(ids)):
        trials+=1;assignment={a:b for b in partition for a in b};rules={};merges=[]
        family_ids={b:'F%02d'%(i+1) for i,b in enumerate(b for b in partition if len(b)>1)}
        for (position,roles),keys in sorted(positions.items(),key=lambda x:repr(x[0])):
            local=defaultdict(list)
            for key in keys:local[assignment[atoms[key]['category'][1]]].append(key)
            for block,kk in sorted(local.items()):
                if len(block)<2 or len(kk)<2:continue
                if len({atoms[k]['category'][1] for k in kk})!=len(kk):continue
                if any(members[a]&members[b] for a,b in combinations(kk,2)):continue
                family=dict(id=family_ids[block],members=list(block),representative=allowed[frozenset(block)],radius_bp=1)
                target=('family',position,roles,block)
                for key in kk:rules[key]=dict(key=target,family=family)
                merges.append(dict(position=position,family=family,source_nodes=len(kk)))
        trial={p:apply_families(vs,rules) for p,vs in views.items()}
        metrics=topology(trial)[3]
        cost=sum(counts[a]*distances[a,allowed[frozenset(b)]] for b in partition for a in b)
        rank=(metrics['C'],max(map(len,partition)),cost,partition)
        if best is None or rank<best[0]:best=(rank,rules,trial,partition,merges,metrics)
    _,rules,trial,partition,merges,metrics=best
    report=dict(partition=partition,families=[dict(id='F%02d'%(i+1),members=list(b),representative=allowed[frozenset(b)],radius_bp=1)
                for i,b in enumerate(b for b in partition if len(b)>1)],
        local_merges=merges,candidate_partitions=trials,metrics=metrics,
        scope='Same position/orientation/role, single-part internal nodes, no co-occurrence on one registered P; direct representative radius <=1 bp.')
    return trial,rules,report
