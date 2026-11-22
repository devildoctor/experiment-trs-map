"""Substitution-only edits on frozen visits, scored on the exported literal graph."""
from collections import Counter, defaultdict
from copy import deepcopy
from fractions import Fraction

from vntrpathstructure import unit_key
from vntrsupernode import CERTAIN


def plan_substitutions(records,tokens,slots,display,visits,policy):
    represented = deepcopy(tokens)
    groups = defaultdict(lambda:defaultdict(list))
    members = {r['id']:set(r['molecules']) for r in records if r['role']!='reference'}
    remaining = {r['id']:int(len(r['sequence'])*Fraction(str(policy['loss_fraction']))) for r in records}
    # The graph writer keys each segment by the aligned slot and literal bases.
    walks = {p:[(slot,t['exact_bases']) for slot,t in zip(slots[p],ts)] for p,ts in tokens.items()}
    nodes = Counter(n for w in walks.values() for n in w)
    edges = Counter(e for w in walks.values() for e in zip(w,w[1:]))
    for p,vs in sorted(visits.items()):
        if p not in members: continue
        for v in vs:
            node = display['nodes'][v['supernode_id']]
            if (node['protected'] or v['structure_assessment']!=CERTAIN or
                v['boundary_assessment']!=CERTAIN or v['position_assessment']!='aligned'): continue
            for i in range(*v['token_interval']):
                t = tokens[p][i]
                group = (node['id'],slots[p][i],unit_key(t),len(t['exact_bases']))
                groups[group][t['exact_bases']].append((p,i))
    candidates = []
    for group,strings in sorted(groups.items()):
        if len(strings)<2: continue
        def rank(s):
            ps = {p for p,i in strings[s]}
            return (-len(set().union(*(members[p] for p in ps))),-len(ps),s)
        target = min(strings,key=rank)
        donor = min(strings[target])
        for source,occurrences in sorted(strings.items()):
            diffs = [i for i,(a,b) in enumerate(zip(source,target)) if a!=b]
            if 0<len(diffs)<=policy['max_unit_substitutions']:
                candidates.append((group,source,target,donor,diffs,occurrences))

    def trial(candidate):
        group,source,target,donor,diffs,occurrences = candidate
        costs = Counter()
        changed = {}
        affected_edges = set()
        for p,i in occurrences:
            if costs[p]+len(diffs)>remaining[p]: continue
            costs[p]+=len(diffs)
            changed[p,i] = (slots[p][i],target)
            if i: affected_edges.add((p,i-1))
            if i+1<len(walks[p]): affected_edges.add((p,i))
        dn,de = Counter(),Counter()
        for (p,i),new in changed.items():
            dn[walks[p][i]]-=1
            dn[new]+=1
        for p,i in affected_edges:
            a,b = walks[p][i:i+2]
            de[a,b]-=1
            de[changed.get((p,i),a),changed.get((p,i+1),b)]+=1
        gain = sum(int(counter[k]>0)-int(counter[k]+v>0)
                   for counter,delta in ((nodes,dn),(edges,de)) for k,v in delta.items())
        return gain,changed,costs,dn,de

    events = defaultdict(list)
    decisions = []
    while candidates:
        best = None
        for index,c in enumerate(candidates):
            gain,changed,costs,dn,de = trial(c)
            if gain<=0: continue
            bp = sum(costs.values())
            order = (-Fraction(gain,bp),-gain,bp,index)
            if best is None or order<best[0]: best=(order,index,gain,changed,costs,dn,de)
        if best is None: break
        _,index,gain,changed,costs,dn,de = best
        group,source,target,(rp,ri),diffs,_ = candidates.pop(index)
        nodes.update(dn)
        edges.update(de)
        batch=len(decisions)+1
        decisions.append(dict(batch=batch,delta_sequence_complexity=gain,edited_bp=sum(costs.values()),
                              source=source,representative=target,representative_occurrence=f'{rp}:u{ri}'))
        for p,cost in costs.items(): remaining[p]-=cost
        for (p,i),node in sorted(changed.items()):
            walks[p][i]=node
            t=represented[p][i]
            t['exact_bases']=target
            for j in diffs:
                pos=t['start']+j
                events[p].append(dict(id=f'{p}:S{pos}',op='S',path=p,occurrence_id=f'{p}:u{i}',
                    supernode_id=group[0],ostart=pos,oend=pos+1,rstart=pos,rend=pos+1,
                    original=source[j],represented=target[j],representative_occurrence=f'{rp}:u{ri}',
                    representative_sequence=target,batch=batch))
    return represented,{p:sorted(es,key=lambda e:e['ostart']) for p,es in events.items()},decisions

