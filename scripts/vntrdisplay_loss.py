"""Opt-in v2 unit correspondence and structure-constrained display loss.

The frozen partition remains evidence, not a genotype. Representatives are
observed strings; every member is compared directly with its final target.
"""
from collections import Counter, defaultdict
from fractions import Fraction
import math

from vntrpathstructure import restore_original
from vntrpathgraph import validate_graph

SELECTION = 'unit-v2-primary-exact-gap-left-lexical'


def alignment_key(token):
    return token['motif'], token.get('orientation', '+')


def display_signature(feature):
    result=[]; tokens=feature['tokens']
    for i,t in enumerate(tokens):
        if t['status']=='unknown':
            result.append(('unknown',t['exact_bases']))
        elif t['status']=='partial' and (i not in (0,len(tokens)-1) or
                                       len(t['template'])-len(t['exact_bases'])>2):
            result.append(('partial',alignment_key(t),t['exact_bases']))
        else:
            result.append(('slot',alignment_key(t)))
    return tuple(result)


def align_units(query, reference):
    """Global M/I/D DP; all secondary scores apply to complete routes.

    I/D describe the observation relative to the reference. Primary ties are
    counted independently of secondary selection. Gap coordinates precede the
    final operation priority, giving left-normalized repeat gaps.
    """
    n,m=len(query),len(reference)
    scores={(0,0,'M'):(0,0,0,(),())}; traces={}; primary={(0,0):(0,1)}
    for i in range(n+1):
        for j in range(m+1):
            if not i and not j: continue
            pc=[]
            transitions=[]
            if i and j:
                a,b=query[i-1],reference[j-1]
                same=a['exact_bases']==b['exact_bases']
                delta=0 if same or alignment_key(a)==alignment_key(b) else 1
                transitions.append((i-1,j-1,'M',delta,len(a['exact_bases']) if same else 0,0 if same else 1))
            if i: transitions.append((i-1,j,'I',1,0,2))
            if j: transitions.append((i,j-1,'D',1,0,3))
            for a,b,state,delta,exact,op in transitions:
                pp,ways=primary[a,b];pc.append((pp+delta,ways))
                for prev in ('M','I','D'):
                    if (a,b,prev) not in scores: continue
                    s=scores[a,b,prev]; opening=state!='M' and prev!=state
                    coords=s[3]+((b,a,op),) if opening else s[3]
                    value=(s[0]+delta,s[1]-exact,s[2]+int(opening),coords,s[4]+(op,))
                    key=(i,j,state)
                    if key not in scores or value<scores[key]:
                        scores[key]=value;traces[key]=(a,b,prev)
            best=min(x[0] for x in pc)
            primary[i,j]=(best,min(2,sum(w for c,w in pc if c==best)))
    key=min(((n,m,s) for s in ('M','I','D') if (n,m,s) in scores),key=scores.get)
    score=scores[key];steps=[]
    while key[:2]!=(0,0):
        a,b,prev=traces[key];i,j,state=key
        steps.append((a if a<i else None,b if b<j else None));key=(a,b,prev)
    return list(reversed(steps)),dict(had_primary_ties=primary[n,m][1]>1,
        selection_policy_version=SELECTION,primary_cost=score[0],exact_paired_bp=-score[1],
        gap_runs=score[2],gap_coordinates=score[3])


def canonical_edits(original, represented):
    """Shortest script, fewest gap runs, leftmost edits, then =/S/I/D.

    Each call is restricted to one corresponding unit, so no alternate script
    can cross a protected unit/module boundary.
    """
    n,m=len(original),len(represented)
    scores={(0,0,'M'):(0,0,(),())};traces={}
    for i in range(n+1):
        for j in range(m+1):
            if not i and not j:continue
            opts=[]
            if i and j:opts.append((i-1,j-1,'M','=' if original[i-1]==represented[j-1] else 'S'))
            if j:opts.append((i,j-1,'I','I'))
            if i:opts.append((i-1,j,'D','D'))
            for a,b,state,op in opts:
                for prev in ('M','I','D'):
                    old=scores.get((a,b,prev))
                    if old is None:continue
                    rank={'=':0,'S':1,'I':2,'D':3}[op]
                    val=(old[0]+(op!='='),old[1]+(state!='M' and state!=prev),
                         old[2]+((a,b,rank),) if op!='=' else old[2],old[3]+(rank,))
                    key=(i,j,state)
                    if key not in scores or val<scores[key]:
                        scores[key]=val;traces[key]=(a,b,prev,op)
    key=min(((n,m,s) for s in ('M','I','D') if (n,m,s) in scores),key=scores.get)
    events=[]
    while key[:2]!=(0,0):
        i,j,_=key;a,b,prev,op=traces[key]
        if op!='=':events.append(dict(op=op,ostart=a,oend=i,rstart=b,rend=j,
                                     original=original[a:i],represented=represented[b:j]))
        key=(a,b,prev)
    return list(reversed(events))


def compare_for_merge(source,target,sf,tf,policy,fraction):
    if display_signature(sf)!=display_signature(tf): return None,'protected_structure'
    budget=math.floor(len(source)*fraction+1e-9)
    if abs(len(source)-len(target))>budget:return None,'budget'
    events=[]
    for i,(a,b) in enumerate(zip(sf['tokens'],tf['tokens'])):
        if a['exact_bases']==b['exact_bases']:continue
        if a['status']=='unknown' or b['status']=='unknown':return None,'unknown_content'
        local=canonical_edits(a['exact_bases'],b['exact_bases'])
        if len(local)>policy['unit_max_edits']:return None,'unit_budget'
        for e in local:
            terminal=None
            if i==0 and e['ostart']==0 and e['rstart']==0 and e['op'] in ('I','D'):terminal='left'
            if i==len(sf['tokens'])-1 and e['oend']==len(a['exact_bases']) and e['rend']==len(b['exact_bases']) and e['op'] in ('I','D'):terminal='right'
            events.append(dict(e,ostart=e['ostart']+a['start'],oend=e['oend']+a['start'],
                               rstart=e['rstart']+b['start'],rend=e['rend']+b['start'],
                               unit=i+1,terminal_completion=terminal))
    if not events or len(events)>budget:return None,'budget'
    for side in ('left','right'):
        if sum(e['terminal_completion']==side for e in events)>policy['terminal_max_edits']:return None,'terminal_budget'
    if restore_original(target,events)!=source:raise ValueError('Residual reconstruction failed')
    return events,'compatible_unit_edits'


def build_unit_graph(records,features,model,backbone='P0001'):
    reference=features[backbone]['tokens']; descriptions={}; walks={};audit={};token_by_key={}
    for record in sorted(records,key=lambda r:r['id']):
        p=record['id'];query=features[p]['tokens']
        steps,info=align_units(query,reference);audit[p]=info
        walk=[];boundary=0;insertion=0
        for qi,rj in steps:
            if qi is not None:
                t=query[qi]
                slot=(boundary,0,insertion) if rj is None else (rj,1,0)
                # Literal bases + comparison position, independent of RT labels.
                key=(slot,t['exact_bases']);descriptions[key]=t['exact_bases']
                token_by_key.setdefault(key,[]).append(t)
                walk.append(key);insertion+=1
            if rj is not None:boundary=rj+1;insertion=0
        walks[p]=walk
    keys=sorted(descriptions);ids={k:i+1 for i,k in enumerate(keys)}
    paths={p:[ids[k] for k in w] for p,w in walks.items()}
    graph=dict(nodes={ids[k]:descriptions[k] for k in keys},paths=paths,
               edges={e for w in paths.values() for e in zip(w,w[1:])},
               coordinates={ids[k]:[list(k[0])] for k in keys},alignments=audit,
               node_tokens={ids[k]:token_by_key[k] for k in keys},backbone_path=backbone)
    validate_graph(graph,{r['id']:r['sequence'] for r in records})
    return graph


def represented_graph(exact,representatives):
    paths={p:exact['paths'][rep] for p,rep in representatives.items()}
    used={n for w in paths.values() for n in w}
    return dict(exact,nodes={n:exact['nodes'][n] for n in used},paths=paths,
                edges={e for w in paths.values() for e in zip(w,w[1:])})


def merge_paths(records,features,model,policy,fraction,exact,display_factory):
    if policy['schema_version']!=2 or policy['model_version']!=model['version']:raise ValueError('Policy mismatch')
    if not 0<=fraction<=.05:raise ValueError('v2 budget must be in [0, .05]')
    originals={r['id']:r['sequence'] for r in records};reps={p:p for p in originals}
    def rank(r):
        f=features[r['id']]
        return (sum(t['status']!='unit' for t in f['tokens']),-len(r['hifi_molecules']),
                -len(r['molecules']),sum(t.get('edit_distance',0) for t in f['tokens']),r['sha256'])
    ordered=sorted((r for r in records if r.get('role')!='reference'),key=rank)
    order={r['id']:i for i,r in enumerate(ordered)};buckets=defaultdict(list)
    for r in ordered:buckets[display_signature(features[r['id']])].append(r)
    candidates={};decisions=[]
    if fraction:
        for bucket in buckets.values():
            for i,target in enumerate(bucket):
                for source in bucket[i+1:]:
                    p,q=source['id'],target['id']
                    events,reason=compare_for_merge(source['sequence'],target['sequence'],features[p],features[q],policy,fraction)
                    if events is not None:candidates[p,q]=events
                    decisions.append(dict(path=p,target_path=q,eligible=events is not None,reason=reason))
    cache={}
    def complexity(mapping):
        active=tuple(sorted(set(mapping.values())))
        if active not in cache:
            minimal={p:p for p in active}
            display=display_factory(represented_graph(exact,minimal),features,minimal,model,compress=True)
            cache[active]=display['metrics']['C']
        return cache[active]
    current=complexity(reps);losses=[];accepted=[];locked=set()
    # Both representative groups and graph-element owner batches can remove a
    # detail shared by several P records. No search over arbitrary subsets.
    while candidates:
        eligible={(p,q):es for (p,q),es in candidates.items() if reps[p]==p and p not in locked and reps[q]==q}
        if not eligible:break
        batches=set();by_target=defaultdict(list);by_source=defaultdict(list)
        for p,q in eligible:by_target[q].append(p);by_source[p].append(q)
        for q,ps in by_target.items():
            batches.add(tuple(sorted((p,q) for p in ps)))
            batches.update(((p,q),) for p in ps)
        owners=defaultdict(set)
        for p in reps:
            for n in exact['paths'][reps[p]]:owners[('N',n)].add(p)
            for e in zip(exact['paths'][reps[p]],exact['paths'][reps[p]][1:]):owners[('E',e)].add(p)
        for member_set in owners.values():
            if not member_set or not member_set.issubset(by_source):continue
            pairs=[]
            for p in sorted(member_set):
                qs=[q for q in by_source[p] if q not in member_set]
                if not qs:break
                q=min(qs,key=lambda q:(len(eligible[p,q]),order[q]));pairs.append((p,q))
            else:batches.add(tuple(pairs))
        best=None
        for batch in sorted(batches):
            sources={p for p,q in batch};targets={q for p,q in batch}
            if sources & targets:continue
            trial=dict(reps);trial.update(batch);after=complexity(trial);gain=current-after
            if gain<=0:continue
            edits=sum(len(eligible[p,q]) for p,q in batch)
            key=(-Fraction(gain,edits),-gain,edits,tuple(order[q] for p,q in batch),batch)
            if best is None or key<best[0]:best=(key,batch,after,gain,edits)
        if best is None:break
        _,batch,after,gain,edits=best
        accepted.append(dict(pairs=batch,C_before=current,C_after=after,delta_C=gain,edits_bp=edits))
        for p,q in batch:
            reps[p]=q;locked.add(q);es=eligible[p,q];counts=Counter(e['op'] for e in es)
            losses.append(dict(path=p,target_path=q,loss_bp=len(es),substitutions=counts['S'],
                insertions=counts['I'],deletions=counts['D'],original_length=len(originals[p]),
                represented_length=len(originals[q]),fraction=len(es)/len(originals[p]),events=es,
                policy_version=policy['version'],protected_structure=True))
        current=after
    return {p:originals[q] for p,q in reps.items()},reps,losses,dict(candidates=decisions,accepted_batches=accepted)
