"""Local, reversible A-E display model. No change to extraction or raw P records."""
from collections import Counter, defaultdict
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
import math

from vntrdisplay_loss import align_units, canonical_edits
from vntrdisplay import graph_metrics
from vntrcore import edit_distance
from vntrpathstructure import restore_original

VERSION='local-v3-AE-20261110'


def signature(t):
    return t['start'],t['end'],t['status'],t['motif'],t.get('orientation','+')


def combine_unknown(tokens):
    result=[]
    for token in tokens:
        t=deepcopy(token);t.pop('cost',None)
        if result and t['status']=='unknown' and result[-1]['status']=='unknown':
            result[-1]['end']=t['end'];result[-1]['exact_bases']+=t['exact_bases']
            result[-1]['annotation_ambiguity']|=t.get('annotation_ambiguity',False)
        else:result.append(t)
    return result


def ambiguity_regions(feature):
    """Exhaust all paths between consecutive mandatory vertices, never a top-k."""
    dag=feature['annotation'].get('optimal_dag')
    if not dag:return []
    outgoing=defaultdict(list)
    for edge in dag['edges']:outgoing[edge['start']].append(edge)
    result=[]
    for start,end in zip(dag['mandatory_vertices'],dag['mandatory_vertices'][1:]):
        ways={end:1}
        for i in range(end-1,start-1,-1):
            ways[i]=sum(ways.get(e['end'],0) for e in outgoing[i] if e['end']<=end)
        if ways.get(start,0)<2:continue
        alternatives=[];stack=[(start,[])]
        while stack:
            i,path=stack.pop()
            if i==end:alternatives.append(combine_unknown(path));continue
            for e in reversed(outgoing[i]):
                if e['end']<=end and ways.get(e['end']):stack.append((e['end'],path+[e]))
        assert len(alternatives)==ways[start]
        alternatives.sort(key=lambda ts:tuple(signature(t) for t in ts))
        layouts={tuple((t['start'],t['end'],t['status'],t.get('orientation')) for t in ts) for ts in alternatives}
        result.append(dict(start=start,end=end,alternatives=alternatives,number=len(alternatives),
            kind='type_only' if len(layouts)==1 else 'boundary_count_or_status',
            complete_unit_counts=sorted({sum(t['status']=='unit' for t in ts) for ts in alternatives})))
    return result


def classify(tokens,regions):
    result=deepcopy(tokens)
    for i,t in enumerate(result):
        overlapping=[(ri,r) for ri,r in enumerate(regions) if r['start']<t['end'] and t['start']<r['end']]
        t['ambiguity_regions']=[ri for ri,r in overlapping]
        t['ambiguity_kind']=('boundary_count_or_status' if any(r['kind']=='boundary_count_or_status' for _,r in overlapping)
            else 'type_only' if overlapping else 'none')
        candidates={t['motif']}
        for _,r in overlapping:
            candidates.update(x['motif'] for alt in r['alternatives'] for x in alt
                              if (x['start'],x['end'],x['status'])==(t['start'],t['end'],t['status']))
        t['candidate_types']=sorted(candidates)
        t['match_kind']=t['status'] if t['status']!='unit' else 'approximate' if t.get('edit_distance',0) else 'exact'
        t['weak_approximate']=t['status']=='unit' and t.get('edit_distance',0)>2
        t['occurrence_index']=i
    return result


def replace_region(tokens,start,end,alternative):
    before=[];after=[]
    for t in tokens:
        if t['end']<=start:before.append(t)
        elif t['start']<start:
            assert t['status']=='unknown'
            before.append(dict(t,end=start,exact_bases=t['exact_bases'][:start-t['start']]))
        if t['start']>=end:after.append(t)
        elif t['end']>end:
            assert t['status']=='unknown'
            after.append(dict(t,start=end,exact_bases=t['exact_bases'][end-t['start']:]))
    result=combine_unknown(before+alternative+after)
    assert ''.join(t['exact_bases'] for t in result)==''.join(t['exact_bases'] for t in tokens)
    return result


def assign_slots(tokens,reference):
    steps,audit=align_units(tokens,reference);slots=[];boundary=0;insertion=0
    for qi,rj in steps:
        if qi is not None:
            slots.append((boundary,0,insertion) if rj is None else (rj,1,0));insertion+=1
        if rj is not None:boundary=rj+1;insertion=0
    # Preserve all primary-optimal partner positions, not merely the chosen trace.
    n,m=len(tokens),len(reference);f=[[0]*(m+1) for _ in range(n+1)]
    def delta(i,j):
        a,b=tokens[i],reference[j]
        return int(a['exact_bases']!=b['exact_bases'] and (a['motif'],a.get('orientation'))!=(b['motif'],b.get('orientation')))
    for i in range(n+1):
        for j in range(m+1):
            if not i:f[i][j]=j
            elif not j:f[i][j]=i
            else:f[i][j]=min(f[i-1][j]+1,f[i][j-1]+1,f[i-1][j-1]+delta(i-1,j-1))
    b=[[0]*(m+1) for _ in range(n+1)]
    for i in range(n,-1,-1):
        for j in range(m,-1,-1):
            if i==n:b[i][j]=m-j
            elif j==m:b[i][j]=n-i
            else:b[i][j]=min(b[i+1][j]+1,b[i][j+1]+1,b[i+1][j+1]+delta(i,j))
    possibilities=[]
    for i in range(n):
        ps=[]
        for j in range(m+1):
            if f[i][j]+1+b[i+1][j]==f[n][m]:ps.append(('gap',j))
            if j<m and f[i][j]+delta(i,j)+b[i+1][j+1]==f[n][m]:ps.append(('unit',j))
        possibilities.append(ps)
    audit.update(selected_steps=steps,selected_slots=slots,primary_partner_options=possibilities)
    return slots,audit


def category(t):
    role=t.get('role','internal');orient=t.get('orientation','+')
    if t['status']!='unit':return t['status'],role,orient
    if t['ambiguity_kind']=='type_only':return 'type_ambiguous',tuple(t['candidate_types']),orient,role
    if t['ambiguity_kind']=='boundary_count_or_status':return 'partition_ambiguous',t['motif'],orient,role
    if t['weak_approximate']:return 'weak_unit',t['motif'],orient,role
    return 'unit',t['motif'],orient,role


def path_view(tokens,slots,reference):
    """A/B/D families, with per-occurrence strings and selected annotation retained."""
    ref_run={};run=0
    for i,t in enumerate(reference):
        if not i or (t['motif'],t['status'],t.get('orientation'))!=(reference[i-1]['motif'],reference[i-1]['status'],reference[i-1].get('orientation')):run+=1
        ref_run[i]=run
    visits=[];i=0;seen=Counter()
    while i<len(tokens):
        t=tokens[i];cat=category(t);j=i+1
        if t['status']=='unit':
            while j<len(tokens) and category(tokens[j])==cat and tokens[j]['status']=='unit':j+=1
        coords=slots[i:j]
        matching=tuple(sorted({ref_run[c[0]] for c in coords if c[1]==1 and c[0]<len(reference)
                              and reference[c[0]]['motif']==t['motif']}))
        pos=('reference_runs',matching) if matching and t['status']=='unit' else ('slot',slots[i])
        if t['status']=='unknown':
            left=(slots[i-1],tokens[i-1]['status'],tokens[i-1]['motif']) if i else ('START',)
            right=(slots[j],tokens[j]['status'],tokens[j]['motif']) if j<len(tokens) else ('END',)
            # At least one full unit long: retain the literal block as a distinct structure.
            extent=('short_unresolved',) if len(t['exact_bases'])<min(len(x['template']) for x in reference if x.get('template')) else ('long_block',t['exact_bases'])
            key=(pos,cat,left,right,extent)
        elif t['status']=='partial':key=(pos,cat,t['motif'],len(t['exact_bases']))
        else:key=(pos,cat)
        # Distinct visits separated by an interruption may not turn into a cycle.
        occurrence=seen[key];seen[key]+=1;key=(key,occurrence)
        visits.append(dict(key=key,count=j-i,variant_vector=[x['exact_bases'] for x in tokens[i:j]],
            tokens=tokens[i:j],occurrence_indices=list(range(i,j)),slots=slots[i:j],
            original_interval=[t['start'],tokens[j-1]['end']],position=pos,category=cat))
        i=j
    return visits


def topology(views):
    nodes={v['key'] for vs in views.values() for v in vs}
    edges={(a['key'],b['key']) for vs in views.values() for a,b in zip(vs,vs[1:])}
    before=defaultdict(set);after=defaultdict(set)
    for a,b in edges:before[b].add(a);after[a].add(b)
    starts={vs[0]['key'] for vs in views.values()};ends={vs[-1]['key'] for vs in views.values()}
    def successor(n):
        if n in ends or len(after[n])!=1:return None
        b=next(iter(after[n]));return b if b not in starts and len(before[b])==1 and b!=n else None
    incoming={b:a for a in nodes if (b:=successor(a)) is not None}
    chains=[];seen=set()
    for head in sorted(nodes,key=lambda n:(n in incoming,repr(n))):
        if head in seen:continue
        chain=[];n=head
        while n is not None and n not in seen:seen.add(n);chain.append(n);n=successor(n)
        chains.append(chain)
    mapping={n:i for i,c in enumerate(chains) for n in c};paths={}
    for p,vs in views.items():
        out=[]
        for v in vs:
            n=mapping[v['key']]
            if not out or out[-1]!=n:out.append(n)
        # Catch accidental cycles or invalid compression immediately.
        assert [k for n in out for k in chains[n]]==[v['key'] for v in vs],p
        paths[p]=out
    es={(a,b) for w in paths.values() for a,b in zip(w,w[1:])}
    return chains,mapping,paths,graph_metrics(range(len(chains)),es)


class ComplexityIndex:
    """Count actual post-chain nodes/edges without rebuilding every path per trial."""
    def __init__(self,views):
        self.views=views;self.nodes=Counter();self.edges=Counter();self.starts=Counter();self.ends=Counter()
        for vs in views.values():self.change(vs,1,self.nodes,self.edges,self.starts,self.ends)

    @staticmethod
    def change(vs,sign,nodes,edges,starts,ends):
        keys=[v['key'] for v in vs]
        for n in set(keys):nodes[n]+=sign
        for e in set(zip(keys,keys[1:])):edges[e]+=sign
        starts[keys[0]]+=sign;ends[keys[-1]]+=sign

    def measure(self,updates):
        ns=self.nodes.copy();es=self.edges.copy();starts=self.starts.copy();ends=self.ends.copy()
        for p,vs in updates.items():
            self.change(self.views[p],-1,ns,es,starts,ends);self.change(vs,1,ns,es,starts,ends)
        active={n for n,c in ns.items() if c};edges=[e for e,c in es.items() if c]
        ind=Counter(b for a,b in edges);out=Counter(a for a,b in edges)
        contracted=sum(a!=b and not ends[a] and not starts[b] and out[a]==1 and ind[b]==1 for a,b in edges)
        n=len(active)-contracted;e=len(edges)-contracted
        return dict(N_d=n,E_d=e,C=n+e)

    def update(self,updates):
        for p,vs in updates.items():
            self.change(self.views[p],-1,self.nodes,self.edges,self.starts,self.ends)
            self.change(vs,1,self.nodes,self.edges,self.starts,self.ends);self.views[p]=vs


def select_interpretations(features):
    """Exhaustive per-region best response against the rest of the full display.

    Coordinate descent reaches a regional fixed point, not a proven joint global
    optimum across every P. The reference annotation is frozen as layout origin.
    """
    regions={p:ambiguity_regions(f) for p,f in features.items()}
    tokens={p:classify(f['tokens'],regions[p]) for p,f in features.items()}
    reference=tokens['P0001'];slots={};alignments={};views={}
    for p,ts in tokens.items():
        slots[p],alignments[p]=assign_slots(ts,reference);views[p]=path_view(ts,slots[p],reference)
    index=ComplexityIndex(views)
    print('AMBIGUITY_REGIONS',sum(len(rs) for rs in regions.values()),'MAX_ALTERNATIVES',
          max((r['number'] for rs in regions.values() for r in rs),default=0),flush=True)
    decisions=[];rounds=0
    while True:
        rounds+=1;changed=False
        print('INTERPRETATION_ROUND',rounds,flush=True)
        for p in sorted(tokens):
            if p=='P0001':continue
            for ri,region in enumerate(regions[p]):
                start,end=region['start'],region['end']
                candidates=[]
                for ai,alt in enumerate(region['alternatives']):
                    trial=classify(replace_region(tokens[p],start,end,alt),regions[p]);ss,aa=assign_slots(trial,reference)
                    vv=path_view(trial,ss,reference);metric=index.measure({p:vv})
                    lexical=tuple(signature(t) for t in trial)
                    candidates.append(((metric['N_d'],metric['E_d'],lexical),ai,trial,ss,aa,vv))
                best=min(candidates,key=lambda x:x[0]);_,ai,trial,ss,aa,vv=best
                old_sig=tuple(signature(t) for t in tokens[p]);new_sig=tuple(signature(t) for t in trial)
                if old_sig!=new_sig:
                    index.update({p:vv});tokens[p],slots[p],alignments[p]=trial,ss,aa;changed=True
                    decisions.append(dict(round=rounds,path=p,region=ri,selected=ai,
                        alternative_count=len(candidates),N_after=best[0][0],E_after=best[0][1]))
                region['selected']=ai
                region['candidate_graph_metrics']=[dict(alternative=x[1],N_d=x[0][0],E_d=x[0][1]) for x in candidates]
        if not changed:break
    assert index.measure({})['C']==topology(views)[3]['C']
    return tokens,slots,alignments,regions,dict(rounds=rounds,decisions=decisions,
        policy='primary-optimal_region_exhaustive_minimum_full_display_N_then_E_then_lexical_fixed_point',
        global_joint_optimum_proven=False,reference_annotation_frozen=True)


def eligible(t,model):
    return (t['status']=='unit' and t['ambiguity_kind']=='none' and not t['weak_approximate']
            and t.get('role','internal')=='internal' and t['motif'] not in model.get('interruptions',[]))


def plan_edits(records,tokens,slots,model,fraction,initial=None):
    """Direct occurrence-to-observed-string substitution; D/ambiguity never edited."""
    if not 0<=fraction<=.10:raise ValueError('v3 fraction outside [0,0.10]')
    by_id={r['id']:r for r in records};reference=tokens['P0001'];originals={p:r['sequence'] for p,r in by_id.items()}
    selected=deepcopy(initial or {});spent=Counter()
    for (p,i),target in selected.items():spent[p]+=edit_distance(tokens[p][i]['exact_bases'],target['exact_bases'])
    local=defaultdict(list)
    for p,ts in sorted(tokens.items()):
        if p=='P0001':continue
        for i,t in enumerate(ts):
            if eligible(t,model):local[(slots[p][i],t.get('orientation','+'),t.get('role','internal'))].append((p,i))
    candidates=[]
    for bucket,members in sorted(local.items(),key=lambda x:repr(x[0])):
        observed=defaultdict(list)
        for p,i in members:observed[(tokens[p][i]['motif'],tokens[p][i]['exact_bases'])].append((p,i))
        targets=sorted(observed,key=lambda key:(min(tokens[p][i].get('edit_distance',0) for p,i in observed[key]),
            -len({m for p,i in observed[key] for m in by_id[p]['hifi_molecules']}),
            -len({m for p,i in observed[key] for m in by_id[p]['molecules']}),key))
        # One frozen observed representative per type at this local position.
        represented_types=set();targets=[key for key in targets if not (key[0] in represented_types or represented_types.add(key[0]))]
        for target_type,target_bases in targets:
            origin=observed[target_type,target_bases][0];target=tokens[origin[0]][origin[1]]
            for source_type in sorted({tokens[p][i]['motif'] for p,i in members}-{target_type}):
                pairs=[]
                for p,i in members:
                    t=tokens[p][i];d=edit_distance(t['exact_bases'],target_bases)
                    if t['motif']==source_type and 0<d<=2:pairs.append((p,i,d))
                if pairs:candidates.append(dict(bucket=bucket,target=target,source_type=source_type,pairs=pairs,
                    target_source=origin,target_rank=targets.index((target_type,target_bases))))
    def represented_tokens(p,mapping):
        ts=deepcopy(tokens[p])
        for i,t in enumerate(ts):
            if (p,i) in mapping:
                target=mapping[p,i]
                t.update(motif=target['motif'],base=target.get('base',target['motif']),template=target['template'],
                    exact_bases=target['exact_bases'],edit_distance=target.get('edit_distance',0),
                    candidate_types=[target['motif']],match_kind=target['match_kind'])
        return ts
    current_tokens={p:represented_tokens(p,selected) for p in tokens}
    views={p:path_view(ts,slots[p],reference) for p,ts in current_tokens.items()}
    index=ComplexityIndex(views)
    current=topology(views)[3]['C'];accepted=[]
    print('LOCAL_CANDIDATES',fraction,len(candidates),flush=True)
    while fraction:
        best=None
        for ci,candidate in enumerate(candidates):
            available=[];additional=Counter()
            for p,i,d in candidate['pairs']:
                if (p,i) in selected:continue
                if spent[p]+additional[p]+d<=math.floor(len(originals[p])*fraction+1e-9):
                    available.append((p,i,d));additional[p]+=d
            if not available:continue
            batches=[available]
            for batch in batches:
                trial=dict(selected);trial.update({(p,i):candidate['target'] for p,i,d in batch})
                ps={p for p,i,d in batch};newtokens={p:represented_tokens(p,trial) for p in ps}
                vv={p:path_view(ts,slots[p],reference) for p,ts in newtokens.items()}
                metrics=index.measure(vv);gain=current-metrics['C']
                if gain<=0:continue
                edits=sum(d for p,i,d in batch)
                rank=(-Fraction(gain,edits),-gain,edits,candidate['target_rank'],ci,tuple(batch))
                if best is None or rank<best[0]:best=(rank,ci,batch,trial,newtokens,vv,metrics,gain)
        if best is None:break
        _,ci,batch,selected,nt,vv,metrics,gain=best
        for p,i,d in batch:spent[p]+=d
        current_tokens.update(nt);index.update(vv)
        accepted.append(dict(candidate=ci,occurrences=batch,target_source=candidates[ci]['target_source'],
            target_type=candidates[ci]['target']['motif'],C_before=current,C_after=metrics['C'],delta_C=gain))
        current=metrics['C']
        if len(accepted)%20==0:print('LOCAL_MERGES',fraction,len(accepted),'C',current,flush=True)
        assert current==topology(views)[3]['C']
    represented={};losses=[]
    for p,ts in current_tokens.items():
        represented[p]=''.join(t['exact_bases'] for t in ts);events=[];rpos=0
        for i,t in enumerate(ts):
            old=tokens[p][i]
            for e in canonical_edits(old['exact_bases'],t['exact_bases']) if (p,i) in selected else []:
                events.append(dict(e,ostart=e['ostart']+old['start'],oend=e['oend']+old['start'],
                    rstart=e['rstart']+rpos,rend=e['rend']+rpos,occurrence=i,
                    original_type=old['motif'],represented_type=t['motif']))
            rpos+=len(t['exact_bases'])
        assert len(events)==spent[p] and len(events)<=math.floor(len(originals[p])*fraction+1e-9)
        assert restore_original(represented[p],events)==originals[p]
        if events:
            counts=Counter(e['op'] for e in events)
            losses.append(dict(path=p,loss_bp=len(events),fraction=len(events)/len(originals[p]),
                substitutions=counts['S'],insertions=counts['I'],deletions=counts['D'],events=events))
    assert represented['P0001']==originals['P0001']
    return current_tokens,views,represented,losses,selected,dict(accepted_batches=accepted,candidates=len(candidates))
