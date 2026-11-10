"""Structure-first path graphs with frozen annotations and reversible display loss.

Repeat counts describe annotations, not genotypes. All motif identities are
protected; lossy edits are limited to explicitly permitted internal positions.
"""
from collections import Counter, defaultdict
from copy import deepcopy
from hashlib import sha256
import json

from vntrpathgraph import digest, validate_graph, compact_graph
from vntrpartition import load_model
from vntrstructure import segment

VERSION = '2.0-structure-first'


def unit_key(t):
    return (t['motif'], t.get('orientation', '+'), t.get('phase', 0),
            t['status'], t.get('role', 'internal'))


def describe(tokens, names):
    """Ordered run/module labels; no periodicity inferred through unknown blocks."""
    labels = []
    for t in tokens:
        if t['status'] != 'unit':
            labels.append(f'{t["status"].upper()}({len(t["exact_bases"])}bp)')
        else:
            labels.append(names.get(t['motif'], t['motif']) +
                          ('~' if t.get('edit_distance', 0) else ''))
    parts = []; i = 0
    while i < len(labels):
        choices = []
        for width in range(1, min(12, (len(labels)-i)//2)+1):
            pattern = labels[i:i+width]
            if any(tokens[j]['status'] != 'unit' for j in range(i, i+width)): continue
            count = 1
            while labels[i+count*width:i+(count+1)*width] == pattern: count += 1
            if count > 1: choices.append((width*count, -width, count, pattern))
        if choices:
            span, neg_width, count, pattern = max(choices)
            text = pattern[0] if len(pattern) == 1 else '('+' '.join(pattern)+')'
            parts.append(text+'^'+str(count)); i += span
        else:
            parts.append(labels[i]); i += 1
    return ' '.join(parts)


def summarize(sequence, tokens, audit, model):
    position = 0
    for t in tokens:
        if t['start'] != position or sequence[t['start']:t['end']] != t['exact_bases']:
            raise ValueError('Partition does not preserve the original sequence')
        position = t['end']
    if position != len(sequence): raise ValueError('Incomplete partition')
    names = model.get('display_labels', {})
    complete = all(t['status'] == 'unit' for t in tokens)
    tied = audit.get('annotation_ambiguity', False)
    runs = []
    for i, t in enumerate(tokens):
        key = unit_key(t)
        if runs and t['status'] == 'unit' and key == runs[-1]['key']:
            runs[-1]['count'] += 1; runs[-1]['end'] = i+1
        else: runs.append(dict(key=key, count=1, start=i, end=i+1))
    interruptions = [dict(unit=i+1, motif=t['motif'], start=t['start'], end=t['end'])
                     for i,t in enumerate(tokens) if t['motif'] in model.get('interruptions', [])]
    signature = [list(unit_key(t)) for t in tokens]
    edited = Counter((t['motif'], t['exact_bases']) for t in tokens if t.get('edit_distance', 0))
    return dict(length=len(sequence), tokens=tokens, annotation=audit, signature=signature,
                expression=describe(tokens, names), runs=[dict(r, key=list(r['key'])) for r in runs],
                complete_units=sum(t['status']=='unit' for t in tokens),
                count_status='unresolved' if tied or not complete else 'model_resolved',
                interruptions=interruptions,
                longest_run=max((r['count'] for r in runs if r['key'][3]=='unit'), default=0),
                unknown_bases=sum(len(t['exact_bases']) for t in tokens if t['status']=='unknown'),
                partial_bases=sum(len(t['exact_bases']) for t in tokens if t['status']=='partial'),
                merge_eligible=complete and not tied and not audit.get('annotation_conflict', False),
                repeated_edited_subtype=any(n>1 for n in edited.values()))


def annotate_paths(records, model_path, locus, imported_rows=None, imported_audit=None,select_representative=False,export_optimal=False):
    model, dictionary = load_model(model_path, locus)
    model_hash = sha256(model_path.read_bytes()).hexdigest()
    imported = defaultdict(list)
    if imported_rows:
        if (imported_audit or {}).get('model_sha256') != model_hash:
            raise ValueError('Imported annotations do not use the requested frozen model')
        for row in imported_rows:
            for f in row['fragments']:
                if 'tokens' in f:
                    imported[(row['sample'], row['platform'], f['id'])].append(f)
    features = {}
    for record in records:
        sequence = record['sequence']; candidates = []
        for o in record['observations']:
            candidates += imported.get((o['sample'], o['platform'], o['fragment']), [])
        if candidates and not select_representative:
            if any(f['sequence'] != sequence for f in candidates): raise ValueError('Annotation sequence mismatch')
            tokens = deepcopy(candidates[0]['tokens']); audit = deepcopy(candidates[0].get('annotation', {}))
            signatures = {json.dumps([(unit_key(t), t['start'], t['end']) for t in f['tokens']]) for f in candidates}
            audit['annotation_conflict'] = len(signatures)>1
            audit['annotation_ambiguity'] = any(f.get('annotation', {}).get('annotation_ambiguity', False) for f in candidates)
            audit['source'] = 'imported_frozen_partition'
        else:
            # A frozen common annotation profile makes the same sequence have the
            # same decomposition regardless of its platform mixture.
            terminal = next((t for t in model.get('terminal_templates', [])
                             if record.get('role')!='reference' and record['samples'] and
                             set(record['samples']).issubset(t['samples']) and
                             sequence.endswith(t['sequence'])), None)
            stop = len(sequence)-len(terminal['sequence']) if terminal else len(sequence)
            tokens, audit = segment((sequence[:stop], dictionary, 'HiFi'), return_details=True, collect_optimal=True,
                                    select_representative=select_representative,export_optimal=export_optimal)
            if terminal:
                token = dict(start=stop, end=len(sequence), exact_bases=sequence[stop:],
                             template=terminal['sequence'], motif=terminal['id'], base=terminal['id'],
                             orientation='+', phase=0, status='unit', role='right_terminal',
                             edit_distance=0, annotation_ambiguity=False)
                tokens.append(token)
                if audit.get('alternative'): audit['alternative'].append(deepcopy(token))
                audit['terminal_scope'] = 'frozen_sample_specific_template_compatibility_not_new_support'
            audit['source'] = 'existing_seed_dp_v2'; audit['annotation_profile'] = 'HiFi_cost_profile'
        features[record['id']] = summarize(sequence, tokens, audit, model)
    return model, model_hash, features


def edit_script(original, represented, allowed=None):
    """Minimum base edits; detect any optimal placement touching forbidden sites.

I adds to represented, D removes from original. Coordinates refer to the
unchanged strings. allowed is in represented/template coordinates.
"""
    n,m = len(original),len(represented)
    cost = [[0]*(m+1) for _ in range(n+1)]
    unsafe = [[False]*(m+1) for _ in range(n+1)]
    def permitted(op, j):
        if allowed is None or op=='=': return True
        return (j-1 in allowed and j in allowed) if op=='D' else j in allowed
    def predecessors(i,j):
        result=[]
        if i and j:
            op='=' if original[i-1]==represented[j-1] else 'S'
            result.append((cost[i-1][j-1]+(op!='='),i-1,j-1,op))
        if i: result.append((cost[i-1][j]+1,i-1,j,'D'))
        if j: result.append((cost[i][j-1]+1,i,j-1,'I'))
        return result
    for i in range(n+1):
        for j in range(m+1):
            if not i and not j: continue
            choices=predecessors(i,j); best=min(x[0] for x in choices); cost[i][j]=best
            unsafe[i][j]=any(unsafe[a][b] or not permitted(op,b)
                             for score,a,b,op in choices if score==best)
    events=[]; i,j=n,m
    while i or j:
        _,a,b,op=next(x for x in predecessors(i,j) if x[0]==cost[i][j])
        if op!='=':
            events.append(dict(op=op, ostart=a, oend=i, rstart=b, rend=j,
                               original=original[a:i], represented=represented[b:j]))
        i,j=a,b
    events.reverse()
    return events, unsafe[n][m]


def restore_original(represented, events):
    parts=[]; cursor=0
    for event in events:
        a,b=event['rstart'],event['rend']
        if a<cursor or represented[a:b]!=event['represented']: raise ValueError('Invalid residual coordinates')
        parts += [represented[cursor:a], event['original']]; cursor=b
    return ''.join(parts)+represented[cursor:]


def compare_for_merge(source, target, sf, tf, policy, budget):
    if sf['signature']!=tf['signature']: return None,'protected_structure'
    if not sf['merge_eligible'] or not tf['merge_eligible']: return None,'unresolved_partition'
    if sf['repeated_edited_subtype'] or tf['repeated_edited_subtype']: return None,'repeated_edited_subtype'
    if abs(len(source)-len(target))>budget: return None,'budget'
    events=[]
    for i,(a,b) in enumerate(zip(sf['tokens'],tf['tokens'])):
        if a['exact_bases']==b['exact_bases']: continue
        if i in (0,len(sf['tokens'])-1) or a.get('role')=='right_terminal': return None,'protected_end'
        if b['exact_bases']!=b['template']: return None,'representative_unit_not_template'
        allowed=set(policy.get('allowed_template_positions',{}).get(a['motif'],[]))
        local,unsafe=edit_script(a['exact_bases'],b['exact_bases'],allowed)
        if unsafe: return None,'protected_position_or_ambiguous_edit'
        for e in local:
            events.append(dict(e, ostart=e['ostart']+a['start'], oend=e['oend']+a['start'],
                               rstart=e['rstart']+b['start'], rend=e['rend']+b['start'], unit=i+1))
    if not events or len(events)>budget: return None,'budget'
    if restore_original(target,events)!=source: raise ValueError('Residual reconstruction failed')
    return events,'allowed_internal_edit'


def merge_paths(records, features, model, policy, budget):
    if not isinstance(budget,int) or not 0<=budget<=2: raise ValueError('Use a 0, 1 or 2 bp whole-path budget')
    if policy.get('schema_version')!=1 or policy.get('model_version')!=model['version']:
        raise ValueError('Structure policy/model mismatch')
    originals={r['id']:r['sequence'] for r in records}; represented=dict(originals)
    representatives={p:p for p in originals}; losses=[]; decisions=[]
    ordered=sorted(records,key=lambda r:(-len(r['hifi_molecules']),-len(r['molecules']),len(r['sequence']),r['sha256']))
    buckets=defaultdict(list)
    for r in ordered: buckets[json.dumps(features[r['id']]['signature'])].append(r)
    for candidates in buckets.values():
        unassigned={r['id'] for r in candidates}
        for target in candidates:
            if target['id'] not in unassigned: continue
            unassigned.remove(target['id'])
            if not budget: continue
            for source in candidates:
                p=source['id']
                if p not in unassigned: continue
                if source.get('role')=='reference' or target.get('role')=='reference': continue
                events,reason=compare_for_merge(source['sequence'],target['sequence'],features[p],features[target['id']],policy,budget)
                decisions.append(dict(path=p,target_path=target['id'],accepted=events is not None,reason=reason))
                if events is None: continue
                represented[p]=target['sequence']; representatives[p]=target['id']; unassigned.remove(p)
                counts=Counter(e['op'] for e in events)
                losses.append(dict(path=p,target_path=target['id'],loss_bp=len(events),
                    substitutions=counts['S'],insertions=counts['I'],deletions=counts['D'],
                    original_length=len(source['sequence']),represented_length=len(target['sequence']),
                    fraction=len(events)/len(source['sequence']),events=events,
                    positions=sorted({min(e['rstart'],len(target['sequence'])-1) for e in events}),
                    policy_version=policy['version'],reason=reason))
    return represented,representatives,losses,decisions


def blocks(feature, model):
    tokens=feature['tokens']; names=model.get('display_labels',{})
    if feature['annotation'].get('annotation_ambiguity') or feature['annotation'].get('annotation_conflict'):
        intervals=feature['annotation'].get('ambiguity_intervals')
        if intervals is None or feature['annotation'].get('annotation_conflict'):
            sequence=''.join(t['exact_bases'] for t in tokens)
            return [dict(key=('unresolved_partition',),sequence=sequence,tokens=tokens,
                         label=f'PARTITION UNRESOLVED ({len(sequence)}bp)',barrier=True)]
        normalized=[]
        for t in tokens:
            unresolved=any(t['start']<b and a<t['end'] for a,b in intervals)
            if unresolved and normalized and normalized[-1]['status']=='unresolved':
                normalized[-1]['end']=t['end'];normalized[-1]['exact_bases']+=t['exact_bases']
            elif unresolved:
                normalized.append(dict(t,motif='UNRESOLVED',status='unresolved',template='',edit_distance=0))
            else:normalized.append(t)
        tokens=normalized
    result=[]
    for t in tokens:
        key=unit_key(t)
        if result and t['status']=='unit' and key==result[-1]['key'] and t.get('role')!='right_terminal':
            result[-1]['tokens'].append(t); result[-1]['sequence']+=t['exact_bases']
        else: result.append(dict(key=key,sequence=t['exact_bases'],tokens=[t],
                                 barrier=t['status']!='unit' or t.get('role')=='right_terminal'))
    for b in result:
        b['label']=(f'PARTITION UNRESOLVED ({len(b["sequence"])}bp)' if b['tokens'][0]['status']=='unresolved'
                    else describe(b['tokens'],names))
    return result


def align_blocks(query, reference):
    """Unit-run edit DP. Count changes stay in whole run nodes, not repeat loops."""
    n,m=len(query),len(reference)
    cost=[[0]*(m+1) for _ in range(n+1)]; ways=[[0]*(m+1) for _ in range(n+1)]; ways[0][0]=1
    def choices(i,j):
        result=[]
        if i and j: result.append((cost[i-1][j-1]+(query[i-1]['key']!=reference[j-1]['key']),i-1,j-1))
        if i: result.append((cost[i-1][j]+1,i-1,j))
        if j: result.append((cost[i][j-1]+1,i,j-1))
        return result
    for i in range(n+1):
        for j in range(m+1):
            if not i and not j: continue
            c=choices(i,j); cost[i][j]=min(x[0] for x in c)
            ways[i][j]=min(2,sum(ways[a][b] for s,a,b in c if s==cost[i][j]))
    steps=[];i,j=n,m
    while i or j:
        _,a,b=next(x for x in choices(i,j) if x[0]==cost[i][j])
        steps.append((a if a<i else None,b if b<j else None));i,j=a,b
    return list(reversed(steps)),ways[n][m]>1


def build_structural_graph(records, features, model, represented, representatives, backbone_path=None):
    if not records: return dict(nodes={},edges=set(),paths={},coordinates={},labels={},alignments={})
    ordered=sorted(records,key=lambda r:(-len(r['hifi_molecules']),-len(r['molecules']),len(r['sequence']),r['sha256']))
    backbone_path=backbone_path or next((r['id'] for r in records if r.get('role')=='reference'),ordered[0]['id'])
    reference=blocks(features[representatives[backbone_path]],model)
    descriptions={}; keyed_paths={}; alignment_audit={}
    def add(slot,b):
        key=(slot,b['key'],b['sequence'],b['label'])
        descriptions[key]=b
        return key
    for record in ordered:
        p=record['id']; query=blocks(features[representatives[p]],model)
        left=0
        def same(a,b):return a['key']==b['key'] and a['sequence']==b['sequence']
        while left<min(len(query),len(reference)) and same(query[left],reference[left]):left+=1
        right=0
        while right<min(len(query)-left,len(reference)-left) and same(query[-1-right],reference[-1-right]):right+=1
        walk=[add((i,1,0,0),query[i]) for i in range(left)]
        qe=len(query)-right; re=len(reference)-right
        qmid=query[left:qe]; rmid=reference[left:re]
        steps,tied=align_blocks(qmid,rmid)
        if tied and qmid:
            ts=[t for b in qmid for t in b['tokens']]; sequence=''.join(b['sequence'] for b in qmid)
            partition_unresolved=any('PARTITION UNRESOLVED' in b['label'] for b in qmid)
            b=dict(key=('correspondence_unresolved',left,re),sequence=sequence,tokens=ts,
                   label=(f'PARTITION UNRESOLVED ({len(sequence)}bp)' if partition_unresolved else
                          'ALIGNMENT ? | '+describe(ts,model.get('display_labels',{}))),barrier=True)
            walk.append(add((left,1,re-left,0),b))
        else:
            j=left; insertion=0
            for qi,rj in steps:
                if qi is not None:
                    slot=(j,0,len(qmid),insertion) if rj is None else (left+rj,1,0,0)
                    walk.append(add(slot,qmid[qi])); insertion+=1
                if rj is not None:j=left+rj+1;insertion=0
        walk += [add((len(reference)-right+i,1,0,0),query[qe+i]) for i in range(right)]
        keyed_paths[p]=walk
        alignment_audit[p]=dict(correspondence_ambiguous=tied,shared_prefix_runs=left,shared_suffix_runs=right,
                               convention='exact_boundary_runs_then_run_type_DP; tied_middle_kept_whole')
    keys=sorted(descriptions,key=lambda k:(k[0],repr(k[1:])))
    ids={k:i+1 for i,k in enumerate(keys)}
    paths={p:[ids[k] for k in walk] for p,walk in keyed_paths.items()}
    raw=dict(nodes={ids[k]:descriptions[k]['sequence'] for k in keys},paths=paths,
             edges={edge for walk in paths.values() for edge in zip(walk,walk[1:])},
             coordinates={ids[k]:[list(k[0])] for k in keys},alignments=alignment_audit,
             barriers={ids[k] for k in keys if descriptions[k]['barrier']})
    validate_graph(raw,represented)
    graph=compact_graph(raw); by_node={ids[k]:descriptions[k] for k in keys}
    graph['labels']={}; graph['structure_mode']=True; graph['node_features']={}
    graph['backbone_path']=backbone_path
    for n,chain in graph['source_nodes'].items():
        bs=[by_node[x] for x in chain]; ts=[t for b in bs for t in b['tokens']]
        label=bs[0]['label'] if len(bs)==1 else describe(ts,model.get('display_labels',{}))
        graph['labels'][n]=label
        graph['node_features'][n]=dict(expression=label,
            unit_count=None if 'PARTITION UNRESOLVED' in label else sum(t['status']=='unit' for t in ts),
            uncertain=any(b['barrier'] and 'UNRESOLVED' in b['label'] or 'ALIGNMENT ?' in b['label'] for b in bs),
            unit_keys=[list(unit_key(t)) for t in ts])
    graph['dictionary_legend']=[model.get('display_labels',{}).get(m['id'],m['id'])+'='+m['sequence'] for m in model['motifs']]
    graph['dictionary_legend'] += [t['id']+'='+t['sequence']+' (right terminal only)' for t in model.get('terminal_templates',[])]
    graph['path_expressions']={p:features[representatives[p]]['expression'] for p in paths}
    validate_graph(graph,represented)
    return graph


def select_subgraph(graph,names):
    paths={p:graph['paths'][p] for p in names}; used={n for w in paths.values() for n in w}
    result=dict(graph,nodes={n:s for n,s in graph['nodes'].items() if n in used},paths=paths,
                edges={e for w in paths.values() for e in zip(w,w[1:])},
                coordinates={n:graph['coordinates'][n] for n in used})
    return result
