"""Ordered repeat blocks on frozen path occurrences; reference-anchored sharing."""
from collections import defaultdict
from hashlib import sha256
import json

from vntrpathstructure import unit_key

CERTAIN = 'invariant_within_model'


def protected(t, interruptions):
    return (t['status'] != 'unit' or t.get('role', 'internal') != 'internal'
            or t['motif'] in interruptions or t.get('structure_assessment') != CERTAIN
            or t.get('boundary_assessment') != CERTAIN)


def select_blocks(tokens, interruptions, width_limit):
    """Minimize visible blocks, then total period width, then leftmost longest."""
    keys = [unit_key(t) for t in tokens]
    size = len(keys)
    barrier = [protected(t, interruptions) for t in tokens]
    best = [(0, 0)] * (size + 1)
    selected = {}
    stop = size
    for i in reversed(range(size)):
        if barrier[i]:
            stop = i
        candidates = [(i+1, 1)]
        for width in range(1, min(width_limit, (stop-i)//2)+1):
            word = keys[i:i+width]
            if any(width % k == 0 and word == word[:k]*(width//k) for k in range(1,width)):
                continue
            end = i+width
            while end+width <= stop and keys[end:end+width] == word:
                end += width
                candidates.append((end,width))
        end,width = min(candidates, key=lambda c:(1+best[c[0]][0],c[1]+best[c[0]][1],-c[0],c[1]))
        best[i] = (1+best[end][0],width+best[end][1])
        selected[i] = end,width
    blocks = []
    i = 0
    while i < size:
        end,width = selected[i]
        count = (end-i)//width
        kind = ('opaque' if tokens[i]['status'] != 'unit' or tokens[i].get('structure_assessment') != CERTAIN
                else 'unit' if count == 1 else 'run' if width == 1 else 'tandem_module')
        blocks.append(dict(start=i,end=end,pattern=tuple(keys[i:i+width]),count=count,
                           kind=kind,barrier=barrier[i],
                           assessed=(tokens[i].get('structure_assessment')==CERTAIN and
                                     tokens[i].get('boundary_assessment')==CERTAIN),
                           literal_anchor=tokens[i]['exact_bases'] if tokens[i]['status']!='unit' else None))
        i = end
    return blocks


def partner_options(query, reference):
    """All partners and gap positions in all maximum-cardinality monotone matches."""
    def prefix(a,b):
        d = [[0]*(len(b)+1) for _ in range(len(a)+1)]
        for i,x in enumerate(a,1):
            for j,y in enumerate(b,1):
                d[i][j] = max(d[i-1][j],d[i][j-1],d[i-1][j-1]+int(x==y))
        return d
    n,m = len(query),len(reference)
    left = prefix(query,reference)
    right = prefix(query[::-1],reference[::-1])
    optimum = left[n][m]
    options = []
    for i,x in enumerate(query):
        matches = [('ref',j) for j,y in enumerate(reference) if x==y and
                   left[i][j]+1+right[n-i-1][m-j-1]==optimum]
        gaps = [('gap',j) for j in range(m+1) if left[i][j]+right[n-i-1][m-j]==optimum]
        options.append(matches+gaps)
    return options


def align_regions(blocks, ref, path, outer_certain):
    # Align protected landmarks first; repeated, ambiguously placed barriers
    # constrain the search only where correspondence is unique. A novel
    # insertion is shareable only at a unique gap with the same ordered context.
    def keys(bs,owner):
        return [(b['pattern'],b['barrier'],b['literal_anchor']) if b['assessed']
                else ('unassessed_or_ambiguous',owner,i) for i,b in enumerate(bs)]
    qkeys,rkeys=keys(blocks,path),keys(ref,'P0001')
    def landmarks(bs):
        return [i for i,b in enumerate(bs) if b['barrier'] and b['assessed']]
    qmarks,rmarks = landmarks(blocks),landmarks(ref)
    local = [('local',path,i) for i in range(len(blocks))]
    if not outer_certain:
        return local, [[] for _ in blocks]
    options = [[] for _ in blocks]
    if path=='P0001':
        return [('ref',i) for i in range(len(blocks))], [[('ref',i)] for i in range(len(blocks))]
    matched = {-1:-1,len(blocks):len(ref)}
    barrier_options = partner_options([qkeys[i] for i in qmarks],[rkeys[i] for i in rmarks])
    for i,opts in zip(qmarks,barrier_options):
        options[i] = [(k,rmarks[j] if k=='ref' or j<len(rmarks) else len(ref)) for k,j in opts]
        if len(opts)==1 and opts[0][0]=='ref':
            j = rmarks[opts[0][1]]
            matched[i] = j
            local[i] = ('ref',j)
    anchors = sorted(matched)
    for a,b in zip(anchors,anchors[1:]):
        lo,hi = matched[a]+1,matched[b]
        query = qkeys[a+1:b]
        opts = partner_options(query,rkeys[lo:hi])
        for offset,op in enumerate(opts):
            i = a+1+offset
            options[i] = [(k,j+lo) for k,j in op]
            if len(op)!=1: continue
            k,j = options[i][0]
            if k=='ref':
                local[i] = (k,j)
            else:
                # The entire adjacent insertion group fixes repeated positions.
                l = r = offset
                while l and opts[l-1]==op: l-=1
                while r+1<len(opts) and opts[r+1]==op: r+=1
                local[i] = ('gap',j,tuple(query[l:r+1]),offset-l)
    return local, options


def make_structure(records,tokens,audits,model,model_hash,policy,policy_hash):
    by_id = {r['id']:r for r in records}
    interruptions = set(model.get('interruptions',[]))
    blocks = {p:select_blocks(ts,interruptions,policy['max_module_width']) for p,ts in tokens.items()}
    nodes,visits = {},{}
    edge_members = defaultdict(list)
    for p,bs in sorted(blocks.items()):
        r = by_id[p]
        outer = (model['boundary_status']=='confirmed' and
                 all(o.get('boundary_status')=='confirmed' and o.get('core_status')=='success'
                     and not o.get('core_ambiguous',True) for o in r['observations']))
        regions,options = align_regions(bs,blocks['P0001'],p,outer)
        visits[p] = []
        for index,(block,region) in enumerate(zip(bs,regions)):
            i,j = block['start'],block['end']
            ts = tokens[p][i:j]
            structure = ts[0].get('structure_assessment','not_assessed')
            boundary = ts[0].get('boundary_assessment','not_assessed') if outer else 'uncertain_outer'
            if structure!=CERTAIN or boundary!=CERTAIN:
                region=('local',p,index)
            identity = [model['locus'],model_hash,policy_hash,region,block['pattern'],block['barrier']]
            digest = sha256(json.dumps(identity,separators=(',',':')).encode()).hexdigest()[:24]
            node_id = 'SN_'+digest
            visit_id = f'{p}:v{index}'
            n = nodes.setdefault(node_id,dict(id=node_id,region_id=region,ordered_pattern=block['pattern'],
                kind=block['kind'],structure_assessment=structure,boundary_assessment=boundary,entry_port=node_id+':L',exit_port=node_id+':R',protected=block['barrier'],
                member_visit_ids=[],observed_repeat_counts=set(),path_ids=set()))
            if block['kind']=='run': n['kind']='run'
            n['member_visit_ids'].append(visit_id)
            n['path_ids'].add(p)
            n['observed_repeat_counts'].add(block['count'])
            homologs = audits[p]['primary_partner_options'][i:j]
            visit = dict(visit_id=visit_id,path_id=p,supernode_id=node_id,token_interval=[i,j],
                repeat_count=block['count'],ordered_occurrence_ids=[f'{p}:u{k}' for k in range(i,j)],
                original_interval=[ts[0]['start'],ts[-1]['end']],represented_interval=[ts[0]['start'],ts[-1]['end']],
                boundary_assessment=boundary,structure_assessment=structure,
                homology_assessment='multiple' if any(len(v)>1 for v in homologs) else 'unique',
                region_partner_options=options[index],position_assessment='path_local' if region[0]=='local' else 'aligned')
            if visits[p]: edge_members[visits[p][-1]['supernode_id'],node_id].append(p)
            visits[p].append(visit)
    def support(paths):
        return len({m for p in paths if by_id[p]['role']!='reference' for m in by_id[p]['molecules']})
    for n in nodes.values():
        n['path_ids']=sorted(n['path_ids'])
        n['path_count']=len(n['path_ids'])
        n['molecule_count']=support(n['path_ids'])
        n['visit_count']=len(n['member_visit_ids'])
        n['observed_repeat_counts']=sorted(n['observed_repeat_counts'])
    edges=[dict(source=a,target=b,path_ids=sorted(set(ps)),path_count=len(set(ps)),
                visit_count=len(ps),molecule_count=support(set(ps))) for (a,b),ps in sorted(edge_members.items())]
    return dict(nodes=nodes,edges=edges,locus=model['locus'],schema_version=1,
                pattern_fields=['motif','orientation','phase','status','role']),visits

