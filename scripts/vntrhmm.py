"""Independent motif-profile HMM with an explicit skip state and one Viterbi trace.

The motif/skip design follows the TRGT idea; this is not a port of TRGT code.
M/I consume observed bases, D consumes template positions, and B joins modules.
Every token retains its observed string; decoding never rewrites bases.
"""
from array import array
from functools import lru_cache
from math import log

from vntrcore import edit_distance, canonical_motif, motif_orientation
from vntrmotifs import base_metrics

VERSION = 'motif-hmm-single-v1'
PARAMETERS = dict(match_emission=.90, mismatch_emission=.10/3,
                  match_step=.90, indel_open=.05, insertion_extend=.25,
                  deletion_extend=.50, module_enter=.50,
                  skip_enter=.50, skip_extend=.50, skip_exit=.50)


@lru_cache(maxsize=16)
def build_model(motifs):
    # Each state is (motif index, kind, template position), with incoming edges.
    states=[(-1,'B',-1),(-1,'X',-1)];incoming=[[],[(0,log(.5)),(1,log(.5))]]
    incoming[0].append((1,log(.5)))
    emitting=[1];silent=[]
    for mi,motif in enumerate(motifs):
        length=len(motif);matches=[];inserts=[];deletes={}
        for kind,dest in [('M',matches),('I',inserts)]:
            for j in range(length):
                dest.append(len(states));states.append((mi,kind,j));incoming.append([])
                emitting.append(dest[-1])
        for j in range(1,length):
            deletes[j]=len(states);states.append((mi,'D',j));incoming.append([])
            silent.append(deletes[j])
        for j,m in enumerate(matches):
            # A leading deletion run can enter at a later template position.
            entry=.5*(.9 if j==0 else .05*.5**(j-1))
            incoming[m].append((0,log(entry)))
            if j:
                incoming[m].extend([(matches[j-1],log(.9)),(inserts[j-1],log(.75))])
                if j>1:incoming[m].append((deletes[j-1],log(.5)))
            incoming[inserts[j]]=[(m,log(.05)),(inserts[j],log(.25))]
            if j:
                incoming[deletes[j]]=[(matches[j-1],log(.05))]
                if j>1:incoming[deletes[j]].append((deletes[j-1],log(.5)))
        incoming[0].extend([(matches[-1],log(.9)),(inserts[-1],log(.75))])
        if length>1:incoming[0].append((deletes[length-1],0.0))
    emissions={}
    for base in 'ACGT':
        emissions[base]=[log(.25) if kind in ('I','X') else
            log(.9 if motifs[mi][j]==base else .1/3) if kind=='M' else 0.0
            for mi,kind,j in states]
    return states,incoming,emitting,silent+[0],emissions


def decode(sequence,motifs,assess=False):
    """Return a single maximum-score trace; fixed state order breaks exact ties."""
    if not sequence:return ([],0.0,[],[0,0]) if assess else ([],0.0)
    if set(sequence)-set('ACGT'):raise ValueError('HMM requires ACGT input')
    states,incoming,emitting,silent,emissions=build_model(tuple(motifs))
    count=len(states);previous=[float('-inf')]*count;previous[0]=0.0
    traces=[array('i',[-1])*count]
    tied_rows = [{}]
    for base in sequence:
        ties = {}
        scores=[float('-inf')]*count;back=array('i',[-1])*count;ems=emissions[base]
        for s in emitting:
            best=float('-inf');parent=-1
            for t,weight in incoming[s]:
                value=previous[t]+weight
                if value>best:best=value;parent=t
            scores[s]=best+ems[s];back[s]=parent
            if assess and parent >= 0:
                alternatives = [t for t, weight in incoming[s] if t != parent
                                and best - (previous[t] + weight) <= 1e-10]
                if alternatives:
                    ties[s] = alternatives
        for s in silent:
            best=float('-inf');parent=-1
            for t,weight in incoming[s]:
                value=scores[t]+weight
                if value>best:best=value;parent=t
            scores[s]=best;back[s]=parent
            if assess and parent >= 0:
                alternatives = [t for t, weight in incoming[s] if t != parent
                                and best - (scores[t] + weight) <= 1e-10]
                if alternatives:
                    ties[s] = alternatives
        if assess: tied_rows.append(ties)
        previous=scores;traces.append(back)
    score=previous[0];i=len(sequence);s=0;trace=[]
    while i or s:
        trace.append((s,i));parent=traces[i][s]
        if parent<0:raise RuntimeError('No complete HMM trace')
        if states[s][1] in ('M','I','X'):i-=1
        s=parent
    trace.reverse()
    spans=[];start=0;mi=None
    for s,end in trace:
        motif,kind,j=states[s]
        if kind=='B':
            if end>start:spans.append((start,end,mi))
            start=end;mi=None
        elif mi is None:mi=motif
    if not assess:
        return spans,score
    # Follow only cells reachable from the optimum terminal. Each emitted base
    # projects to a motif and a token-start bit, independent of internal M/I/D.
    options = [set() for _ in sequence]
    reachable = {(len(sequence), 0)}
    work = list(reachable)
    while work:
        end, state = work.pop()
        if (end, state) == (0, 0):
            continue
        motif, kind, _ = states[state]
        emits = kind in ('M', 'I', 'X')
        for parent in [traces[end][state], *tied_rows[end].get(state, ())]:
            if emits:
                options[end-1].add((motif, parent == 0))
            cell = (end-int(emits), parent)
            if cell not in reachable:
                reachable.add(cell)
                work.append(cell)
    # Count complete motif tokens along complete traces (SKIP is excluded).
    rank = {state:i for i,state in enumerate(emitting + silent)}
    ranges = {(0,0):(0,0)}
    for end,state in sorted(reachable, key=lambda c:(c[0], rank[c[1]])):
        if (end,state) == (0,0):
            continue
        kind = states[state][1]
        parents = [traces[end][state], *tied_rows[end].get(state, ())]
        bounds = [ranges[end-int(kind in ('M','I','X')), t] for t in parents]
        increments = [int(kind == 'B' and states[t][0] >= 0) for t in parents]
        ranges[end,state] = (min(b[0]+v for b,v in zip(bounds,increments)),
                             max(b[1]+v for b,v in zip(bounds,increments)))
    return spans,score,options,list(ranges[len(sequence),0])


def segment_hmm(task,return_details=False,assess_structure=False,**unused):
    sequence,dictionary,_=task
    records=sorted((m for m in dictionary if m.get('level','base')=='base'),key=lambda m:m['id'])
    motifs=[m['sequence'] for m in records]
    if not motifs:raise ValueError('HMM needs a frozen motif dictionary')
    if set(sequence)-set('ACGT'):
        if set(sequence)-set('ACGTRYSWKMBDHVN'):
            raise ValueError('Invalid DNA alphabet')
        opaque=dict(start=0,end=len(sequence),exact_bases=sequence,status='unknown',
                    motif='U',base='U',template='',orientation='?',phase=0,edit_distance=0,
                    annotation_ambiguity=False,boundary_assessment='not_assessed',
                    structure_assessment='not_assessed')
        audit=dict(source=VERSION,decoder='opaque_iupac',ambiguity_scope='not_assessed',
                   annotation_ambiguity=False,alternative=[],metrics=base_metrics([opaque]))
        return ([opaque],audit) if return_details else [opaque]
    decoded=decode(sequence,motifs,assess=assess_structure)
    spans,score=decoded[:2];tokens=[]
    for start,end,mi in spans:
        bases=sequence[start:end]
        if mi==-1:
            if tokens and tokens[-1]['status']=='unknown':
                tokens[-1]['end']=end;tokens[-1]['exact_bases']+=bases;continue
            t=dict(status='unknown',motif='U',base='U',template='',orientation='?',
                   family='UNKNOWN',phase=0,canonical_orientation='?',edit_distance=0)
        else:
            record=records[mi];template=record['sequence'];family=canonical_motif(template)
            phase,orientation=motif_orientation(template,family)
            t=dict(status='unit',motif=record['id'],base=record['base'],template=template,
                   orientation='+',family=family,phase=phase,canonical_orientation=orientation,
                   edit_distance=edit_distance(template,bases))
        t.update(start=start,end=end,exact_bases=bases,annotation_ambiguity=False,
                 phase_policy='fixed_template_phase',phase_reference='canonical_family',
                 orientation_reference='catalog_template')
        tokens.append(t)
    assert ''.join(t['exact_bases'] for t in tokens)==sequence
    audit=dict(source=VERSION,decoder='single_viterbi',viterbi_score=score,
               parameters=PARAMETERS,annotation_ambiguity=False,structure_assessment='not_assessed',
               ambiguity_scope='not_enumerated_single_trace',alternative=[],
               metrics=base_metrics(tokens))
    if assess_structure:
        options, count_range = decoded[2:]
        for t in tokens:
            start,end=t['start'],t['end']
            # Both outer cuts must occur in every complete optimal projection.
            certain = all(b for _,b in options[start]) and (
                end == len(options) or all(b for _,b in options[end]))
            invariant = certain and all(len(opts) == 1 for opts in options[start:end])
            t['boundary_assessment'] = 'invariant_within_model' if certain else 'ambiguous'
            t['structure_assessment'] = 'invariant_within_model' if invariant else 'ambiguous'
            t['annotation_ambiguity'] = not invariant
        audit.update(ambiguity_scope='optimal_boundary_and_motif_projections',
                     structure_assessment='ambiguous' if any(t['annotation_ambiguity'] for t in tokens) else 'invariant_within_model',
                     tie_tolerance=1e-10, complete_motif_count_range=count_range,
                     annotation_ambiguity=any(t['annotation_ambiguity'] for t in tokens),
                     alternative_base_projections=[
                         dict(offset=i, choices=[[records[m]['id'] if m>=0 else 'SKIP',b]
                                                for m,b in sorted(opts)])
                         for i,opts in enumerate(options) if len(opts)>1])
    return (tokens,audit) if return_details else tokens
