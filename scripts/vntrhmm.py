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


def decode(sequence,motifs):
    """Return a single maximum-score trace; fixed state order breaks exact ties."""
    if not sequence:return [],0.0
    if set(sequence)-set('ACGT'):raise ValueError('HMM requires ACGT input')
    states,incoming,emitting,silent,emissions=build_model(tuple(motifs))
    count=len(states);previous=[float('-inf')]*count;previous[0]=0.0
    traces=[array('i',[-1])*count]
    for base in sequence:
        scores=[float('-inf')]*count;back=array('i',[-1])*count;ems=emissions[base]
        for s in emitting:
            best=float('-inf');parent=-1
            for t,weight in incoming[s]:
                value=previous[t]+weight
                if value>best:best=value;parent=t
            scores[s]=best+ems[s];back[s]=parent
        for s in silent:
            best=float('-inf');parent=-1
            for t,weight in incoming[s]:
                value=scores[t]+weight
                if value>best:best=value;parent=t
            scores[s]=best;back[s]=parent
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
    return spans,score


def segment_hmm(task,return_details=False,**unused):
    sequence,dictionary,_=task
    records=sorted((m for m in dictionary if m.get('level','base')=='base'),key=lambda m:m['id'])
    motifs=[m['sequence'] for m in records]
    if not motifs:raise ValueError('HMM needs a frozen motif dictionary')
    spans,score=decode(sequence,motifs);tokens=[]
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
               parameters=PARAMETERS,annotation_ambiguity=False,
               ambiguity_scope='not_enumerated_single_trace',alternative=[],
               metrics=base_metrics(tokens))
    return (tokens,audit) if return_details else tokens
