"""Frozen, sample-balanced multi-scale flank discovery inside the supplied window."""
import argparse
from collections import defaultdict
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from vntrflankalign import (VERSION, SCORING, NORMALIZATION, canonical_hash, digest,
    cigar_segments, flank_search, all_exact, score_flank, traceback)

WINDOWS = dict(W1=[100,10], W2=[100,50,10], W3=[100,25,10], W4=[100,50,25,10],
               W5=[150,10], W6=[150,50,10], W7=[150,25,10], W8=[150,50,25,10])
SAMPLES = ['LCL5','LCL6','LCL7','LCL8']
FETCH = dict(ATXN2=[111598450,111599519], D1S80=[2458776,2460240])
MODEL_FILE = Path(__file__).resolve().parents[1]/'catalogs/vclust_flank_models.3469834.json'


def load_snapshot(path):
    path = Path(path)
    meta = json.loads((path.parent/'summary.json').read_text(encoding='utf-8'))['metadata']
    sources = {x['bam']:(x['sample'], x['platform']) for x in meta['bam_audits']}
    rows = []
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            row = json.loads(line)
            row['sample'], row['platform'] = sources[row['bam']]
            rows.append(row)
    return rows


def quality(record, minimum):
    return not record['flag'] & (4|256|512|1024|2048) and minimum <= record['mapq'] < 255


def discovery_records(records):
    groups = defaultdict(list)
    for r in records:
        if r['sample'] in SAMPLES and r['platform']=='HiFi' and quality(r,50):
            groups[r['molecule_key']].append(r)
    selected, conflicts = [], []
    for key, group in sorted(groups.items()):
        signatures = {(r['chrom'],r['pos'],r['cigar'],r['query_sequence']) for r in group}
        if len(signatures) != 1:
            conflicts.append(key)
        else:
            selected.append(min(group,key=lambda r:(r['record_id'],r['bam'])))
    return selected, conflicts


def build_profile(records, fetch_interval, dna):
    fl, fr = fetch_interval
    chosen, conflicts = discovery_records(records)
    cov = {s:np.zeros(fr-fl,dtype=np.int32) for s in SAMPLES}
    alt = {s:np.zeros(fr-fl,dtype=np.float64) for s in SAMPLES}
    for record in chosen:
        c, a, query = cov[record['sample']], alt[record['sample']], record['query_sequence'].upper()
        for r0,r1,q0,q1,op in cigar_segments(record):
            begin,end = max(r0,fl),min(r1,fr)
            if op in 'M=XD' and begin<end:
                x,y=begin-fl,end-fl
                c[x:y] += 1
                if op in 'XD':
                    a[x:y] += 1
                elif op=='M':
                    for r in range(begin,end):
                        a[r-fl] += int(query[q0+r-r0] != dna[r-fl] or dna[r-fl] not in 'ACGT')
            elif op in 'IS' and fl<=r0<fr:
                a[r0-fl] += q1-q0 if op=='I' else 1
    depths = {s:float(cov[s].mean()) for s in SAMPLES}
    values = sum(alt[s]/max(depths[s],1.0) for s in SAMPLES)/len(SAMPLES)
    # The first interval has an inclusive upper endpoint, exactly as vclust.
    bins=np.where(values<=.10,0,np.where(values<.25,1,np.where(values<.75,2,
          np.where(values<1.5,3,np.where(values<5,4,5))))).astype(np.int8)
    model=json.loads(MODEL_FILE.read_text(encoding='utf-8'))
    return dict(records=chosen,conflicts=conflicts,coverage=cov,alts=alt,depths=depths,
                values=values,bins=bins,fetch_interval=list(fetch_interval),dna=dna,
                model=model,model_sha256=sha256(MODEL_FILE.read_bytes()).hexdigest())


def assess_window(profile, boundary, side, width):
    fl,fr=profile['fetch_interval']
    start,end=(boundary-width,boundary) if side=='left' else (boundary,boundary+width)
    if not fl<=start<end<=fr:
        return dict(status='outside_fetch', probability=None, minimum_coverage={})
    x,y=start-fl,end-fl
    coverage={s:int(profile['coverage'][s][x:y].min()) for s in SAMPLES}
    if any(b not in 'ACGT' for b in profile['dna'][x:y]):
        return dict(status='reference_unknown',probability=None,minimum_coverage=coverage)
    vals=profile['bins'][x:y]
    if side=='left':vals=vals[::-1]
    model=profile['model'];prior=model['prior_vc']
    ll_ref=math.log(1-prior)+sum(math.log(model['MODEL_REF'][6*i+int(v)]) for i,v in enumerate(vals))
    ll_vc=math.log(prior)+sum(math.log(model['MODEL_VC'][6*i+int(v)]) for i,v in enumerate(vals))
    difference=ll_vc-ll_ref
    probability=(math.exp(-difference)/(1+math.exp(-difference)) if difference>=0
                 else 1/(1+math.exp(difference)))
    return dict(status='evaluable' if min(coverage.values())>=5 else 'coverage_insufficient',
                probability=probability,minimum_coverage=coverage)


def anchor_support(profile,boundary,side):
    fl,fr=profile['fetch_interval']
    interval=(boundary-250,boundary) if side=='left' else (boundary,boundary+250)
    if not fl<=interval[0]<interval[1]<=fr:
        return dict(status='outside_fetch',counts={s:0 for s in SAMPLES},members={s:[] for s in SAMPLES})
    pattern=profile['dna'][interval[0]-fl:interval[1]-fl]
    members={s:[] for s in SAMPLES};evidence=[];incomplete={s:[] for s in SAMPLES}
    for record in profile['records']:
        search=flank_search(record,interval,side)
        sequence=record['query_sequence'].upper()
        endpoints=all_exact(pattern,sequence,search,side) if search['status']=='ready' else []
        matrices=None
        if endpoints:
            scored=dict(search,endpoints=endpoints,match_threshold=175)
        else:
            scored,matrices=score_flank(pattern,sequence,search)
            endpoints=scored['endpoints']
        if scored.get('status')=='search_incomplete':
            incomplete[record['sample']].append(record['molecule_key'])
        if endpoints:
            representative=min(endpoints,key=lambda e:(e[1],-e[3],abs(e[0]-search['original_hint']),e[0]))
            members[record['sample']].append(record['molecule_key'])
            evidence.append(dict(molecule=record['molecule_key'],record_id=record['record_id'],sample=record['sample'],
                endpoint=representative,trace=traceback(pattern,sequence,scored,representative,matrices)))
    counts={s:len(set(members[s])) for s in SAMPLES}
    # A lower bound of five completed successes is enough. Otherwise unfinished
    # searches that could cross the threshold prevent declaring this point failed.
    status='pass' if min(counts.values())>=5 else 'search_incomplete' if any(counts[s]<5<=counts[s]+len(incomplete[s]) for s in SAMPLES) else 'anchor_support_insufficient'
    return dict(status=status,counts=counts,members=members,evidence=evidence,incomplete=incomplete)


def base_model(locus,fetch_interval,dna,anchor_length=250):
    return dict(schema_version=1,extraction_version=VERSION,locus=locus.name,chrom=locus.chrom,
                core_interval=[locus.start,locus.end],fetch_interval=list(fetch_interval),
                reference_sha256=digest(dna),anchor_length=anchor_length,analysis_interval=None,
                scoring=SCORING,normalization_policy=NORMALIZATION,search_padding=120,max_dp_cells=2000000,
                purpose='ablation' if anchor_length==100 else 'core_only',diagnostic_anchor_length=100,
                coordinate_system='0-based half-open BAM-stored query / reference forward',
                discovery_samples=SAMPLES,discovery_min_mapq=50,minimum_position_molecules=5,
                minimum_anchor_molecules=5,minimum_match_fraction=.7,stable_probability=.5,
                molecule_key='sample,namespace,QNAME')


def discover(locus,profile,groups=None,support_lookup=None):
    groups=groups or WINDOWS
    fl,fr=profile['fetch_interval'];models={};audits=[];score_cache={};anchor_cache={}
    for name,widths in groups.items():
        model=base_model(locus,[fl,fr],profile['dna'])
        model.update(purpose='multiscale_context',window_group=name,windows=widths,
                     probability_model_sha256=profile['model_sha256'],depths=profile['depths'],
                     discovery_conflicts=profile['conflicts'],sides={})
        for side,origin,cap in [('left',locus.start,locus.start-fl-250),('right',locus.end,fr-locus.end-250)]:
            outcome=dict(status='no_candidate_within_cap',cap=cap,boundary=None)
            for distance in range(max(-1,cap)+1):
                boundary=origin-distance if side=='left' else origin+distance
                key=(side,boundary)
                if key not in score_cache:
                    score_cache[key]={str(w):assess_window(profile,boundary,side,w) for w in (10,25,50,100,150)}
                scores=score_cache[key]
                row=dict(group=name,side=side,distance=distance,boundary=boundary,scores=scores,
                         anchor_status='not_evaluated',anchor_counts={},status='window_failed')
                if all(scores[str(w)]['status']=='evaluable' and scores[str(w)]['probability']>=.5 for w in widths):
                    if key not in anchor_cache:
                        anchor_cache[key]=(support_lookup(boundary,side) if support_lookup else anchor_support(profile,boundary,side))
                    support=anchor_cache[key]
                    row.update(anchor_status=support['status'],anchor_counts=support['counts'],status=support['status'])
                    if support['status']=='pass':
                        outcome=dict(status='determined',cap=cap,boundary=boundary,distance=distance,
                                     support=support,window_scores={str(w):scores[str(w)] for w in widths})
                    elif support['status']=='search_incomplete':
                        outcome=dict(status='search_incomplete',cap=cap,boundary=None,distance=distance)
                audits.append(row)
                if outcome['status'] in ('determined','search_incomplete'):
                    break
            model['sides'][side]=outcome
        if all(model['sides'][s]['status']=='determined' for s in ('left','right')):
            model['analysis_interval']=[model['sides'][s]['boundary'] for s in ('left','right')]
            model['support_intersections']={s:sorted(set(model['sides']['left']['support']['members'][s]) & set(model['sides']['right']['support']['members'][s])) for s in SAMPLES}
        models[name]=model
    return models,audits


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--catalog',type=Path,required=True)
    parser.add_argument('--locus',required=True,choices=sorted(FETCH))
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    from vntrgraph import read_catalog,IndexedFastaReference
    from vntrreads import write_table
    loci=read_catalog(args.catalog)
    locus=next(l for l in loci if l.name==args.locus)
    reference=IndexedFastaReference(args.reference)
    fl,fr=FETCH[locus.name];dna=reference.fetch(locus.chrom,fl,fr).upper()
    profile=build_profile(load_snapshot(args.records),[fl,fr],dna)
    models,audits=discover(locus,profile)
    args.out.mkdir(parents=True,exist_ok=False)
    for group,model in models.items():
        model['input_records_sha256']=sha256(args.records.read_bytes()).hexdigest()
        (args.out/f'{group}.json').write_text(json.dumps(model,indent=2)+'\n')
    write_table(args.out/'boundary_candidates.tsv',audits,['group','side','distance','boundary','scores','anchor_status','anchor_counts','status'])
    print(json.dumps({g:m['analysis_interval'] for g,m in models.items()}))


if __name__=='__main__':main()
