"""Vamos-style effective motif selection for the small frozen locus dictionaries."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from vntrcore import edit_distance


def select_motifs(model,tokens,q=.1):
    """Enumerate the same objective/constraints as vamos effMotifILP.py.

    Counts are unit occurrences across registered P, including REF, not allele
    frequencies. Terminal-only templates are outside the compression problem.
    """
    if not 0<=q<=1:raise ValueError('q must be in [0,1]')
    sequences={m['id']:m['sequence'] for m in model['motifs']}
    counts=Counter(t['motif'] for ts in tokens.values() for t in ts
                   if t['status']=='unit' and t['motif'] in sequences)
    ids=sorted(sequences,key=lambda m:(-counts[m],m))
    distances={a:{b:edit_distance(sequences[a],sequences[b]) for b in ids} for a in ids}
    # The upstream pipeline sorts motifs by descending occurrence count. Its
    # percentile array contains each pair distance twice, weighted by count[i].
    histogram=Counter()
    for i,a in enumerate(ids):
        for b in ids[i+1:]:histogram[distances[a][b]]+=2*counts[a]
    total=sum(histogram.values());quantile=Fraction(0)
    if total:
        rank=Fraction(str(q))*(total-1);lo=rank.numerator//rank.denominator;hi=lo+(rank.denominator!=1)
        def at(index):
            offset=0
            for distance,weight in sorted(histogram.items()):
                offset+=weight
                if index<offset:return distance
        quantile=at(lo)+(rank-lo)*(at(hi)-at(lo))
    removed=int(sum(counts.values())*Fraction(str(q)));delta=int(quantile*removed)
    best=None;examined=0
    for size in range(1,len(ids)+1):
        for subset in combinations(ids,size):
            examined+=1;mapping={}
            for a in ids:
                choices=[b for b in subset if counts[b]>=counts[a]]
                if not choices:break
                mapping[a]=min(choices,key=lambda b:(distances[a][b],-counts[b],b))
            if len(mapping)!=len(ids):continue
            cost=sum(counts[a]*distances[a][b] for a,b in mapping.items())
            if cost>delta:continue
            kept=sorted(set(mapping.values()))
            rank=(delta*len(kept)+cost,len(kept),cost,kept)
            if best is None or rank<best[0]:best=(rank,kept,mapping,cost)
    _,kept,mapping,cost=best
    return dict(q=q,original_ids=sorted(ids),selected_ids=kept,mapping=mapping,
        counts={m:counts[m] for m in ids},distances=distances,delta=delta,
        pair_distance_quantile=float(quantile),weighted_template_cost=cost,
        objective=delta*len(kept)+cost,candidate_subsets=examined,
        algorithm='exact_subset_enumeration_vamos_objective',
        count_scope='Each registered P once, including reference; terminal-only templates excluded.',
        source='https://github.com/ChaissonLab/vamos/blob/3dd3b429922399e9f14e9bfebc0ea3056b18fd0a/buildTRs/lib/effMotifILP.py',
        cost_scope='Frequency-weighted template distance, distinct from actual per-P base edits.')
