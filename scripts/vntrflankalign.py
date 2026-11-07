"""Complete affine flank placement with explicit local ambiguity and provenance.

Coordinates index BAM-stored SEQ; a reverse flag never reverses that SEQ again.
The reference template is consumed completely, query ends are free. This is a
bounded local locator, not a genome-wide uniqueness or genotype assessment.
"""
from hashlib import sha256
import json
import math
from pathlib import Path
import re

import numpy as np
from numba import njit

VERSION = 'flank-v2.1'
SCORING = dict(mismatch=2, gap_open=5, gap_extend=1)
NORMALIZATION = 'cost,matches,cigar_distance,query_coordinates-v1'
FIELDS = ['extraction_version', 'boundary_model_sha256', 'extraction_status',
          'core_status', 'analysis_status', 'core_outer_relation', 'core_ambiguous',
          'analysis_ambiguous', 'threshold_path_dependent', 'normalization_reason']
OPS = re.compile(r'(\d+)([MIDNSHP=X])')
INF = 100000000


def digest(value):
    return sha256(value.encode('ascii')).hexdigest()


def canonical_hash(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                             ensure_ascii=True).encode()).hexdigest()


@njit(cache=True)
def affine_matrix(pattern, text):
    """Gotoh states M/D/I, including min/max matches among every optimal state.

    Opposite adjacent gaps are dominated by a mismatch with this fixed scoring.
    Final M/D states exclude terminal query insertions (free query suffix).
    """
    m, n = len(pattern), len(text)
    cost = np.full((3, m+1, n+1), INF, dtype=np.int32)
    low = np.full((3, m+1, n+1), INF, dtype=np.int32)
    high = np.full((3, m+1, n+1), -INF, dtype=np.int32)
    prev = np.full((3, m+1, n+1), 255, dtype=np.uint8)
    for j in range(n+1):
        cost[0, 0, j] = low[0, 0, j] = high[0, 0, j] = 0
    for i in range(1, m+1):
        cost[1, i, 0] = 5+i
        low[1, i, 0] = high[1, i, 0] = 0
        prev[1, i, 0] = 0 if i == 1 else 1
    for i in range(1, m+1):
        for j in range(1, n+1):
            base = pattern[i-1]
            match = int(base == text[j-1] and
                        (base == 65 or base == 67 or base == 71 or base == 84))
            for state in range(3):
                best, lo, hi, parent = INF, INF, -INF, 255
                for source in range(3):
                    if state == 0:
                        pi, pj = i-1, j-1
                        penalty, extra = (0 if match else 2), match
                    elif state == 1:
                        if source == 2:
                            continue
                        pi, pj = i-1, j
                        penalty, extra = (1 if source == 1 else 6), 0
                    else:
                        if source == 1:
                            continue
                        pi, pj = i, j-1
                        penalty, extra = (1 if source == 2 else 6), 0
                    value = cost[source, pi, pj] + penalty
                    a, b = low[source, pi, pj]+extra, high[source, pi, pj]+extra
                    if value < best:
                        best, lo, hi, parent = value, a, b, source
                    elif value == best:
                        lo = min(lo, a)
                        if b > hi:
                            hi, parent = b, source
                cost[state, i, j], low[state, i, j] = best, lo
                high[state, i, j], prev[state, i, j] = hi, parent
    return cost, low, high, prev


def encode(sequence):
    return np.frombuffer(sequence.encode('ascii'), dtype=np.uint8)


def cigar_segments(record):
    cigar, sequence = record['cigar'], record['query_sequence']
    ops = [(int(n), op) for n, op in OPS.findall(cigar)]
    if not ops or ''.join(f'{n}{op}' for n, op in ops) != cigar or any(n <= 0 for n, _ in ops):
        raise ValueError('invalid_cigar')
    if sum(n for n, op in ops if op in 'MIS=X') != len(sequence):
        raise ValueError('cigar_sequence_length_mismatch')
    r, q, result = record['pos'], 0, []
    for n, op in ops:
        rn, qn = r+(n if op in 'MDN=X' else 0), q+(n if op in 'MIS=X' else 0)
        result.append((r, rn, q, qn, op))
        r, q = rn, qn
    return result


def cigar_hint(segments, coordinate):
    # Fixed original hint; at an insertion the first reference-consuming
    # operation touching the boundary wins, matching the archived diagnostic.
    for r0, r1, q0, q1, op in segments:
        if op in 'M=X' and r0 <= coordinate <= r1:
            return q0+coordinate-r0
        if op == 'D' and r0 <= coordinate <= r1:
            return q0
        if op == 'N' and r0 < coordinate < r1:
            return None
    return None


def flank_search(record, template_interval, side, padding=120, shift=0):
    sequence = record['query_sequence'].upper()
    segments = cigar_segments(record)
    start, end = template_interval
    if any(op == 'N' and r0 < end and r1 > start for r0, r1, _, _, op in segments):
        return dict(status='reference_skip_in_template', side=side)
    a, b = cigar_hint(segments, start), cigar_hint(segments, end)
    if a is None or b is None:
        return dict(status='template_hint_unavailable', side=side)
    requested = [min(a, b)+shift-padding, max(a, b)+shift+padding]
    actual = [max(0, requested[0]), min(len(sequence), requested[1])]
    if actual[0] >= actual[1]:
        return dict(status='empty_search_interval', side=side)
    return dict(status='ready', side=side, reference_interval=list(template_interval),
                original_hint=b if side == 'left' else a,
                requested_search_interval=requested, search_interval=actual,
                physical_read_edge_clipped=actual != requested)


def all_exact(template, sequence, search, side):
    if any(b not in 'ACGT' for b in template):
        return []
    a, b = search['search_interval']
    endpoints, start = [], sequence.find(template, a, b)
    while start >= 0:
        endpoints.append([start+len(template) if side == 'left' else start,
                          0, len(template), len(template)])
        start = sequence.find(template, start+1, b)
    return endpoints


def score_flank(template, sequence, search, max_cells=2000000):
    if search['status'] != 'ready':
        return dict(search, complete=False, endpoints=[]), None
    if any(b not in 'ACGT' for b in template):
        return dict(search, status='reference_unknown', complete=False, endpoints=[]), None
    a, b = search['search_interval']
    if (len(template)+1)*(b-a+1) > max_cells:
        return dict(search, status='search_incomplete', reason='max_dp_cells',
                    max_dp_cells=max_cells, complete=False, endpoints=[]), None
    text = sequence[a:b]
    p, t = (template, text) if search['side'] == 'left' else (template[::-1], text[::-1])
    matrices = affine_matrix(encode(p), encode(t))
    cost, low, high, _ = matrices
    need, endpoints = math.ceil(.7*len(template)), []
    for j in range(1, len(t)+1):
        best = min(int(cost[0, -1, j]), int(cost[1, -1, j]))
        states = [s for s in (0, 1) if cost[s, -1, j] == best]
        lo = min(int(low[s, -1, j]) for s in states)
        hi = max(int(high[s, -1, j]) for s in states)
        if hi >= need:
            cut = a+j if search['side'] == 'left' else b-j
            endpoints.append([cut, best, lo, hi])
    endpoints.sort()
    return dict(search, status='scored', complete=True, endpoints=endpoints,
                endpoint_columns=['cut', 'cost', 'min_matches', 'max_matches'],
                match_threshold=need, anchor_length=len(template)), matrices


def pair_endpoints(left, right, hints):
    """All minimum-cost valid pairs, without a quadratic scan of all endpoints."""
    left, right = sorted(left), sorted(right)
    pos, active, left_cost, best, pairs = 0, [], INF, INF, []
    for r in right:
        while pos < len(left) and left[pos][0] < r[0]:
            l = left[pos]
            if l[1] < left_cost:
                left_cost, active = l[1], [l]
            elif l[1] == left_cost:
                active.append(l)
            pos += 1
        if not active:
            continue
        total = left_cost+r[1]
        if total < best:
            best, pairs = total, []
        if total == best:
            for l in active:
                pairs.append(dict(qstart=l[0], qend=r[0], cost=total,
                                  matches=l[3]+r[3], left=l, right=r,
                                  cigar_distance=abs(l[0]-hints[0])+abs(r[0]-hints[1])))
    pairs.sort(key=lambda p: (p['cost'], -p['matches'], p['cigar_distance'], p['qstart'], p['qend']))
    return pairs


def traceback(template, sequence, evidence, endpoint, matrices):
    a, b = evidence['search_interval']
    side, cut = evidence['side'], endpoint[0]
    if endpoint[1] == 0:
        start, end = (cut-len(template), cut) if side == 'left' else (cut, cut+len(template))
        result = dict(query_interval=[start, end], cigar=f'{len(template)}=',
                      cost=0, matches=len(template))
    else:
        cost, _, high, prev = matrices
        i, j = len(template), cut-a if side == 'left' else b-cut
        last = j
        state = min((0, 1), key=lambda s: (int(cost[s, i, j]), -int(high[s, i, j]), s))
        p, t = (template, sequence[a:b]) if side == 'left' else (template[::-1], sequence[a:b][::-1])
        ops = []
        while i:
            source = int(prev[state, i, j])
            if state == 0:
                ops.append('=' if p[i-1] == t[j-1] and p[i-1] in 'ACGT' else 'X')
                i, j = i-1, j-1
            elif state == 1:
                ops.append('D'); i -= 1
            else:
                ops.append('I'); j -= 1
            state = source
        ops.reverse()
        if side == 'right':
            ops.reverse()
        packed = []
        for op in ops:
            if packed and packed[-1][1] == op:
                packed[-1][0] += 1
            else:
                packed.append([1, op])
        result = dict(query_interval=[a+j, a+last] if side == 'left' else [b-last, b-j],
                      cigar=''.join(f'{n}{op}' for n, op in packed),
                      cost=endpoint[1], matches=endpoint[3])
    result['reference_interval'] = evidence['reference_interval']
    verify_alignment(template, sequence, result)
    return result


def verify_alignment(template, sequence, path):
    start, end = path['query_interval']
    query = sequence[start:end]
    i = j = score = matches = 0
    for n, op in OPS.findall(path['cigar']):
        n = int(n)
        if op in '=X':
            for _ in range(n):
                equal = template[i] == query[j] and template[i] in 'ACGT'
                if equal != (op == '='):
                    raise ValueError('flank_trace_base_mismatch')
                score += 0 if equal else 2
                matches += int(equal); i += 1; j += 1
        elif op == 'I':
            j += n; score += 5+n
        elif op == 'D':
            i += n; score += 5+n
        else:
            raise ValueError('invalid_flank_trace_op')
    if (i, j, score, matches) != (len(template), len(query), path['cost'], path['matches']):
        raise ValueError('flank_trace_not_reconstructable')


def locate_interval(record, interval, reference_fetch, anchor_length=250, padding=120,
                    shifts=(0, 0), max_cells=2000000, diagnostic=False):
    sequence = record['query_sequence'].upper()
    start, end = interval
    templates = [(start-anchor_length, start), (end, end+anchor_length)]
    evidence, patterns, exact = [], [], []
    for index, side in enumerate(('left', 'right')):
        side_padding = padding[index] if isinstance(padding, (tuple, list)) else padding
        search = flank_search(record, templates[index], side, side_padding, shifts[index])
        template = reference_fetch(*templates[index])
        patterns.append(template)
        if len(template) != anchor_length:
            search = dict(search, status='reference_template_unavailable')
        evidence.append(search)
        exact.append(all_exact(template, sequence, search, side) if search['status'] == 'ready' else [])
    hints = [x.get('original_hint', 0) for x in evidence]
    pairs = pair_endpoints(exact[0], exact[1], hints)
    matrices = [None, None]
    if pairs:
        evidence = [dict(e, status='scored', complete=True, endpoints=x,
                         endpoint_columns=['cut', 'cost', 'min_matches', 'max_matches'],
                         anchor_length=anchor_length, match_threshold=math.ceil(.7*anchor_length),
                         enumerated='all_zero_cost_endpoints') for e, x in zip(evidence, exact)]
    else:
        for i in range(2):
            evidence[i], matrices[i] = score_flank(patterns[i], sequence, evidence[i], max_cells)
        if all(e['complete'] for e in evidence):
            pairs = pair_endpoints(evidence[0]['endpoints'], evidence[1]['endpoints'], hints)
    complete = all(e.get('complete', False) for e in evidence)
    result = dict(reference_interval=list(interval), anchor_length=anchor_length,
                  status='success' if pairs else 'no_qualified_pair' if complete else
                         'search_incomplete' if any(e['status'] == 'search_incomplete' for e in evidence) else 'unavailable',
                  complete=complete, left=evidence[0], right=evidence[1], pairs=pairs,
                  ambiguous=len(pairs)>1, distinct_sequence_count=len({sequence[p['qstart']:p['qend']] for p in pairs}),
                  threshold_path_dependent=any(p[s][2] < math.ceil(.7*anchor_length) for p in pairs for s in ('left', 'right')),
                  normalization_policy=NORMALIZATION)
    if pairs:
        representative = pairs[0]
        result['representative'] = {k: representative[k] for k in ('qstart', 'qend', 'cost', 'matches', 'cigar_distance')}
        reason = 'only_minimum_cost_pair'
        remaining = pairs
        for key, reverse in [('matches', True), ('cigar_distance', False)]:
            if len(remaining) <= 1:
                break
            val = (max if reverse else min)(p[key] for p in remaining)
            new = [p for p in remaining if p[key] == val]
            if len(new) < len(remaining):
                reason = key
            remaining = new
        if len(remaining)>1:
            reason = 'query_coordinates'
        result['normalization_reason'] = reason
        for i, side in enumerate(('left', 'right')):
            used = {tuple(p[side]) for p in pairs}
            evidence[i]['optimal_traces'] = {str(e[0]): traceback(patterns[i], sequence, evidence[i], e, matrices[i]) for e in sorted(used)}
    return result


def core_eligible(fragment):
    evidence = fragment.get('core_anchor_evidence', {})
    representative = evidence.get('representative', {})
    return (fragment.get('extraction_version') == VERSION and fragment.get('core_status') == 'success'
            and evidence.get('complete') and evidence.get('status') == 'success'
            and (fragment.get('qstart'), fragment.get('qend')) == (representative.get('qstart'), representative.get('qend'))
            and fragment.get('core_outer_relation') != 'conflict'
            and bool(fragment.get('sequence'))
            and fragment['qend']-fragment['qstart'] == len(fragment['sequence']))


def validate_model(model, locus, reference):
    required = {'schema_version', 'extraction_version', 'locus', 'chrom', 'core_interval',
                'fetch_interval', 'reference_sha256', 'anchor_length', 'analysis_interval',
                'scoring', 'normalization_policy', 'search_padding', 'max_dp_cells'}
    if not required <= model.keys() or model['schema_version'] != 1 or model['extraction_version'] != VERSION:
        raise ValueError('invalid_flank_boundary_model_schema')
    if (model['locus'], model['chrom'], model['core_interval']) != (locus.name, locus.chrom, [locus.start, locus.end]):
        raise ValueError('flank_boundary_model_locus_mismatch')
    if model['scoring'] != SCORING or model['normalization_policy'] != NORMALIZATION:
        raise ValueError('unsupported_flank_scoring_or_normalization')
    if (type(model['search_padding']) is not int or model['search_padding'] < 0
            or type(model['max_dp_cells']) is not int or model['max_dp_cells'] < 1):
        raise ValueError('invalid_flank_search_limits')
    for key in ('fetch_interval', 'core_interval', 'analysis_interval'):
        interval = model[key]
        if key == 'analysis_interval' and interval is None:
            continue
        if (not isinstance(interval, list) or len(interval) != 2
                or any(type(v) is not int for v in interval) or not 0 <= interval[0] < interval[1]):
            raise ValueError('invalid_flank_reference_interval')
    length = model['anchor_length']
    if length != 250 and not (length == 100 and model.get('purpose') == 'ablation'):
        raise ValueError('unsupported_flank_anchor_length')
    fl, fr = model['fetch_interval']
    dna = reference.fetch(locus.chrom, fl, fr).upper()
    if len(dna) != fr-fl or digest(dna) != model['reference_sha256']:
        raise ValueError('flank_reference_hash_mismatch')
    for interval in (model['core_interval'], model['analysis_interval']):
        if interval is not None and not fl <= interval[0]-length < interval[0] < interval[1] < interval[1]+length <= fr:
            raise ValueError('flank_template_outside_fetch_interval')
    outer = model['analysis_interval']
    if outer and not outer[0] <= locus.start < locus.end <= outer[1]:
        raise ValueError('flank_analysis_interval_does_not_contain_core')
    return dna


class FlankContext:
    def __init__(self, model, locus, reference):
        self.model = json.loads(json.dumps(model))
        self.dna = validate_model(model, locus, reference)
        self.start, self.end = model['fetch_interval']
        self.model_hash = canonical_hash(model)

    def fetch(self, start, end):
        return self.dna[start-self.start:end-self.start] if self.start <= start <= end <= self.end else ''

    def locate(self, record, interval, length=None, **kwargs):
        return locate_interval(record, interval, self.fetch, length or self.model['anchor_length'],
                               padding=kwargs.pop('padding', self.model['search_padding']),
                               max_cells=self.model['max_dp_cells'], **kwargs)

    def extract(self, meta, sequence, query, reference_interval, boundary_status):
        raw = dict(meta, query_sequence=sequence)
        core = self.locate(raw, self.model['core_interval'])
        outer_interval = self.model['analysis_interval']
        outer = (dict(core) if outer_interval == self.model['core_interval'] else self.locate(raw, outer_interval)) if outer_interval else dict(status='not_frozen', pairs=[], ambiguous=False)
        relation, status = 'not_checked', 'display_only'
        if core['status'] == 'success':
            a, b = core['representative']['qstart'], core['representative']['qend']
            query, status = (a, b), 'core_only'
            if outer['status'] == 'success':
                oa, ob = outer['representative']['qstart'], outer['representative']['qend']
                relation = 'consistent' if oa <= a < b <= ob else 'conflict'
                status = 'both_consistent' if relation == 'consistent' else 'core_outer_conflict'
                outer['any_compatible_pair'] = any(o['qstart'] <= c['qstart'] < c['qend'] <= o['qend'] for o in outer['pairs'] for c in core['pairs'])
        elif outer['status'] == 'success':
            status = 'analysis_only'
        result = dict(meta, status='displayed' if query and query[0]<query[1] else 'no_observed_repeat_bases',
                      reference_interval=reference_interval, boundary_status=boundary_status,
                      extraction_version=VERSION, boundary_model_sha256=self.model_hash,
                      extraction_status=status, core_status=core['status'], analysis_status=outer['status'],
                      core_outer_relation=relation, core_ambiguous=core['ambiguous'],
                      analysis_ambiguous=outer['ambiguous'], threshold_path_dependent=core.get('threshold_path_dependent',False),
                      normalization_reason=core.get('normalization_reason'), core_anchor_evidence=core,
                      analysis_anchor_evidence=outer, anchor='both' if core['status']=='success' else 'cigar_only',
                      anchor_scores=[])
        if query and query[0]<query[1]:
            result.update(qstart=query[0], qend=query[1], sequence=sequence[query[0]:query[1]])
        if outer['status']=='success':
            oa, ob = outer['representative']['qstart'], outer['representative']['qend']
            result.update(analysis_query_interval=[oa, ob], analysis_sequence=sequence[oa:ob])
            if relation == 'consistent':
                result.update(left_context=sequence[oa:query[0]], right_context=sequence[query[1]:ob])
                assert result['left_context']+result['sequence']+result['right_context']==result['analysis_sequence']
        diagnostic = self.locate(raw, self.model['core_interval'], 100)
        result['diagnostic_100bp'] = diagnostic
        result['diagnostic_100bp_agrees'] = (diagnostic.get('representative',{}).get('qstart'),diagnostic.get('representative',{}).get('qend')) == (core.get('representative',{}).get('qstart'),core.get('representative',{}).get('qend')) if diagnostic['status']==core['status']=='success' else None
        return result


def load_context(path, locus, reference):
    return FlankContext(json.loads(Path(path).read_text(encoding='utf-8')), locus, reference)


def write_flank_outputs(destination, rows, paths=()):
    from html import escape
    from vntrreads import write_table
    membership={(o['row_id'],o['record_id']):p['id'] for p in paths for o in p['observations']}
    contexts, links, audit = {}, [], []
    for row in rows:
        for record in row['records']:
            identity=dict(row_id=row['id'],record_id=record['record_id'],molecule=row['molecule_key'],
                          sample=row['sample'],platform=row['platform'],read_id=row['read_id'])
            audit.append(dict(identity,**{k:v for k,v in record.items() if k not in identity and k!='analysis_sequence'}))
            sequence=record.get('analysis_sequence')
            if sequence is not None:
                key=digest(sequence);contexts[key]=sequence
                links.append(dict(identity,context_sha256=key,context_length=len(sequence),
                    core_path=membership.get((row['id'],record['record_id']),''),
                    analysis_query_interval=record['analysis_query_interval'],
                    **{k:record.get(k) for k in FIELDS}))
    (destination/'extraction_audit.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=True)+'\n' for r in audit),encoding='utf-8')
    (destination/'analysis_context.fa').write_text(''.join(f'>{key}\n{sequence}\n' for key,sequence in sorted(contexts.items())),encoding='ascii')
    columns=['row_id','record_id','molecule','sample','platform','read_id','context_sha256','context_length',
             'core_path','analysis_query_interval']+FIELDS
    write_table(destination/'context_members.tsv',links,columns)
    table=''.join('<tr>'+''.join('<td>'+escape(str(r[k]))+'</td>' for k in
         ['sample','platform','read_id','core_path','context_length','extraction_status','core_ambiguous','analysis_ambiguous'])+'</tr>' for r in links)
    html='''<!doctype html><meta charset="utf-8"><title>外围观测</title>
<style>body{font-family:Arial,"Microsoft YaHei",sans-serif;margin:28px;color:#263445}td,th{padding:7px;border-bottom:1px solid #ddd}p{line-height:1.7}</style>
<h1>外围原始观测与核心路径</h1><p>外围与核心分别定位；上下文差异不增加核心 P。规范代表仍可能具有切点歧义。</p>
<p><a href="analysis_context.fa">外围序列</a> · <a href="context_members.tsv">完整关联表</a> · <a href="extraction_audit.jsonl">全部定位证据</a> · <a href="index.html">核心结构图</a></p>
<table><tr><th>样本</th><th>平台</th><th>read</th><th>核心P</th><th>外围bp</th><th>状态</th><th>核心歧义</th><th>外围歧义</th></tr>'''+table+'</table>'
    (destination/'contexts.html').write_text(html,encoding='utf-8')
