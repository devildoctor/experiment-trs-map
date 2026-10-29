"""Shared motif evidence and exact-slice annotation rules.

Exact primitive periods describe sequence, not a unique biological repeat unit.
Catalog units and observed compound patterns remain separately reportable.
"""
from functools import lru_cache
from dataclasses import dataclass


# 求字符串的最短精确周期；这是序列性质，不自动等于生物学重复单元定义。
@lru_cache(maxsize=4096)
def primitive_unit(sequence):
    for width in range(1, len(sequence) + 1):
        if len(sequence) % width == 0 and sequence == sequence[:width] * (len(sequence) // width):
            return sequence[:width]
    return sequence


# 枚举循环旋转和反向互补用于家族比较；分解时仍需遵守具体模板相位。
@lru_cache(maxsize=4096)
def forms(sequence):
    reverse = sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return tuple(sorted({s[i:] + s[:i] for s in (sequence, reverse) for i in range(len(s))}))


def canonical(sequence):
    return forms(sequence)[0]


def definition(sequence):
    primitive = primitive_unit(sequence)
    return dict(primitive_sequence=primitive, repeat_factor=len(sequence) // len(primitive),
                level='compound' if primitive != sequence else 'base',
                family=canonical(sequence), phase_policy='fixed_template_phase')


# 从局部精确串联重复中提出模板；重叠命中按跨度优先保留，未穷举所有可能分解。
def exact_runs(sequence, max_period=100, min_period=1, min_copies=3, min_span=12):
    """Local primitive exact repeats; coordinates are 0-based, half-open."""
    hits = []
    for width in range(min_period, min(max_period, len(sequence) // min_copies) + 1):
        i = 0
        while i + width * min_copies <= len(sequence):
            unit = sequence[i:i + width]
            if set(unit) - set('ACGT') or primitive_unit(unit) != unit:
                i += 1
                continue
            end = i + width
            while sequence[end:end + width] == unit:
                end += width
            if end - i >= max(min_span, min_copies * width):
                hits.append(dict(start=i, end=end, sequence=unit, count=(end-i)//width, period=width))
                i = end
            else:
                i += 1
    selected = []
    for hit in sorted(hits, key=lambda h: (-(h['end']-h['start']), h['period'], h['start'])):
        if not any(hit['start'] < h['end'] and h['start'] < hit['end'] for h in selected):
            selected.append(hit)
    return sorted(selected, key=lambda h: h['start'])


def exact_occurrences(sequence, motif):
    hits = set()
    for word in forms(motif):
        start = sequence.find(word)
        while start >= 0:
            hits.add((start, start + len(word)))
            start = sequence.find(word, start + 1)
    return sorted(hits)


# 单个子串命中不够：须相邻同家族重复，或被两个直接相邻的单元夹住。
def context_support(sequence, motifs):
    """Exact observations count only in tandem or between two adjacent units.

    A variant can occur once per molecule. A lone substring is not tandem evidence.
    Matching a long candidate never supplies context to its own subinterval.
    """
    motifs = tuple(dict.fromkeys(motifs))
    hits = {m: exact_occurrences(sequence, m) for m in motifs}
    starts, ends = {}, {}
    for motif, positions in hits.items():
        for start, end in positions:
            starts.setdefault(start, []).append((end, motif))
            ends.setdefault(end, []).append((start, motif))
    supported = {}
    for motif, positions in hits.items():
        family = canonical(motif)
        accepted = []
        for start, end in positions:
            left = ends.get(start, [])
            right = starts.get(end, [])
            tandem = any(canonical(m) == family for _, m in left + right)
            if tandem or (left and right):
                accepted.append((start, end))
        supported[motif] = dict(occurrences=positions, context_occurrences=accepted)
    return supported


@dataclass(eq=False)
class PathState:
    score: float
    previous: object
    step: tuple


# 切片有向无环图中每个终点只保留两个不同前缀，用于最佳分解及一个备选方案。
class TwoPaths:
    """Two distinct paths in the slice DAG; scores are costs, not probabilities."""
    def __init__(self, length):
        self.paths = [[] for _ in range(length + 1)]
        self.paths[0] = [PathState(0., None, ())]
        self.cost = [float('inf')] * (length + 1)
        self.cost[0] = 0.

    def offer(self, start, end, edge_cost, step):
        for prefix in self.paths[start]:
            choices = self.paths[end]
            if any(p.previous is prefix and p.step == step for p in choices):
                continue
            candidate = PathState(prefix.score + edge_cost, prefix, step)
            choices.append(candidate)
            choices.sort(key=lambda p: p.score)
            del choices[2:]
            self.cost[end] = choices[0].score

    def trace(self, rank=0):
        if rank >= len(self.paths[-1]):
            return []
        node = self.paths[-1][rank]
        steps = []
        while node.previous is not None:
            steps.append(node.step)
            node = node.previous
        return list(reversed(steps))

    # 分数差为描述性代价差；发现并列即标歧义，但不会枚举全部并列方案。
    def audit(self):
        best = self.paths[-1][0].score
        second = self.paths[-1][1].score if len(self.paths[-1]) > 1 else None
        return dict(primary_score=best, alternative_score=second,
                    score_gap=None if second is None else max(0., second-best),
                    annotation_ambiguity=second is not None and abs(second-best) < 1e-8,
                    score_definition='descriptive_cost_not_probability',
                    phase_policy='fixed_template_phase', alternative_limit=1)


# 完整、编辑匹配、部分及未知碱基分别统计；full_unit_fraction 不代表序列正确率。
def base_metrics(tokens):
    result = dict(full_units=0, exact_full_bases=0, edited_full_bases=0,
                  partial_bases=0, unknown_bases=0)
    for token in tokens:
        t = token if isinstance(token, dict) else vars(token)
        width = len(t['exact_bases'])
        if t['status'] in {'full', 'unit'}:
            result['full_units'] += 1
            edited = t.get('edit_distance', 0) > 0 or any(c in t.get('edit_operations', '') for c in 'XID')
            result['edited_full_bases' if edited else 'exact_full_bases'] += width
        elif t['status'] == 'partial':
            result['partial_bases'] += width
        else:
            result['unknown_bases'] += width
    total = sum(result[k] for k in ('exact_full_bases', 'edited_full_bases', 'partial_bases', 'unknown_bases'))
    result['full_unit_fraction'] = (result['exact_full_bases'] + result['edited_full_bases']) / max(1, total)
    return result


# 目录复合单元和实际观察到的复合排列分别报告；不同目录模式可能重叠，不能相加。
def compound_annotations(tokens, catalog=()):
    """Catalog blocks must align to complete base tokens; each pattern has its own count.

    Different catalog patterns may overlap, so their counts must never be summed.
    Unlisted exact harmonics are not enumerated as additional biological units.
    """
    result = observed_compounds(tokens)
    rows = [t if isinstance(t, dict) else vars(t) for t in tokens]
    for motif in dict.fromkeys(catalog):
        matches = []
        i = 0
        while i < len(rows):
            seq, j = '', i
            while j < len(rows) and rows[j]['status'] in {'full', 'unit'} and len(seq) < len(motif):
                seq += rows[j]['exact_bases']; j += 1
            if j-i >= 2 and seq in forms(motif):
                matches.append((i,j,seq)); i = j
            else:
                i += 1
        for start,end,seq in matches:
            result.append(dict(start_token=start,end_token=end,start=rows[start]['start'],end=rows[end-1]['end'],
                               pattern=[t.get('motif',t.get('family')) for t in rows[start:end]],
                               count=1,units_per_block=end-start,base_unit_count=end-start,
                               exact_bases=seq,catalog_motif=motif,basis='catalog_block_on_base_tokens'))
    return result


# 在完整 token 序列上找重复的混合模式，遇到 UNKNOWN 或 partial 不跨越连接。
def observed_compounds(tokens, max_period=12):
    """Mixed patterns of whole observed base tokens; never cross UNKNOWN/partial."""
    rows = [t if isinstance(t, dict) else vars(t) for t in tokens]
    labels = [(t.get('motif', t.get('family')), t.get('orientation', '+'), t['exact_bases'])
              if t['status'] in {'full', 'unit'} else None for t in rows]
    result, i = [], 0
    while i < len(rows):
        best = None
        for width in range(2, min(max_period, (len(rows)-i)//2) + 1):
            pattern = labels[i:i+width]
            if None in pattern or len(set(pattern)) < 2:
                continue
            if any(width % p == 0 and pattern == pattern[:p] * (width//p) for p in range(1, width)):
                continue
            count = 1
            while labels[i+count*width:i+(count+1)*width] == pattern:
                count += 1
            if count >= 2 and (best is None or count*width > best['base_unit_count']):
                end = i + count*width
                best = dict(start_token=i, end_token=end, start=rows[i]['start'], end=rows[end-1]['end'],
                            pattern=[x[0] + ('r' if x[1] == '-' else '') for x in pattern],
                            count=count, units_per_block=width, base_unit_count=count*width,
                            exact_bases=''.join(t['exact_bases'] for t in rows[i:end]),
                            basis='observed_base_token_pattern')
        if best:
            result.append(best); i = best['end_token']
        else:
            i += 1
    return result
