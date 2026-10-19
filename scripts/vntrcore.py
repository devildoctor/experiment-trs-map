"""Sequence and evidence contracts for the VNTR point graph (no external packages).

All scores are descriptive costs, never calibrated probabilities. Candidate
generation and joint read assignment are separate stages. Graph traversal is
not an allele caller: only saved full paths are reportable.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace, asdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json


from vntrmotifs import primitive_unit, definition, forms, exact_runs, context_support, TwoPaths, base_metrics

VERSION = "2.1-phase1"
DNA = set("ACGT")


def revcomp(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def compact(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


@dataclass(frozen=True)
class Read:
    name: str
    sample: str
    platform: str
    sequence: str
    quality: str = ""
    namespace: str = ""
    query_sequence: str = ""
    qL: int = 0
    qR: int = 0
    anchor_scores: tuple = ()
    boundary_status: str = "candidate_window"
    records: tuple = ()
    platforms: tuple = ()
    read_group: str = ""
    library: str = ""
    hp: str = ""
    original_query_interval: tuple = ()

    @property
    def molecule_key(self):
        return compact([self.sample, self.namespace, self.name])

    @property
    def observed_platforms(self):
        return self.platforms or (self.platform,)


@dataclass
class Token:
    family: str
    exact_bases: str
    start: int
    end: int
    phase: int = 0
    orientation: str = "+"
    edit_operations: str = ""
    status: str = "full"
    annotation_ambiguity: bool = False
    level: str = "base"
    phase_policy: str = "fixed_template_phase"

    def __iter__(self):
        return iter((self.family, self.exact_bases))


@dataclass
class Allele:
    name: str
    sample: str
    sequence: str
    support: dict = field(default_factory=dict)
    read_ids: tuple = ()
    identity: float = 0.0
    tokens: list = field(default_factory=list)
    uncertain: bool = False
    dbg_k: int = 0
    graph_support: float = 0.0
    unitigs: int = 0
    seed_source: tuple = ()
    seed_molecules: tuple = ()
    call_status: str = "candidate"
    ambiguous_ids: tuple = ()
    consensus_residual: float = 0.0
    flags: list = field(default_factory=list)
    blocks: list = field(default_factory=list)
    platform_consistency: dict = field(default_factory=dict)
    annotation: dict = field(default_factory=dict)
    compound_details: list = field(default_factory=list)


def deduplicate_reads(reads):
    """Never concatenate SEQ; conflicting full observations remain unresolved.

Namespace is an explicit manifest value, otherwise BAM LB, otherwise empty.
Without provenance, matching names within a sample are conservatively treated
as the same molecule across files/platforms. Original record IDs are retained.
"""
    grouped = defaultdict(list)
    for read in reads:
        grouped[read.molecule_key].append(read)
    result, audit = [], []
    for key, group in sorted(grouped.items()):
        sequences = {r.sequence for r in group}
        record_ids = tuple(sorted({x for r in group for x in r.records}))
        platforms = tuple(sorted({p for r in group for p in r.observed_platforms}))
        if len(sequences) != 1:
            audit.append(dict(molecule_key=key, status="conflicting_sequences",
                              observations=len(group), records=record_ids, platforms=platforms,
                              sequence_hashes=sorted(sha256(s.encode()).hexdigest() for s in sequences)))
            continue
        best = min(group, key=lambda r: (r.platform != "HiFi", -sum(r.anchor_scores),
                                         -sum(ord(q) - 33 for q in r.quality), r.records))
        result.append(replace(best, records=record_ids, platforms=platforms))
        if len(group) > 1:
            audit.append(dict(molecule_key=key, status="deduplicated", observations=len(group),
                              records=record_ids, platforms=platforms))
    return result, audit


@lru_cache(maxsize=2048)
def edit_distance(a, b):
    """Exact end-to-end Levenshtein distance, Myers bit vectors, all lengths.

Python integers have arbitrary width. Equal lengths never imply distance 0.
The algorithm uses O(min(n,m)) bits and no length-only fallback.
"""
    if a == b:
        return 0
    if len(a) > len(b):
        a, b = b, a
    if not a:
        return len(b)
    peq = defaultdict(int)
    for i, char in enumerate(a):
        peq[char] |= 1 << i
    mask = (1 << len(a)) - 1
    top = 1 << (len(a) - 1)
    pv, mv, score = mask, 0, len(a)
    for char in b:
        eq = peq[char]
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | ~(xh | pv)
        mh = pv & xh
        score += bool(ph & top) - bool(mh & top)
        ph = (ph << 1) | 1
        mh <<= 1
        pv = (mh | ~(xv | ph)) & mask
        mv = ph & xv
    return score


def similarity(a, b, minimum=.8):
    return 1.0 - edit_distance(a, b) / max(1, len(a), len(b))


class AlignmentUnresolved(ValueError):
    pass


def align_to_seed(seed, read, band=8, max_cells=6_000_000):
    """Adaptive banded exact traceback, checked against exact distance.

Trace bytes are bounded. Resource exhaustion is explicit, never replaced by
length difference. M consumes both, D seed only, I read only.
"""
    if seed == read:
        return [("M", i, i) for i in range(len(seed))]
    n, m = len(seed), len(read)
    exact = edit_distance(seed, read)
    width = max(abs(n - m), min(band, 16), 1)
    while True:
        cells = sum(min(m, i + width) - max(0, i - width) + 1 for i in range(n + 1))
        if cells > max_cells:
            raise AlignmentUnresolved(f"traceback_resource_limit:{n}x{m}:band={width}")
        previous = {j: j for j in range(min(m, width) + 1)}
        trace = [(0, bytearray([2] * len(previous)))]
        for i in range(1, n + 1):
            lo, hi = max(0, i - width), min(m, i + width)
            current, row = {}, bytearray(hi - lo + 1)
            for j in range(lo, hi + 1):
                choices = [(previous.get(j - 1, n + m + 1) + (seed[i - 1] != read[j - 1]), 0)
                           if j else (n + m + 1, 0),
                           (previous.get(j, n + m + 1) + 1, 1),
                           (current.get(j - 1, n + m + 1) + 1, 2)]
                current[j], row[j - lo] = min(choices)
            trace.append((lo, row))
            previous = current
        if previous.get(m) == exact:
            break
        width = min(max(n, m), max(width * 2, exact))
    operations, i, j = [], n, m
    while i or j:
        lo, row = trace[i]
        op = row[j - lo]
        if op == 0:
            i, j = i - 1, j - 1
            operations.append(("M", i, j))
        elif op == 1:
            i -= 1
            operations.append(("D", i, j))
        else:
            j -= 1
            operations.append(("I", i, j))
    return operations[::-1]


def alignment_counts(seed, read):
    counts = dict(substitution=0, insertion=0, deletion=0, homopolymer_indel=0,
                  indel_ambiguity=False)
    operations = align_to_seed(seed, read)
    for op, i, j in operations:
        if op == "M":
            counts["substitution"] += seed[i] != read[j]
        else:
            seq, pos = (read, j) if op == "I" else (seed, i)
            homo = ((pos > 0 and seq[pos - 1] == seq[pos]) or
                    (pos + 1 < len(seq) and seq[pos + 1] == seq[pos]))
            counts["homopolymer_indel" if homo else "insertion" if op == "I" else "deletion"] += 1
            # Repeat indels can have equivalent placements; this flag is not a
            # calibrated probability or a complete enumeration of alignments.
            counts["indel_ambiguity"] |= homo
    counts["edit_distance"] = edit_distance(seed, read)
    return counts


def make_consensus(reads):
    """Deterministic real-read medoid plus one vote per independent molecule.

HiFi can provide the initial backbone when present. No reference correction,
random sampling, POA claim, or uncalibrated Q weighting is used.
"""
    reads, _ = deduplicate_reads(reads)
    if not reads:
        raise ValueError("consensus_without_independent_molecules")
    reliable = [r for r in reads if r.platform == "HiFi"] or reads
    choices = sorted({r.sequence for r in reliable})
    frequencies = Counter(r.sequence for r in reads)
    seed = min(choices, key=lambda s: (sum(edit_distance(s, t) * n for t, n in frequencies.items()), s))
    columns = [Counter() for _ in seed]
    insertions = [Counter() for _ in range(len(seed) + 1)]
    for read in reads:
        inserted = defaultdict(list)
        for op, i, j in align_to_seed(seed, read.sequence):
            if op == "M":
                columns[i][read.sequence[j]] += 1
            elif op == "D":
                columns[i]["-"] += 1
            else:
                inserted[i].append(read.sequence[j])
        for i in range(len(insertions)):
            insertions[i]["".join(inserted.get(i, []))] += 1
    result = []
    for i in range(len(seed) + 1):
        insertion, count = min(insertions[i].items(), key=lambda kv: (-kv[1], kv[0]))
        if count > len(reads) / 2:
            result.append(insertion)
        if i < len(seed):
            base, count = min(columns[i].items(), key=lambda kv: (-kv[1], kv[0]))
            if count <= len(reads) / 2:
                base = seed[i]
            if base != "-":
                result.append(base)
    return "".join(result)


def _candidate_groups(source, minimum):
    exact = defaultdict(list)
    for read in source:
        exact[read.sequence].append(read)
    clusters = [list(exact[s]) for s in sorted(exact) if len(exact[s]) >= minimum]
    stable = {g[0].sequence for g in clusters}
    residual = sorted((r for r in source if r.sequence not in stable),
                      key=lambda r: (r.platform != "HiFi", len(r.sequence), r.molecule_key))
    for read in residual:
        threshold = .96 if read.platform == "HiFi" else .82
        scores = sorted(((similarity(read.sequence, g[0].sequence), i)
                         for i, g in enumerate(clusters)), reverse=True)
        # Seed ties are preserved as separate residual groups, never alternating.
        if scores and scores[0][0] >= threshold and (len(scores) == 1 or scores[0][0] > scores[1][0]):
            clusters[scores[0][1]].append(read)
        else:
            clusters.append([read])
    return clusters


def split_stable_variants(group, minimum):
    """Split repeated residual substitution/indel patterns before consensus.

    A split needs independent support on both sides. Alignment resource failures
    retain the whole group and an explicit flag. No variant is corrected to the
    reference; ungrouped residual molecules remain available for joint scoring.
    """
    if len(group) < 2 * minimum or len({r.sequence for r in group}) == 1:
        return [(group, [])]
    seed = min({r.sequence for r in group},
               key=lambda s: (sum(edit_distance(s,r.sequence) for r in group), s))
    columns = defaultdict(lambda: defaultdict(list))
    try:
        for read in group:
            insertions = defaultdict(list)
            for op,i,j in align_to_seed(seed,read.sequence):
                if op == "M":
                    columns[(i,"base")][read.sequence[j]].append(read)
                elif op == "D":
                    columns[(i,"base")]["-"].append(read)
                else:
                    insertions[i].append(read.sequence[j])
            for i in range(len(seed)+1):
                columns[(i,"insertion")]["".join(insertions.get(i,[]))].append(read)
    except AlignmentUnresolved as error:
        return [(group,[str(error)])]
    choices = []
    for site,variants in columns.items():
        stable = {v:rs for v,rs in variants.items() if len(rs) >= minimum}
        if len(stable) >= 2:
            choices.append((min(len(rs) for rs in stable.values()), site, stable))
    if not choices:
        return [(group,[])]
    _,_,split = max(choices,key=lambda x:(x[0],x[1]))
    used, result = set(), []
    for variant,reads in sorted(split.items()):
        used.update(r.molecule_key for r in reads)
        for subgroup,flags in split_stable_variants(reads,minimum):
            result.append((subgroup,flags+["stable_residual_split"]))
    residual = [r for r in group if r.molecule_key not in used]
    if residual:
        result.append((residual,["residual_minor_pattern"]))
    return result


def generate_candidates(reads, minimum):
    by_sample, candidates = defaultdict(list), []
    for read in reads:
        by_sample[read.sample].append(read)
    for sample, source in sorted(by_sample.items()):
        groups = [item for group in _candidate_groups(source,minimum)
                  for item in split_stable_variants(group,minimum)]
        for group,split_flags in groups:
            if len(group) < minimum:
                continue
            flags = list(split_flags)
            try:
                sequence = make_consensus(group)
            except AlignmentUnresolved as error:
                sequence = min(group, key=lambda r: (r.platform != "HiFi", r.molecule_key)).sequence
                flags.append(str(error))
            if not sequence:
                continue
            same = next((a for a in candidates if a.sample == sample and a.sequence == sequence), None)
            ids = {r.molecule_key for r in group}
            platforms = {p for r in group for p in r.observed_platforms}
            if same:
                same.seed_molecules = tuple(sorted(set(same.seed_molecules) | ids))
                same.seed_source = tuple(sorted(set(same.seed_source) | platforms))
                same.flags.extend(flags)
            else:
                candidates.append(Allele("", sample, sequence, seed_source=tuple(sorted(platforms)),
                                         seed_molecules=tuple(sorted(ids)), flags=flags))
    for sample in sorted(by_sample):
        for i, allele in enumerate(sorted((a for a in candidates if a.sample == sample),
                                          key=lambda a: (len(a.sequence), a.sequence)), 1):
            allele.name = f"{sample}_C{i}"
    return candidates


def assign_reads(reads, candidates, min_margin=.5):
    assignments = []
    for read in reads:
        scores = {}
        threshold = .96 if read.platform == "HiFi" else .82
        for candidate in candidates:
            if candidate.sample != read.sample:
                continue
            distance = edit_distance(candidate.sequence, read.sequence)
            identity = 1 - distance / max(1, len(candidate.sequence), len(read.sequence))
            if identity < threshold:
                scores[candidate.name] = dict(status="poor_fit", identity=identity, edit_distance=distance)
                continue
            try:
                stats = alignment_counts(candidate.sequence, read.sequence)
                weight = 1.0 if read.platform == "HiFi" else .5
                cost = stats["substitution"] + stats["insertion"] + stats["deletion"] + weight * stats["homopolymer_indel"]
                scores[candidate.name] = dict(stats, status="scored", cost=cost, identity=identity)
            except AlignmentUnresolved as error:
                scores[candidate.name] = dict(status="unresolved", reason=str(error), identity=identity,
                                               edit_distance=distance)
        ranked = sorted((s["cost"], name) for name, s in scores.items() if s["status"] == "scored")
        unresolved = any(s["status"] == "unresolved" for s in scores.values())
        gap = ranked[1][0] - ranked[0][0] if len(ranked) > 1 else None
        near = [name for cost, name in ranked if cost - ranked[0][0] < min_margin] if ranked else []
        if unresolved:
            state, chosen = "unresolved", ""
        elif not ranked:
            state, chosen = "unexplained", ""
        elif len(near) > 1:
            state, chosen = "ambiguous", ""
        else:
            state, chosen = "assigned", ranked[0][1]
        assignments.append(dict(molecule_key=read.molecule_key, sample=read.sample,
                                platform=read.platform, platforms=read.observed_platforms,
                                status=state, candidate=chosen, best=ranked[0][1] if ranked else "",
                                score_gap=gap, near_candidates=near, scores=scores))
    return assignments


def call_candidates(reads, minimum=3, ploidy=2, min_margin=.5):
    reads, _ = deduplicate_reads(reads)
    candidates = generate_candidates(reads, minimum)
    assignments = assign_reads(reads, candidates, min_margin)
    # Repeated residual structures from every platform get another opportunity
    # to form candidates after consensus changes the original seed groups.
    residual_ids = {x["molecule_key"] for x in assignments if x["status"] == "unexplained"}
    extra = generate_candidates([r for r in reads if r.molecule_key in residual_ids], minimum)
    for allele in extra:
        if not any(a.sample == allele.sample and a.sequence == allele.sequence for a in candidates):
            allele.name = f"{allele.sample}_C{1 + sum(a.sample == allele.sample for a in candidates)}"
            allele.flags.append("residual_discovery")
            candidates.append(allele)
    if extra:
        assignments = assign_reads(reads, candidates, min_margin)
    lookup = {r.molecule_key: r for r in reads}
    for allele in candidates:
        assigned = [x for x in assignments if x["candidate"] == allele.name]
        allele.read_ids = tuple(sorted(x["molecule_key"] for x in assigned))
        allele.ambiguous_ids = tuple(sorted(x["molecule_key"] for x in assignments
                                            if x["status"] == "ambiguous" and allele.name in x["near_candidates"]))
        support = Counter(p for key in allele.read_ids for p in lookup[key].observed_platforms)
        allele.support = dict(support)
        identities = [x["scores"][allele.name]["identity"] for x in assigned]
        allele.identity = sum(identities) / max(1, len(identities))
        allele.consensus_residual = 1 - allele.identity
        for platform in sorted(support):
            platform_reads = [lookup[key] for key in allele.read_ids if platform in lookup[key].observed_platforms]
            record = dict(molecules=len(platform_reads),status="insufficient_support")
            if len(platform_reads) >= minimum:
                try:
                    platform_sequence = make_consensus(platform_reads)
                    distance = edit_distance(allele.sequence,platform_sequence)
                    record.update(status="concordant" if distance == 0 else "discordant",
                                  edit_distance=distance,sequence_sha256=sha256(platform_sequence.encode()).hexdigest())
                    if distance:
                        allele.flags.append("platform_consensus_conflict")
                except AlignmentUnresolved as error:
                    record.update(status="unresolved",reason=str(error))
                    allele.flags.append("platform_consensus_unresolved")
            allele.platform_consistency[platform] = record
        if any(x["scores"][allele.name].get("indel_ambiguity") for x in assigned):
            allele.flags.append("indel_placement_ambiguity")
        allele.call_status = "supported_sequence" if len(allele.read_ids) >= minimum else "candidate_only"
    genotypes = []
    for sample in sorted({r.sample for r in reads}):
        source = [x for x in assignments if x["sample"] == sample]
        supported = [a for a in candidates if a.sample == sample and len(a.read_ids) >= minimum]
        combinations_scored = []
        if len(supported) <= 12:
            for count in range(1, min(ploidy, len(supported)) + 1):
                for group in combinations(supported, count):
                    names = tuple(a.name for a in group)
                    cost = 0.0
                    for assignment in source:
                        read = lookup[assignment["molecule_key"]]
                        reject = max(1.0, len(read.sequence) * (.04 if read.platform == "HiFi" else .18))
                        values = [assignment["scores"][n]["cost"] for n in names
                                  if assignment["scores"].get(n, {}).get("status") == "scored"]
                        cost += min([reject] + values)
                    combinations_scored.append(dict(candidates=names, cost=cost))
        combinations_scored.sort(key=lambda x: (x["cost"], len(x["candidates"]), x["candidates"]))
        state = ("no_supported_sequence" if not supported else
                 "ploidy_conflict" if len(supported) > ploidy else
                 "one_observed_sequence_not_homozygous" if len(supported) == 1 else "supported_sequence_set_unphased")
        for allele in supported:
            allele.uncertain = len(supported) > ploidy or any(x not in {"stable_residual_split","residual_discovery"} for x in allele.flags)
            if len(supported) > ploidy:
                allele.call_status = "candidate_ploidy_conflict"
        genotypes.append(dict(sample=sample, ploidy=ploidy, status=state,
                              supported_candidates=[a.name for a in supported],
                              assignments=dict(Counter(x["status"] for x in source)),
                              combinations=combinations_scored,
                              combination_status="descriptive_not_probability" if len(supported) <= 12 else "too_many_candidates",
                              phase="locus_local_only"))
    return candidates, assignments, genotypes


def cluster_reads(reads, minimum):
    return [a for a in call_candidates(reads, minimum)[0] if len(a.read_ids) >= minimum]


def canonical_motif(sequence):
    return min(s[i:] + s[:i] for s in (sequence, revcomp(sequence)) for i in range(len(s)))


def motif_orientation(template, family):
    for orientation, seq in (("+", family), ("-", revcomp(family))):
        for phase in range(len(seq)):
            if seq[phase:] + seq[:phase] == template:
                return phase, orientation
    return 0, "+"


def discover_motifs(sequences, supplied=()):
    """Propose units and compound blocks; evidence validation is separate."""
    supplied = tuple(dict.fromkeys(tuple(supplied) + tuple(primitive_unit(m) for m in supplied)))
    scores = Counter({m: 1e6 for m in supplied})
    for sequence in sequences:
        for hit in exact_runs(sequence, max_period=200, min_period=2):
            scores[hit["sequence"]] += hit["count"]
        if len(sequence) < 6:
            continue
        windows = [sequence] if len(sequence) <= 600 else [sequence[i:i+600] for i in range(0, len(sequence)-90, 300)]
        for window in windows:
            periods = []
            for period in range(2, min(200, len(window)//3) + 1):
                agreement = sum(a == b for a,b in zip(window, window[period:])) / (len(window)-period)
                if agreement >= .70:
                    periods.append((period, agreement))
            if not periods:
                continue
            best = max(a for _, a in periods)
            selected = sorted(set(sorted(periods,key=lambda x:(-x[1],x[0]))[:8] +
                                  [(p,a) for p,a in periods if a >= best-.26][:8]))
            for period, agreement in selected:
                units = Counter(window[i:i+period] for i in range(0,len(window)-period+1,period))
                for unit,count in units.items():
                    if count < 2 or set(unit)-DNA:
                        continue
                    # Skip trivial harmonics; retain nontrivial composite units.
                    if any(len(unit)%d == 0 and unit == unit[:d]*(len(unit)//d) for d in range(1,len(unit))):
                        continue
                    scores[unit] += count * agreement
    selected = list(dict.fromkeys(supplied))
    families = {canonical_motif(m) for m in selected}
    for motif, _ in sorted(scores.items(), key=lambda x: (-x[1],len(x[0]),x[0])):
        family = canonical_motif(motif)
        if family not in families:
            selected.append(motif)
            families.add(family)
        if len(selected) >= max(32, len(supplied)):
            break
    return tuple(selected)


def motif_dictionary(reference, alleles, reads, supplied, minimum):
    proposed = list(discover_motifs([reference] + [a.sequence for a in alleles], supplied))
    # Nominate observed interruptions from consensus slices, then require raw-molecule evidence.
    known = {canonical_motif(m) for m in proposed}
    for allele in alleles:
        tokens, _ = decompose(allele.sequence, proposed)
        for i, token in enumerate(tokens):
            if token.status != 'full' or not 2 <= len(token.exact_bases) <= 200:
                continue
            if i == 0 or i+1 == len(tokens) or any(t.status != 'full' for t in (tokens[i-1], tokens[i+1])):
                continue
            unit = primitive_unit(token.exact_bases)
            if set(unit)-DNA or canonical_motif(unit) in known:
                continue
            proposed.append(unit); known.add(canonical_motif(unit))
    catalog_templates = set(supplied) | {primitive_unit(m) for m in supplied}
    evidence = [(r, context_support(r.sequence, proposed)) for r in reads]
    records = []
    for motif in proposed:
        locations = [dict(molecule_key=r.molecule_key, start=a, end=b)
                     for r, hits in evidence for a,b in hits[motif]['context_occurrences']]
        ids = sorted({p['molecule_key'] for p in locations})
        substring_ids = {r.molecule_key for r,hits in evidence if hits[motif]['occurrences']}
        validated = len(ids) >= minimum
        source = 'catalog' if motif in supplied else 'catalog_primitive' if motif in catalog_templates else 'discovered'
        record = dict(motif=motif, source=source, molecule_ids=ids, read_support=len(ids),
                      substring_read_support=len(substring_ids), support_positions=locations,
                      validated=validated, usable=validated or motif in catalog_templates,
                      evidence_status='supported' if validated else 'observed_below_threshold' if ids else 'not_observed_in_repeat_context')
        record.update(definition(motif)); records.append(record)
    usable = [r for r in records if r['usable']]
    for record in usable:
        unit = record['motif']
        shorter = {v for r in usable if len(r['motif']) < len(unit) for v in forms(r['motif'])}
        reachable = {0}
        for i in range(len(unit)):
            if i in reachable:
                for motif in shorter:
                    if unit.startswith(motif,i): reachable.add(i+len(motif))
        if len(unit) in reachable: record['level'] = 'compound'
    return records


def cigar_edit(template, exact):
    ops = []
    for op,i,j in align_to_seed(template,exact):
        mark = ("=" if template[i] == exact[j] else "X") if op == "M" else op
        if ops and ops[-1][1] == mark:
            ops[-1][0] += 1
        else:
            ops.append([1,mark])
    return "".join(f"{n}{op}" for n,op in ops)


def decompose(sequence, motifs, level="base", return_details=False):
    """Two best exact-slice paths under explicit fixed template phase.

    Partial bases never count as fully explained copies. Primitive powers are
    excluded from the base layer; the sequence is retained in every path.
    """
    n = len(sequence)
    templates = {primitive_unit(m) if level == 'base' else m for m in motifs}
    templates = tuple(sorted(templates | {revcomp(m) for m in templates}, key=lambda m:(len(m),m)))
    dp = TwoPaths(n)
    for i in range(n):
        weak_terminal = any(len(m)>1 and ((i==0 and sequence[i]==m[-1]) or (i==n-1 and sequence[i]==m[0])) for m in templates)
        dp.offer(i,i+1,.8 if weak_terminal else 2.,(i,i+1,'UNKNOWN',sequence[i:i+1],'unknown'))
        for motif in templates:
            length = len(motif)
            allowed = 0 if length == 1 else 1 if length <= 6 else min(3,max(1,int(length*.10)))
            for width in range(max(1 if length == 1 else 2,length-allowed),length+allowed+1):
                if i+width > n: continue
                exact = sequence[i:i+width]
                distance = edit_distance(motif,exact)
                if distance <= allowed:
                    terminal_partial = width < length and ((i+width == n and motif.startswith(exact)) or (i == 0 and motif.endswith(exact)))
                    if not terminal_partial:
                        dp.offer(i,i+width,.2+distance*1.5,(i,i+width,motif,exact,'full'))
            if i == 0:
                for width in range(2,min(length,n+1)):
                    if sequence[:width] == motif[-width:]:
                        dp.offer(i,width,.8,(i,width,motif,sequence[:width],'partial'))
            remaining = n-i
            if 2 <= remaining < length and sequence[i:] == motif[:remaining]:
                dp.offer(i,n,.8,(i,n,motif,sequence[i:],'partial'))
    audit = dp.audit()
    primary, alternate = dp.trace(), dp.trace(1)
    alternate_set = set(alternate)
    def convert(steps):
        tokens = []
        for start,end,motif,exact,status in steps:
            family = 'UNKNOWN' if status == 'unknown' else canonical_motif(motif)
            phase, orientation = (0,'+') if status == 'unknown' else motif_orientation(motif,family)
            edits = '' if status == 'unknown' else cigar_edit(motif,exact)
            tied = audit['annotation_ambiguity'] and (start,end,motif,exact,status) not in alternate_set
            if status == 'unknown' and tokens and tokens[-1].status == 'unknown':
                tokens[-1].exact_bases += exact; tokens[-1].end = end
                tokens[-1].annotation_ambiguity |= tied
            else:
                tokens.append(Token(family,exact,start,end,phase,orientation,edits,status,tied,level))
        assert ''.join(t.exact_bases for t in tokens) == sequence
        assert all(t.exact_bases == sequence[t.start:t.end] for t in tokens)
        return tokens
    tokens = convert(primary)
    audit['metrics'] = base_metrics(tokens)
    audit['alternative'] = [asdict(t) for t in convert(alternate)] if alternate else []
    audit['alternative_metrics'] = base_metrics(audit['alternative']) if alternate else None
    explained = audit['metrics']['full_unit_fraction'] if n else 0.
    return (tokens, explained, audit) if return_details else (tokens, explained)


def token_counts(tokens):
    return dict(full=dict(Counter(t.family for t in tokens if t.status == "full")),
                partial=dict(Counter(t.family for t in tokens if t.status == "partial")),
                unknown_bases=sum(len(t.exact_bases) for t in tokens if t.status == "unknown"),
                partial_bases=sum(len(t.exact_bases) for t in tokens if t.status == "partial"))


def build_graph(alleles, left, right):
    """Context nodes shared only by identical complete token paths.

Different positions get different nodes, even when bases match. This avoids
unobserved recombinant paths; display compression is handled separately.
"""
    nodes = {("flank","LEFT",left):"L",("flank","RIGHT",right):"R"}
    edges, paths = set(), []
    for allele in alleles:
        path_key = sha256(compact([(t.family,t.exact_bases,t.status) for t in allele.tokens]).encode()).hexdigest()
        ids = ["L"]
        for index,token in enumerate(allele.tokens):
            key = (f"{path_key}:{index}",token.family,token.exact_bases)
            if key not in nodes:
                nodes[key] = f"N{len(nodes)-1}"
            ids.append(nodes[key])
        ids.append("R")
        edges.update(zip(ids,ids[1:]))
        paths.append((allele,ids))
    sequence_by_node = {node:key[2] for key,node in nodes.items()}
    for allele,ids in paths:
        assert "".join(sequence_by_node[x] for x in ids) == left+allele.sequence+right
    return nodes,edges,paths


def graph_evidence(nodes, edges, paths):
    node_reads, edge_reads = defaultdict(set), defaultdict(set)
    for allele,ids in paths:
        for node in set(ids):
            node_reads[node].update(allele.read_ids)
        for edge in set(zip(ids,ids[1:])):
            edge_reads[edge].update(allele.read_ids)
    return dict(nodes=[dict(id=node,context=key[0],family=key[1],exact_bases=key[2],
                            molecule_ids=sorted(node_reads[node]),read_support=len(node_reads[node]))
                       for key,node in nodes.items()],
                edges=[dict(start=a,end=b,orientation="++",overlap="0M",molecule_ids=sorted(edge_reads[a,b]),
                            read_support=len(edge_reads[a,b])) for a,b in sorted(edges)],
                paths=[dict(id=a.name,sample=a.sample,sequence=a.sequence,
                            nodes=ids,orientation="+",molecule_ids=a.read_ids,read_support=len(a.read_ids),
                            platform_support=a.support,ambiguous_reads=len(a.ambiguous_ids),
                            occurrence_count=dict(Counter(ids)),motif_copy_count=token_counts(a.tokens),
                            count_scope="observed_sequence_annotation",
                            definitive_copy_count=None if a.uncertain else token_counts(a.tokens)["full"],
                            uncertain=a.uncertain,flags=a.flags,
                            call_status=a.call_status,flanks="reference_context_25bp") for a,ids in paths],
                contract="Report only saved full paths; support is unique molecules, not traversals.")


def dbg_evidence(allele, reads, destination=None):
    """Directed exact-string k-mers; evidence only, never prune/invent paths."""
    selected = [r for r in reads if r.molecule_key in allele.read_ids]
    k = min(31,max(1,len(allele.sequence)//2))
    edges = defaultdict(set)
    ordered = []
    for read in selected:
        path = [read.sequence[i:i+k] for i in range(max(0,len(read.sequence)-k+1))]
        ordered.append(dict(molecule_key=read.molecule_key,platforms=read.observed_platforms,path=path))
        for edge in set(zip(path,path[1:])):
            edges[edge].add(read.molecule_key)
    sequence = allele.sequence
    walked = [sequence[i:i+k] for i in range(max(0,len(sequence)-k+1))]
    path_edges = list(zip(walked,walked[1:]))
    supported = sum(bool(edges.get(e)) for e in path_edges)/max(1,len(path_edges))
    if not path_edges:
        supported = float(any(r.sequence == sequence for r in selected))
    incoming, outgoing = defaultdict(set),defaultdict(set)
    for a,b in edges:
        outgoing[a].add(b)
        incoming[b].add(a)
    branching = sum(len(incoming[n]) != 1 or len(outgoing[n]) != 1 for n in set(incoming)|set(outgoing))
    if destination is not None:
        destination.write_text(compact(dict(k=k,encoding="oriented_exact_string",overlap=k-1,
                                            edges=[dict(start=a,end=b,molecule_ids=sorted(ids)) for (a,b),ids in sorted(edges.items())],
                                            candidate_path=walked,read_paths=ordered)),encoding="utf-8")
    return k,supported,max(1,branching)
