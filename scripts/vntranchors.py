"""Reference flank grading and conservative, molecule-aware BAM extraction.

The supplied per-locus BAM is the entire input universe. Stored SEQ is never
stitched across records. All genomic coordinates are zero-based, half-open.
"""
from collections import Counter, defaultdict
from bisect import bisect_left
from dataclasses import dataclass
from math import ceil, floor
from pathlib import Path
import csv
import json
import re
import subprocess

from vntrcore import DNA, Read, compact, revcomp

LENGTHS = (150, 250, 350, 500)
IDENTITY = {'HiFi': .90, 'ONT': .80, 'ONTUL': .80}
OBSERVED_FRACTION = .90
LOCALIZED_FRACTION = .80
POSITION_TOLERANCE = 50
MAX_CANDIDATES = 100


def dependencies():
    try:
        import edlib
        import ahocorasick
    except ImportError as exc:
        raise RuntimeError('Anchor mode requires edlib==1.3.9.post1 and pyahocorasick==2.2.0') from exc
    return edlib, ahocorasick


def canonical(kmer):
    return min(kmer, revcomp(kmer))


def nonoverlapping(positions):
    chosen = []
    for pos in sorted(positions):
        if not chosen or pos >= chosen[-1] + 31:
            chosen.append(pos)
    return chosen


def load_equivalences(path, reference):
    """Only explicit, documented, base-identical coordinate blocks may collapse.

    An _alt suffix alone never proves equivalence. The table must be derived
    from an authoritative placement; the local sequence verifies every base.
    """
    if path is None:
        return []
    with Path(path).open(encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    for row in rows:
        for key in ('alt_start','alt_end','primary_start','primary_end'):
            row[key] = int(row[key])
        if not row.get('evidence_url') or row['strand'] not in ('+','-'):
            raise ValueError('Reference equivalence requires a placement evidence URL and strand')
        alt=reference.fetch(row['alt_contig'],row['alt_start'],row['alt_end'])
        primary=reference.fetch(row['primary_contig'],row['primary_start'],row['primary_end'])
        if not alt or (alt if row['strand']=='+' else revcomp(alt)) != primary:
            raise ValueError('Reference equivalence block is not base-identical: '+str(row))
    for i,row in enumerate(rows):
        for other in rows[i+1:]:
            if row['alt_contig']==other['alt_contig'] and max(row['alt_start'],other['alt_start']) < min(row['alt_end'],other['alt_end']):
                raise ValueError('Overlapping reference equivalence blocks')
    return rows


def normalized_position(chrom, pos, blocks, length=1):
    for block in blocks:
        if chrom==block['alt_contig'] and block['alt_start']<=pos and pos+length<=block['alt_end']:
            if block['strand']=='+':
                target=block['primary_start']+pos-block['alt_start']
            else:
                target=block['primary_end']-(pos-block['alt_start'])-length
            return block['primary_contig'],target,block['strand']=='-'
    return chrom,pos,False


def reference_profiles(path, loci, reference, lengths=LENGTHS, equivalences=()):
    """Count every canonical candidate 31-mer over every reference contig.

    Counts and contig counts are complete. First 100 occurrence positions per
    kmer are retained; truncation is explicit (only repeated kmers can truncate).
    A 30-base overlap preserves cross-chunk matches without double counting.
    """
    _, aho = dependencies()
    profiles, words = {}, set()
    for locus in loci:
        for side in ('L', 'R'):
            templates = []
            for length in lengths:
                start = locus.start-length if side == 'L' else locus.end
                end = locus.start if side == 'L' else locus.end+length
                if start < 0 or end > reference.contig_length(locus.chrom):
                    continue
                sequence = reference.fetch(locus.chrom, start, end)
                kmers = [(i, canonical(sequence[i:i+31])) for i in range(length-30)
                         if set(sequence[i:i+31]) <= DNA]
                words.update(word for _, word in kmers)
                templates.append(dict(length=length, start=start, end=end, sequence=sequence, kmers=kmers))
            profiles[locus.name+':'+side] = templates
    automaton = aho.Automaton()
    for word in words:
        automaton.add_word(word, word)
        automaton.add_word(revcomp(word), word)
    automaton.make_automaton()
    counts, contigs, positions = Counter(), defaultdict(Counter), defaultdict(list)
    with Path(path).open(encoding='ascii') as handle:
        chrom, parts, size, offset = None, [], 0, 0
        def scan(chunk):
            if not words:
                return
            for end, word in automaton.iter(chunk):
                counts[word] += 1
                contigs[word][chrom] += 1
                if len(positions[word]) < 100:
                    positions[word].append([chrom, offset+end-30])
        for line in handle:
            if line.startswith('>'):
                scan(''.join(parts))
                chrom, parts, size, offset = line[1:].split()[0], [], 0, 0
            else:
                part = line.strip().upper()
                parts.append(part)
                size += len(part)
                if size >= 4_000_000:
                    chunk = ''.join(parts)
                    scan(chunk)
                    offset += len(chunk)-30
                    parts, size = [chunk[-30:]], 30
        scan(''.join(parts))
    for templates in profiles.values():
        for template in templates:
            kmers=template.pop('kmers')
            raw_unique = [i for i,word in kmers if counts[word]==1]
            unique = [i for i,word in kmers if counts[word]<=100 and
                      len({normalized_position(c,p,equivalences,31)[:2] for c,p in positions[word]})==1]
            independent = nonoverlapping(unique)
            template.update(unique_starts=unique, unique_count=len(unique),
                            raw_unique_starts=raw_unique,raw_unique_count=len(raw_unique),
                            reference_equivalences=list(equivalences),
                            uniqueness_scope='explicit_locus_equivalence' if equivalences else 'all_reference_contigs',
                            nonoverlap_starts=independent, nonoverlap_count=len(independent),
                            grade='preferred' if len(independent)>=2 else 'limited' if unique else 'nonunique')
    probes = {word:dict(count=counts[word], contigs=dict(contigs[word]),
                       positions=positions[word], positions_truncated=counts[word]>100,
                       locus_count=len({normalized_position(c,p,equivalences,31)[:2] for c,p in positions[word]}) if counts[word]<=100 else None)
              for word in sorted(words)}
    return profiles, probes


def endpoint_costs(pattern, text):
    """Exact infix Levenshtein cost at every target endpoint, with Myers bits.

    No exact seed requirement. Free target prefix means the shifted Ph has a
    zero low bit (global alignment would insert a one). Tested against full DP.
    """
    size = len(pattern)
    if not size:
        return [0]*(len(text)+1)
    masks = defaultdict(int)
    for i, char in enumerate(pattern):
        masks[char] |= 1 << i
    mask, high = (1 << size)-1, 1 << (size-1)
    pv, mv, score = mask, 0, size
    result = [score]
    for char in text:
        eq = masks[char]
        xv = eq | mv
        xh = (((eq & pv)+pv) ^ pv) | eq
        ph = mv | ~(xh | pv)
        mh = pv & xh
        score += bool(ph & high) - bool(mh & high)
        ph = (ph << 1) & mask
        mh = (mh << 1) & mask
        pv = (mh | ~(xv | ph)) & mask
        mv = ph & xv
        result.append(score)
    return result


class AnchorOverflow(RuntimeError):
    pass


def flank_candidates(flank, sequence, minimum_identity):
    """Enumerate local endpoint minima, including approximate alternatives.

    Adjacent equal-cost endpoints are all kept. A limit produces unresolved,
    never a truncated candidate list that appears uniquely mapped.
    """
    edlib, _ = dependencies()
    length = len(flank)
    budget = floor(length*(1-minimum_identity)/minimum_identity + 1e-8)
    costs = endpoint_costs(flank, sequence)
    ends, i = [], 1
    while i < len(costs):
        j = i
        while j+1 < len(costs) and costs[j+1] == costs[i]:
            j += 1
        if costs[i] <= budget and costs[i] <= costs[i-1] and (j+1 == len(costs) or costs[j] < costs[j+1]):
            ends.extend(range(i, j+1))
        i = j+1
    if len(ends) > MAX_CANDIDATES:
        raise AnchorOverflow('too_many_anchor_endpoints')
    hits = []
    for end in ends:
        lo = max(0, end-length-budget)
        reverse = edlib.align(flank[::-1], sequence[lo:end][::-1], mode='SHW', task='locations', k=budget)
        if reverse['editDistance'] != costs[end]:
            raise RuntimeError('Infix endpoint/traceback disagreement')
        for _, reverse_end in reverse['locations']:
            start = end-reverse_end-1
            if start >= end:
                continue
            aligned = edlib.align(flank, sequence[start:end], mode='NW', task='path', k=budget)
            nice = edlib.getNiceAlignment(aligned, flank, sequence[start:end])
            a, b = nice['query_aligned'], nice['target_aligned']
            matches = sum(x == y and x in DNA for x, y in zip(a, b))
            observed = sum(x != '-' and y != '-' for x, y in zip(a, b))
            identity = matches/len(a)
            if identity+1e-12 < minimum_identity or observed/length < OBSERVED_FRACTION:
                continue
            mapping, ref_i, query_i = {}, 0, start
            for x, y in zip(a, b):
                if x != '-' and y != '-':
                    mapping[ref_i] = query_i
                ref_i += x != '-'
                query_i += y != '-'
            hits.append(dict(start=start, end=end, edits=aligned['editDistance'], identity=identity,
                             observed_fraction=observed/length, mapping=mapping,
                             inner_left_gap=(a[-1]=='-' or b[-1]=='-'),
                             inner_right_gap=(a[0]=='-' or b[0]=='-')))
            if len(hits) > MAX_CANDIDATES:
                raise AnchorOverflow('too_many_anchor_locations')
    return hits


@dataclass
class Record:
    name: str
    flag: int
    chrom: str
    pos: int
    mapq: int
    cigar: str
    sequence: str
    quality: str
    tags: dict
    namespace: str
    library: str
    source: str
    ops: list
    hard_left: int
    hard_right: int
    full_length: int
    platform: str = ''

    def original_base(self, q):
        q += self.hard_left
        return self.full_length-1-q if self.flag & 16 else q

    def original_interval(self, a, b):
        a, b = a+self.hard_left, b+self.hard_left
        return (self.full_length-b, self.full_length-a) if self.flag & 16 else (a, b)

    def original_sequence(self):
        # Only for comparison across records. Stored SEQ used for extraction is unchanged.
        return revcomp(self.sequence) if self.flag & 16 else self.sequence

    def original_start(self):
        return self.hard_right if self.flag & 16 else self.hard_left


def mapped_bases(record, start, end):
    ref, query = record.pos, 0
    for size, op in record.ops:
        if op in 'M=X':
            for pos in range(max(ref, start), min(ref+size, end)):
                yield pos, query+pos-ref
            ref += size
            query += size
        elif op in 'DN':
            ref += size
        elif op in 'IS':
            query += size


def parse_records(lines, sample, reference, locus, source, namespace, platform=''):
    records, audit, groups, header_length = [], Counter(), {}, None
    for index, line in enumerate(lines, 1):
        if line.startswith('@'):
            fields = line.rstrip().split('\t')
            tags = dict(x.split(':', 1) for x in fields[1:] if ':' in x)
            if fields[0] == '@SQ' and tags.get('SN') == locus.chrom:
                header_length = int(tags['LN'])
            if fields[0] == '@RG':
                groups[tags.get('ID', '')] = tags
            continue
        fields = line.rstrip('\n').split('\t')
        audit['input_records'] += 1
        if len(fields) < 11:
            audit['invalid_sam'] += 1
            continue
        name, flag, chrom, pos, mq, cigar = fields[:6]
        flag, pos, mq = int(flag), int(pos)-1, int(mq)
        tags = {p[0]:p[2] for item in fields[11:] if len(p:=item.split(':',2)) == 3}
        library = groups.get(tags.get('RG',''), {}).get('LB','')
        ns = namespace if namespace is not None else library
        sequence, quality = fields[9].upper(), fields[10]
        ops = [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]
        if cigar != '*' and ''.join(str(n)+op for n,op in ops) != cigar:
            audit['invalid_cigar'] += 1
            continue
        if sequence != '*' and cigar != '*' and sum(n for n,op in ops if op in 'MIS=X') != len(sequence):
            audit['invalid_cigar_query_length'] += 1
            continue
        if quality != '*' and (sequence == '*' or len(quality) != len(sequence)):
            audit['invalid_quality_length'] += 1
            continue
        hleft = ops[0][0] if ops and ops[0][1] == 'H' else 0
        hright = ops[-1][0] if ops and ops[-1][1] == 'H' else 0
        stored_length = len(sequence) if sequence != '*' else sum(n for n,op in ops if op in 'MIS=X')
        records.append(Record(name,flag,chrom,pos,mq,cigar,sequence,quality,tags,ns,library,
                              f'{source}:record{index}:flag{flag}:RG={tags.get("RG", "")}',
                              ops,hleft,hright,hleft+stored_length+hright,platform))
        for key, mask in [('unmapped',4),('secondary',256),('supplementary',2048),('qcfail',512),('duplicate',1024)]:
            if flag & mask:
                audit[key+'_retained'] += 1
        if mq < 20:
            audit['low_mapq_retained'] += 1
        if mq == 255:
            audit['mapq_unknown_retained'] += 1
    if header_length != reference.contig_length(locus.chrom):
        raise ValueError(f'BAM/FASTA contig length mismatch: {source}')
    return records, audit


def consistent_sequences(group):
    available = [r for r in group if r.sequence != '*']
    # Secondary representations can be trimmed inconsistently by upstream tools;
    # a contradictory identity/length is unresolved, not permission to stitch.
    if len({r.full_length for r in available}) > 1:
        return False
    for i, a in enumerate(available):
        for b in available[i+1:]:
            lo, hi = max(a.original_start(), b.original_start()), min(a.original_start()+len(a.sequence), b.original_start()+len(b.sequence))
            if lo < hi:
                sa = a.original_sequence()[lo-a.original_start():hi-a.original_start()]
                sb = b.original_sequence()[lo-b.original_start():hi-b.original_start()]
                if any(x != y and x in DNA and y in DNA for x,y in zip(sa,sb)):
                    return False
    return True


def localization(provider, hit, template, group, locus, min_mapq):
    matched, sources = set(), set()
    blocks=template.get('reference_equivalences',[])
    for record in group:
        if record.flag & (4|256|512|1024) or record.mapq < min_mapq or record.mapq == 255:
            continue
        if record.full_length != provider.full_length:
            continue
        intervals=[(template['start'],template['end'])] if record.chrom==locus.chrom else []
        intervals += [(b['alt_start'],b['alt_end']) for b in blocks
                      if b['alt_contig']==record.chrom and b['primary_contig']==locus.chrom
                      and max(b['primary_start'],template['start'])<min(b['primary_end'],template['end'])]
        for start,end in intervals:
            for refpos, q in mapped_bases(record,start,end):
                chrom,refpos,reverse=normalized_position(record.chrom,refpos,blocks)
                if chrom!=locus.chrom or (bool(record.flag&16)^reverse)!=bool(provider.flag&16):
                    continue
                offset = refpos-template['start']
                candidate = hit['mapping'].get(offset)
                if candidate is not None and abs(record.original_base(q)-provider.original_base(candidate)) <= POSITION_TOLERANCE:
                    matched.add(offset)
                    sources.add(record.source)
    return len(matched)/template['length'], sorted(sources)


def competing_mapping(provider, hit, template, group, locus):
    """Other mappings of the same query flank are ambiguity, even at low MAPQ.

    An alignment elsewhere on an unrelated portion of an ultralong read does
    not veto this locus. Compare original-query coordinates, not just SA tags.
    """
    expected = {provider.original_base(q): template['start']+i for i,q in hit['mapping'].items()}
    for record in group:
        if record.flag & (4|512|1024) or record.full_length != provider.full_length:
            continue
        bad = 0
        targets = {(record.full_length-1-raw if record.flag & 16 else raw)-record.hard_left:target
                   for raw,target in expected.items()}
        query_positions = sorted(targets)
        ref, query = record.pos, 0
        for size, op in record.ops:
            if op in 'M=X':
                begin=bisect_left(query_positions,query)
                end=bisect_left(query_positions,query+size,begin)
                for q in query_positions[begin:end]:
                    chrom,mapped,reverse=normalized_position(record.chrom,ref+q-query,template.get('reference_equivalences',[]))
                    if chrom != locus.chrom or (bool(record.flag&16)^reverse) != bool(provider.flag&16) or abs(mapped-targets[q]) > POSITION_TOLERANCE:
                        bad += 1
                ref += size
                query += size
            elif op in 'DN':
                ref += size
            elif op in 'IS':
                query += size
        if bad/template['length'] >= LOCALIZED_FRACTION:
            return record.source
    return ''


def choose_pair(left, right):
    pairs = [(a['edits']+b['edits'], a, b) for a in left for b in right if a['end'] <= b['start']]
    if not pairs:
        return None, 'missing_or_wrong_order'
    pairs.sort(key=lambda x:(x[0],x[1]['end'],x[2]['start']))
    best = pairs[0]
    qL, qR = best[1]['end'], best[2]['start']
    for score, a, b in pairs[1:]:
        if (a['end'], b['start']) == (qL, qR):
            continue
        if score == best[0]:
            return None, 'tied_cutpoints'
        different = abs(a['start']-best[1]['start']) > 10 or abs(b['start']-best[2]['start']) > 10
        gap = max(3, ceil(.01*(best[1]['length']+best[2]['length'])))
        if different and score-best[0] < gap:
            return None, 'competing_anchor_pair'
    return (best[1], best[2]), ''


def extract_records(lines, sample, platform, locus, reference, profiles, min_mapq,
                    source='', namespace=None):
    records, audit = parse_records(lines,sample,reference,locus,source,namespace,platform)
    return extract_parsed(records,audit,sample,locus,profiles,min_mapq)


def extract_parsed(records, audit, sample, locus, profiles, min_mapq):
    groups = defaultdict(list)
    for r in records:
        groups[(r.namespace,r.name)].append(r)
    reads, decisions = [], []
    for (ns,name), group in groups.items():
        row = dict(molecule_key=compact([sample,ns,name]),sample=sample,
                   platforms=sorted({r.platform for r in group}),
                   record_count=len(group),sources=[r.source for r in group],status='',observations=[])
        if not consistent_sequences(group):
            row['status'] = 'inconsistent_record_sequences_or_lengths'
            audit[row['status']] += 1
            decisions.append(row)
            continue
        observations = []
        for provider in group:
            platform = provider.platform
            info = dict(record=provider.source,status='',template_pairs=[])
            row['observations'].append(info)
            if provider.sequence == '*' or provider.flag & (512|1024):
                info['status'] = 'no_sequence_or_qc_excluded'
                continue
            if provider.flag & 4:
                info['status'] = 'no_target_orientation'
                continue
            if set(provider.sequence)-set('ACGTN'):
                info['status'] = 'unsupported_query_alphabet'
                continue
            hits = {'L':[], 'R':[]}
            try:
                for side in ('L','R'):
                    for template in profiles[locus.name+':'+side]:
                        found = flank_candidates(template['sequence'],provider.sequence,IDENTITY[platform])
                        for hit in found:
                            hit.update(length=template['length'],grade=template['grade'])
                            fraction, sources = localization(provider,hit,template,group,locus,min_mapq)
                            hit.update(localized_fraction=fraction,localization_sources=sources,
                                       competitor=competing_mapping(provider,hit,template,group,locus))
                            covered_unique = [i for i in template['unique_starts']
                                              if sum(j in hit['mapping'] for j in range(i,i+31)) >= 28]
                            hit['covered_unique_nonoverlap'] = len(nonoverlapping(covered_unique))
                            hits[side].append(hit)
            except AnchorOverflow as exc:
                info['status'] = str(exc)
                continue
            accepted_pairs, unresolved = [], []
            # Compare candidates within each template pair; unlike raw edit sums,
            # scores across different template lengths are not comparable.
            for ll in sorted({h['length'] for h in hits['L']}):
                for rr in sorted({h['length'] for h in hits['R']}):
                    pair, reason = choose_pair([h for h in hits['L'] if h['length']==ll],
                                               [h for h in hits['R'] if h['length']==rr])
                    if pair is None:
                        unresolved.append(reason)
                        continue
                    a,b = pair
                    detail = dict(left_length=ll,right_length=rr,qL=a['end'],qR=b['start'],
                                  identity=[a['identity'],b['identity']],grades=[a['grade'],b['grade']],
                                  localized_fraction=[a['localized_fraction'],b['localized_fraction']],
                                  covered_unique_nonoverlap=[a['covered_unique_nonoverlap'],b['covered_unique_nonoverlap']],
                                  localization_sources=sorted(set(a['localization_sources']+b['localization_sources'])),
                                  competitor=a['competitor'] or b['competitor'])
                    info['template_pairs'].append(detail)
                    if detail['competitor']:
                        detail['status'] = 'competing_genomic_mapping'
                    elif a['inner_left_gap'] or b['inner_right_gap']:
                        detail['status'] = 'anchor_inner_boundary_gap'
                    elif min(a['localized_fraction'],b['localized_fraction']) < LOCALIZED_FRACTION:
                        detail['status'] = 'sequence_only_unverified'
                    elif min(a['covered_unique_nonoverlap'],b['covered_unique_nonoverlap']) < 2:
                        detail['status'] = 'limited_reference_discrimination'
                    else:
                        detail['status'] = 'eligible_pair'
                        accepted_pairs.append((a,b,detail))
            info['unresolved_pairs'] = sorted(set(unresolved))
            if not accepted_pairs:
                info['status'] = next((d['status'] for d in info['template_pairs']),
                                      unresolved[0] if unresolved else 'missing_or_ambiguous_anchor')
                continue
            cuts = {(a['end'],b['start']) for a,b,_ in accepted_pairs}
            # All localized, sufficiently discriminating template lengths must
            # agree; a longer flank never silently overrides a shorter one.
            if len(cuts) != 1:
                info['status'] = 'template_boundary_conflict'
                continue
            a,b,detail = min(accepted_pairs,key=lambda x:(x[0]['length']+x[1]['length'],x[0]['length']))
            qL,qR = a['end'],b['start']
            observed = provider.sequence[qL:qR]
            info.update(qL=qL,qR=qR,sequence=observed,original_interval=provider.original_interval(qL,qR))
            if not observed:
                info['status'] = 'zero_length_candidate'
                continue
            if set(observed)-DNA:
                info['status'] = 'ambiguous_repeat_bases'
                continue
            info['status'] = 'complete_support'
            read = Read(name,sample,platform,observed,
                        provider.quality[qL:qR] if provider.quality!='*' else '',ns,
                        provider.sequence,qL,qR,(a['identity'],b['identity']),locus.boundary_status,
                        tuple(sorted({provider.source,*detail['localization_sources']})),(platform,),
                        provider.tags.get('RG',''),provider.library,provider.tags.get('HP',''),provider.original_interval(qL,qR))
            observations.append(read)
        statuses = [i['status'] for i in row['observations']]
        conflicts = {'template_boundary_conflict','tied_cutpoints','competing_genomic_mapping','competing_anchor_pair'}
        if any(s in conflicts for s in statuses):
            row['status'] = 'unresolved_molecule_conflict'
        elif observations:
            if len({(r.sequence,r.original_query_interval) for r in observations}) != 1:
                row['status'] = 'conflicting_complete_observations'
            else:
                row['status'] = 'complete_support'
                reads.extend(observations)
                audit['spanning_alignment'] += len(observations)
        else:
            row['status'] = statuses[0] if statuses else 'no_usable_record'
        audit['molecule_'+row['status']] += 1
        for info in row['observations']:
            audit[info['status']] += 1
        decisions.append(row)
    audit['spanning_molecules_before_global_dedup'] = len({r.molecule_key for r in reads})
    return reads, audit, decisions


def extract_bam(bam,sample,platform,locus,reference,samtools,profiles,min_mapq,namespace=None):
    # Entire supplied BAM, no region argument or MAPQ/FLAG prefilter.
    command = [str(samtools),'view','-h',str(bam)]
    result = subprocess.run(command,capture_output=True,text=True)
    if result.returncode:
        raise RuntimeError(f'samtools view failed: {bam}: {result.stderr.strip()}')
    return extract_records(result.stdout.splitlines(),sample,platform,locus,reference,profiles,
                           min_mapq,str(bam),namespace)


def extract_bams(bams,locus,reference,samtools,profiles,min_mapq,namespaces):
    """Combine records across supplied files before molecule evidence decisions."""
    by_sample, per_bam, failures = defaultdict(list), [], []
    for sample,platform,bam in bams:
        result = subprocess.run([str(samtools),'view','-h',str(bam)],capture_output=True,text=True)
        if result.returncode:
            failures.append(f'samtools view failed: {bam}: {result.stderr.strip()}')
            continue
        try:
            records,audit = parse_records(result.stdout.splitlines(),sample,reference,locus,str(bam),
                                          namespaces.get(str(bam.resolve())),platform)
        except ValueError as exc:
            failures.append(str(exc))
            continue
        by_sample[sample].extend(records)
        per_bam.append(dict(bam=str(bam),sample=sample,platform=platform,reasons=dict(audit)))
    reads,details = [],[]
    for sample,records in by_sample.items():
        found,_,decisions=extract_parsed(records,Counter(),sample,locus,profiles,min_mapq)
        reads.extend(found)
        details.extend(decisions)
    for row in per_bam:
        prefix=row['bam']+':record'
        row['spanning_observations']=sum(any(s.startswith(prefix) for s in r.records) for r in reads)
        row['molecule_decisions']=dict(Counter(d['status'] for d in details if any(s.startswith(prefix) for s in d['sources'])))
    return reads,per_bam,details,failures
