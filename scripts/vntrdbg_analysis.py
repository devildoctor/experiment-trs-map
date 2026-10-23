#!/usr/bin/env python3
"""Conservative, evidence-preserving D1S80 branch and 16-bp unit analysis.

An exact substring in one observed record is a joint witness for its ordered
branch interval. It is not an allele call. Only a two-flank-anchored observation
identical to the entire candidate is a full-length exact witness.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
from pathlib import Path

from vntrcore import decompose, revcomp
from vntrdbg import raw_graph

VERSION = 'd1s80_joint_units_v1'
UNIT_BP = 16


def read_tsv(path):
    csv.field_size_limit(64 * 1024 * 1024)
    with Path(path).open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, fields):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', extrasaction='raise')
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False, sort_keys=True)
                             if isinstance(value, (list, dict)) else value
                             for key, value in row.items()})


def occurrences(haystack, needle):
    """All starts, including overlapping occurrences."""
    found = []
    start = 0
    while True:
        pos = haystack.find(needle, start)
        if pos < 0:
            return found
        found.append(pos)
        start = pos + 1


def branch_queries(records, k, minimum=3):
    """Build reproducible base-coordinate queries, independent of unitig IDs.

    Low-support observed edges do not create branch definitions. Reference and
    candidate edges remain in the topology as explicitly labeled hypotheses.
    Neighboring branch-edge positions form one block; only nonoverlapping blocks
    are paired, so one long bubble is not mistaken for two independent choices.
    """
    raw = raw_graph(records, k)
    sample_by_molecule = {key: value['sample'] for key, value in raw['molecules'].items()}
    hypothetical = set()
    for record in records:
        if record['kind'] in {'reference', 'candidate'}:
            path = raw['paths'].get(record['id'], [])
            hypothetical.update(zip(path, path[1:]))
    kept = set(hypothetical)
    for edge, keys in raw['edge_support'].items():
        if max(Counter(sample_by_molecule[key] for key in keys).values(), default=0) >= minimum:
            kept.add(edge)
    incoming, outgoing = defaultdict(set), defaultdict(set)
    for left, right in kept:
        outgoing[left].add(right)
        incoming[right].add(left)
    queries = []
    for record in records:
        if record['kind'] != 'candidate':
            continue
        sequence = record['sequence']
        path = raw['paths'].get(record['id'], [])
        positions = [i for i, (left, right) in enumerate(zip(path, path[1:]))
                     if len(outgoing[left]) > 1 or len(incoming[right]) > 1]
        blocks = []
        for pos in positions:
            if blocks and pos <= blocks[-1][1] + 1:
                blocks[-1] = (blocks[-1][0], pos)
            else:
                blocks.append((pos, pos))
        # A k+1 base signature spans [edge position, edge position+k+1).
        # Merge overlapping signatures before defining independent branches.
        merged = []
        for first, last in blocks:
            if merged and first < merged[-1][1] + k + 1:
                merged[-1] = (merged[-1][0], last)
            else:
                merged.append((first, last))
        candidate_id = record['id'].removeprefix('CAND_')
        for index, (first, last) in enumerate(merged, 1):
            end = last + k + 1
            signature = sequence[first:end]
            queries.append(dict(query_id=f'{candidate_id}.k{k}.B{index:03}',
                                candidate_id=candidate_id, sample=record['sample'],
                                kind='branch_block', start0=first, end0=end,
                                branch_edge_count=last-first+1, exact_bases=signature,
                                candidate_occurrences=len(occurrences(sequence, signature))))
        for index, (left, right) in enumerate(zip(merged, merged[1:]), 1):
            start, end = left[0], right[1] + k + 1
            signature = sequence[start:end]
            queries.append(dict(query_id=f'{candidate_id}.k{k}.J{index:03}',
                                candidate_id=candidate_id, sample=record['sample'],
                                kind='adjacent_branch_pair', start0=start, end0=end,
                                branch_edge_count=(left[1]-left[0]+1)+(right[1]-right[0]+1),
                                exact_bases=signature,
                                candidate_occurrences=len(occurrences(sequence, signature))))
        queries.append(dict(query_id=f'{candidate_id}.k{k}.FULL', candidate_id=candidate_id,
                            sample=record['sample'], kind='full_candidate_exact', start0=0,
                            end0=len(sequence), branch_edge_count=len(positions),
                            exact_bases=sequence, candidate_occurrences=1))
    return queries


def joint_support(records, queries):
    """Count each independent molecule once; retain every exact OBS occurrence.

    Nonmatching observations are uninformative: this version does not infer
    conflicts from sequencing errors, approximate alignment, or absent coverage.
    """
    observed = [record for record in records if record['kind'] == 'observed']
    by_sample = defaultdict(lambda: defaultdict(list))
    for record in observed:
        by_sample[record['sample']][record['molecule_key']].append(record)
    summaries, evidence = [], []
    for query in queries:
        molecules = by_sample[query['sample']]
        support = 0
        assigned = 0
        for key, paths in sorted(molecules.items()):
            hits = []
            for path in paths:
                if query['kind'] == 'full_candidate_exact':
                    positions = [0] if (path['coverage'] == 'two_flank_anchored'
                                        and path['sequence'] == query['exact_bases']) else []
                else:
                    positions = occurrences(path['sequence'], query['exact_bases'])
                hits.extend(dict(path_id=path['id'], start0=pos,
                                 end0=pos+len(query['exact_bases']), coverage=path['coverage'])
                            for pos in positions)
            if hits:
                support += 1
                assigned += key in query.get('assigned_molecule_keys', set())
            evidence.append(dict(query_id=query['query_id'], candidate_id=query['candidate_id'],
                                 sample=query['sample'], molecule_key=key,
                                 decision='exact_support' if hits else 'uninformative',
                                 occurrence_count=len(hits), occurrences_json=hits,
                                 assigned_to_candidate=key in query.get('assigned_molecule_keys', set())))
        summaries.append(dict(query_id=query['query_id'], candidate_id=query['candidate_id'],
                              sample=query['sample'], kind=query['kind'], k=query['k'],
                              start0=query['start0'], end0=query['end0'],
                              query_bp=len(query['exact_bases']),
                              query_sha256=sha256(query['exact_bases'].encode()).hexdigest(),
                              candidate_occurrences=query['candidate_occurrences'],
                              branch_edge_count=query['branch_edge_count'],
                              tested_molecules=len(molecules), exact_support_molecules=support,
                              assigned_exact_support_molecules=assigned,
                              uninformative_molecules=len(molecules)-support,
                              conflict_molecules='NA', support_fraction='NA',
                              interpretation='exact_ordered_interval_lower_bound'))
    return summaries, evidence


def template_for(token):
    if token['status'] == 'unknown':
        return ''
    family = token['family']
    oriented = family if token['orientation'] == '+' else revcomp(family)
    phase = token['phase']
    return oriented[phase:] + oriented[:phase]


def primitive_units(candidates, motifs, legacy_tokens):
    templates = sorted({row['motif'] for row in motifs
                        if row['usable'] == 'True' and row['level'] == 'base'
                        and len(row['motif']) == UNIT_BP})
    if not templates:
        raise ValueError('No usable 16-bp D1S80 templates')
    legacy = defaultdict(list)
    for row in legacy_tokens:
        legacy[row['allele']].append(row)
    units, summaries, counts = [], [], []
    for candidate in candidates:
        name = candidate['id'].removeprefix('CAND_')
        sequence = candidate['sequence']
        primary, _, audit = decompose(sequence, tuple(templates), return_details=True)
        schemes = [('primary', primary)]
        if audit['alternative']:
            schemes.append(('alternative', audit['alternative']))
        for scheme, tokens in schemes:
            tokens = [dict(vars(item)) if not isinstance(item, dict) else dict(item)
                      for item in tokens]
            if ''.join(token['exact_bases'] for token in tokens) != sequence:
                raise AssertionError(f'{name} {scheme}: units do not reconstruct sequence')
            position = 0
            for index, token in enumerate(tokens, 1):
                if token['start'] != position or token['end'] != position+len(token['exact_bases']):
                    raise AssertionError(f'{name} {scheme}: coordinates do not partition sequence')
                position = token['end']
                parents = [i for i, old in enumerate(legacy[name], 1)
                           if int(old['start']) < token['end'] and int(old['end']) > token['start']]
                units.append(dict(candidate_id=name, sample=candidate['sample'],
                                  decomposition=scheme, unit_index=index,
                                  start0=token['start'], end0=token['end'],
                                  status=token['status'], family=token['family'],
                                  template_bases=template_for(token),
                                  exact_bases=token['exact_bases'],
                                  orientation=token['orientation'], phase_offset=token['phase'],
                                  edit_operations=token['edit_operations'],
                                  annotation_ambiguity=token['annotation_ambiguity'],
                                  legacy_token_indexes_json=parents))
            if scheme == 'primary':
                full = sum(token['status'] == 'full' for token in tokens)
                partial = sum(token['status'] == 'partial' for token in tokens)
                unknown = sum(len(token['exact_bases']) for token in tokens
                              if token['status'] == 'unknown')
                alternate = audit['alternative']
                alt_count = sum(item['status'] == 'full' for item in alternate) if alternate else None
                count_status = ('unresolved' if unknown or partial else
                                'ambiguous' if audit['annotation_ambiguity'] and alt_count != full else
                                'resolved_count' if not audit['annotation_ambiguity'] or alt_count == full
                                else 'ambiguous')
                summaries.append(dict(candidate_id=name, sample=candidate['sample'],
                                      sequence_bp=len(sequence), nominal_period_bp=UNIT_BP,
                                      length_equiv_cn=len(sequence)/UNIT_BP,
                                      full_16bp_units=full, partial_units=partial,
                                      unknown_bp=unknown, primitive_cn=full if count_status == 'resolved_count' else 'NA',
                                      primitive_cn_status=count_status,
                                      annotation_ambiguity=audit['annotation_ambiguity'],
                                      alternative_full_units=alt_count if alt_count is not None else 'NA',
                                      legacy_full_tokens=sum(row['status'] == 'full' for row in legacy[name]),
                                      method_version=VERSION))
        grouped = defaultdict(lambda: [0, 0])
        for row in units:
            if row['candidate_id'] == name and row['decomposition'] == 'primary':
                key = (row['family'], row['status'])
                grouped[key][0] += 1
                grouped[key][1] += len(row['exact_bases'])
        for (family, status), (segments, bp) in sorted(grouped.items()):
            counts.append(dict(candidate_id=name, sample=candidate['sample'],
                               decomposition='primary', family=family, status=status,
                               segment_count=segments, observed_bp=bp))
    return units, summaries, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dbg-dir', type=Path, required=True)
    parser.add_argument('--pointgraph-dir', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--k', type=int, nargs='+', default=[21, 33])
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    inputs = [args.dbg_dir/'input_paths.json', args.pointgraph_dir/'motifs.tsv',
              args.pointgraph_dir/'tokens.tsv']
    records = json.loads(inputs[0].read_text(encoding='utf-8'))
    if not records or any(set(item['sequence'])-set('ACGT') for item in records):
        raise ValueError('Expected nonempty exact ACGT input paths')
    candidates = [item for item in records if item['kind'] == 'candidate']
    motifs, legacy = read_tsv(inputs[1]), read_tsv(inputs[2])
    units, unit_summary, motif_counts = primitive_units(candidates, motifs, legacy)
    all_queries = []
    for k in args.k:
        for query in branch_queries(records, k):
            query['k'] = k
            candidate = next(item for item in candidates
                             if item['id'] == 'CAND_' + query['candidate_id'])
            query['assigned_molecule_keys'] = set(candidate.get('candidate_molecule_ids', []))
            all_queries.append(query)
    joint, evidence = joint_support(records, all_queries)
    write_tsv(args.out/'joint_queries.tsv',
              ({key: value for key, value in query.items() if key != 'assigned_molecule_keys'}
               for query in all_queries),
              ['query_id','candidate_id','sample','kind','k','start0','end0',
               'branch_edge_count','exact_bases','candidate_occurrences'])
    write_tsv(args.out/'joint_support.tsv', joint, list(joint[0]))
    write_tsv(args.out/'joint_molecules.tsv', evidence, list(evidence[0]))
    write_tsv(args.out/'primitive_units.tsv', units, list(units[0]))
    write_tsv(args.out/'unit_summary.tsv', unit_summary, list(unit_summary[0]))
    write_tsv(args.out/'motif_counts.tsv', motif_counts, list(motif_counts[0]))
    manifest = dict(method_version=VERSION, source_inputs={str(path): sha256(path.read_bytes()).hexdigest()
                                                          for path in inputs}, k=args.k,
                    candidate_count=len(candidates), observed_path_count=sum(x['kind']=='observed' for x in records),
                    joint_query_count=len(joint), unit_row_count=len(units),
                    caveats=['exact sequence witnesses are a lower bound',
                             'nonmatching observations are uninformative, not conflicts',
                             'local-only observations never support full-candidate queries',
                             'candidate assignment is not used as a support filter',
                             'unit decomposition does not change candidate sequences or genotypes'])
    (args.out/'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: manifest[key] for key in ('method_version','candidate_count',
                                                     'observed_path_count','joint_query_count','unit_row_count')}))


if __name__ == '__main__':
    main()
