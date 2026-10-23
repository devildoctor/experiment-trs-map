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

VERSION = 'd1s80_joint_units_v2'
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
                                branch_edge_count=sum(first <= p <= last for p in positions),
                                branch_positions_json=[p for p in positions if first <= p <= last],
                                exact_bases=signature,
                                candidate_occurrences=len(occurrences(sequence, signature))))
        for index, (left, right) in enumerate(zip(merged, merged[1:]), 1):
            start, end = left[0], right[1] + k + 1
            signature = sequence[start:end]
            queries.append(dict(query_id=f'{candidate_id}.k{k}.J{index:03}',
                                candidate_id=candidate_id, sample=record['sample'],
                                kind='adjacent_branch_pair', start0=start, end0=end,
                                branch_edge_count=sum(start <= p <= right[1] for p in positions),
                                branch_positions_json=[p for p in positions if start <= p <= right[1]],
                                exact_bases=signature,
                                candidate_occurrences=len(occurrences(sequence, signature))))
        queries.append(dict(query_id=f'{candidate_id}.k{k}.FULL', candidate_id=candidate_id,
                            sample=record['sample'], kind='full_candidate_exact', start0=0,
                            end0=len(sequence), branch_edge_count=len(positions),
                            branch_positions_json=positions,
                            exact_bases=sequence, candidate_occurrences=1))
    return queries


def anchored_interval(sequence, left, right):
    """A uniquely placed, ordered pair; do not select among repeat occurrences."""
    starts, ends = occurrences(sequence, left), occurrences(sequence, right)
    if len(starts) != 1 or len(ends) != 1 or starts[0] + len(left) > ends[0]:
        return None
    return starts[0], ends[0] + len(right)


def define_events(records, queries):
    """Graph intervals with explicit exact states, not normalized SNV events.

    Alternative interval states require the same unique flanks in another
    candidate from the same sample. These states define tests, never votes.
    """
    candidates = {r['id'].removeprefix('CAND_'): r for r in records if r['kind'] == 'candidate'}
    events = []
    for query in queries:
        candidate = candidates[query['candidate_id']]
        full = query['kind'] == 'full_candidate_exact'
        left, right = query['exact_bases'][:query['k']], query['exact_bases'][-query['k']:]
        interval = anchored_interval(candidate['sequence'], left, right)
        localized = interval == (query['start0'], query['end0'])
        states = defaultdict(list)
        states[query['exact_bases']].append(query['candidate_id'])
        for name, other in candidates.items():
            if name == query['candidate_id'] or other['sample'] != query['sample']:
                continue
            bounds = (0, len(other['sequence'])) if full else (
                anchored_interval(other['sequence'], left, right) if localized else None)
            if bounds:
                states[other['sequence'][bounds[0]:bounds[1]]].append(name)
        event_ids = [query['query_id']]
        if query['kind'] == 'adjacent_branch_pair':
            index = int(query['query_id'].rsplit('J', 1)[1])
            prefix = query['query_id'].rsplit('J', 1)[0]
            event_ids = [f'{prefix}B{index:03}', f'{prefix}B{index+1:03}']
        events.append(dict(query_id=query['query_id'], event_id=query['query_id'],
                           candidate_id=query['candidate_id'], sample=query['sample'],
                           kind=query['kind'], start0=query['start0'], end0=query['end0'],
                           target_state_id=sha256(query['exact_bases'].encode()).hexdigest(),
                           left_anchor='' if full else left, right_anchor='' if full else right,
                           localization_status='full_locus_boundaries' if full else
                           'unique_context' if localized else 'unresolved_context',
                           states_json=[dict(state_id=sha256(seq.encode()).hexdigest(),
                                             exact_bases=seq, candidate_ids=sorted(names))
                                        for seq, names in sorted(states.items())],
                           ordered_event_ids_json=event_ids))
    return events


def assess_observation(path, query, event, candidate_sequence):
    result = dict(path_id=path['id'], assessment='uninformative',
                  reason='missing_or_repeated_context', start0=None, end0=None, state_id=None)
    if query['kind'] == 'full_candidate_exact':
        if path['coverage'] != 'two_flank_anchored':
            result['reason'] = 'local_only_cannot_resolve_full_locus'
            return result
        bounds = (0, len(path['sequence']))
    elif path['coverage'] == 'two_flank_anchored' and path['sequence'] == candidate_sequence:
        result.update(assessment='support', reason='exact_full_target',
                      start0=query['start0'], end0=query['end0'], state_id=event['target_state_id'])
        return result
    elif event['localization_status'] == 'unique_context':
        bounds = anchored_interval(path['sequence'], event['left_anchor'], event['right_anchor'])
    else:
        return result
    if bounds is None:
        return result
    start, end = bounds
    result.update(start0=start, end0=end, reason='unmodeled_sequence_or_error')
    for state in event['states_json']:
        if path['sequence'][start:end] == state['exact_bases']:
            is_target = state['state_id'] == event['target_state_id']
            result.update(assessment='support' if is_target else 'conflict',
                          reason='exact_target_state' if is_target else 'exact_competitor_state',
                          state_id=state['state_id'])
            break
    return result


def joint_support(records, queries, events=None):
    """One vote per molecule; classify explicit states and retain raw witnesses.

    Unknown mismatches do not become conflicts. Mixed resolved states from one
    molecule are uninformative. Separate fragments are never joined.
    """
    if events is None:
        events = define_events(records, queries)
    event_by_id = {event['query_id']: event for event in events}
    candidates = {r['id'].removeprefix('CAND_'): r for r in records if r['kind'] == 'candidate'}
    observed = [record for record in records if record['kind'] == 'observed']
    by_sample = defaultdict(lambda: defaultdict(list))
    for record in observed:
        by_sample[record['sample']][record['molecule_key']].append(record)
    summaries, evidence = [], []
    for query in queries:
        molecules = by_sample[query['sample']]
        exact = 0
        assigned = 0
        decisions = Counter()
        event = event_by_id[query['query_id']]
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
                exact += 1
                assigned += key in query.get('assigned_molecule_keys', set())
            observations = [assess_observation(path, query, event,
                            candidates[query['candidate_id']]['sequence']) for path in paths]
            resolved = {o['assessment'] for o in observations if o['assessment'] != 'uninformative'}
            decision = next(iter(resolved)) if len(resolved) == 1 else 'uninformative'
            reason = ('mixed_states_in_one_molecule' if len(resolved) > 1 else
                      'resolved_exact_state' if resolved else 'no_resolved_state')
            decisions[decision] += 1
            evidence.append(dict(query_id=query['query_id'], candidate_id=query['candidate_id'],
                                 sample=query['sample'], molecule_key=key,
                                 decision=decision, reason=reason,
                                 occurrence_count=len(hits), occurrences_json=hits,
                                 observations_json=observations,
                                 assigned_to_candidate=key in query.get('assigned_molecule_keys', set())))
        informative = decisions['support'] + decisions['conflict']
        summaries.append(dict(query_id=query['query_id'], candidate_id=query['candidate_id'],
                              sample=query['sample'], kind=query['kind'], k=query['k'],
                              start0=query['start0'], end0=query['end0'],
                              query_bp=len(query['exact_bases']),
                              query_sha256=sha256(query['exact_bases'].encode()).hexdigest(),
                              candidate_occurrences=query['candidate_occurrences'],
                              branch_edge_count=query['branch_edge_count'],
                              tested_molecules=len(molecules), exact_support_molecules=exact,
                              assigned_exact_support_molecules=assigned,
                              informative_molecules=informative,
                              support_molecules=decisions['support'],
                              conflict_molecules=decisions['conflict'],
                              uninformative_molecules=decisions['uninformative'],
                              support_fraction=decisions['support']/informative if informative else 'NA',
                              status='computed' if informative else 'no_informative_molecules',
                              interpretation='exact_defined_states_not_allele_frequency'))
    return summaries, evidence


def template_for(token):
    if token['status'] == 'unknown':
        return ''
    family = token['family']
    oriented = family if token['orientation'] == '+' else revcomp(family)
    phase = token['phase']
    return oriented[phase:] + oriented[:phase]


def primitive_units(records, motifs, legacy_tokens):
    """Annotate every saved sequence; fragment counts never become locus CN.

    The core decomposer retains two paths, not every tied annotation. A tied
    result therefore remains ambiguous even when the two displayed counts agree.
    """
    from vntrmotifs import primitive_unit

    templates = sorted({row['motif'] for row in motifs
                        if row['usable'] in ('True', True) and row['level'] == 'base'
                        and len(row['motif']) == UNIT_BP
                        and len(primitive_unit(row['motif'])) == UNIT_BP})
    if not templates:
        raise ValueError('No usable 16-bp D1S80 templates')
    legacy = defaultdict(list)
    for row in legacy_tokens:
        legacy[row['allele']].append(row)
    units, summaries, counts = [], [], []
    for record in records:
        path_id = record['id']
        is_candidate = record['kind'] == 'candidate'
        name = path_id.removeprefix('CAND_') if is_candidate else 'NA'
        sequence = record['sequence']
        sequence_hash = sha256(sequence.encode()).hexdigest()
        if record.get('sequence_sha256', sequence_hash) != sequence_hash:
            raise ValueError(f'{path_id}: sequence hash does not match input')
        coverage = record.get('coverage', 'unknown')
        fragment = record['kind'] == 'observed' and coverage != 'two_flank_anchored'
        count_scope = ('observed_fragment' if fragment else
                       'observed_full_interval' if record['kind'] == 'observed' else
                       'candidate_sequence' if is_candidate else 'reference_catalog_interval')
        metadata = dict(path_id=path_id, kind=record['kind'], candidate_id=name,
                        sample=record.get('sample', ''), coverage=coverage,
                        molecule_key=record.get('molecule_key', 'NA'),
                        platforms_json=record.get('platforms', []),
                        sequence_sha256=sequence_hash, count_scope=count_scope)
        old_tokens = legacy.get(name, []) if is_candidate else []
        primary, _, audit = decompose(sequence, tuple(templates), return_details=True)
        schemes = {'primary': [dict(vars(item)) for item in primary]}
        if audit['alternative']:
            schemes['alternative'] = [dict(item) for item in audit['alternative']]
        def signature(token):
            return tuple(token[key] for key in ('start', 'end', 'status', 'family',
                                                'orientation', 'phase', 'exact_bases',
                                                'edit_operations'))
        signatures = {scheme: {signature(token) for token in tokens}
                      for scheme, tokens in schemes.items()}
        path_units = []
        for scheme, tokens in schemes.items():
            if ''.join(token['exact_bases'] for token in tokens) != sequence:
                raise AssertionError(f'{path_id} {scheme}: units do not reconstruct sequence')
            position = 0
            for index, token in enumerate(tokens, 1):
                if token['start'] != position or token['end'] != position+len(token['exact_bases']):
                    raise AssertionError(f'{path_id} {scheme}: coordinates do not partition sequence')
                if token['exact_bases'] != sequence[token['start']:token['end']]:
                    raise AssertionError(f'{path_id} {scheme}: unit differs from source slice')
                position = token['end']
                parents = [i for i, old in enumerate(old_tokens, 1)
                           if int(old['start']) < token['end'] and int(old['end']) > token['start']]
                other = 'alternative' if scheme == 'primary' else 'primary'
                token_ambiguous = (audit['annotation_ambiguity']
                                   and signature(token) not in signatures.get(other, set()))
                path_units.append(dict(**metadata, decomposition=scheme, unit_index=index,
                                       start0=token['start'], end0=token['end'],
                                       status=token['status'], family=token['family'],
                                       template_bases=template_for(token),
                                       exact_bases=token['exact_bases'],
                                       observed_unit_bp=len(token['exact_bases']),
                                       orientation=token['orientation'], phase_offset=token['phase'],
                                       phase_policy=audit['phase_policy'],
                                       edit_operations=token['edit_operations'],
                                       annotation_ambiguity=token_ambiguous,
                                       legacy_token_indexes_json=parents if is_candidate else 'NA'))
            if scheme == 'primary':
                full = sum(token['status'] == 'full' for token in tokens)
                partial = sum(token['status'] == 'partial' for token in tokens)
                unknown = sum(len(token['exact_bases']) for token in tokens
                              if token['status'] == 'unknown')
                alternate = audit['alternative']
                alt_count = sum(item['status'] == 'full' for item in alternate) if alternate else None
                annotation_status = ('unresolved' if unknown or partial else
                                     'ambiguous' if audit['annotation_ambiguity'] else 'resolved_count')
                count_status = 'fragment_only' if fragment else annotation_status
                summaries.append(dict(**metadata,
                                      sequence_bp=len(sequence), nominal_period_bp=UNIT_BP,
                                      length_equiv_cn=len(sequence)/UNIT_BP,
                                      full_16bp_units=full, partial_units=partial,
                                      unknown_bp=unknown, primitive_cn=full if count_status == 'resolved_count' else 'NA',
                                      primitive_cn_status=count_status,
                                      annotation_count_status=annotation_status,
                                      annotation_ambiguity=audit['annotation_ambiguity'],
                                      alternative_full_units=alt_count if alt_count is not None else 'NA',
                                      primary_score=audit['primary_score'],
                                      alternative_score=audit['alternative_score'] if alternate else 'NA',
                                      score_gap=audit['score_gap'] if alternate else 'NA',
                                      alternative_is_tied=bool(alternate) and audit['annotation_ambiguity'],
                                      alternative_limit=audit['alternative_limit'],
                                      score_definition=audit['score_definition'],
                                      phase_policy=audit['phase_policy'],
                                      legacy_full_tokens=sum(row['status'] == 'full' for row in old_tokens)
                                      if is_candidate else 'NA',
                                      method_version=VERSION))
        units.extend(path_units)
        grouped = defaultdict(lambda: [0, 0])
        for row in path_units:
            if row['decomposition'] == 'primary':
                key = (row['family'], row['status'])
                grouped[key][0] += 1
                grouped[key][1] += len(row['exact_bases'])
        for (family, status), (segments, bp) in sorted(grouped.items()):
            counts.append(dict(**metadata,
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
    if len(set(args.k)) != len(args.k) or any(k < 2 for k in args.k):
        parser.error('--k must contain distinct integers >= 2')
    inputs = [args.dbg_dir/'input_paths.json', args.pointgraph_dir/'motifs.tsv',
              args.pointgraph_dir/'tokens.tsv', args.dbg_dir/'input_summary.json']
    source_summary = json.loads(inputs[3].read_text(encoding='utf-8'))
    if source_summary['locus'] != 'D1S80':
        raise ValueError('This 16-bp analysis is explicitly restricted to D1S80')
    for path in inputs[1:3]:
        if sha256(path.read_bytes()).hexdigest() != source_summary['source_evidence_hashes'][path.name]:
            raise ValueError(f'Pointgraph source differs from frozen dBG provenance: {path.name}')
    records = json.loads(inputs[0].read_text(encoding='utf-8'))
    if not records or any(not item['sequence'] or set(item['sequence'])-set('ACGT') for item in records):
        raise ValueError('Expected nonempty exact ACGT input paths')
    if len({item['id'] for item in records}) != len(records):
        raise ValueError('Duplicate path IDs')
    if Counter(item['kind'] for item in records) != source_summary['path_types']:
        raise ValueError('Path counts differ from frozen input summary')
    for record in records:
        if sha256(record['sequence'].encode()).hexdigest() != record['sequence_sha256']:
            raise ValueError(f'Sequence hash mismatch: {record["id"]}')
    candidates = [item for item in records if item['kind'] == 'candidate']
    if not candidates or not any(item['kind'] == 'observed' for item in records):
        raise ValueError('Candidate and observed paths are required')
    args.out.mkdir(parents=True, exist_ok=False)
    motifs, legacy = read_tsv(inputs[1]), read_tsv(inputs[2])
    units, unit_summary, motif_counts = primitive_units(records, motifs, legacy)
    all_queries = []
    for k in args.k:
        for query in branch_queries(records, k):
            query['k'] = k
            candidate = next(item for item in candidates
                             if item['id'] == 'CAND_' + query['candidate_id'])
            query['assigned_molecule_keys'] = set(candidate.get('candidate_molecule_ids', []))
            all_queries.append(query)
    events = define_events(records, all_queries)
    joint, evidence = joint_support(records, all_queries, events)
    write_tsv(args.out/'joint_queries.tsv',
              ({key: value for key, value in query.items() if key != 'assigned_molecule_keys'}
               for query in all_queries),
              ['query_id','candidate_id','sample','kind','k','start0','end0',
               'branch_edge_count','branch_positions_json','exact_bases','candidate_occurrences'])
    write_tsv(args.out/'events.tsv', events, list(events[0]))
    write_tsv(args.out/'joint_support.tsv', joint, list(joint[0]))
    write_tsv(args.out/'joint_molecules.tsv', evidence, list(evidence[0]))
    molecule_events = [dict(query_id=row['query_id'], event_id=row['query_id'],
                           sample=row['sample'], molecule_key=row['molecule_key'], **observation)
                       for row in evidence for observation in row['observations_json']]
    write_tsv(args.out/'molecule_events.tsv', molecule_events, list(molecule_events[0]))
    write_tsv(args.out/'primitive_units.tsv', units, list(units[0]))
    write_tsv(args.out/'unit_summary.tsv', unit_summary, list(unit_summary[0]))
    write_tsv(args.out/'motif_counts.tsv', motif_counts, list(motif_counts[0]))
    manifest = dict(method_version=VERSION, source_inputs={str(path): sha256(path.read_bytes()).hexdigest()
                                                          for path in inputs}, k=args.k,
                    candidate_count=len(candidates), observed_path_count=sum(x['kind']=='observed' for x in records),
                    decomposed_path_count=len(unit_summary), path_types=dict(Counter(x['kind'] for x in records)),
                    joint_query_count=len(joint), unit_row_count=len(units), locus_id='D1S80',
                    caveats=['exact sequence witnesses are a lower bound',
                             'only exact defined alternate states count as conflicts; other mismatches are uninformative',
                             'events are graph intervals, not normalized biological variants',
                             'support fraction is conditional on exact modeled states, not an allele frequency',
                             'same-molecule mixed states are uninformative; fragments are never joined',
                             'local-only observations never support full-candidate queries',
                             'local-only unit counts do not establish whole-locus copy number',
                             'one second-best decomposition is retained, not all optimal paths',
                             'candidate assignment is not used as a support filter',
                             'unit decomposition does not change candidate sequences or genotypes'])
    (args.out/'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: manifest[key] for key in ('method_version','candidate_count',
                                                     'observed_path_count','decomposed_path_count',
                                                     'joint_query_count','unit_row_count')}))


if __name__ == '__main__':
    main()
