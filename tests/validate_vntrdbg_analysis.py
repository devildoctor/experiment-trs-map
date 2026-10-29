"""Check exported evidence against source sequences without importing its producer."""
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import sys


def check(ok, message):
    if not ok:
        raise AssertionError(message)


def rows(path):
    csv.field_size_limit(64 * 1024 * 1024)
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def indexed(items, field, label):
    result = {row[field]: row for row in items}
    check(len(result) == len(items), f'{label}: duplicate {field}')
    return result


def starts(sequence, needle):
    return [i for i in range(len(sequence)-len(needle)+1)
            if sequence[i:i+len(needle)] == needle] if needle else []


def truth(value):
    check(value in ('True', 'False'), f'invalid boolean: {value!r}')
    return value == 'True'


def metadata(row, source):
    for field in ('kind', 'sample', 'coverage'):
        check(row[field] == source.get(field, '' if field == 'sample' else 'unknown'),
              f'{source["id"]}: wrong {field}')
    name = source['id'].removeprefix('CAND_') if source['kind'] == 'candidate' else 'NA'
    check(row['candidate_id'] == name, f'{source["id"]}: wrong candidate id')
    check(row['molecule_key'] == source.get('molecule_key', 'NA'), 'wrong molecule metadata')
    check(json.loads(row['platforms_json']) == source.get('platforms', []), 'wrong platform metadata')
    check(row['sequence_sha256'] == sha256(source['sequence'].encode()).hexdigest(), 'wrong sequence hash')
    fragment = source['kind'] == 'observed' and source['coverage'] != 'two_flank_anchored'
    scope = ('observed_fragment' if fragment else 'observed_full_interval' if source['kind'] == 'observed'
             else 'candidate_sequence' if source['kind'] == 'candidate' else 'reference_catalog_interval')
    check(row['count_scope'] == scope, 'wrong unit count scope')


def expected_queries(records, ks, minimum=3):
    """Rebuild support from k+1 strings, then union overlapping signature intervals."""
    expected = {}
    for k in ks:
        support, keep = defaultdict(lambda: defaultdict(set)), set()
        for record in records:
            sequence = record['sequence']
            edges = {sequence[i:i+k+1] for i in range(len(sequence)-k)}
            if record['kind'] in ('candidate', 'reference'):
                keep.update(edges)
            elif record['kind'] == 'observed':
                for edge in edges:
                    support[edge][record['sample']].add(record['molecule_key'])
        keep.update(edge for edge, samples in support.items() if any(len(v) >= minimum for v in samples.values()))
        outgoing, incoming = Counter(e[:-1] for e in keep), Counter(e[1:] for e in keep)
        for record in records:
            if record['kind'] != 'candidate':
                continue
            sequence, name = record['sequence'], record['id'].removeprefix('CAND_')
            branch = [i for i in range(len(sequence)-k)
                      if outgoing[sequence[i:i+k]] > 1 or incoming[sequence[i+1:i+k+1]] > 1]
            blocks = []
            for pos in branch:
                if blocks and pos < blocks[-1][1]:
                    blocks[-1] = blocks[-1][0], pos+k+1
                else:
                    blocks.append((pos, pos+k+1))
            definitions = [(f'B{i:03}', 'branch_block', a, b) for i, (a, b) in enumerate(blocks, 1)]
            definitions += [(f'J{i:03}', 'adjacent_branch_pair', left[0], right[1])
                            for i, (left, right) in enumerate(zip(blocks, blocks[1:]), 1)]
            definitions.append(('FULL', 'full_candidate_exact', 0, len(sequence)))
            for suffix, kind, a, b in definitions:
                expected[f'{name}.k{k}.{suffix}'] = dict(candidate_id=name, sample=record['sample'],
                    kind=kind, k=k, start0=a, end0=b, exact_bases=sequence[a:b],
                    branch_positions=[p for p in branch if a <= p and p+k+1 <= b])
    return expected


def anchor_interval(sequence, left, right):
    a, b = starts(sequence, left), starts(sequence, right)
    if len(a) == len(b) == 1 and a[0]+len(left) <= b[0]:
        return a[0], b[0]+len(right)
    return None


def assessment(source, query, event, states, candidate):
    sequence, target = source['sequence'], event['target_state_id']
    if query['kind'] == 'full_candidate_exact':
        interval = (0, len(sequence)) if source['coverage'] == 'two_flank_anchored' else None
    elif source['coverage'] == 'two_flank_anchored' and sequence == candidate['sequence']:
        return 'support', int(query['start0']), int(query['end0']), target, True
    elif event['localization_status'] == 'unique_context':
        interval = anchor_interval(sequence, event['left_anchor'], event['right_anchor'])
    else:
        interval = None
    if interval is not None:
        a, b = interval
        found = [s['state_id'] for s in states.values() if s['exact_bases'] == sequence[a:b]]
        if found:
            return ('support' if found[0] == target else 'conflict'), a, b, found[0], False
        return 'uninformative', a, b, None, False
    return 'uninformative', None, None, None, False


def validate_joint(records, manifest, out):
    candidates = {r['id'].removeprefix('CAND_'): r for r in records if r['kind'] == 'candidate'}
    source_by_id = {r['id']: r for r in records}
    molecules = defaultdict(lambda: defaultdict(list))
    for source in records:
        if source['kind'] == 'observed':
            molecules[source['sample']][source['molecule_key']].append(source)
    queries = indexed(rows(out/'joint_queries.tsv'), 'query_id', 'queries')
    summaries = indexed(rows(out/'joint_support.tsv'), 'query_id', 'support')
    events = indexed(rows(out/'events.tsv'), 'query_id', 'events')
    event_ids = indexed(list(events.values()), 'event_id', 'event ids')
    # 旧结果未记录门槛时沿用历史默认 3；新结果必须按实际参数重建分支。
    minimum = manifest.get('branch_min_support', 3)
    check(type(minimum) is int and minimum >= 1, 'invalid branch support threshold')
    expected = expected_queries(records, manifest['k'], minimum)
    check(set(queries) == set(summaries) == set(events) == set(expected), 'query set is incomplete')
    evidence = rows(out/'joint_molecules.tsv')
    flattened = [dict(query_id=r['query_id'], event_id=r['query_id'], sample=r['sample'],
                      molecule_key=r['molecule_key'], **o)
                 for r in evidence for o in json.loads(r['observations_json'])]
    actual_flat = rows(out/'molecule_events.tsv')
    serialized = [{k: '' if v is None else str(v) for k, v in r.items()} for r in flattened]
    check(actual_flat == serialized, 'molecule_events differs from nested observation evidence')
    by_query = defaultdict(dict)
    for row in evidence:
        qid, key = row['query_id'], row['molecule_key']
        check(qid in queries and key not in by_query[qid], 'unexpected/repeated molecule query')
        by_query[qid][key] = row
    for qid, query in queries.items():
        candidate = candidates[query['candidate_id']]
        for field in ('candidate_id', 'sample', 'kind', 'exact_bases'):
            check(query[field] == expected[qid][field], f'{qid}: wrong {field}')
        for field in ('k', 'start0', 'end0'):
            check(int(query[field]) == expected[qid][field], f'{qid}: wrong {field}')
        branch = expected[qid]['branch_positions']
        check(json.loads(query['branch_positions_json']) == branch, f'{qid}: wrong branch positions')
        check(int(query['branch_edge_count']) == len(branch), f'{qid}: nonbranch positions counted')
        check(int(query['candidate_occurrences']) == len(starts(candidate['sequence'], query['exact_bases'])),
              f'{qid}: wrong candidate occurrences')
        event = events[qid]
        for field in ('candidate_id', 'sample', 'kind', 'start0', 'end0'):
            check(event[field] == query[field], f'{qid}: inconsistent event {field}')
        state_list = json.loads(event['states_json'])
        states = indexed(state_list, 'state_id', qid+' states')
        check(len({s['exact_bases'] for s in state_list}) == len(states), f'{qid}: duplicate sequence states')
        check(event['target_state_id'] in states, f'{qid}: absent target state')
        check(states[event['target_state_id']]['exact_bases'] == query['exact_bases'], f'{qid}: wrong target state')
        status = event['localization_status']
        check(status in ('full_locus_boundaries', 'unique_context', 'unresolved_context'), f'{qid}: unknown context')
        if query['kind'] == 'full_candidate_exact':
            check(status == 'full_locus_boundaries', f'{qid}: full state lacks locus boundary')
        elif status == 'unique_context':
            check(anchor_interval(candidate['sequence'], event['left_anchor'], event['right_anchor']) ==
                  (int(query['start0']), int(query['end0'])), f'{qid}: ambiguous target anchors')
        expected_states = defaultdict(set)
        expected_states[query['exact_bases']].add(query['candidate_id'])
        for name, other in candidates.items():
            if other['sample'] != query['sample'] or name == query['candidate_id']:
                continue
            bounds = (0, len(other['sequence'])) if query['kind'] == 'full_candidate_exact' else (
                anchor_interval(other['sequence'], event['left_anchor'], event['right_anchor'])
                if status == 'unique_context' else None)
            if bounds is not None:
                expected_states[other['sequence'][bounds[0]:bounds[1]]].add(name)
        check({s['exact_bases']: set(s['candidate_ids']) for s in state_list} == dict(expected_states),
              f'{qid}: missing or invented candidate state')
        for state in state_list:
            check(state['state_id'] == sha256(state['exact_bases'].encode()).hexdigest(), f'{qid}: state hash mismatch')
            check(bool(state['candidate_ids']) and len(set(state['candidate_ids'])) == len(state['candidate_ids']),
                  f'{qid}: state candidate sources absent/repeated')
            for name in state['candidate_ids']:
                name = name.removeprefix('CAND_')
                other = candidates[name]
                check(other['sample'] == query['sample'], f'{qid}: state crosses samples')
                if query['kind'] == 'full_candidate_exact':
                    sequence = other['sequence']
                elif name == query['candidate_id']:
                    sequence = other['sequence'][int(query['start0']):int(query['end0'])]
                else:
                    interval = anchor_interval(other['sequence'], event['left_anchor'], event['right_anchor'])
                    check(status == 'unique_context' and interval is not None, f'{qid}: alternative lacks unique context')
                    sequence = other['sequence'][interval[0]:interval[1]]
                check(sequence == state['exact_bases'], f'{qid}: state changed source sequence')
        ordered = json.loads(event['ordered_event_ids_json'])
        expected_order = [qid]
        if query['kind'] == 'adjacent_branch_pair':
            prefix, number = qid.rsplit('J', 1)
            expected_order = [f'{prefix}B{int(number):03}', f'{prefix}B{int(number)+1:03}']
        check(ordered == expected_order, f'{qid}: wrong constituent event order')
        check(len(ordered) == len(set(ordered)), f'{qid}: duplicate ordered event')
        previous = -1
        for eid in ordered:
            check(eid in event_ids, f'{qid}: nonexistent ordered event')
            child = queries[event_ids[eid]['query_id']]
            check(child['candidate_id'] == query['candidate_id'] and child['k'] == query['k'], f'{qid}: event context mismatch')
            check(int(child['start0']) >= previous, f'{qid}: reversed event order')
            previous = int(child['start0'])
        sample_molecules = molecules[query['sample']]
        check(set(by_query[qid]) == set(sample_molecules), f'{qid}: omitted/extra sample molecules')
        decisions, exact_count, assigned_count = Counter(), 0, 0
        assigned = set(candidate.get('candidate_molecule_ids', []))
        for key, paths in sample_molecules.items():
            row = by_query[qid][key]
            check(row['candidate_id'] == query['candidate_id'] and row['sample'] == query['sample'], 'wrong evidence metadata')
            check(truth(row['assigned_to_candidate']) == (key in assigned), f'{qid}/{key}: wrong assignment')
            hits, expected_hits = json.loads(row['occurrences_json']), []
            for source in paths:
                positions = (starts(source['sequence'], query['exact_bases']) if query['kind'] != 'full_candidate_exact'
                             else [0] if source['coverage'] == 'two_flank_anchored'
                             and source['sequence'] == query['exact_bases'] else [])
                expected_hits.extend((source['id'], pos, pos+len(query['exact_bases']), source['coverage']) for pos in positions)
            actual_hits = [(h['path_id'], h['start0'], h['end0'], h['coverage']) for h in hits]
            check(sorted(actual_hits) == sorted(expected_hits), f'{qid}/{key}: missing/false exact occurrence')
            check(int(row['occurrence_count']) == len(hits), f'{qid}/{key}: wrong hit count')
            exact_count += bool(hits)
            assigned_count += bool(hits) and key in assigned
            observed = indexed(json.loads(row['observations_json']), 'path_id', qid+'/'+key+' observations')
            check(set(observed) == {p['id'] for p in paths}, f'{qid}/{key}: missing path assessment')
            informative = set()
            for path_id, observed_state in observed.items():
                decision, a, b, state, full_target = assessment(source_by_id[path_id], query, event, states, candidate)
                check(observed_state['assessment'] == decision, f'{qid}/{path_id}: false assessment')
                check(bool(observed_state['reason']), f'{qid}/{path_id}: reason missing')
                check((observed_state['start0'], observed_state['end0'], observed_state['state_id']) == (a, b, state),
                      f'{qid}/{path_id}: wrong state slice')
                if decision != 'uninformative':
                    if full_target:
                        check(observed_state['reason'] == 'exact_full_target', f'{qid}/{path_id}: full target not labeled')
                    informative.add(decision)
                else:
                    check(observed_state.get('state_id') in (None, '', 'NA'), f'{qid}/{path_id}: uninformative state assigned')
            decision = next(iter(informative)) if len(informative) == 1 else 'uninformative'
            check(row['decision'] == decision and bool(row['reason']), f'{qid}/{key}: wrong molecule decision')
            reason = ('mixed_states_in_one_molecule' if len(informative) > 1 else
                      'resolved_exact_state' if informative else 'no_resolved_state')
            check(row['reason'] == reason, f'{qid}/{key}: wrong molecule decision reason')
            decisions[decision] += 1
        summary = summaries[qid]
        for field in ('candidate_id', 'sample', 'kind', 'k', 'start0', 'end0', 'candidate_occurrences', 'branch_edge_count'):
            check(summary[field] == query[field], f'{qid}: inconsistent summary {field}')
        check(int(summary['query_bp']) == len(query['exact_bases']), f'{qid}: wrong query length')
        check(summary['query_sha256'] == sha256(query['exact_bases'].encode()).hexdigest(), f'{qid}: wrong query hash')
        informative_count = decisions['support']+decisions['conflict']
        totals = dict(tested_molecules=len(sample_molecules), exact_support_molecules=exact_count,
                      assigned_exact_support_molecules=assigned_count, support_molecules=decisions['support'],
                      conflict_molecules=decisions['conflict'], informative_molecules=informative_count,
                      uninformative_molecules=decisions['uninformative'])
        for field, count in totals.items():
            check(int(summary[field]) == count, f'{qid}: wrong {field}')
        check(math.isclose(float(summary['support_fraction']), decisions['support']/informative_count,
                           abs_tol=1e-9, rel_tol=1e-9) if informative_count else summary['support_fraction'] == 'NA',
              f'{qid}: wrong informative support fraction')
    return len(queries), len(evidence)


def validate_edit(row):
    template, sequence, cigar = row['template_bases'], row['exact_bases'], row['edit_operations']
    if row['status'] == 'unknown':
        check(template == cigar == '' and row['family'] == 'UNKNOWN', 'unknown unit has template')
        return
    check(row['status'] in ('full', 'partial') and len(template) == 16, 'invalid primitive template')
    family = row['family']
    check(len(family) == 16 and row['orientation'] in ('+', '-'), 'invalid family/orientation')
    oriented = family if row['orientation'] == '+' else family.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    phase = int(row['phase_offset'])
    check(0 <= phase < 16 and oriented[phase:]+oriented[:phase] == template, 'wrong phase/template')
    operations = re.findall(r'([1-9][0-9]*)([=XID])', cigar)
    check(bool(operations) and ''.join(n+op for n, op in operations) == cigar, 'malformed CIGAR')
    a = b = 0
    for count, operation in operations:
        n = int(count)
        if operation in ('=', 'X'):
            check(a+n <= len(template) and b+n <= len(sequence), 'CIGAR alignment overrun')
            check(all((x == y) == (operation == '=') for x, y in zip(template[a:a+n], sequence[b:b+n])),
                  'CIGAR base comparison incorrect')
            a += n
            b += n
        elif operation == 'I':
            b += n
        else:
            a += n
        check(a <= len(template) and b <= len(sequence), 'CIGAR consumption overrun')
    check(a == len(template) and b == len(sequence), 'CIGAR must consume full template and observed slice')


def validate_units(records, out):
    sources = {r['id']: r for r in records}
    units = rows(out/'primitive_units.tsv')
    summaries = indexed(rows(out/'unit_summary.tsv'), 'path_id', 'unit summaries')
    check(set(summaries) == set(sources), 'unit summaries omit source paths')
    grouped = defaultdict(lambda: defaultdict(list))
    for row in units:
        check(row['path_id'] in sources and row['decomposition'] in ('primary', 'alternative'), 'unknown unit path/scheme')
        metadata(row, sources[row['path_id']])
        validate_edit(row)
        grouped[row['path_id']][row['decomposition']].append(row)
    check(set(grouped) == set(sources), 'unit decomposition omits source paths')
    expected_counts = {}
    for path_id, source in sources.items():
        summary, sequence = summaries[path_id], source['sequence']
        metadata(summary, source)
        check('primary' in grouped[path_id], f'{path_id}: primary scheme absent')
        for scheme, values in grouped[path_id].items():
            values.sort(key=lambda r: int(r['unit_index']))
            check([int(r['unit_index']) for r in values] == list(range(1, len(values)+1)), 'wrong unit indexes')
            pos = 0
            for row in values:
                a, b = int(row['start0']), int(row['end0'])
                check(a == pos and a < b <= len(sequence), f'{path_id}/{scheme}: noncontiguous units')
                check(row['exact_bases'] == sequence[a:b] and int(row['observed_unit_bp']) == b-a, 'unit changed source slice')
                if row['status'] == 'partial':
                    check((a == 0 or b == len(sequence)) and b-a < 16, 'internal/full-length partial unit')
                pos = b
            check(pos == len(sequence), f'{path_id}/{scheme}: incomplete sequence')
        primary, alternate = grouped[path_id]['primary'], grouped[path_id].get('alternative', [])
        full = sum(r['status'] == 'full' for r in primary)
        partial = sum(r['status'] == 'partial' for r in primary)
        unknown = sum(len(r['exact_bases']) for r in primary if r['status'] == 'unknown')
        for field, value in dict(sequence_bp=len(sequence), nominal_period_bp=16,
                                 full_16bp_units=full, partial_units=partial, unknown_bp=unknown).items():
            check(int(summary[field]) == value, f'{path_id}: wrong {field}')
        check(math.isclose(float(summary['length_equiv_cn']), len(sequence)/16), 'wrong length quotient')
        ambiguous = truth(summary['annotation_ambiguity'])
        status = 'unresolved' if unknown or partial else 'ambiguous' if ambiguous else 'resolved_count'
        count_status = 'fragment_only' if summary['count_scope'] == 'observed_fragment' else status
        withheld = summary.get('cn_policy') == 'withheld_terminal_validation'
        reported_status = 'withheld' if withheld and count_status == 'resolved_count' else count_status
        check(summary['annotation_count_status'] == status and summary['primitive_cn_status'] == reported_status,
              f'{path_id}: overclaimed annotation/count status')
        check(summary['primitive_cn'] == (str(full) if count_status == 'resolved_count' and not withheld else 'NA'),
              'overclaimed primitive CN')
        check(summary['alternative_full_units'] == (str(sum(r['status'] == 'full' for r in alternate)) if alternate else 'NA'),
              'wrong alternative unit count')
        check(int(summary['alternative_limit']) == 1 and summary['score_definition'] == 'descriptive_cost_not_probability',
              'incorrect alternative limit/score interpretation')
        check(summary['phase_policy'] == 'fixed_template_phase', 'wrong phase policy')
        if alternate:
            gap = float(summary['alternative_score'])-float(summary['primary_score'])
            check(gap >= -1e-8 and math.isclose(float(summary['score_gap']), max(0., gap), abs_tol=1e-8), 'wrong score gap')
            check(ambiguous == (abs(gap) < 1e-8) == truth(summary['alternative_is_tied']), 'wrong tie status')
        else:
            check(not ambiguous and not truth(summary['alternative_is_tied']), 'absent alternative tied')
        fields = ('start0', 'end0', 'status', 'family', 'orientation', 'phase_offset', 'exact_bases', 'edit_operations')
        signatures = {scheme: {tuple(r[f] for f in fields) for r in values} for scheme, values in grouped[path_id].items()}
        for scheme, values in grouped[path_id].items():
            other = 'alternative' if scheme == 'primary' else 'primary'
            for row in values:
                check(truth(row['annotation_ambiguity']) == (ambiguous and tuple(row[f] for f in fields) not in signatures.get(other, set())),
                      'wrong per-unit ambiguity flag')
        for row in primary:
            key = path_id, row['family'], row['status']
            count, bp = expected_counts.get(key, (0, 0))
            expected_counts[key] = count+1, bp+len(row['exact_bases'])
    counts = {}
    for row in rows(out/'motif_counts.tsv'):
        key = row['path_id'], row['family'], row['status']
        check(key not in counts and key[0] in sources and row['decomposition'] == 'primary', 'wrong/duplicate motif count')
        metadata(row, sources[key[0]])
        counts[key] = int(row['segment_count']), int(row['observed_bp'])
    check(counts == expected_counts, 'motif counts missing or incorrect')
    return len(units)


def main(out, input_file):
    records = json.loads(input_file.read_text(encoding='utf-8'))
    check(bool(records), 'empty input')
    indexed(records, 'id', 'input paths')
    for record in records:
        check(bool(record['sequence']) and set(record['sequence']) <= set('ACGT'), 'invalid source sequence')
        if 'sequence_sha256' in record:
            check(record['sequence_sha256'] == sha256(record['sequence'].encode()).hexdigest(), 'input hash mismatch')
    manifest = json.loads((out/'manifest.json').read_text(encoding='utf-8'))
    terminal = manifest.get('terminal_units')
    for summary in rows(out/'unit_summary.tsv'):
        if terminal:
            check(terminal.get('enabled') is True and terminal.get('cn_policy') == 'withheld_terminal_validation',
                  'terminal manifest lacks explicit CN withholding')
            check(summary.get('cn_policy') == 'withheld_terminal_validation' and summary['primitive_cn'] == 'NA',
                  'terminal run must withhold all CN')
        else:
            check('cn_policy' not in summary, 'default run unexpectedly overrides CN policy')
    ks = manifest['k']
    check(bool(ks) and len(set(ks)) == len(ks) and all(isinstance(k, int) and k >= 2 for k in ks), 'invalid k values')
    digest, matched = sha256(input_file.read_bytes()).hexdigest(), False
    for name, expected_hash in manifest['source_inputs'].items():
        path = Path(name)
        check(path.is_file(), f'missing manifest source: {name}')
        check(sha256(path.read_bytes()).hexdigest() == expected_hash, f'changed manifest source: {name}')
        if path.name == input_file.name:
            check(expected_hash == digest, 'validator input differs from run input')
            matched = True
    check(matched, 'manifest omits input path hash')
    queries, molecule_rows = validate_joint(records, manifest, out)
    unit_rows = validate_units(records, out)
    check(manifest['candidate_count'] == sum(r['kind'] == 'candidate' for r in records), 'wrong manifest candidate count')
    check(manifest['observed_path_count'] == sum(r['kind'] == 'observed' for r in records), 'wrong manifest observed count')
    check(manifest['decomposed_path_count'] == len(records) and
          manifest['path_types'] == dict(Counter(r['kind'] for r in records)), 'wrong manifest path coverage')
    check(manifest['joint_query_count'] == queries and manifest['unit_row_count'] == unit_rows, 'wrong manifest output counts')
    print(json.dumps(dict(status='PASS', queries=queries, molecule_query_rows=molecule_rows,
                          unit_rows=unit_rows, input_paths=len(records),
                          branch_min_support=manifest.get('branch_min_support', 3),
                          branch_min_support_source='recorded manifest' if 'branch_min_support' in manifest else 'inferred legacy default',
                          candidates=sum(r['kind'] == 'candidate' for r in records))))


if __name__ == '__main__':
    main(Path(sys.argv[1]), Path(sys.argv[2]))
