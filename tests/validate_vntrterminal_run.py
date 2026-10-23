"""Independently validate and compare the fixed 199-path D1S80 terminal run.

Run on the compute node: NEW_OUT BASELINE_OUT INPUT_PATHS POINTGRAPH_DIR.
Raw sequences and molecule identities are read in place; reports contain only
candidate-level counts, coordinates, costs and aggregate validation evidence.
"""
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

from validate_vntrdbg_analysis import check, indexed, main as validate_analysis, rows, truth


UNIT = 16
JOINT_FILES = ('joint_queries.tsv', 'events.tsv', 'joint_support.tsv',
               'joint_molecules.tsv', 'molecule_events.tsv')


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def period(sequence):
    return next(width for width in range(1, len(sequence)+1)
                if len(sequence) % width == 0 and sequence == sequence[:width]*(len(sequence)//width))


def distance(a, b):
    """Plain row-wise Levenshtein, independent of the producer's bit vectors."""
    previous = list(range(len(b)+1))
    for i, x in enumerate(a, 1):
        current = [i]
        for j, y in enumerate(b, 1):
            current.append(min(previous[j]+1, current[j-1]+1, previous[j-1]+(x != y)))
        previous = current
    return previous[-1]


def fasta_slice(path, chrom, start, end):
    entries = [line.split('\t') for line in Path(str(path)+'.fai').read_text().splitlines()]
    found = [entry for entry in entries if entry[0] == chrom]
    check(len(found) == 1, 'reference chromosome missing or duplicate')
    length, offset, bases, width = map(int, found[0][1:5])
    check(0 <= start < end <= length and 0 < bases <= width, 'invalid reference interval')
    chunks = []
    with path.open('rb') as handle:
        while start < end:
            take = min(end-start, bases-start % bases)
            handle.seek(offset+(start//bases)*width+start % bases)
            chunks.append(handle.read(take))
            start += take
    result = b''.join(chunks).decode('ascii').upper()
    check(set(result) <= set('ACGT'), 'ambiguous reference context')
    return result


def reference_context(manifest, records, pointgraph):
    terminal = manifest['terminal_units']
    params = {row['parameter']: row['value'] for row in rows(pointgraph/'parameters.tsv')}
    check(params['boundary_shift'] == '0', 'original boundary shift must be zero')
    catalogs = [Path(name) for name in manifest['source_inputs'] if name.endswith('.bed')]
    check(len(catalogs) == 1 and digest(catalogs[0]) == params['catalog_sha256'], 'catalog changed')
    locus = []
    for line in catalogs[0].read_text(encoding='utf-8').splitlines():
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        info = dict(item.split('=', 1) for item in fields[3].split(';') if '=' in item)
        if info.get('ID') == 'D1S80':
            locus.append((fields[0], int(info['REPEAT_START']), int(info['REPEAT_END']), info['BOUNDARY']))
    check(len(locus) == 1 and locus[0][3] == 'confirmed', 'original D1S80 interval missing')
    chrom, start, end, _ = locus[0]
    check((terminal['chrom'], terminal['repeat_start0'], terminal['repeat_end0']) == (chrom, start, end),
          'manifest boundary differs from frozen catalog')
    reference = Path(terminal['reference_path'])
    stat = reference.stat()
    check(stat.st_size == int(params['reference_size']) == terminal['reference_size'] and
          stat.st_mtime_ns == int(params['reference_mtime_ns']) == terminal['reference_mtime_ns'],
          'reference identity changed')
    check(digest(Path(str(reference)+'.fai')) == params['reference_fai_sha256'], 'reference index changed')
    repeat, outer = fasta_slice(reference, chrom, start, end), fasta_slice(reference, chrom, end, end+UNIT)
    refs = [record['sequence'] for record in records if record['kind'] == 'reference']
    check(refs == [repeat], 'input reference differs from unchanged catalog interval')
    check(sha256(repeat.encode()).hexdigest() == terminal['reference_repeat_sha256'] and
          sha256(outer.encode()).hexdigest() == terminal['right_context_sha256'], 'reference context hash changed')
    return repeat[-UNIT:], outer


def independently_assess(records, pointgraph, templates, outer):
    observations, stored, decisions = defaultdict(list), defaultdict(list), defaultdict(list)
    for record in records:
        if record['kind'] == 'observed':
            observations[record['sample'], record['molecule_key']].append(record)
    for row in rows(pointgraph/'molecules.tsv'):
        stored[row['sample'], row['molecule_key']].append(row)
    for line in (pointgraph/'anchor_decisions.jsonl').read_text(encoding='utf-8').splitlines():
        if line:
            row = json.loads(line)
            decisions[row['molecule_key']].append(row)
    assessed = {}
    for key, paths in observations.items():
        full = [path for path in paths if path['coverage'] == 'two_flank_anchored']
        result = dict(paths=paths, platforms=set(), valid=False, eligible=False,
                      terminal='', inner=None, outer=False)
        assessed[key] = result
        if not full:
            continue
        result['platforms'] = {platform for path in full for platform in path.get('platforms', [])}
        first = full[0]
        signatures = {(path['sequence'], tuple(path.get('stored_query_interval', [])),
                       tuple(path.get('original_query_interval', [])), tuple(sorted(path.get('source_records', []))))
                      for path in full}
        anchor = decisions[key[1]]
        coherent = (len(signatures) == 1 and bool(stored[key]) and bool(anchor)
                    and all(path['sequence_sha256'] == sha256(path['sequence'].encode()).hexdigest()
                            and path.get('anchor_status') == 'complete_support' for path in full)
                    and all(row.get('sample') == key[0] and row.get('status') == 'complete_support'
                            and set(first.get('source_records', [])).issubset(row.get('sources', [])) for row in anchor))
        if not coherent:
            continue
        slices = []
        for raw in stored[key]:
            a, b, whole = int(raw['qL']), int(raw['qR']), raw['query_sequence']
            valid = (0 <= a < b <= len(whole) and whole[a:b] == raw['sequence'] == first['sequence']
                     and [a, b] == first.get('stored_query_interval')
                     and json.loads(raw['original_query_interval']) == first.get('original_query_interval')
                     and sorted(json.loads(raw['records'])) == sorted(first.get('source_records', [])))
            if not valid:
                break
            slices.append((a, b, whole[b:b+UNIT], raw['sequence']))
        if len(slices) != len(stored[key]) or len(set(slices)) != 1:
            continue
        a, b, after, sequence = slices[0]
        result.update(valid=True, qL=a, qR=b)
        if len(sequence) < 2*UNIT:
            continue
        inner = min(distance(sequence[-2*UNIT:-UNIT], template) for template in templates)
        result.update(terminal=sequence[-UNIT:], inner=inner, outer=after == outer,
                      eligible=inner <= 1 and after == outer)
    return assessed


def validate_terminal_evidence(records, manifest, out, pointgraph):
    terminal = manifest['terminal_units']
    for name in ('molecules.tsv', 'anchor_decisions.jsonl', 'parameters.tsv', 'motifs.tsv'):
        expected = [value for path, value in manifest['source_inputs'].items() if Path(path).name == name]
        check(expected == [digest(pointgraph/name)], 'pointgraph directory differs from manifest source')
    check(terminal['candidate_votes'] is False and terminal['side'] == 'right', 'unexpected terminal scope')
    check(terminal['phase_policy'] == 'exact_oriented_template_no_rotation_pooling', 'terminal phase scope changed')
    minimum = terminal['min_support']
    check(isinstance(minimum, int) and minimum >= 3, 'insufficient molecule threshold')
    ref_terminal, outer = reference_context(manifest, records, pointgraph)
    ordinary = {row['motif'] for row in rows(pointgraph/'motifs.tsv') if row['usable'] == 'True'
                and row['level'] == 'base' and len(row['motif']) == UNIT and period(row['motif']) == UNIT}
    ordinary |= {reverse_complement(template) for template in ordinary}
    check(bool(ordinary), 'ordinary primitive dictionary missing')
    assessed = independently_assess(records, pointgraph, ordinary, outer)
    nominations = {ref_terminal} | {item['terminal'] for item in assessed.values() if item['terminal']}
    nominations -= ordinary
    samples = {record['sample'] for record in records if record['kind'] != 'reference'}
    candidates = rows(out/'terminal_unit_candidates.tsv')
    by_template_sample = {(row['template_bases'], row['sample_id']): row for row in candidates}
    check(len(by_template_sample) == len(candidates), 'duplicate template/sample row')
    check(set(by_template_sample) == {(template, sample) for template in nominations for sample in samples},
          'missing or invented terminal nominations')
    evidence = rows(out/'terminal_molecules.tsv')
    by_evidence = {(row['terminal_id'], row['sample_id'], row['molecule_key']): row for row in evidence}
    check(len(by_evidence) == len(evidence), 'duplicate terminal/molecule evidence')
    expected_evidence, accepted, common = set(), defaultdict(set), {}
    for (template, sample), row in by_template_sample.items():
        identifier = 'D1S80_R_'+sha256(template.encode()).hexdigest()[:16]
        check(row['terminal_id'] == identifier and row['side'] == 'right' and row['role'] == 'terminal'
              and row['locus_id'] == 'D1S80' and int(row['template_length']) == UNIT, 'incorrect terminal identity')
        forms = [sequence[offset:]+sequence[:offset] for sequence in (template, reverse_complement(template))
                 for offset in range(UNIT)]
        check(row['family'] == min(forms) and row['orientation'] in ('+', '-'), 'incorrect terminal family/orientation')
        oriented = row['family'] if row['orientation'] == '+' else reverse_complement(row['family'])
        phase = int(row['phase'])
        check(0 <= phase < UNIT and oriented[phase:]+oriented[:phase] == template, 'incorrect terminal template phase')
        check(truth(row['reference_terminal']) == (template == ref_terminal), 'wrong reference-terminal flag')
        expected_coordinates = (str(terminal['repeat_end0']-UNIT), str(terminal['repeat_end0'])) if template == ref_terminal else ('NA', 'NA')
        check((row['reference_start0'], row['reference_end0']) == expected_coordinates, 'wrong terminal reference coordinates')
        check(int(row['min_support']) == minimum and row['dictionary_version'] == terminal['dictionary_version'],
              'terminal threshold/version mismatch')
        support = {key for key, item in assessed.items() if key[0] == sample and item['eligible']
                   and item['terminal'] == template and period(template) == UNIT}
        platforms = Counter(platform for key in support for platform in assessed[key]['platforms'])
        check(int(row['support_molecules']) == len(support) and json.loads(row['platform_counts_json']) == dict(platforms),
              'terminal support/platform totals do not match raw independent molecules')
        proposed = {key for key, item in assessed.items() if key[0] == sample and item['terminal'] == template}
        check(int(row['nominated_observed_molecules']) == len(proposed), 'wrong nomination molecule count')
        source = (['reference'] if template == ref_terminal else [])
        if any(item['terminal'] == template for item in assessed.values()):
            source.append('complete_observed')
        check(json.loads(row['source_json']) == sorted(source), 'nomination uses wrong source types')
        status = ('rejected_nonprimitive' if period(template) != UNIT else
                  'accepted' if len(support) >= minimum else 'insufficient_independent_support')
        check(row['status'] == status, 'terminal acceptance disagrees with independent threshold')
        if status == 'accepted':
            accepted[sample].add(identifier)
        if template == ref_terminal:
            common[sample] = dict(support_molecules=len(support), platform_counts=dict(sorted(platforms.items())))
        for key, item in assessed.items():
            if key[0] != sample:
                continue
            eid = identifier, sample, key[1]
            expected_evidence.add(eid)
            check(eid in by_evidence, 'molecule evidence missing')
            witness = by_evidence[eid]
            check(witness['decision'] == ('support' if key in support else 'uninformative'), 'incorrect molecule support')
            check(set(json.loads(witness['path_ids_json'])) == {path['id'] for path in item['paths']},
                  'molecule path provenance incomplete')
            check(truth(witness['boundary_evidence_valid']) == item['valid'], 'invalid boundary validation flag')
            if item['valid']:
                check((int(witness['qL']), int(witness['qR'])) == (item['qL'], item['qR']), 'changed raw query interval')
            if item['terminal']:
                check(witness['terminal_bases'] == item['terminal'] and int(witness['inner_distance']) == item['inner']
                      and truth(witness['outer_exact']) == item['outer'], 'incorrect terminal/context assessment')
            if key in support:
                check(set(json.loads(witness['platforms_json'])) == item['platforms'], 'wrong supporting platform provenance')
                check(witness['validation_reason'] == 'exact_terminal_with_context', 'support lacks exact context reason')
    check(set(by_evidence) == expected_evidence, 'extra molecule evidence or candidate self-vote')
    check({sample: sorted(ids) for sample, ids in accepted.items()} == terminal['accepted_by_sample'], 'accepted dictionary mismatch')
    check(len(nominations) == terminal['proposed_template_count'] and
          sum(map(len, accepted.values())) == terminal['accepted_template_sample_count'], 'manifest terminal totals mismatch')
    check({sample: item['support_molecules'] for sample, item in common.items()} ==
          {'LCL5': 14, 'LCL6': 16, 'LCL7': 17, 'LCL8': 12}, 'common terminal support differs from diagnostic evidence')
    return dict(common_terminal_support=common, common_terminal_total=sum(x['support_molecules'] for x in common.values()),
                proposed_templates=len(nominations), accepted_template_samples=sum(map(len, accepted.values())),
                tested_molecules=len(assessed), evidence_rows=len(evidence)), by_template_sample


def compare_paths(records, out, baseline):
    new = indexed(rows(out/'unit_summary.tsv'), 'path_id', 'new summaries')
    old = indexed(rows(baseline/'unit_summary.tsv'), 'path_id', 'old summaries')
    grouped = []
    for directory in (baseline, out):
        groups = defaultdict(list)
        for row in rows(directory/'primitive_units.tsv'):
            if row['decomposition'] == 'primary':
                groups[row['path_id']].append(row)
        for values in groups.values():
            values.sort(key=lambda row: int(row['unit_index']))
        grouped.append(groups)
    fields = ('start0', 'end0', 'status', 'family', 'orientation', 'phase_offset', 'template_bases',
              'exact_bases', 'edit_operations')
    comparison = []
    for record in records:
        pid, sequence = record['id'], record['sequence']
        check(old[pid]['sequence_sha256'] == new[pid]['sequence_sha256'] == record['sequence_sha256'], 'source sequence changed')
        check(new[pid]['primitive_cn'] == 'NA', 'CN was released')
        if record['kind'] == 'observed':
            continue
        before, after = grouped[0][pid], grouped[1][pid]
        prefix = [[tuple(row[field] for field in fields) for row in values if int(row['end0']) <= len(sequence)-32]
                  for values in (before, after)]
        check(prefix[0] == prefix[1], 'annotation before final 32 bases changed')
        check(int(old[pid]['unknown_bp']) == 15 and int(new[pid]['unknown_bp']) == 0, 'expected terminal unknown recovery failed')
        check(int(new[pid]['full_16bp_units']) == int(old[pid]['full_16bp_units'])+1 and
              int(new[pid]['partial_units']) == 0, 'expected full-unit recovery failed')
        tail = after[-2:]
        check([(int(row['start0']), int(row['end0'])) for row in tail] ==
              [(len(sequence)-32, len(sequence)-16), (len(sequence)-16, len(sequence))], 'wrong final unit intervals')
        check(all(row['status'] == 'full' and row['edit_operations'] == '16=' for row in tail), 'tail not two exact full units')
        comparison.append(dict(path_id=pid, kind=record['kind'], sequence_bp=len(sequence),
            old_full_units=old[pid]['full_16bp_units'], new_full_units=new[pid]['full_16bp_units'],
            old_unknown_bp=old[pid]['unknown_bp'], new_unknown_bp=new[pid]['unknown_bp'],
            old_partial_units=old[pid]['partial_units'], new_partial_units=new[pid]['partial_units'],
            old_score=old[pid]['primary_score'], new_score=new[pid]['primary_score'],
            new_penultimate_start0=tail[0]['start0'], new_last_start0=tail[1]['start0'], new_last_end0=tail[1]['end0'],
            prefix_unchanged=True, sequence_unchanged=True, primitive_cn=new[pid]['primitive_cn'],
            annotation_ambiguity=new[pid]['annotation_ambiguity']))
    check(len(comparison) == 9 and sum(row['kind'] == 'candidate' for row in comparison) == 8, 'comparison must cover eight candidates and reference')
    c10 = next(row for row in comparison if row['path_id'] == 'CAND_LCL5_C10')
    check(int(c10['sequence_bp']) == 432 and int(c10['old_full_units']) == 26 and int(c10['new_full_units']) == 27
          and math.isclose(float(c10['old_score']), 40.) and math.isclose(float(c10['new_score']), 9.9), 'C10 diagnostic expectation failed')
    return comparison, c10


def main(out, baseline, input_file, pointgraph):
    validate_analysis(out, input_file)
    validate_analysis(baseline, input_file)
    records = json.loads(input_file.read_text(encoding='utf-8'))
    check(len(records) == 199, 'this comparator expects the frozen 199-path D1S80 input')
    manifest = json.loads((out/'manifest.json').read_text(encoding='utf-8'))
    previous = json.loads((baseline/'manifest.json').read_text(encoding='utf-8'))
    check('terminal_units' not in previous and manifest['terminal_units']['enabled'] is True, 'wrong before/after run modes')
    check(all(manifest['source_inputs'].get(name) == value for name, value in previous['source_inputs'].items()),
          'baseline source provenance differs')
    unchanged_joint = {}
    for name in JOINT_FILES:
        check(digest(out/name) == digest(baseline/name), f'joint output changed: {name}')
        unchanged_joint[name] = digest(out/name)
    evidence, candidates = validate_terminal_evidence(records, manifest, out, pointgraph)
    for row in rows(out/'primitive_units.tsv'):
        identifier = row['terminal_template_id']
        if identifier == 'NA':
            continue
        source = next(record for record in records if record['id'] == row['path_id'])
        check((source['kind'] != 'observed' or source['coverage'] == 'two_flank_anchored')
              and int(row['end0']) == len(source['sequence']),
              'terminal match outside complete right boundary')
        matching = [candidate for (template, sample), candidate in candidates.items()
                    if template == row['template_bases'] and candidate['terminal_id'] == identifier
                    and candidate['status'] == 'accepted' and (source['kind'] == 'reference' or source['sample'] == sample)]
        check(bool(matching), 'terminal match lacks accepted own-sample template')
    comparison, c10 = compare_paths(records, out, baseline)
    with (out/'comparison.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparison[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(comparison)
    result = dict(status='PASS', validator='independent_fasta_and_plain_levenshtein', input_paths=len(records),
                  candidate_count=8, unchanged_joint_file_hashes=unchanged_joint,
                  input_paths_sha256=digest(input_file), all_cn_withheld=True,
                  catalog_and_reference_identity_unchanged=True, original_sequences_preserved=True,
                  candidate_and_reference_prefixes_preserved=True, **evidence, LCL5_C10=c10)
    (out/'verification.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    check(len(sys.argv) == 5, 'Usage: validate_vntrterminal_run.py NEW_OUT BASELINE_OUT INPUT_PATHS POINTGRAPH_DIR')
    main(*(Path(arg) for arg in sys.argv[1:]))
