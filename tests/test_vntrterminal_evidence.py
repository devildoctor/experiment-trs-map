"""Independent-molecule and provenance contracts for right-terminal templates."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from vntrterminal import validate_terminals


MOTIF = 'CCTGGTGGTCCTCCTT'
TERMINAL = 'AGTCGATACGTTAGCA'
RIGHT = 'TACGACCTAGGATCGA'
LEFT = 'AACCGGTT'


def observation(sample, name, terminal=TERMINAL, inner=MOTIF, outer=RIGHT):
    sequence = MOTIF + inner + terminal
    key = json.dumps([sample, 'library', name], separators=(',', ':'))
    sources = [f'synthetic/{sample}/{name}']
    interval = [len(LEFT), len(LEFT) + len(sequence)]
    record = dict(id=f'OBS_{sample}_{name}', kind='observed', sample=sample,
                  molecule_key=key, sequence=sequence,
                  sequence_sha256=sha256(sequence.encode()).hexdigest(),
                  coverage='two_flank_anchored', anchor_status='complete_support',
                  stored_query_interval=interval, original_query_interval=interval,
                  source_records=sources, platforms=['HiFi'])
    molecule = dict(sample=sample, molecule_key=key, sequence=sequence,
                    query_sequence=LEFT + sequence + outer,
                    qL=str(interval[0]), qR=str(interval[1]),
                    original_query_interval=json.dumps(interval),
                    records=json.dumps(sources))
    decision = dict(sample=sample, molecule_key=key, status='complete_support',
                    sources=sources)
    return record, molecule, decision


def dataset(samples=None, terminal=TERMINAL):
    sequence = MOTIF * 2 + terminal
    records = [dict(id='REF', kind='reference', sample='', sequence=sequence,
                    sequence_sha256=sha256(sequence.encode()).hexdigest(),
                    coverage='catalog_interval', platforms=[])]
    molecules, decisions = [], []
    for sample, count in (samples or {'S': 3}).items():
        for index in range(count):
            record, molecule, decision = observation(sample, f'r{index}', terminal)
            records.append(record)
            molecules.append(molecule)
            decisions.append(decision)
    return records, molecules, decisions


def evaluate(data):
    return validate_terminals(*data, [MOTIF], RIGHT, 'synthetic', 100, 148)


def candidate(rows, sample='S', terminal=TERMINAL):
    return next(row for row in rows
                if row['sample_id'] == sample and row['template_bases'] == terminal)


class TerminalEvidenceTests(unittest.TestCase):
    def test_three_independent_molecules_accept_without_mutating_inputs(self):
        data = dataset()
        before = deepcopy(data)
        rows, evidence, accepted = evaluate(data)
        row = candidate(rows)
        self.assertEqual(row['support_molecules'], 3)
        self.assertEqual(row['platform_counts_json'], {'HiFi': 3})
        self.assertEqual(row['status'], 'accepted')
        self.assertEqual(accepted['S'][TERMINAL], row['terminal_id'])
        self.assertEqual(sum(item['decision'] == 'support' for item in evidence), 3)
        self.assertEqual(data, before)

    def test_duplicate_observations_records_and_decisions_do_not_add_votes(self):
        data = dataset({'S': 2})
        records, molecules, decisions = data
        for index in range(3):
            duplicate = deepcopy(records[1])
            duplicate['id'] += f'_duplicate_{index}'
            records.append(duplicate)
            molecules.append(deepcopy(molecules[0]))
            decisions.append(deepcopy(decisions[0]))
        rows, _, accepted = evaluate(data)
        self.assertEqual(candidate(rows)['support_molecules'], 2)
        self.assertEqual(accepted, {})

    def test_reference_and_candidate_never_vote_or_candidate_nominate(self):
        data = dataset({'S': 2})
        sequence = MOTIF * 2 + TERMINAL
        other_terminal = TERMINAL[1:] + TERMINAL[:1]
        for index, terminal in enumerate([TERMINAL, other_terminal]):
            data[0].append(dict(id=f'CAND_C{index}', kind='candidate', sample='S',
                                sequence=sequence[:-16] + terminal,
                                coverage='candidate_sequence', platforms=[]))
        rows, _, accepted = evaluate(data)
        self.assertEqual(candidate(rows)['support_molecules'], 2)
        self.assertNotIn(other_terminal, {row['template_bases'] for row in rows})
        self.assertEqual(accepted, {})

    def test_local_observation_neither_votes_nor_nominates(self):
        data = dataset({'S': 2})
        local, _, _ = observation('S', 'local')
        local['coverage'] = 'alignment_local_only'
        data[0].append(local)
        rows, evidence, accepted = evaluate(data)
        self.assertEqual(candidate(rows)['support_molecules'], 2)
        self.assertEqual(accepted, {})
        self.assertEqual(next(row for row in evidence
                              if row['molecule_key'] == local['molecule_key'])['validation_reason'],
                         'local_only')

    def test_support_threshold_is_per_sample(self):
        rows, _, accepted = evaluate(dataset({'S': 3, 'T': 2}))
        self.assertEqual(candidate(rows, 'S')['status'], 'accepted')
        self.assertEqual(candidate(rows, 'T')['support_molecules'], 2)
        self.assertEqual(set(accepted), {'S'})

    def test_source_slice_boundary_and_record_provenance_must_agree(self):
        for field in ('records', 'qR', 'query_sequence', 'original_query_interval'):
            with self.subTest(field=field):
                data = dataset()
                row = data[1][0]
                if field == 'records':
                    row[field] = json.dumps(['unrelated/source'])
                elif field == 'qR':
                    row[field] = str(int(row[field]) - 1)
                elif field == 'query_sequence':
                    row[field] = row[field][:8] + 'A' + row[field][9:]
                else:
                    row[field] = json.dumps([9, 57])
                rows, evidence, accepted = evaluate(data)
                self.assertEqual(candidate(rows)['support_molecules'], 2)
                self.assertEqual(accepted, {})
                self.assertEqual(next(item for item in evidence
                                      if item['molecule_key'] == row['molecule_key'])['decision'],
                                 'uninformative')

    def test_missing_conflicting_or_noncomplete_anchor_decisions_do_not_support(self):
        for status in (None, 'unresolved_molecule_conflict',
                       'conflicting_complete_observations',
                       'inconsistent_record_sequences_or_lengths', 'unknown_status'):
            with self.subTest(status=status):
                data = dataset()
                if status is None:
                    data[2].pop(0)
                else:
                    data[2][0]['status'] = status
                rows, _, accepted = evaluate(data)
                self.assertEqual(candidate(rows)['support_molecules'], 2)
                self.assertEqual(accepted, {})

    def test_anchor_decision_sample_and_sources_must_match(self):
        for field, value in [('sample', 'OTHER'), ('sources', ['unrelated/source'])]:
            with self.subTest(field=field):
                data = dataset()
                data[2][0][field] = value
                rows, _, accepted = evaluate(data)
                self.assertEqual(candidate(rows)['support_molecules'], 2)
                self.assertEqual(accepted, {})

    def test_absent_source_provenance_cannot_validate_a_molecule(self):
        data = dataset()
        data[0][1]['source_records'] = []
        data[1][0]['records'] = '[]'
        data[2][0]['sources'] = []
        rows, _, accepted = evaluate(data)
        self.assertEqual(candidate(rows)['support_molecules'], 2)
        self.assertEqual(accepted, {})

    def test_frozen_anchor_status_and_hash_must_agree(self):
        for field, value in [('anchor_status', 'missing_or_ambiguous_anchor'),
                             ('sequence_sha256', '0' * 64)]:
            with self.subTest(field=field):
                data = dataset()
                data[0][1][field] = value
                rows, _, accepted = evaluate(data)
                self.assertEqual(candidate(rows)['support_molecules'], 2)
                self.assertEqual(accepted, {})

    def test_exact_outer_context_and_original_inner_template_are_required(self):
        for kwargs in ({'outer': 'A' * 16}, {'inner': 'A' * 16},
                       {'inner': TERMINAL}, {'outer': RIGHT[:-1]}):
            with self.subTest(kwargs=kwargs):
                data = dataset({'S': 2})
                record, molecule, decision = observation('S', 'bad_context', **kwargs)
                data[0].append(record)
                data[1].append(molecule)
                data[2].append(decision)
                rows, _, accepted = evaluate(data)
                self.assertEqual(candidate(rows)['support_molecules'], 2)
                self.assertEqual(accepted, {})

    def test_cyclic_phase_states_do_not_pool_support(self):
        data = dataset({'S': 2})
        rotated = TERMINAL[1:] + TERMINAL[:1]
        record, molecule, decision = observation('S', 'rotated', rotated)
        data[0].append(record)
        data[1].append(molecule)
        data[2].append(decision)
        rows, _, accepted = evaluate(data)
        self.assertEqual(candidate(rows)['support_molecules'], 2)
        self.assertEqual(candidate(rows, terminal=rotated)['support_molecules'], 1)
        self.assertEqual(candidate(rows)['family'], candidate(rows, terminal=rotated)['family'])
        self.assertEqual(accepted, {})

    def test_conflicting_full_observations_exclude_entire_molecule(self):
        data = dataset()
        other = deepcopy(data[0][1])
        other['id'] += '_conflict'
        other['sequence'] = 'A' + other['sequence'][1:]
        other['sequence_sha256'] = sha256(other['sequence'].encode()).hexdigest()
        data[0].append(other)
        rows, _, accepted = evaluate(data)
        self.assertEqual(candidate(rows)['support_molecules'], 2)
        self.assertEqual(accepted, {})

    def test_local_metadata_cannot_supply_platform_support_for_full_observation(self):
        data = dataset()
        local = deepcopy(data[0][1])
        local.update(id='OBS_local_same_molecule', coverage='alignment_local_only',
                     platforms=['ONT'], source_records=['local/source'])
        data[0].append(local)
        rows, _, _ = evaluate(data)
        self.assertEqual(candidate(rows)['platform_counts_json'], {'HiFi': 3})

    def test_nonprimitive_terminal_remains_rejected(self):
        terminal = 'ACGT' * 4
        rows, _, accepted = evaluate(dataset(terminal=terminal))
        self.assertEqual(candidate(rows, terminal=terminal)['status'], 'rejected_nonprimitive')
        self.assertEqual(accepted, {})

    def test_ambiguous_terminal_bases_never_form_an_accepted_template(self):
        data = dataset({'S': 1})
        del data[0][1:]
        data[1].clear()
        data[2].clear()
        for index in range(3):
            record, molecule, decision = observation('S', f'ambiguous{index}',
                                                      'N' + TERMINAL[1:])
            data[0].append(record)
            data[1].append(molecule)
            data[2].append(decision)
        _, _, accepted = evaluate(data)
        self.assertEqual(accepted, {})

    def test_invalid_reference_span_and_support_threshold_raise(self):
        data = dataset()
        with self.assertRaises(ValueError):
            validate_terminals(*data, [MOTIF], RIGHT, 'synthetic', 100, 149)
        with self.assertRaises(ValueError):
            validate_terminals(*data, [MOTIF], RIGHT, 'synthetic', 100, 148, min_support=2)


if __name__ == '__main__':
    unittest.main()
