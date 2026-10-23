"""Sequence and scope contracts for D1S80 primitive-unit annotations."""
from hashlib import sha256
from pathlib import Path
import re
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from vntrdbg_analysis import primitive_units


MOTIF = 'CCTGGTGGTCCTCCTT'


def motif_rows(*motifs):
    return [dict(motif=motif, usable='True', level='base') for motif in motifs]


def record(identifier, sequence, kind='observed', coverage='two_flank_anchored'):
    result = dict(id=identifier, sequence=sequence, kind=kind, sample='S',
                  coverage=coverage, platforms=['HiFi'],
                  sequence_sha256=sha256(sequence.encode()).hexdigest())
    if kind == 'observed':
        result['molecule_key'] = 'molecule_' + identifier
    return result


class PrimitivePathTests(unittest.TestCase):
    def assert_partition(self, records, units):
        for source in records:
            source_units = [unit for unit in units if unit['path_id'] == source['id']]
            self.assertTrue(source_units)
            for scheme in {unit['decomposition'] for unit in source_units}:
                selected = [unit for unit in source_units if unit['decomposition'] == scheme]
                self.assertEqual([unit['unit_index'] for unit in selected],
                                 list(range(1, len(selected) + 1)))
                position = 0
                for unit in selected:
                    self.assertEqual(unit['start0'], position)
                    self.assertEqual(unit['exact_bases'],
                                     source['sequence'][unit['start0']:unit['end0']])
                    self.assertEqual(unit['observed_unit_bp'], unit['end0'] - unit['start0'])
                    position = unit['end0']
                self.assertEqual(position, len(source['sequence']))
                self.assertEqual(''.join(unit['exact_bases'] for unit in selected), source['sequence'])

    def test_all_path_kinds_and_short_fragments_preserve_source(self):
        records = [record('REF', MOTIF * 2, 'reference', 'catalog_interval'),
                   record('CAND_C1', MOTIF * 2, 'candidate', 'candidate_sequence'),
                   record('OBS_full', MOTIF * 2),
                   record('OBS_18', MOTIF + MOTIF[:2], coverage='alignment_local_only'),
                   record('OBS_22', MOTIF + MOTIF[:6], coverage='alignment_local_only')]
        legacy = [dict(allele='C1', start='0', end='32', status='full')]
        units, summaries, counts = primitive_units(records, motif_rows(MOTIF), legacy)
        self.assertEqual({row['path_id'] for row in summaries}, {row['id'] for row in records})
        self.assertEqual(len(summaries), len(records))
        self.assertEqual({row['path_id'] for row in counts}, {row['id'] for row in records})
        self.assert_partition(records, units)
        for source, summary in zip(records, summaries):
            self.assertEqual(summary['kind'], source['kind'])
            self.assertEqual(summary['coverage'], source['coverage'])
            self.assertEqual(summary['sequence_sha256'], source['sequence_sha256'])
            self.assertEqual(summary['platforms_json'], ['HiFi'])
            self.assertEqual(summary['molecule_key'], source.get('molecule_key', 'NA'))
            selected = [row for row in units if row['path_id'] == source['id']]
            if source['kind'] == 'candidate':
                self.assertEqual(summary['candidate_id'], 'C1')
                self.assertEqual(summary['legacy_full_tokens'], 1)
                self.assertTrue(all(row['legacy_token_indexes_json'] == [1] for row in selected))
            else:
                self.assertEqual(summary['candidate_id'], 'NA')
                self.assertEqual(summary['legacy_full_tokens'], 'NA')
                self.assertTrue(all(row['legacy_token_indexes_json'] == 'NA' for row in selected))
            if source['coverage'] == 'alignment_local_only':
                self.assertEqual(summary['count_scope'], 'observed_fragment')
                self.assertEqual(summary['primitive_cn'], 'NA')
                self.assertEqual(summary['primitive_cn_status'], 'fragment_only')

    def test_complete_fragment_does_not_become_locus_copy_number(self):
        records = [record('OBS_fragment', MOTIF, coverage='alignment_local_only')]
        _, summaries, _ = primitive_units(records, motif_rows(MOTIF), [])
        summary = summaries[0]
        self.assertEqual(summary['full_16bp_units'], 1)
        self.assertEqual(summary['annotation_count_status'], 'resolved_count')
        self.assertEqual(summary['primitive_cn'], 'NA')
        self.assertEqual(summary['primitive_cn_status'], 'fragment_only')

    def test_second_best_alternative_is_not_a_tied_annotation(self):
        records = [record('OBS_unique', MOTIF)]
        units, summaries, _ = primitive_units(records, motif_rows(MOTIF), [])
        summary = summaries[0]
        self.assertTrue(any(unit['decomposition'] == 'alternative' for unit in units))
        self.assertGreater(summary['score_gap'], 0)
        self.assertFalse(summary['alternative_is_tied'])
        self.assertFalse(summary['annotation_ambiguity'])
        self.assertEqual(summary['alternative_limit'], 1)
        self.assertEqual(summary['phase_policy'], 'fixed_template_phase')

    def test_equal_count_tie_stays_ambiguous_in_both_schemes(self):
        records = [record('OBS_tie', 'G' + MOTIF[1:])]
        units, summaries, _ = primitive_units(records, motif_rows(MOTIF, 'A' + MOTIF[1:]), [])
        summary = summaries[0]
        self.assertTrue(summary['alternative_is_tied'])
        self.assertEqual(summary['full_16bp_units'], summary['alternative_full_units'])
        self.assertEqual(summary['primitive_cn'], 'NA')
        self.assertEqual(summary['primitive_cn_status'], 'ambiguous')
        self.assertTrue(all(unit['annotation_ambiguity'] for unit in units))
        self.assert_partition(records, units)

    def test_substitution_and_indel_cigars_reconstruct_unit_slices(self):
        examples = [('X', 'A' + MOTIF[1:]),
                    ('I', MOTIF[:8] + 'A' + MOTIF[8:]),
                    ('D', MOTIF[:8] + MOTIF[9:])]
        for expected, sequence in examples:
            with self.subTest(edit=expected):
                records = [record('OBS_edit', sequence)]
                units, _, _ = primitive_units(records, motif_rows(MOTIF), [])
                primary = [unit for unit in units if unit['decomposition'] == 'primary']
                self.assertTrue(any(expected in unit['edit_operations'] for unit in primary))
                self.assert_partition(records, units)
                for unit in units:
                    if unit['status'] == 'unknown':
                        continue
                    operations = re.findall(r'(\d+)([=XID])', unit['edit_operations'])
                    self.assertEqual(''.join(n + op for n, op in operations), unit['edit_operations'])
                    template, exact = unit['template_bases'], unit['exact_bases']
                    template_pos = exact_pos = 0
                    for raw_length, operation in operations:
                        length = int(raw_length)
                        if operation in '=X':
                            for offset in range(length):
                                equal = template[template_pos + offset] == exact[exact_pos + offset]
                                self.assertEqual(equal, operation == '=')
                        if operation in '=XD':
                            template_pos += length
                        if operation in '=XI':
                            exact_pos += length
                    self.assertEqual(template_pos, len(template))
                    self.assertEqual(exact_pos, len(exact))

    def test_harmonic_templates_and_incorrect_hash_are_rejected(self):
        with self.assertRaisesRegex(ValueError, '16-bp'):
            primitive_units([record('OBS', MOTIF)], motif_rows('ACGT' * 4), [])
        source = record('OBS', MOTIF)
        source['sequence_sha256'] = 'incorrect'
        with self.assertRaisesRegex(ValueError, 'hash'):
            primitive_units([source], motif_rows(MOTIF), [])


if __name__ == '__main__':
    unittest.main()
