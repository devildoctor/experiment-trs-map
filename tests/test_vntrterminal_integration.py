"""Terminal evidence is opt-in, sample-specific and never publishes CN."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrdbg_analysis import primitive_units
from test_vntrdbg_units import MOTIF, motif_rows, record

TERMINAL = 'ACGATTCGAGCTAGCA'


class TerminalAnnotationTests(unittest.TestCase):
    def test_sample_scope_local_exclusion_and_cn_withholding(self):
        sequence = MOTIF*2+TERMINAL
        sources = [record('REF', sequence, 'reference', 'catalog_interval'),
                   record('CAND_S', sequence, 'candidate', 'candidate_sequence'),
                   record('CAND_T', sequence, 'candidate', 'candidate_sequence'),
                   record('OBS_full', sequence),
                   record('OBS_local', sequence, coverage='alignment_local_only')]
        sources[2]['sample'] = 'T'
        old_units, old, _ = primitive_units(sources, motif_rows(MOTIF), [])
        units, summaries, _ = primitive_units(sources, motif_rows(MOTIF), [], {'S': {TERMINAL: 'R1'}})
        for source, before, after in zip(sources, old, summaries):
            self.assertEqual(after['sequence_sha256'], before['sequence_sha256'])
            self.assertEqual(after['primitive_cn'], 'NA')
            self.assertEqual(after['cn_policy'], 'withheld_terminal_validation')
            selected = [r for r in units if r['path_id'] == source['id'] and r['decomposition'] == 'primary']
            self.assertEqual(''.join(r['exact_bases'] for r in selected), sequence)
            if source['id'] in ('CAND_T', 'OBS_local'):
                self.assertEqual(after['available_terminal_ids_json'], [])
                self.assertEqual(after['unknown_bp'], before['unknown_bp'])
            else:
                self.assertEqual(after['unknown_bp'], 0)
                self.assertEqual(after['annotated_full_units'], 3)
                self.assertEqual(selected[-1]['terminal_template_id'], 'R1')
                self.assertEqual(selected[-1]['end0'], len(sequence))

    def test_unique_count_is_still_withheld_in_opt_in_mode(self):
        source = record('CAND_unique', MOTIF*2, 'candidate', 'candidate_sequence')
        _, baseline, _ = primitive_units([source], motif_rows(MOTIF), [])
        _, opt_in, _ = primitive_units([source], motif_rows(MOTIF), [], {})
        self.assertEqual(baseline[0]['primitive_cn'], 2)
        self.assertEqual(opt_in[0]['primitive_cn'], 'NA')
        self.assertEqual(opt_in[0]['annotation_count_status'], 'resolved_count')
        self.assertEqual(opt_in[0]['primitive_cn_status'], 'withheld')


if __name__ == '__main__':
    unittest.main()
