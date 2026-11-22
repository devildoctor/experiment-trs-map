"""Synthetic contracts for optional, oriented right-terminal templates."""
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from vntrcore import canonical_motif, decompose, revcomp


MOTIF = 'CCTGGTGGTCCTCCTT'
TERMINAL = 'ACGTTACCGATGGCAT'


def serialized(result):
    tokens, explained, audit = result
    return [list(map(asdict, tokens)), explained, audit]


class RightTerminalTests(unittest.TestCase):
    def decompose_terminal(self, sequence, motifs=(MOTIF,)):
        result = decompose(sequence, motifs, return_details=True,
                           right_terminal_motifs=(TERMINAL,))
        tokens, _, audit = result
        for rows in ([asdict(token) for token in tokens], audit['alternative']):
            if not rows:
                continue
            self.assertEqual(''.join(row['exact_bases'] for row in rows), sequence)
            for row in rows:
                self.assertEqual(row['exact_bases'], sequence[row['start']:row['end']])
                if row['family'] == canonical_motif(TERMINAL):
                    self.assertEqual(row['end'], len(sequence))
        return result

    def test_default_results_match_prechange_snapshot(self):
        # Complete tokens and both-path audits from base commit 2081705.
        sequences = ['', MOTIF, MOTIF * 2, MOTIF + MOTIF[:15],
                     MOTIF[:8] + 'A' + MOTIF[8:], 'N' + MOTIF + 'NN', 'CAGCAGCAG']
        # 逐字段审计确认：仅第 6 个样例的两个备选 token 歧义标记由 False 修正为 True。
        corrected_ambiguity = {(5, 0), (5, 1)}
        outputs = []
        for case_index, sequence in enumerate(sequences):
            result = decompose(sequence, [MOTIF, 'CAG'], return_details=True)
            self.assertEqual(result, decompose(sequence, [MOTIF, 'CAG'],
                             return_details=True, right_terminal_motifs=()))
            snapshot = serialized(result)
            for token_index, token in enumerate(snapshot[2]['alternative']):
                self.assertIs(token['annotation_ambiguity'],
                              (case_index, token_index) in corrected_ambiguity)
            # 单独验证修正值后，仅将该字段投影回旧值，保留其余全部字段的旧快照契约。
            snapshot[2] = dict(snapshot[2], alternative=[dict(token, annotation_ambiguity=False)
                               for token in snapshot[2]['alternative']])
            outputs.append(snapshot)
        digest = sha256(json.dumps(outputs, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        self.assertEqual(digest, '0b9bcd108bea5ebff1510b199404a6cf6b4cdb41302e73fc8094841d48566f26')

    def test_missing_terminal_is_explained_at_right_end(self):
        sequence = MOTIF + TERMINAL
        original = decompose(sequence, [MOTIF], return_details=True)
        self.assertGreater(original[2]['metrics']['unknown_bases'], 0)
        tokens, _, audit = self.decompose_terminal(sequence)
        self.assertEqual([(token.start, token.end, token.edit_operations) for token in tokens],
                         [(0, 16, '16='), (16, 32, '16=')])
        self.assertEqual(audit['metrics']['full_units'], 2)
        self.assertEqual(audit['metrics']['unknown_bases'], 0)

    def test_internal_and_left_partial_matches_are_not_offered(self):
        for sequence in (TERMINAL + MOTIF, TERMINAL[-8:] + MOTIF):
            with self.subTest(sequence=sequence):
                result = self.decompose_terminal(sequence)
                self.assertEqual(result, decompose(sequence, [MOTIF], return_details=True))

    def test_phase_and_orientation_are_not_expanded(self):
        for sequence in (TERMINAL[5:] + TERMINAL[:5], revcomp(TERMINAL)):
            with self.subTest(sequence=sequence):
                _, _, audit = self.decompose_terminal(sequence, motifs=())
                self.assertEqual(audit['metrics']['full_units'], 0)
                self.assertGreater(audit['metrics']['unknown_bases'], 0)
        _, _, audit = self.decompose_terminal(TERMINAL, motifs=())
        self.assertEqual(audit['metrics']['full_units'], 1)
        self.assertEqual(audit['metrics']['exact_full_bases'], 16)

    def test_true_fifteen_base_prefix_remains_partial(self):
        tokens, _, audit = self.decompose_terminal(MOTIF + TERMINAL[:15])
        self.assertEqual((tokens[-1].start, tokens[-1].end, tokens[-1].status),
                         (16, 31, 'partial'))
        self.assertEqual(audit['metrics']['full_units'], 1)
        self.assertEqual(audit['metrics']['partial_bases'], 15)

    def test_seventeen_base_insertion_is_preserved(self):
        inserted = TERMINAL[:8] + 'A' + TERMINAL[8:]
        tokens, _, audit = self.decompose_terminal(MOTIF + inserted)
        self.assertEqual((tokens[-1].start, tokens[-1].end, tokens[-1].status),
                         (16, 33, 'full'))
        self.assertEqual(tokens[-1].exact_bases, inserted)
        self.assertEqual(tokens[-1].edit_operations, '8=1I8=')
        self.assertEqual(audit['metrics']['edited_full_bases'], 17)

    def test_terminal_template_does_not_change_unknown_costs(self):
        for sequence in ('A', 'TNNN', 'NNNA'):
            with self.subTest(sequence=sequence):
                result = self.decompose_terminal(sequence, motifs=())
                self.assertEqual(result, decompose(sequence, [], return_details=True))
                self.assertEqual(result[2]['primary_score'], len(sequence) * 2.)

    def test_ordinary_duplicate_retains_internal_scope(self):
        sequence = TERMINAL + MOTIF
        original = decompose(sequence, [MOTIF, TERMINAL], return_details=True)
        with_terminal = decompose(sequence, [MOTIF, TERMINAL], return_details=True,
                                  right_terminal_motifs=(TERMINAL,))
        self.assertEqual(original, with_terminal)
        self.assertEqual(with_terminal[2]['metrics']['full_units'], 2)


if __name__ == '__main__':
    unittest.main()
