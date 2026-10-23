import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from vntrdbg_analysis import branch_queries, joint_support, occurrences, primitive_units


def item(identifier, sequence, kind='observed', key=None, coverage='two_flank_anchored'):
    result = dict(id=identifier, sequence=sequence, kind=kind, sample='S',
                  coverage=coverage, platforms=['HiFi'])
    if kind == 'observed':
        result['molecule_key'] = key or identifier
    return result


class JointSupportTests(unittest.TestCase):
    def test_exact_interval_requires_one_contiguous_observation(self):
        records = [item('o1', 'AAACCC', key='split', coverage='alignment_local_only'),
                   item('o2', 'GGGTTT', key='split', coverage='alignment_local_only'),
                   item('o3', 'AAACCCGGGTTT', key='whole'),
                   item('o4', 'AAACCCGGGTTT', key='whole')]
        query = dict(query_id='Q', candidate_id='C', sample='S', kind='adjacent_branch_pair',
                     exact_bases='CCCGGG', start0=3, end0=9, k=3,
                     candidate_occurrences=1, branch_edge_count=2)
        summary, evidence = joint_support(records, [query])
        self.assertEqual(summary[0]['tested_molecules'], 2)
        self.assertEqual(summary[0]['exact_support_molecules'], 1)
        self.assertEqual(next(row for row in evidence if row['molecule_key']=='split')['decision'],
                         'uninformative')
        self.assertEqual(next(row for row in evidence if row['molecule_key']=='whole')['occurrence_count'], 2)

    def test_full_candidate_rejects_local_and_longer_observation(self):
        records = [item('a', 'AACCGG', key='local', coverage='alignment_local_only'),
                   item('b', 'TAACCGG', key='long'), item('c', 'AACCGG', key='exact')]
        query = dict(query_id='F', candidate_id='C', sample='S', kind='full_candidate_exact',
                     exact_bases='AACCGG', start0=0, end0=6, k=3,
                     candidate_occurrences=1, branch_edge_count=2)
        summary, _ = joint_support(records, [query])
        self.assertEqual(summary[0]['exact_support_molecules'], 1)
        self.assertEqual(summary[0]['conflict_molecules'], 'NA')

    def test_candidate_derived_branch_pairs_are_ordered(self):
        records = [item('CAND_C1', 'AAACCCGATCGTACGATCGGGTTT', kind='candidate'),
                   item('CAND_C2', 'AAAGCCGATCGTACGATCGGATTT', kind='candidate')]
        queries = branch_queries(records, 3, 1)
        pairs = [row for row in queries if row['candidate_id']=='C1'
                 and row['kind']=='adjacent_branch_pair']
        self.assertTrue(pairs)
        for row in pairs:
            self.assertEqual(row['exact_bases'], records[0]['sequence'][row['start0']:row['end0']])
        self.assertEqual(occurrences('AAAA', 'AA'), [0, 1, 2])


class PrimitiveUnitTests(unittest.TestCase):
    def test_48_bp_template_is_three_16_bp_units(self):
        motif = 'CCTGGTGGTCCTCCTT'
        sequence = motif * 3
        candidates = [item('CAND_C1', sequence, kind='candidate')]
        motifs = [dict(motif=motif, usable='True', level='base'),
                  dict(motif=sequence, usable='True', level='base')]
        legacy = [dict(allele='C1', start='0', end='48', status='full')]
        units, summary, counts = primitive_units(candidates, motifs, legacy)
        primary = [row for row in units if row['decomposition']=='primary']
        self.assertEqual(''.join(row['exact_bases'] for row in primary), sequence)
        self.assertEqual(len(primary), 3)
        self.assertEqual(summary[0]['full_16bp_units'], 3)
        self.assertEqual(summary[0]['legacy_full_tokens'], 1)
        self.assertEqual(sum(row['segment_count'] for row in counts), 3)
        self.assertTrue(all(row['legacy_token_indexes_json']==[1] for row in primary))


if __name__ == '__main__':
    unittest.main()
