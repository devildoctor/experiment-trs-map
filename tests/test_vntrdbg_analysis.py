import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from vntrdbg_analysis import branch_queries, joint_support, occurrences, primitive_units
from vntrdbg_analysis import define_events, anchored_interval


def item(identifier, sequence, kind='observed', key=None, coverage='two_flank_anchored'):
    result = dict(id=identifier, sequence=sequence, kind=kind, sample='S',
                  coverage=coverage, platforms=['HiFi'])
    if kind == 'observed':
        result['molecule_key'] = key or identifier
    return result


class JointSupportTests(unittest.TestCase):
    def test_exact_interval_requires_one_contiguous_observation(self):
        records = [item('CAND_C', 'AAACCCGGGTTT', kind='candidate'),
                   item('o1', 'AAACCC', key='split', coverage='alignment_local_only'),
                   item('o2', 'GGGTTT', key='split', coverage='alignment_local_only'),
                   item('o3', 'AAACCCGGGTTT', key='whole'),
                   item('o4', 'AAACCCGGGTTT', key='whole')]
        query = dict(query_id='C.k3.J001', candidate_id='C', sample='S', kind='adjacent_branch_pair',
                     exact_bases='CCCGGG', start0=3, end0=9, k=3,
                     candidate_occurrences=1, branch_edge_count=2)
        summary, evidence = joint_support(records, [query])
        self.assertEqual(summary[0]['tested_molecules'], 2)
        self.assertEqual(summary[0]['exact_support_molecules'], 1)
        self.assertEqual(next(row for row in evidence if row['molecule_key']=='split')['decision'],
                         'uninformative')
        self.assertEqual(next(row for row in evidence if row['molecule_key']=='whole')['occurrence_count'], 2)

    def test_full_candidate_rejects_local_and_longer_observation(self):
        records = [item('CAND_C', 'AACCGG', kind='candidate'),
                   item('a', 'AACCGG', key='local', coverage='alignment_local_only'),
                   item('b', 'TAACCGG', key='long'), item('c', 'AACCGG', key='exact')]
        query = dict(query_id='F', candidate_id='C', sample='S', kind='full_candidate_exact',
                     exact_bases='AACCGG', start0=0, end0=6, k=3,
                     candidate_occurrences=1, branch_edge_count=2)
        summary, _ = joint_support(records, [query])
        self.assertEqual(summary[0]['exact_support_molecules'], 1)
        self.assertEqual(summary[0]['conflict_molecules'], 0)
        self.assertEqual(summary[0]['support_molecules'], 1)
        self.assertEqual(summary[0]['uninformative_molecules'], 2)

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

    def test_merged_branch_counts_real_positions_not_span(self):
        from vntrdbg import raw_graph
        records = [item('CAND_C1', 'AAACCCGATCGTACGATCGGGTTT', kind='candidate'),
                   item('CAND_C2', 'AAAGCCGATCGTACGATCGGATTT', kind='candidate')]
        raw = raw_graph(records, 3)
        edges = raw['edges']
        query_rows = branch_queries(records, 3)
        saw_gap = False
        for q in query_rows:
            path = raw['paths']['CAND_'+q['candidate_id']]
            expected = [i for i in range(q['start0'], max(q['start0'], q['end0']-3))
                        if sum(a==path[i] for a,b in edges)>1
                        or sum(b==path[i+1] for a,b in edges)>1]
            self.assertEqual(q['branch_positions_json'], expected)
            self.assertEqual(q['branch_edge_count'], len(expected))
            saw_gap |= bool(expected and len(expected)<expected[-1]-expected[0]+1)
        self.assertTrue(saw_gap)

    def test_explicit_competitors_mixed_molecule_and_sample_boundary(self):
        target, alt = 'AAACCCGGG', 'AAATTTGGG'
        records = [item('CAND_C', target, kind='candidate'),
                   item('CAND_A', alt, kind='candidate'),
                   item('CAND_EQ', target, kind='candidate'),
                   item('own', target), item('alt', alt), item('error', 'AAACCTGGG'),
                   item('mix1', target, key='mixed'), item('mix2', alt, key='mixed'),
                   item('foreign', 'AAAGGGCCC', kind='candidate')]
        records[-1]['sample'] = 'OTHER'
        records.append(item('foreign_observation', 'AAAGGGCCC'))
        query = dict(query_id='F', candidate_id='C', sample='S', kind='full_candidate_exact',
                     exact_bases=target, start0=0, end0=9, k=3,
                     candidate_occurrences=1, branch_edge_count=2)
        events = define_events(records, [query])
        self.assertEqual(len(events[0]['states_json']), 2)
        summary, evidence = joint_support(records, [query], events)
        self.assertEqual([summary[0][k] for k in ('support_molecules','conflict_molecules',
                                                'uninformative_molecules')], [1, 1, 3])
        self.assertEqual(summary[0]['support_fraction'], .5)
        self.assertEqual(next(r for r in evidence if r['molecule_key']=='mixed')['reason'],
                         'mixed_states_in_one_molecule')

    def test_interval_alternative_needs_unique_ordered_context(self):
        target, alt = 'AAACCCGGG', 'AAATTTGGG'
        records = [item('CAND_C', target, kind='candidate'), item('CAND_A', alt, kind='candidate'),
                   item('target', target, coverage='alignment_local_only'),
                   item('alternative', alt, coverage='alignment_local_only'),
                   item('repeat', target+target, coverage='alignment_local_only'),
                   item('reversed', 'GGGCCCAAA', coverage='alignment_local_only'),
                   item('unknown', 'AAACCTGGG', coverage='alignment_local_only')]
        query = dict(query_id='C.k3.B001', candidate_id='C', sample='S', kind='branch_block',
                     exact_bases=target, start0=0, end0=9, k=3,
                     candidate_occurrences=1, branch_edge_count=2)
        summary, evidence = joint_support(records, [query])
        self.assertEqual([summary[0][k] for k in ('support_molecules','conflict_molecules',
                                                'uninformative_molecules')], [1, 1, 3])
        repeated = next(r for r in evidence if r['molecule_key']=='repeat')
        self.assertEqual(repeated['occurrence_count'], 2)
        self.assertEqual(repeated['decision'], 'uninformative')
        self.assertIsNone(anchored_interval('GGGCCCAAA', 'AAA', 'GGG'))

    def test_repeated_candidate_context_remains_unresolved(self):
        target = 'AAACCCGGG'
        records = [item('CAND_C', target+target, kind='candidate'),
                   item('local', target, coverage='alignment_local_only')]
        query = dict(query_id='C.k3.B001', candidate_id='C', sample='S', kind='branch_block',
                     exact_bases=target, start0=0, end0=9, k=3,
                     candidate_occurrences=2, branch_edge_count=2)
        summary, _ = joint_support(records, [query])
        self.assertEqual(summary[0]['exact_support_molecules'], 1)
        self.assertEqual(summary[0]['informative_molecules'], 0)
        self.assertEqual(summary[0]['support_fraction'], 'NA')


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
