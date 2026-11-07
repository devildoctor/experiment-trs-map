import itertools
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))

from vntrflankalign import (affine_matrix, encode, pair_endpoints, locate_interval,
                           flank_search, score_flank, traceback, verify_alignment, cigar_segments, cigar_hint)


def exhaustive(pattern, text):
    """Independent enumeration, no DP recurrence or pruning by the global best."""
    # A full reference deletion can end at any free query prefix; it has no
    # matches and will fail the anchor threshold, but participates in the DP.
    answers = {j: (5+len(pattern), 0, 0) for j in range(1,len(text)+1)}
    for begin in range(len(text)):
        for end in range(begin+1, len(text)+1):
            sub = text[begin:end]
            def visit(i, j, last, cost, matches, first_query, last_query):
                if i == len(pattern) and j == len(sub):
                    if first_query != 'M' or last_query != 'M':
                        return
                    old = answers.get(end, (10**9, 10**9, -1))
                    if cost < old[0]:
                        answers[end] = (cost, matches, matches)
                    elif cost == old[0]:
                        answers[end] = (cost, min(matches, old[1]), max(matches, old[2]))
                    return
                if i < len(pattern) and j < len(sub):
                    equal = pattern[i] == sub[j]
                    visit(i+1, j+1, 'M', cost+(0 if equal else 2), matches+int(equal), first_query or 'M', 'M')
                if i < len(pattern):
                    visit(i+1, j, 'D', cost+(1 if last == 'D' else 6), matches, first_query, last_query)
                if j < len(sub):
                    visit(i, j+1, 'I', cost+(1 if last == 'I' else 6), matches, first_query or 'I', 'I')
            visit(0, 0, '', 0, 0, '', '')
    return answers


class AffineTests(unittest.TestCase):
    def test_all_endpoints_against_exhaustive(self):
        cases = 0
        patterns = [''.join(p) for n in range(1,4) for p in itertools.product('AC', repeat=n)]
        texts = [''.join(t) for n in range(1,5) for t in itertools.product('AC', repeat=n)]
        for p in patterns:
            for t in texts:
                costs, low, high, _ = affine_matrix(encode(p), encode(t))
                expected = exhaustive(p, t)
                for j in range(1, len(t)+1):
                    best = min(int(costs[s,-1,j]) for s in (0,1))
                    states = [s for s in (0,1) if costs[s,-1,j] == best]
                    got = best, min(int(low[s,-1,j]) for s in states), max(int(high[s,-1,j]) for s in states)
                    self.assertEqual(got, expected[j], (p,t,j))
                cases += 1
        self.assertEqual(cases, 420)

    def test_joint_suboptimal_endpoint_and_order_independence(self):
        left = [[20,0,250,250],[8,6,249,249]]
        right = [[12,0,250,250],[25,9,247,247]]
        got = pair_endpoints(left,right,[8,12])
        self.assertEqual((got[0]['qstart'],got[0]['qend'],got[0]['cost']),(8,12,6))
        self.assertEqual(got,pair_endpoints(left[::-1],right[::-1],[8,12]))

    def test_normalization_preserves_cost_ties(self):
        pairs=pair_endpoints([[1,6,247,247],[2,6,249,249]],[[9,0,250,250]],[1,9])
        self.assertEqual(len(pairs),2)
        self.assertEqual(pairs[0]['qstart'],2)

    def test_tied_paths_can_cross_match_threshold(self):
        pattern,query='CCACAACAAC','CACCCAAAACAAACA'
        search=dict(status='ready',side='left',reference_interval=[0,10],search_interval=[0,15])
        evidence,matrices=score_flank(pattern,query,search)
        endpoint=next(e for e in evidence['endpoints'] if e[0]==14)
        self.assertEqual(endpoint,[14,8,6,9])
        self.assertLess(endpoint[2],evidence['match_threshold'])
        self.assertGreaterEqual(endpoint[3],evidence['match_threshold'])
        path=traceback(pattern,query,evidence,endpoint,matrices)
        self.assertEqual((path['matches'],path['cost']),(9,8))

    def test_exact_multi_hits_all_retained(self):
        ref='ACGACGTTTTGCA'
        record=dict(pos=0,cigar=f'{len(ref)}=',query_sequence=ref)
        result=locate_interval(record,[6,10],lambda a,b:ref[a:b] if a>=0 else '',3,padding=30)
        self.assertEqual(result['status'],'success')
        self.assertEqual({p['qstart'] for p in result['pairs']},{3,6})
        self.assertTrue(result['ambiguous'])
        self.assertEqual(result['representative']['qstart'],6)

    def test_affine_trace_and_reverse_side(self):
        ref='ACGTCGATTTTTTGGCATGCA'
        query=ref[:17]+'A'+ref[17:]
        record=dict(pos=0,cigar=f'17=1I{len(ref)-17}=',query_sequence=query)
        result=locate_interval(record,[8,13],lambda a,b:ref[a:b] if a>=0 else '',7,padding=10)
        self.assertEqual(result['status'],'success')
        for side in ('left','right'):
            for path in result[side]['optimal_traces'].values():
                a,b=path['reference_interval']
                verify_alignment(ref[a:b],query,path)
        result2=locate_interval(dict(record,flag=16),[8,13],lambda a,b:ref[a:b] if a>=0 else '',7,padding=10)
        self.assertEqual(result,result2)

    def test_unfinished_is_not_unique(self):
        ref='ACGTCGATTTTTTGGCATGCA'
        query=ref.replace('C','A')
        result=locate_interval(dict(pos=0,cigar=f'{len(query)}M',query_sequence=query),[8,13],
                               lambda a,b:ref[a:b] if a>=0 else '',7,padding=10,max_cells=2)
        self.assertEqual(result['status'],'search_incomplete')
        self.assertNotIn('representative',result)

    def test_cigar_deletion_skip_and_query_consumption(self):
        record=dict(pos=10,cigar='2S3=2D4=3N2=1H',query_sequence='AACGTTTTTAA')
        segments=cigar_segments(record)
        self.assertEqual(cigar_hint(segments,14),5)
        self.assertIsNone(cigar_hint(segments,20))
        self.assertEqual(flank_search(record,[18,23],'right')['status'],'reference_skip_in_template')
        with self.assertRaisesRegex(ValueError,'length_mismatch'):
            cigar_segments(dict(record,query_sequence='A'))


if __name__=='__main__':
    unittest.main()
