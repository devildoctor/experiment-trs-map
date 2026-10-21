import random
import sys
from pathlib import Path
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrdbg import raw_graph,compact_graph,supported_records
from vntrdbg_inputs import intervals
from vntranchors import Record


def record(identifier,seq,molecule=None,sample='S',kind='observed'):
    return dict(id=identifier,sequence=seq,kind=kind,molecule_key=molecule or identifier,
                sample=sample,platforms=['HiFi'],coverage='test')


class TestDBG(unittest.TestCase):
    def graph(self,records,k):
        return compact_graph(records,raw_graph(records,k),k,'TEST','all_evidence')

    def test_no_overlap_only_edges(self):
        rr=[record('p','AAC'),record('q','ACG')]
        self.assertFalse(raw_graph(rr,3)['edges'])

    def test_endpoint_splits_unitig(self):
        graph=self.graph([record('a','AACGT'),record('b','ACGT')],3)
        self.assertGreaterEqual(len(graph['nodes']),2)
        self.assertEqual(len(graph['paths']),2)

    def test_cycle_and_repeated_visits(self):
        graph=self.graph([record('a','ACACACAC'),record('b','ACAC')],3)
        self.assertTrue(any(e['source']==e['target'] for e in graph['edges']) or len(graph['edges'])>=2)
        self.assertTrue(all(n['molecule_support']==2 for n in graph['nodes'].values()))

    def test_duplicate_molecule_one_vote(self):
        graph=self.graph([record('a','AACGT','same'),record('b','AACGT','same')],3)
        self.assertTrue(all(n['molecule_support']==1 for n in graph['nodes'].values()))

    def test_reference_candidate_not_observed_support(self):
        graph=self.graph([record('ref','AACGT',kind='reference'),record('c','AACGT',kind='candidate')],3)
        self.assertTrue(all(n['molecule_support']==0 for n in graph['nodes'].values()))

    def test_filter_does_not_pool_samples(self):
        rr=[record('a','AACGT',sample='A'),record('b','AACGT',sample='B')]
        kept,_,edges=supported_records(rr,raw_graph(rr,3),3,2)
        self.assertFalse(kept)
        self.assertFalse(edges)

    def test_filter_splits_without_bridging(self):
        rr=[record('a','AAACCCGGG'),record('b','AAACC'),record('c','CCGGG')]
        kept,_,_=supported_records(rr,raw_graph(rr,3),3,2)
        pieces=[r for r in kept if r.get('parent_id')=='a']
        self.assertGreater(len(pieces),1)
        self.assertTrue(all(not r['filter_complete'] for r in pieces))
        self.graph(kept,3)

    def test_orientation_not_canonicalized(self):
        graph=self.graph([record('a','AAAC'),record('b','GTTT')],3)
        self.assertEqual(len(graph['nodes']),2)

    def test_random_paths_exact(self):
        rng=random.Random(21)
        for _ in range(25):
            rr=[record(f'r{i}',''.join(rng.choices('ACGT',k=rng.randint(10,70)))) for i in range(8)]
            self.graph(rr,rng.choice([3,5,7]))

    def test_cigar_skip_splits_and_internal_insertions_preserved(self):
        r=Record('x',0,'chrT',10,60,'3M2I3M5N4M','A'*12,'',{},'','','source',[(3,'M'),(2,'I'),(3,'M'),(5,'N'),(4,'M')],0,0,12,'HiFi')
        self.assertEqual(intervals(r,10,25),[(0,8),(8,12)])


if __name__=='__main__':
    unittest.main()
