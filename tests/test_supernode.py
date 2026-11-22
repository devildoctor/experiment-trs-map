from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrhmm import decode,segment_hmm
from vntrsupernode import select_blocks,align_regions,partner_options,CERTAIN


def units(word):
    return [dict(status='unit',motif=x,start=i,end=i+1,exact_bases=x,orientation='+',
                 phase=0,boundary_assessment=CERTAIN,structure_assessment=CERTAIN) for i,x in enumerate(word)]


class SupernodeTests(unittest.TestCase):
    def parse(self,word):
        return select_blocks(units(word),{'X'},12)

    def test_patterns_preserve_phase_and_order(self):
        for word,expected in [('AAABB',[('A',3),('B',2)]),('ABABAB',[('AB',3)]),
                ('AAABBB',[('A',3),('B',3)]),('ABABABA',[('AB',3),('A',1)]),
                ('AABAABAAB',[('AAB',3)]),('AAAXAA',[('A',3),('X',1),('A',2)])]:
            blocks=self.parse(word)
            self.assertEqual([(''.join(k[0] for k in b['pattern']),b['count']) for b in blocks],expected)
        a=self.parse('ABAB'); b=self.parse('BABA')
        self.assertNotEqual(a[0]['pattern'],b[0]['pattern'])

    def test_correspondence_ignores_counts_but_preserves_positions(self):
        ref=self.parse('AAAAAAAXAA')
        q=self.parse('AAAAAAAAAXAAAA')
        self.assertEqual(align_regions(q,ref,'P0002',True)[0],[('ref',0),('ref',1),('ref',2)])
        self.assertEqual(partner_options(['A'],['A','A']),[[('ref',0),('ref',1)]])
        self.assertTrue(all(r[0]=='local' for r in align_regions(q,ref,'P0002',False)[0]))

    def test_all_barriers_split_blocks(self):
        for change in [dict(status='partial'),dict(status='unknown'),dict(role='right_terminal'),
                       dict(structure_assessment='ambiguous'),dict(boundary_assessment='not_assessed')]:
            ts=units('AAAAA');ts[2].update(change)
            self.assertEqual([(b['start'],b['end']) for b in select_blocks(ts,[],12)],[(0,2),(2,3),(3,5)])

    def test_ambiguous_chosen_type_is_never_a_shared_anchor(self):
        ref=units('AXB');query=units('AXB')
        for ts in (ref,query):ts[1]['structure_assessment']='ambiguous'
        bs=[select_blocks(ts,[],12) for ts in (ref,query)]
        regions,opts=align_regions(bs[1],bs[0],'P0002',True)
        self.assertNotIn(('ref',1),opts[1])
        self.assertEqual((regions[0],regions[2]),(('ref',0),('ref',2)))

    def test_optimal_projections_with_seven_or_eight_units(self):
        # Six A units and either two C units or one CC. Paths tie exactly.
        inf=float('-inf')
        states=[(-1,'B',-1),(0,'M',0),(1,'M',0),(2,'M',0),(2,'M',1)]
        model=(states,[[(1,0),(2,0),(4,0)],[(0,0)],[(0,0)],[(0,0)],[(3,0)]],
               [1,2,3,4],[0],{'A':[0,0,inf,inf,inf],'C':[0,inf,0,0,0]})
        with patch('vntrhmm.build_model',return_value=model):
            spans,_,opts,counts=decode('AAAAAACC',['A','C','CC'],assess=True)
        self.assertEqual(counts,[7,8])
        self.assertEqual(opts[5],{(0,True)})
        self.assertEqual(opts[7],{(1,True),(2,False)})

    def test_same_projection_multiple_states_and_unassessed_trace(self):
        states=[(-1,'B',-1),(0,'M',0),(0,'M',1)]
        model=(states,[[(1,0),(2,0)],[(0,0)],[(0,0)]],[1,2],[0],{'A':[0,0,0]})
        with patch('vntrhmm.build_model',return_value=model):
            self.assertEqual(decode('A',['A'],True)[2:],([{(0,True)}],[1,1]))
        t,a=segment_hmm(('AAAA',[dict(id='A',base='A',sequence='AA')],'HiFi'),return_details=True)
        self.assertEqual(a['ambiguity_scope'],'not_enumerated_single_trace')
        self.assertTrue(all(x.get('structure_assessment','not_assessed')=='not_assessed' for x in t))
