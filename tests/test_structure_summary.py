import copy
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from test_structure import row
from vntrstructure import seed_dictionary,add_reference_units,segment,learn_variants,summarize
from vntrsummary import structural_summary,exact_runs,local_candidates


def high(name,seq,dictionary,**kwargs):
    r=row(name,seq,dictionary,**kwargs)
    for f in r['fragments']:f['anchor']='both';f['anchor_scores']=[1,1]
    return r


class SummaryTests(unittest.TestCase):
    def setUp(self):self.d=add_reference_units(seed_dictionary(['CTG']),'CTG'*9+'TTG'+'CTG'*13)
    def test_sequential_names_and_distinct_reference_interruption(self):
        self.assertEqual([(m['id'],m['sequence']) for m in self.d],[('M1','CTG'),('M2','TTG')])
        rows=[high(str(i),'CTG'*5+'ATG'+'CTG'*5,self.d) for i in range(3)]
        dic,_=learn_variants(rows,self.d)
        self.assertEqual([m['id'] for m in dic],[f'M{i+1}' for i in range(len(dic))])
        self.assertTrue(all(m['base']==m['id'] for m in dic))
    def test_lengths_group_but_interruptions_split(self):
        rows=[high(str(i),'CTG'*(i+4)+'TTG'+'CTG'*4,self.d) for i in range(3)]
        rows += [high('double'+str(i),'CTG'*4+'TTG'+'CTG'*4+'TTG'+'CTG'*4,self.d) for i in range(3)]
        s=structural_summary(rows)
        self.assertEqual(len(s['classes']),2)
        self.assertEqual(s['audit']['supported_classes'],2)
        self.assertEqual(sorted(len(g['subtypes']) for g in s['classes']),[1,3])
        for g in s['classes']:
            for st in g['subtypes']:
                self.assertIn(st['representative']['sequence'],[r['fragments'][0]['sequence'] for r in rows])
    def test_unknown_and_partial_not_joined(self):
        rows=[high('unknown','CTG'*4+'NNN'+'TTG'+'CTG'*4,self.d),high('partial','CTG'*4+'TTG',self.d)]
        rows[1]['fragments'][0]['coverage']='left_partial'
        s=structural_summary(rows)
        self.assertFalse(s['classes']);self.assertEqual(s['audit']['unassigned_fragments'],2)
        self.assertFalse(any(m['kind']=='tandem' for m in s['modules']))
    def test_same_molecule_not_multiple_support(self):
        a=high('same','CTG'*4+'TTG'+'CTG'*4,self.d,copies=2);b=copy.deepcopy(a);b['id']='other_bam'
        s=structural_summary([a,b]);g=s['classes'][0]
        self.assertEqual(g['support'],1);self.assertEqual(g['hifi_support'],1);self.assertEqual(g['status'],'candidate')
        self.assertTrue(all(m['support']==1 for m in s['modules']))
    def test_molecule_conflict_kept_unassigned(self):
        a=high('same','CTG'*4+'TTG'+'CTG'*4,self.d);b=high('same','CTG'*7+'TTG'+'CTG'*4,self.d,bam='copy.bam')
        s=structural_summary([a,b]);self.assertFalse(s['classes']);self.assertEqual(s['audit']['conflicting_molecules'],1)
    def test_shared_module_and_tandem_separate(self):
        rows=[high(str(i),'CTG'*8+'TTG'+'CTG'*4+'TTG'+'CTG'*8,self.d) for i in range(3)]
        s=structural_summary(rows)
        self.assertTrue(any(m['kind']=='shared_block' and m['structure']=='M2×1 → M1×4 → M2×1' and m['support']==3 for m in s['modules']))
        self.assertTrue(any(m['kind']=='tandem' and m['status']=='supported' for m in s['modules']))
        for m in s['modules']:
            for x in m['instances']:
                f=next(r['fragments'][0] for r in rows if r['id']==x['row_id'])
                self.assertEqual(x['sequence'],f['sequence'][x['start']:x['end']])
    def test_reference_independent_period_after_offset(self):
        hits=exact_runs('N'*17+'GATTACA'*4+'N'*8)
        self.assertEqual([(h['start'],h['sequence'],h['count']) for h in hits],[(17,'GATTACA',4)])
        dic=add_reference_units(seed_dictionary(['CTG']),'CTG'*4+'N'*8+'GATTACA'*4)
        self.assertTrue(any(m['sequence']=='GATTACA' and m['source']=='reference_period' for m in dic))
    def test_ont_does_not_establish_high_confidence_class(self):
        s=structural_summary([high(str(i),'CTG'*4+'TTG'+'CTG'*4,self.d,platform='ONT') for i in range(5)])
        self.assertEqual(s['classes'][0]['status'],'candidate');self.assertEqual(s['classes'][0]['hifi_support'],0)
    def test_longer_composite_than_old_six_token_limit(self):
        dic=seed_dictionary(['ACGT','AACCTG']);r=high('r',('ACGT'*7+'AACCTG')*3,dic)
        hits=local_candidates(r['fragments'][0]['tokens'])
        self.assertTrue(any(h['kind']=='tandem' and len(h['pattern'])==8 and h['count']==3 for h in hits))
    def test_family_no_longer_erases_interruptions(self):
        r=high('r','CTG'*8+'TTG'+'CTG'*4+'TTG'+'CTG'*8,self.d);summarize([r],self.d)
        self.assertEqual(r['fragments'][0]['family_structure'],'M1×8 → M2×1 → M1×4 → M2×1 → M1×8')
    def test_single_base_reference_repeat_reconstructs(self):
        tokens=segment(('A'*15+'N',seed_dictionary(['A']),'HiFi'))
        self.assertEqual(''.join(t['exact_bases'] for t in tokens),'A'*15+'N')
        self.assertEqual(sum(t['status']=='unit' for t in tokens),15)


if __name__=='__main__':unittest.main()
