import copy
import json
from pathlib import Path
import random
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrflankbounds import build_profile, assess_window, base_model, discover, SAMPLES
from vntrflankalign import FlankContext, core_eligible
from vntrreads import extract_record
from vntrpathgraph import select_paths
from vntrsummary import strong_boundary


class ProfileTests(unittest.TestCase):
    def records(self):
        return [dict(sample=s,platform='HiFi',flag=0,mapq=60,chrom='chr1',pos=0,
                     cigar='2S4=2D2=1I2=3N2=1H',query_sequence='AA'+'AAAA'+'AA'+'C'+'AA'+'AA',
                     molecule_key=json.dumps([s,'',str(i)]),record_id=s+str(i),bam=s)
                for s in SAMPLES for i in range(5)]

    def test_event_coverage_and_fixed_normalization(self):
        rows=self.records();p=build_profile(rows,[0,15],'A'*15)
        for sample in SAMPLES:
            self.assertEqual(list(p['coverage'][sample]),[5]*10+[0]*3+[5]*2)
            self.assertEqual(p['alts'][sample][0],5)  # one softclip event, not two bases
            self.assertEqual(p['alts'][sample][4],5)  # deletion is observed
            self.assertEqual(p['alts'][sample][8],5)  # insertion at reference cursor
            self.assertEqual(p['depths'][sample],4)
        q=build_profile(rows+rows,[0,15],'A'*15)
        self.assertEqual(p['depths'],q['depths'])
        self.assertEqual(assess_window(p,0,'right',15)['status'],'coverage_insufficient')

    def test_no_missing_sample_average_and_conflicting_molecule(self):
        rows=self.records();bad=copy.deepcopy(rows[0]);bad['query_sequence']='C'+bad['query_sequence'][1:]
        p=build_profile(rows+[bad],[0,15],'A'*15)
        self.assertEqual(len(p['conflicts']),1)
        self.assertEqual(int(p['coverage']['LCL5'][1]),4)
        p=build_profile([r for r in rows if r['sample']!='LCL8'],[0,15],'A'*15)
        self.assertEqual(assess_window(p,0,'right',4)['status'],'coverage_insufficient')

    def test_terminal_same_point_recheck_and_cap(self):
        locus=SimpleNamespace(name='X',chrom='chr1',start=252,end=260)
        profile=dict(fetch_interval=[0,512],dna='A'*512,model_sha256='x',depths={},conflicts=[])
        def score(p,b,side,w):
            distance=252-b if side=='left' else b-260
            passed=(distance==2 or (distance==0 and w>=50) or (distance==1 and w<50))
            return dict(status='evaluable',probability=float(passed),minimum_coverage={s:5 for s in SAMPLES})
        def support(p,b,side):
            return dict(status='pass',counts={s:5 for s in SAMPLES},members={s:[side+str(i) for i in range(5)] for s in SAMPLES})
        with patch('vntrflankbounds.assess_window',score),patch('vntrflankbounds.anchor_support',support):
            models,audit=discover(locus,profile,{'X':[100,10]})
        self.assertEqual(models['X']['analysis_interval'],[250,262])
        self.assertTrue(all(not v for v in models['X']['support_intersections'].values()))
        self.assertEqual(len(audit),6)

    def test_incomplete_first_candidate_cannot_be_skipped(self):
        locus=SimpleNamespace(name='X',chrom='chr1',start=252,end=260)
        profile=dict(fetch_interval=[0,512],dna='A'*512,model_sha256='x',depths={},conflicts=[])
        with patch('vntrflankbounds.assess_window',return_value=dict(status='evaluable',probability=1)), \
             patch('vntrflankbounds.anchor_support',return_value=dict(status='search_incomplete',counts={})):
            models,audit=discover(locus,profile,{'X':[100,10]})
        self.assertIsNone(models['X']['analysis_interval'])
        self.assertEqual(len(audit),2)


class IntegrationTests(unittest.TestCase):
    def setUp(self):
        rng=random.Random(61)
        self.dna=''.join(rng.choice('ACGT') for _ in range(1100))
        self.locus=SimpleNamespace(name='X',chrom='chr1',start=500,end=570,boundary_status='catalog')
        self.reference=SimpleNamespace(fetch=lambda chrom,a,b:self.dna[a:b])
        self.model=base_model(self.locus,[0,1100],self.dna)
        self.fields=['r','0','chr1','1','60','1100=','*','0','0',self.dna,'*']

    def fragment(self,model):
        context=FlankContext(model,self.locus,self.reference)
        fragment=extract_record(self.fields,self.locus,'','',context)
        fragment.update(id='F1',record_id='rec1')
        return fragment

    def row(self,fragment):
        return dict(id='R1',read_id='r',sample='LCL5',namespace='',platform='HiFi',bam='x',fragments=[fragment])

    def test_core_only_and_outer_exact_reconstruction(self):
        core=self.fragment(self.model)
        self.assertEqual(core['sequence'],self.dna[500:570])
        self.assertEqual(core['extraction_status'],'core_only')
        model=dict(self.model,analysis_interval=[480,600]);both=self.fragment(model)
        self.assertEqual(both['extraction_status'],'both_consistent')
        self.assertEqual(both['core_anchor_evidence']['pairs'],core['core_anchor_evidence']['pairs'])
        self.assertEqual(both['analysis_sequence'],both['left_context']+both['sequence']+both['right_context'])
        paths,excluded=select_paths([self.row(both)])
        self.assertEqual(len(paths),1);self.assertFalse(excluded)
        self.assertTrue(strong_boundary(self.row(both),both))
        both['core_ambiguous']=True
        self.assertEqual(len(select_paths([self.row(both)])[0]),1)
        self.assertFalse(strong_boundary(self.row(both),both))

    def test_conflict_and_old_scores_cannot_bypass(self):
        fragment=self.fragment(self.model)
        fragment.update(core_outer_relation='conflict',extraction_status='core_outer_conflict',anchor_scores=[1,1])
        self.assertFalse(core_eligible(fragment))
        self.assertFalse(select_paths([self.row(fragment)])[0])
        fragment.update(core_outer_relation='not_checked',core_status='unavailable',extraction_status='display_only')
        self.assertFalse(select_paths([self.row(fragment)])[0])

    def test_model_identity_and_template_extent(self):
        with self.assertRaisesRegex(ValueError,'hash_mismatch'):
            self.fragment(dict(self.model,reference_sha256='bad'))
        with self.assertRaisesRegex(ValueError,'outside_fetch'):
            self.fragment(dict(self.model,analysis_interval=[100,570]))
        with self.assertRaisesRegex(ValueError,'normalization'):
            self.fragment(dict(self.model,normalization_policy='other'))

    def test_invalid_limits_and_slice_length_are_rejected(self):
        for updates in ({'search_padding':-1},{'max_dp_cells':0},{'max_dp_cells':True},
                        {'analysis_interval':[570,500]},{'fetch_interval':[0.0,1100]}):
            with self.assertRaises(ValueError):
                self.fragment(dict(self.model,**updates))
        fragment=self.fragment(self.model);fragment['sequence']=fragment['sequence'][:-1]
        self.assertFalse(core_eligible(fragment))

    def test_core_and_outer_same_target_do_not_alias(self):
        fragment=self.fragment(dict(self.model,analysis_interval=[500,570]))
        core,outer=fragment['core_anchor_evidence'],fragment['analysis_anchor_evidence']
        self.assertEqual(core['pairs'],outer['pairs'])
        self.assertEqual(core['representative'],outer['representative'])
        self.assertNotIn('any_compatible_pair',core)
        self.assertTrue(outer['any_compatible_pair'])

    def test_outer_incomplete_does_not_invalidate_completed_core(self):
        model=dict(self.model,analysis_interval=[480,600])
        context=FlankContext(model,self.locus,self.reference)
        full=self.fragment(model)
        core=full['core_anchor_evidence']
        with patch.object(context,'locate',side_effect=[core,dict(status='search_incomplete',pairs=[],ambiguous=False),core]):
            fragment=extract_record(self.fields,self.locus,'','',context)
        self.assertEqual(fragment['extraction_status'],'core_only')
        self.assertTrue(core_eligible(fragment))
        with patch.object(context,'locate',side_effect=[dict(status='search_incomplete',pairs=[],ambiguous=False),core,core]):
            fragment=extract_record(self.fields,self.locus,'','',context)
        self.assertEqual(fragment['extraction_status'],'analysis_only')
        self.assertFalse(core_eligible(fragment))


if __name__=='__main__':unittest.main()
