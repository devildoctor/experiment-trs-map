"""B2 default, reference identity, evidence admission and failure accounting."""
import copy
from pathlib import Path
import random
import sys
from types import SimpleNamespace
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrflankalign import build_context, FlankContext, core_eligible
from vntrreads import extract_record
from vntrpathgraph import select_paths
from vntrsummary import strong_boundary

class B2Tests(unittest.TestCase):
    def setUp(self):
        rng=random.Random(61)
        self.dna=''.join(rng.choice('ACGT') for _ in range(1100))
        self.locus=SimpleNamespace(name='X',chrom='chr1',start=500,end=570,boundary_status='confirmed')
        self.reference=SimpleNamespace(fetch=lambda chrom,a,b:self.dna[a:b],contig_length=lambda c:len(self.dna))
        self.context=build_context(self.locus,self.reference)
        self.fields=['r','0','chr1','1','60','1100=','*','0','0',self.dna,'*']

    def fragment(self):
        f=extract_record(self.fields,self.locus,self.context)
        f.update(id='F1',record_id='rec1')
        return f

    def row(self,f):
        return dict(id='R1',read_id='r',sample='S',namespace='',platform='HiFi',bam='x',fragments=[f])

    def test_exact_default_and_no_other_strategy(self):
        f=self.fragment()
        self.assertEqual((f['qstart'],f['qend'],f['sequence']),(500,570,self.dna[500:570]))
        self.assertTrue(core_eligible(f))
        self.assertEqual(self.context.model['anchor_length'],250)
        self.assertNotIn('diagnostic_100bp',f)
        self.assertNotIn('analysis_anchor_evidence',f)
        self.assertEqual(len(select_paths([self.row(f)])[0]),1)

    def test_ambiguity_retained_but_not_strong(self):
        f=self.fragment();self.assertTrue(strong_boundary(self.row(f),f))
        f['core_ambiguous']=True
        self.assertEqual(len(select_paths([self.row(f)])[0]),1)
        self.assertFalse(strong_boundary(self.row(f),f))

    def test_no_old_score_or_version_bypass(self):
        for change in [dict(core_status='unavailable'),dict(extraction_version='flank-v2.1'),
                       dict(core_anchor_evidence={}),dict(sequence=self.dna[500:569])]:
            f=self.fragment();f.update(change,anchor_scores=[1,1])
            self.assertFalse(core_eligible(f));self.assertFalse(select_paths([self.row(f)])[0])

    def test_model_rejects_ablation_outer_and_mismatch(self):
        for changes in [dict(anchor_length=100,purpose='ablation'),dict(analysis_interval=[480,590]),
                        dict(reference_sha256='bad'),dict(normalization_policy='other'),
                        dict(search_padding=-1),dict(max_dp_cells=0),dict(core_interval=[500,571])]:
            with self.assertRaises(ValueError):FlankContext(dict(self.context.model,**changes),self.locus,self.reference)

    def test_partial_preserved_but_not_complete(self):
        self.fields[3]='251';self.fields[5]='280=';self.fields[9]=self.dna[250:530]
        f=self.fragment()
        self.assertEqual(f['sequence'],self.dna[500:530])
        self.assertEqual(f['anchor'],'left_only')
        self.assertFalse(core_eligible(f))

    def test_missing_reference_template_is_unavailable(self):
        self.locus.start=100;self.locus.end=170
        self.context=build_context(self.locus,self.reference)
        f=self.fragment()
        self.assertEqual(f['core_anchor_evidence']['left']['status'],'reference_template_unavailable')
        self.assertFalse(core_eligible(f))

if __name__=='__main__':unittest.main()
