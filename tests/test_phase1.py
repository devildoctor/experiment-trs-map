"""Phase 1 biological counterexamples and lossless-output contracts."""
import copy
from dataclasses import asdict
from pathlib import Path
import random
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrcore import Read,Allele,Token,decompose,motif_dictionary,canonical_motif,revcomp,discover_motifs
from vntrmotifs import exact_runs,context_support,compound_annotations,base_metrics
from vntrstructure import seed_dictionary,segment,learn_variants,annotate_rows
from vntrreads import annotate


def reads(sequence,n=3):
    return [Read('r'+str(i),'S','HiFi',sequence) for i in range(n)]


def structure_rows(sequence,dictionary,n=3):
    return [dict(id='r'+str(i),read_id='r'+str(i),sample='S',platform='HiFi',bam='synthetic.bam',
                 record_count=1,records=[],status='displayed',fragments=[dict(id='r'+str(i)+'F1',sequence=sequence,
                 qstart=0,qend=len(sequence),mapq=60,alignment_type='primary',flags=[],coverage='spans_candidate_interval',
                 anchor='both',tokens=segment((sequence,dictionary,'HiFi')))]) for i in range(n)]


class Phase1Tests(unittest.TestCase):
    def test_01_primitive_and_catalog_compound(self):
        seq='ACGT'*30
        records=motif_dictionary(seq,[Allele('a','S',seq)],reads(seq),('ACGTACGT',),3)
        base=[m['motif'] for m in records if m['usable'] and m['level']=='base']
        ts,_,audit=decompose(seq,base,return_details=True)
        self.assertEqual(len(ts),30)
        compounds=compound_annotations(ts,['ACGTACGT','ACGTACGTACGT'])
        self.assertEqual(sum(x['count'] for x in compounds if x.get('catalog_motif')=='ACGTACGT'),15)
        self.assertEqual(sum(x['count'] for x in compounds if x.get('catalog_motif')=='ACGTACGTACGT'),10)
        self.assertEqual(len(segment((seq,seed_dictionary(['ACGTACGT']),'HiFi'))),30)
        self.assertEqual(audit['metrics']['exact_full_bases'],120)

    def test_02_local_mixed_periods(self):
        seq='ACGT'*25+'GATAA'*20
        self.assertEqual([(h['period'],h['count']) for h in exact_runs(seq,max_period=200,min_period=2)],[(4,25),(5,20)])
        records=motif_dictionary(seq,[Allele('a','S',seq)],reads(seq),(),3)
        base=[m['motif'] for m in records if m['usable'] and m['level']=='base']
        ts,_,_=decompose(seq,base,return_details=True)
        self.assertEqual(len(ts),45)
        self.assertEqual([len(t.exact_bases) for t in ts],[4]*25+[5]*20)

    def test_03_lone_substring_is_not_repeat_evidence(self):
        motif='AACCGTACGATGAACC';seq='N'*20+motif+'N'*20
        record=next(m for m in motif_dictionary('',[],reads(seq),(motif,),3) if m['motif']==motif)
        self.assertEqual(record['substring_read_support'],3)
        self.assertEqual(record['read_support'],0)
        self.assertFalse(record['validated']);self.assertTrue(record['usable'])

    def test_04_single_interruption_per_molecule(self):
        a='AACCGTACGATGAACC';b='AACCGTTCGATGAACC';seq=a+a+b+a+a
        records=motif_dictionary(a*5,[Allele('a','S',seq)],reads(seq),(a,),3)
        candidate=next(m for m in records if canonical_motif(m['motif'])==canonical_motif(b))
        self.assertEqual(candidate['read_support'],3);self.assertTrue(candidate['validated'])
        dic=seed_dictionary([a]);learned,_=learn_variants(structure_rows(seq,dic),dic)
        self.assertIn(canonical_motif(b),{canonical_motif(m['sequence']) for m in learned})

    def test_05_zero_support_catalog_not_validated(self):
        record=motif_dictionary('',[],[],('AACCGT',),3)[0]
        self.assertFalse(record['validated']);self.assertEqual(record['read_support'],0)
        self.assertTrue(record['usable']);self.assertEqual(record['evidence_status'],'not_observed_in_repeat_context')

    def test_06_phase_policy_and_raw_sequence(self):
        motif='AACCGT';shifted=motif[2:]+motif[:2];seq=shifted*5
        self.assertEqual(canonical_motif(motif),canonical_motif(shifted))
        ts,_,audit=decompose(seq,[motif],return_details=True)
        self.assertEqual(''.join(t.exact_bases for t in ts),seq)
        self.assertEqual(audit['phase_policy'],'fixed_template_phase')
        self.assertGreater(audit['metrics']['partial_bases']+audit['metrics']['unknown_bases'],0)
        self.assertEqual(''.join(t['exact_bases'] for t in audit['alternative']),seq)

    def test_07_mononucleotide_is_one_base_each(self):
        ts,explained=decompose('AAAA',['A']);self.assertEqual(len(ts),4);self.assertEqual(explained,1)
        self.assertEqual(len(segment(('AAAA',seed_dictionary(['A']),'HiFi'))),4)
        self.assertEqual(len(annotate(('AAAA',('A',),'HiFi'))),4)

    def test_08_one_base_partial_never_explained(self):
        ts,explained,audit=decompose('A',['ACGT'],return_details=True)
        self.assertEqual(explained,0);self.assertEqual(audit['metrics']['unknown_bases'],1)
        for tokens in [segment(('A',seed_dictionary(['ACGT']),'HiFi')),annotate(('A',('ACGT',),'HiFi'))]:
            self.assertEqual(base_metrics(tokens)['unknown_bases'],1)
        ts,explained,audit=decompose('AC',['ACGT'],return_details=True)
        self.assertEqual(explained,0);self.assertEqual(audit['metrics']['partial_bases'],2)
        seq='ACGT'*5+'A'
        for tokens in [decompose(seq,['ACGT'])[0],segment((seq,seed_dictionary(['ACGT']),'HiFi')),annotate((seq,('ACGT',),'HiFi'))]:
            self.assertEqual(base_metrics(tokens)['full_units'],5)
            self.assertEqual(base_metrics(tokens)['unknown_bases'],1)

    def test_09_harmonics_not_learned_as_60bp_unit(self):
        seq='ACGT'*30;dic=seed_dictionary(['TTTTTTTT'])
        rows=structure_rows(seq,dic);learned,audit=learn_variants(rows,dic)
        acgt=[m for m in learned if canonical_motif(m['sequence'])==canonical_motif('ACGT')]
        self.assertEqual(len(acgt),1)
        self.assertFalse(any(len(m['sequence'])>8 for m in learned))
        annotate_rows(rows,learned,1)
        self.assertEqual([len(r['fragments'][0]['tokens']) for r in rows],[30]*3)

    def test_10_duplicate_molecules_do_not_inflate_support(self):
        seq='ACGT'*8;rs=reads(seq,1)*3
        record=next(m for m in motif_dictionary(seq,[Allele('a','S',seq)],rs,('ACGT',),3) if m['motif']=='ACGT')
        self.assertEqual(record['read_support'],1);self.assertFalse(record['validated'])

    def test_11_compounds_do_not_cross_unknown_or_partial(self):
        tokens=[Token('ACGT','ACGT',0,4),Token('UNKNOWN','NN',4,6,status='unknown'),Token('ACGT','ACGT',6,10)]
        self.assertEqual(compound_annotations(tokens,['ACGTACGT']),[])
        tokens[1].status='partial';self.assertEqual(compound_annotations(tokens,['ACGTACGT']),[])

    def test_12_primary_and_alternative_lossless_and_scored(self):
        for seq in ['ACGT'*5,'ACGTACATACGT','NNACGTACGTNN','CGTACGTAC','AAAA','']:
            ts,_,audit=decompose(seq,['ACGT','ACAT','A'],return_details=True)
            for path in [list(map(asdict,ts)),audit['alternative']]:
                if path:
                    self.assertEqual(''.join(t['exact_bases'] for t in path),seq)
                    self.assertTrue(all(t['exact_bases']==seq[t['start']:t['end']] for t in path))
            if audit['alternative_score'] is not None:self.assertGreaterEqual(audit['score_gap'],0)
            dic=seed_dictionary(['ACGT','ACAT','A']);ts,audit=segment((seq,dic,'HiFi'),True)
            self.assertEqual(''.join(t['exact_bases'] for t in ts),seq)
            if audit['alternative']:self.assertEqual(''.join(t['exact_bases'] for t in audit['alternative']),seq)

    def test_13_reverse_motif_and_compound_alignment(self):
        seq=revcomp('AACCGTAC')*4;dic=seed_dictionary(['AACCGTAC'])
        ts=segment((seq,dic,'HiFi'));self.assertEqual(len(ts),4)
        self.assertTrue(all(t['orientation']=='-' for t in ts))
        self.assertEqual(len(compound_annotations(ts,['AACCGTAC'*2])),2)
        self.assertEqual(len(decompose(seq,['AACCGTAC'])[0]),4)

    def test_14_equal_cost_ambiguity_retained(self):
        seq='ACGT';ts,_,audit=decompose(seq,['ACAT','ACCT'],return_details=True)
        self.assertTrue(audit['annotation_ambiguity']);self.assertAlmostEqual(audit['score_gap'],0)
        self.assertNotEqual(ts[0].family,audit['alternative'][0]['family'])
        ts,audit=segment((seq,seed_dictionary(['ACAT','ACCT']),'HiFi'),True)
        self.assertTrue(audit['annotation_ambiguity']);self.assertAlmostEqual(audit['score_gap'],0)


if __name__=='__main__':unittest.main()
