import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrstructure import seed_dictionary,add_reference_units,segment,learn_variants,summarize,compound_blocks,emit
from vntrcore import revcomp
from vntrgraph import Locus,read_catalog


def row(name,sequence,dictionary,platform='HiFi',bam='input.bam',copies=1):
    f=dict(id=name+'_F1',sequence=sequence,qstart=10,qend=10+len(sequence),coverage='spans_candidate_interval',
           anchor='left_only',mapq=60,alignment_type='primary',flags=[],tokens=segment((sequence,dictionary,platform)))
    r=dict(id=name,read_id=name,sample='S1',platform=platform,bam=bam,fragments=[],record_count=copies,records=[],status='displayed')
    for i in range(copies):
        g=copy.deepcopy(f);g['id']=name+f'_F{i+1}';r['fragments'].append(g)
    return r


class StructureTests(unittest.TestCase):
    def test_no_junk_swallowing(self):
        motif='CCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGAAGAATCCATCCCAT'
        seq='N'*12+motif+'N'*12
        ts=segment((seq,seed_dictionary([motif]),'HiFi'))
        u=[t for t in ts if t['status']=='unit']
        self.assertEqual([(t['start'],t['end'],t['edit_distance']) for t in u],[(12,66,0)])
    def test_exact_repeat_reuses_node_and_support_not_copy_count(self):
        dic=seed_dictionary(['ACGT'])
        rows=[row('r1','ACGT'*4,dic,copies=2),row('r2','ACGT'*3,dic)]
        g=summarize(rows,dic)
        self.assertEqual(len(g['exact_nodes']),1)
        self.assertEqual(g['exact_nodes'][0]['support'],2)
        self.assertEqual(g['exact_nodes'][0]['occurrences'],11)
        self.assertEqual(g['edges'][0]['support'],2)
        self.assertEqual(g['edges'][0]['source'],g['edges'][0]['target'])
    def test_cross_bam_duplicates_not_independent_support(self):
        dic=seed_dictionary(['ACGT']);g=summarize([row('same','ACGT'*2,dic),row('same','ACGT'*2,dic,bam='other.bam')],dic)
        self.assertEqual(g['structure_nodes'][0]['support'],1)
    def test_reverse_motif_retained(self):
        motif='AACCGTACGATG'
        ts=segment((revcomp(motif)*3,seed_dictionary([motif]),'HiFi'))
        self.assertTrue(all(t['orientation']=='-' for t in ts));self.assertEqual(len(ts),3)
    def test_unknown_not_discarded(self):
        dic=seed_dictionary(['ACGT']);seq='ACGT'*3+'NNNNNNNN'+'ACGT'*2
        r=row('r',seq,dic);g=summarize([r],dic)
        self.assertEqual(''.join(t['exact_bases'] for t in r['fragments'][0]['tokens']),seq)
        self.assertIn('unresolved_internal_segment',r['fragments'][0]['features'])
        self.assertTrue(any(n['id'].startswith('U') for n in g['structure_nodes']))
    def test_compound_is_observed_not_invented(self):
        dic=seed_dictionary(['AACCTG','GGTCAA']);r=row('r',('AACCTG'+'GGTCAA')*3,dic)
        g=summarize([r],dic);cs=r['fragments'][0]['compounds']
        self.assertEqual(cs[0]['count'],3);self.assertEqual(cs[0]['pattern'],['M1','M2'])
    def test_learning_needs_distinct_hifi_reads(self):
        motif='AACCGTACGATGAACC';variant='AACCGTTCGATGAACC';dic=seed_dictionary([motif])
        rows=[row('r'+str(i),variant*2,dic) for i in range(3)]
        learned,audit=learn_variants(rows,dic)
        self.assertEqual(len(learned),2);self.assertEqual(learned[1]['support'],3)
        self.assertEqual(audit[0]['status'],'accepted')
        self.assertEqual(len(learn_variants([row('same',variant*5,dic,copies=4)],dic)[0]),1)
    def test_reference_interruption_is_added(self):
        dic=add_reference_units(seed_dictionary(['CTG']),'CTG'*6+'TTG'+'CTG'*5)
        self.assertIn('TTG',{m['sequence'] for m in dic})
        self.assertTrue(any(m['source']=='reference_unit' for m in dic))
    def test_short_interruption_three_reads(self):
        dic=seed_dictionary(['CTG'])
        rows=[row(f'r{i}','CTG'*4+'TTG'+'CTG'*4,dic) for i in range(3)]
        learned,audit=learn_variants(rows,dic)
        self.assertIn('TTG',{m['sequence'] for m in learned})
        self.assertTrue(any(c['status']=='accepted' for c in audit))
    def test_one_bp_terminal_not_spurious_partial(self):
        dic=seed_dictionary(['AACCGTACGATG'])
        ts=segment(('AACCGTACGATG'+'A',dic,'HiFi'))
        self.assertEqual(ts[-1]['status'],'unknown')
    def test_output_reconstructs_and_all_rows_kept(self):
        dic=seed_dictionary(['ACGT']);rows=[row('r1','ACGT'*4,dic),row('r2','ACGT'*3,dic)]
        rows.append(dict(id='r3',read_id='r3',sample='S',platform='ONT',bam='b',fragments=[],records=[],record_count=1,status='flank_only'))
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);a=emit(out,Locus('X','chr1',10,30,('ACGT',)),rows,dic,dict(read_rows=3,rows_with_repeat_sequence=2,motif_candidates=[]))
            self.assertEqual(a['shared_structure_nodes'],1)
            d=json.loads((out/'X.structures.json').read_text(encoding='utf-8'))
            self.assertEqual(len(d['rows']),3)
            self.assertIn('RC:i:2',(out/'X.reads.gfa').read_text(encoding='utf-8'))
            self.assertNotIn('__DATA__',(out/'X.structure.html').read_text(encoding='utf-8'))
    def test_catalog_seeds(self):
        loci=read_catalog(Path(__file__).resolve().parents[1]/'catalogs'/'grch38_structure_seeds_il1rn_per3_abca7_v1.bed')
        self.assertEqual({l.name:tuple(map(len,l.motifs)) for l in loci},{'IL1RN':(86,),'PER3':(54,54),'ABCA7':(25,26)})

    def test_same_read_across_platform_labels_not_double_counted(self):
        dic=seed_dictionary(['ACGT']);g=summarize([row('same','ACGT'*2,dic),row('same','ACGT'*2,dic,platform='ONT')],dic)
        self.assertEqual(g['structure_nodes'][0]['support'],1)


if __name__=='__main__':unittest.main()
