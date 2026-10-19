"""Acceptance cases from the source-review document, plus failure boundaries."""
import importlib.util
import random
import sys
import tempfile
from dataclasses import replace
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import vntrcore as core
import vntrgraph as app


def naive_distance(a,b):
    row = list(range(len(b)+1))
    for i,x in enumerate(a,1):
        new = [i]
        for j,y in enumerate(b,1):
            new.append(min(row[j]+1,new[-1]+1,row[j-1]+(x!=y)))
        row = new
    return row[-1]


class DistanceTests(unittest.TestCase):
    def test_same_length_counterexample(self):
        self.assertEqual(core.edit_distance('ACGT'*30,'ACAT'*30),30)
        self.assertLess(core.similarity('ACGT'*30,'ACAT'*30),1)

    def test_bitvectors_against_independent_matrix(self):
        rng = random.Random(141)
        for _ in range(100):
            a = ''.join(rng.choices('ACGT',k=rng.randrange(45)))
            b = ''.join(rng.choices('ACGT',k=rng.randrange(45)))
            self.assertEqual(core.edit_distance(a,b),naive_distance(a,b),(a,b))

    def test_traceback_reconstructs_and_has_exact_distance(self):
        rng = random.Random(18)
        for _ in range(50):
            a = ''.join(rng.choices('ACGT',k=rng.randrange(1,45)))
            b = ''.join(rng.choices('ACGT',k=rng.randrange(1,45)))
            ops = core.align_to_seed(a,b,band=1)
            self.assertEqual(''.join(b[j] for op,i,j in ops if op!='D'),b)
            self.assertEqual(''.join(a[i] for op,i,j in ops if op!='I'),a)
            self.assertEqual(sum(op!='M' or a[i]!=b[j] for op,i,j in ops),naive_distance(a,b))

    def test_long_sequence_not_length_only(self):
        self.assertEqual(core.edit_distance('ACGT'*8000,'ACAT'*8000),8000)

    def test_traceback_limit_explicit(self):
        with self.assertRaises(core.AlignmentUnresolved):
            core.align_to_seed('A'*100,'C'*100,max_cells=100)


class ExtractionTests(unittest.TestCase):
    def setUp(self):
        rng = random.Random(51)
        self.left = ''.join(rng.choices('ACGT',k=150))
        self.right = ''.join(rng.choices('ACGT',k=150))
        self.repeat = 'ACGT'*4
        self.seq = self.left+self.repeat+self.right
        class Reference:
            def __init__(ref,seq): ref.seq=seq
            def fetch(ref,chrom,start,end): return ref.seq[start:end]
            def contig_length(ref,chrom): return len(ref.seq)
        self.ref = Reference(self.seq)
        self.locus = app.Locus('T','chrT',150,166,('ACGT',))
        self.header = [f'@SQ\tSN:chrT\tLN:{len(self.seq)}']

    def record(self,name='read1',flag=0,cigar='316M',seq=None,pos=1,rg=''):
        seq = self.seq if seq is None else seq
        return f'{name}\t{flag}\tchrT\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t'+('I'*len(seq))+(f'\tRG:Z:{rg}' if rg else '')

    def extract(self,records,platform='HiFi',header=None):
        profiles={}
        for side,start,seq in [('L',0,self.left),('R',166,self.right)]:
            profiles['T:'+side]=[dict(length=150,start=start,end=start+150,sequence=seq,
                                     grade='preferred',unique_starts=list(range(120)))]
        return app.extract_sam_records((header or self.header)+records,'S',platform,self.locus,self.ref,150,20,'test.bam',profiles=profiles)

    def test_reverse_flag_does_not_flip_stored_seq(self):
        reads,_ = self.extract([self.record('f',0),self.record('r',16)])
        self.assertEqual([r.sequence for r in reads],[self.repeat,self.repeat])
        self.assertEqual(reads[1].query_sequence,self.seq)

    def test_supplementary_counts_once(self):
        reads,_ = self.extract([self.record(flag=0),self.record(flag=2048)])
        unique,audit = core.deduplicate_reads(reads)
        self.assertEqual(len(unique),1)
        self.assertEqual(unique[0].sequence,self.repeat)
        self.assertEqual(len(unique[0].records),2)
        self.assertEqual(audit[0]['observations'],2)

    def test_hard_clipped_records_not_concatenated(self):
        reads,audit = self.extract([self.record(flag=2048,cigar='160M156H',seq=self.seq[:160]),
                                   self.record(flag=2048,cigar='160H156M',seq=self.seq[160:],pos=161)])
        self.assertEqual(reads,[])
        self.assertEqual(audit['missing_or_ambiguous_anchor'],2)

    def test_original_query_interval_records_hard_clips(self):
        reads,_ = self.extract([self.record(flag=16,cigar='10H316M20H')])
        self.assertEqual(reads[0].original_query_interval,(170,186))
        self.assertEqual(reads[0].sequence,self.repeat)

    def test_library_identity(self):
        header = self.header+['@RG\tID:r1\tLB:lib1','@RG\tID:r2\tLB:lib2']
        reads,_ = self.extract([self.record(rg='r1'),self.record(rg='r2')],header=header)
        self.assertEqual(len(core.deduplicate_reads(reads)[0]),2)

    def test_cross_platform_duplicate_is_one_molecule(self):
        reads1,_ = self.extract([self.record()],platform='ONT')
        reads2,_ = self.extract([self.record()],platform='ONTUL')
        unique,_ = core.deduplicate_reads(reads1+reads2)
        self.assertEqual(len(unique),1)
        self.assertEqual(set(unique[0].observed_platforms),{'ONT','ONTUL'})

    def test_same_molecule_conflict_not_double_support(self):
        a = core.Read('r','S','ONT','ACGT')
        b = core.Read('r','S','ONTUL','ACAT')
        unique,audit = core.deduplicate_reads([a,b])
        self.assertEqual(unique,[])
        self.assertEqual(audit[0]['status'],'conflicting_sequences')

    def test_duplicate_exact_flank_is_ambiguous(self):
        self.assertIsNone(app.flank_hit('ACGT','ACGTACGT',0,.9))

    def test_large_boundary_insertion_does_not_hide_ontul_anchor(self):
        expansion='ACGT'*200
        seq=self.left+expansion+self.right
        reads,_=self.extract([self.record(cigar='166M784I150M',seq=seq)],platform='ONTUL')
        self.assertEqual(len(reads),1)
        self.assertEqual(reads[0].sequence,expansion)

    def test_window_and_confirmed_boundaries_separate(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as folder:
            p = Path(folder)/'loci.bed'
            p.write_text('chr1\t100\t200\tID=T;REPEAT_START=110;REPEAT_END=190;BOUNDARY=confirmed\n')
            locus = app.read_catalog(p)[0]
            self.assertEqual((locus.search_start,locus.search_end,locus.start,locus.end),(100,200,110,190))


class CandidateTests(unittest.TestCase):
    def reads(self,seq,platform='HiFi',prefix='a',count=3):
        return [core.Read(f'{prefix}{i}','S',platform,seq) for i in range(count)]

    def test_ontul_only_expansion_can_seed(self):
        reads = self.reads('ACGT'*10)+self.reads('ACGT'*60,'ONTUL','u')
        alleles,assignments,genotypes = core.call_candidates(reads)
        self.assertEqual({a.sequence for a in alleles},{'ACGT'*10,'ACGT'*60})
        self.assertTrue(all(len(a.read_ids)==3 for a in alleles))

    def test_joint_support_crosses_platform_threshold(self):
        reads = self.reads('ACGT'*10,count=2)+self.reads('ACGT'*10,'ONT','o',1)
        alleles = core.cluster_reads(reads,3)
        self.assertEqual(len(alleles),1)
        self.assertEqual(alleles[0].support,{'HiFi':2,'ONT':1})

    def test_ambiguous_assignment_not_alternated(self):
        a = 'ACGT'*30
        b = a[:40]+'T'+a[41:80]+'T'+a[81:]
        mid = a[:40]+'T'+a[41:]
        reads = self.reads(a)+self.reads(b,prefix='b')+[core.Read('tie','S','HiFi',mid)]
        alleles,assignments,_ = core.call_candidates(reads)
        tied = next(x for x in assignments if x['molecule_key'].endswith(',"tie"]'))
        self.assertEqual(tied['status'],'ambiguous')
        self.assertEqual(tied['score_gap'],0)
        self.assertEqual(sorted(len(a.read_ids) for a in alleles),[3,3])

    def test_three_candidates_retained_and_not_forced_diploid(self):
        reads = sum((self.reads('ACGT'*n,prefix=str(n)) for n in (10,20,30)),[])
        alleles,_,genotypes = core.call_candidates(reads)
        self.assertEqual(len(alleles),3)
        self.assertEqual(genotypes[0]['status'],'ploidy_conflict')
        self.assertTrue(all(a.uncertain for a in alleles))

    def test_one_observed_sequence_not_homozygous(self):
        _,_,genotypes = core.call_candidates(self.reads('ACGT'*10))
        self.assertEqual(genotypes[0]['status'],'one_observed_sequence_not_homozygous')

    def test_stable_variant_with_distinct_read_errors_splits(self):
        sequence = 'ACGT'*50
        reads = []
        for group,base in enumerate(('A','T')):
            for i in range(3):
                s = list(sequence)
                s[40] = base
                s[100+i*4] = 'C'
                reads.append(core.Read(f'{group}-{i}','S','HiFi',''.join(s)))
        alleles,_,_ = core.call_candidates(reads)
        supported = [a for a in alleles if len(a.read_ids)>=3]
        self.assertEqual(len(supported),2)
        self.assertEqual({a.sequence[40] for a in supported},{'A','T'})

    def test_dedup_counts_not_platform_sum(self):
        reads = self.reads('ACGT'*10,'ONT')
        reads += [replace(r,platform='ONTUL') for r in reads]
        allele = core.cluster_reads(reads,3)[0]
        self.assertEqual(len(allele.read_ids),3)
        self.assertEqual(sum(allele.support.values()),6)


class MotifGraphTests(unittest.TestCase):
    def allele(self,name,sequence,motifs):
        tokens,_ = core.decompose(sequence,motifs)
        return core.Allele(name,'S',sequence,{'HiFi':1},('["S","","r"]',),1,tokens)

    def test_variant_keeps_exact_bases(self):
        allele = self.allele('A','ACGTACATACGT',('ACGT',))
        self.assertEqual(''.join(t.exact_bases for t in allele.tokens),allele.sequence)
        self.assertIn('ACAT',[t.exact_bases for t in allele.tokens])
        self.assertTrue(any('X' in t.edit_operations for t in allele.tokens))

    def test_partial_and_mixed_units_count_separately(self):
        allele = self.allele('A','ACGTACATAC',('ACGT','ACAT'))
        counts = core.token_counts(allele.tokens)
        self.assertEqual(sum(counts['full'].values()),2)
        self.assertEqual(sum(counts['partial'].values()),1)

    def test_compound_and_base_dictionary_retained(self):
        sequence = 'ACGTACAT'*10
        reads = [core.Read(str(i),'S','HiFi',sequence) for i in range(3)]
        alleles = core.cluster_reads(reads,3)
        dictionary = core.motif_dictionary(sequence,alleles,reads,('ACGT','ACAT','ACGTACAT'),3)
        self.assertEqual(next(x for x in dictionary if x['motif']=='ACGTACAT')['level'],'compound')
        self.assertTrue(all(x['validated'] for x in dictionary))

    def test_unknown_new_unit_requires_molecule_evidence(self):
        dictionary = core.motif_dictionary('GATTACA'*10,[],[],(),3)
        self.assertTrue(dictionary)
        self.assertTrue(all(not x['validated'] for x in dictionary))

    def test_context_nodes_prevent_recombination(self):
        a = self.allele('A','ACGTTTAA',('ACGT','TTAA','GGCC','ATAT'))
        b = self.allele('B','GGCCACGTATAT',('ACGT','TTAA','GGCC','ATAT'))
        nodes,edges,paths = core.build_graph([a,b],'LEFT','RIGHT')
        adjacency = {}
        for x,y in edges: adjacency.setdefault(x,[]).append(y)
        possible = []
        def walk(path):
            if path[-1]=='R': possible.append(path); return
            for nxt in adjacency.get(path[-1],[]): walk(path+[nxt])
        walk(['L'])
        self.assertEqual({tuple(p) for p in possible},{tuple(ids) for _,ids in paths})

    def test_one_read_ten_units(self):
        allele = self.allele('A','ACGT'*10,('ACGT',))
        nodes,edges,paths = core.build_graph([allele],'TT','AA')
        graph = core.graph_evidence(nodes,edges,paths)
        self.assertEqual(graph['paths'][0]['read_support'],1)
        self.assertEqual(sum(graph['paths'][0]['motif_copy_count']['full'].values()),10)
        self.assertTrue(all(n['read_support']==1 for n in graph['nodes']))

    def test_gfa_svg_roundtrip_same_observed_bases(self):
        allele = self.allele('A','ACGT'*10+'ACATAC',('ACGT','ACAT'))
        with tempfile.TemporaryDirectory(dir=ROOT) as folder:
            folder=Path(folder)
            nodes,edges,paths = core.build_graph([allele],'TT','AA')
            app.write_gfa(folder/'T.gfa',nodes,edges,paths)
            app.write_svg(folder/'T.svg','T',[allele])
            seqs,observed = {},[]
            for line in (folder/'T.gfa').read_text().splitlines():
                p=line.split('\t')
                if p[0]=='S': seqs[p[1]]=p[2]
                if p[0]=='P': observed.append(''.join(seqs[x[:-1]] for x in p[2].split(',')))
            self.assertEqual(observed,['TT'+allele.sequence+'AA'])
            svg=ET.parse(folder/'T.svg')
            rebuilt=''.join(e.attrib['data-exact']*int(e.attrib['data-count']) for e in svg.iter() if 'data-exact' in e.attrib)
            self.assertEqual(rebuilt,allele.sequence)

    def test_over_30000_bases_not_silently_skipped(self):
        sequence='ACGT'*7501
        tokens,ratio=core.decompose(sequence,('ACGT',))
        self.assertEqual(len(tokens),7501)
        self.assertEqual(ratio,1)


if __name__=='__main__': unittest.main()
