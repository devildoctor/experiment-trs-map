"""Acceptance tests for input scope, anchor ambiguity, and physical read evidence."""
from pathlib import Path
from unittest.mock import patch
import random
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import vntranchors as a
import vntrgraph as app
from vntrcore import deduplicate_reads, revcomp


def dp_costs(pattern,text):
    row = [0]*(len(text)+1)
    for i,x in enumerate(pattern,1):
        new = [i]
        for j,y in enumerate(text,1):
            new.append(min(row[j]+1,new[-1]+1,row[j-1]+(x!=y)))
        row = new
    return row


class AnchorTests(unittest.TestCase):
    def setUp(self):
        rng=random.Random(51)
        self.left=''.join(rng.choices('ACGT',k=500))
        self.right=''.join(rng.choices('ACGT',k=500))
        self.repeat='ACGT'*4
        self.seq=self.left+self.repeat+self.right
        class Reference:
            def fetch(ref,chrom,start,end):return self.seq[start:end]
            def contig_length(ref,chrom):return len(self.seq)
        self.ref=Reference()
        self.locus=app.Locus('T','chrT',500,516,('ACGT',))
        self.profiles={}
        for side in ('L','R'):
            ts=[]
            for length in a.LENGTHS:
                start=500-length if side=='L' else 516
                ts.append(dict(length=length,start=start,end=start+length,sequence=self.seq[start:start+length],
                               grade='preferred',unique_starts=list(range(length-30))))
            self.profiles['T:'+side]=ts

    def rec(self,seq=None,name='r',flag=0,pos=1,mapq=60,cigar=None,chrom='chrT'):
        seq=self.seq if seq is None else seq
        cigar=cigar or f'{len(seq)}M'
        return f'{name}\t{flag}\t{chrom}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t'+('I'*len(seq) if seq!='*' else '*')

    def extract(self,records,platform='HiFi'):
        return a.extract_records([f'@SQ\tSN:chrT\tLN:{len(self.seq)}']+records,'S',platform,
                                 self.locus,self.ref,self.profiles,20,'test')

    def test_myers_infix_matches_independent_dp(self):
        rng=random.Random(72)
        for _ in range(100):
            pattern=''.join(rng.choices('ACGT',k=rng.randrange(1,40)))
            text=''.join(rng.choices('ACGTN',k=rng.randrange(1,100)))
            self.assertEqual(a.endpoint_costs(pattern,text),dp_costs(pattern,text))

    def test_complete_and_reverse_stored_sequence(self):
        reads,_,_=self.extract([self.rec(name='f'),self.rec(name='r',flag=16)])
        self.assertEqual([r.sequence for r in reads],[self.repeat]*2)
        self.assertEqual(reads[1].query_sequence,self.seq)

    def test_low_mapq_and_unknown_retained_but_not_support(self):
        for mq in (0,19,255):
            reads,audit,decisions=self.extract([self.rec(mapq=mq)])
            self.assertFalse(reads)
            self.assertEqual(decisions[0]['status'],'sequence_only_unverified')
            self.assertEqual(audit['input_records'],1)

    def test_secondary_is_not_sole_location_evidence(self):
        reads,audit,decisions=self.extract([self.rec(flag=256)])
        self.assertFalse(reads)
        self.assertEqual(audit['secondary_retained'],1)

    def test_softclipped_sequence_alone_does_not_support(self):
        sequence=self.seq+'ACGTACGTAC'
        reads,_,decisions=self.extract([self.rec(seq=sequence,pos=501,cigar='1016S10M')])
        self.assertFalse(reads)
        self.assertEqual(decisions[0]['status'],'sequence_only_unverified')

    def test_supplementary_validates_full_low_mapq_provider_once(self):
        records=[self.rec(mapq=5),
                 self.rec(seq=self.left,flag=2048,cigar='500M516H'),
                 self.rec(seq=self.right,flag=2048,pos=517,cigar='516H500M')]
        reads,_,_=self.extract(records)
        unique,_=deduplicate_reads(reads)
        self.assertEqual(len(unique),1)
        self.assertEqual(unique[0].sequence,self.repeat)
        self.assertEqual(len(unique[0].records),3)

    def test_records_from_different_supplied_files_are_combined(self):
        header=[f'@SQ\tSN:chrT\tLN:{len(self.seq)}']
        primary,_=a.parse_records(header+[self.rec(mapq=5)],'S',self.ref,self.locus,'one.bam',None,'ONT')
        supplements,_=a.parse_records(header+[
            self.rec(seq=self.left,flag=2048,cigar='500M516H'),
            self.rec(seq=self.right,flag=2048,pos=517,cigar='516H500M')],
            'S',self.ref,self.locus,'two.bam',None,'ONTUL')
        reads,_,_=a.extract_parsed(primary+supplements,a.Counter(),'S',self.locus,self.profiles,20)
        self.assertEqual(len(deduplicate_reads(reads)[0]),1)
        self.assertEqual(reads[0].sequence,self.repeat)

    def test_outside_mapping_of_unrelated_read_part_does_not_veto(self):
        extension='G'*600
        primary=self.rec(seq=self.seq+extension,cigar='1016M600S')
        other=self.rec(seq=extension,flag=2048,chrom='chrOther',cigar='1016H600M')
        reads,_,_=self.extract([primary,other])
        self.assertTrue(reads)

    def test_reverse_hardclips_map_to_original_coordinates(self):
        reads,_,_=self.extract([self.rec(flag=16,cigar='10H1016M20H')])
        self.assertEqual(reads[0].original_query_interval,(520,536))

    def test_zero_length_is_reported_not_invented(self):
        reads,_,decisions=self.extract([self.rec(seq=self.left+self.right,cigar='500M16D500M')])
        self.assertFalse(reads)
        self.assertEqual(decisions[0]['status'],'zero_length_candidate')

    def test_hardclips_are_never_joined(self):
        reads,_,_=self.extract([self.rec(seq=self.left,flag=2048,cigar='500M516H'),
                                self.rec(seq=self.right,flag=2048,pos=517,cigar='516H500M')])
        self.assertFalse(reads)

    def test_other_locus_secondary_competes_even_low_mapq(self):
        reads,_,decisions=self.extract([self.rec(),self.rec(flag=256,mapq=0,chrom='chrOther')])
        self.assertFalse(reads)
        self.assertEqual(decisions[0]['status'],'unresolved_molecule_conflict')

    def test_large_expansion_does_not_require_wider_bam_region(self):
        repeat='ACGT'*1250
        reads,_,_=self.extract([self.rec(seq=self.left+repeat+self.right,cigar='500M4984I516M')])
        self.assertEqual(reads[0].sequence,repeat)

    def test_short_template_retained_for_short_flanks(self):
        seq=self.left[-150:]+self.repeat+self.right[:150]
        reads,_,_=self.extract([self.rec(seq=seq,pos=351,cigar='316M')])
        self.assertEqual(reads[0].sequence,self.repeat)

    def test_regular_errors_destroying_fixed_seeds_still_match(self):
        flank=self.left[-150:]
        mutated=list(flank)
        for i in range(10,150,18):mutated[i]=next(c for c in 'ACGT' if c!=mutated[i])
        hits=a.flank_candidates(flank,''.join(mutated),.9)
        self.assertTrue(any(h['start']==0 and h['end']==150 for h in hits))

    def test_flank_insertion_exceeding_old_eight_base_limit(self):
        flank=self.left[-150:]
        seq=flank[:70]+'G'*20+flank[70:]
        hits=a.flank_candidates(flank,seq,.8)
        self.assertTrue(any(h['start']==0 and h['end']==170 for h in hits))

    def test_exact_hit_does_not_hide_nearby_alternative(self):
        flank=self.left[-150:]
        near=('A' if flank[0]!='A' else 'C')+flank[1:]
        hits=a.flank_candidates(flank,near+'T'*60+flank,.9)
        self.assertTrue(any(h['start']<=1 for h in hits))
        self.assertTrue(any(h['start']==210 for h in hits))

    def test_whole_bam_command_has_no_region_or_filters(self):
        class Result:
            returncode=0
            stdout=f'@SQ\tSN:chrT\tLN:{len(self.seq)}\n'+self.rec()
        with patch('vntranchors.subprocess.run',return_value=Result()) as run:
            a.extract_bam(Path('given.bam'),'S','HiFi',self.locus,self.ref,'samtools',self.profiles,20)
        self.assertEqual(run.call_args.args[0],['samtools','view','-h','given.bam'])

    def test_unknown_repeat_sequence_not_claimed_as_complete(self):
        reads,_,decisions=self.extract([self.rec(seq=self.left+'N'*16+self.right)])
        self.assertFalse(reads)
        self.assertEqual(decisions[0]['status'],'ambiguous_repeat_bases')

    def test_full_reference_profile_including_reverse_complement(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'ref.fa'
            # One extra copy of a reference probe must count as repeated even in reverse orientation.
            path.write_text('>chrT\n'+self.seq+'\n>other\n'+revcomp(self.left[-31:])+'\n')
            profiles,probes=a.reference_profiles(path,[self.locus],self.ref,[150])
            word=a.canonical(self.left[-31:])
            self.assertEqual(probes[word]['count'],2)
            self.assertEqual(profiles['T:L'][0]['grade'],'preferred')
            self.assertNotIn(119,profiles['T:L'][0]['unique_starts'])

    def test_overlapping_unique_windows_not_independent(self):
        self.assertEqual(a.nonoverlapping(range(31)),[0])
        self.assertEqual(a.nonoverlapping(range(32)),[0,31])

    def test_documented_alt_collapses_only_declared_blocks(self):
        blocks=[dict(alt_contig='chrAlt',alt_start=0,alt_end=500,
                     primary_contig='chrT',primary_start=0,primary_end=500,strand='+')]
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'ref.fa'
            path.write_text('>chrT\n'+self.seq+'\n>chrAlt\n'+self.seq+'\n')
            raw,_=a.reference_profiles(path,[self.locus],self.ref,[150])
            collapsed,probes=a.reference_profiles(path,[self.locus],self.ref,[150],blocks)
            self.assertEqual(raw['T:L'][0]['grade'],'nonunique')
            self.assertEqual(collapsed['T:L'][0]['grade'],'preferred')
            self.assertEqual(collapsed['T:L'][0]['raw_unique_count'],0)
            self.assertEqual(collapsed['T:R'][0]['grade'],'nonunique')
            self.assertEqual(probes[a.canonical(self.left[-31:])]['count'],2)
            self.assertEqual(probes[a.canonical(self.left[-31:])]['locus_count'],1)
            path.write_text(path.read_text()+'>true_paralog\n'+self.seq+'\n')
            third,_=a.reference_profiles(path,[self.locus],self.ref,[150],blocks)
            self.assertEqual(third['T:L'][0]['grade'],'nonunique')

    def test_equivalent_alt_secondary_is_not_a_competing_locus(self):
        blocks=[dict(alt_contig='chrOther',alt_start=x,alt_end=y,primary_contig='chrT',
                     primary_start=x,primary_end=y,strand='+') for x,y in [(0,500),(516,1016)]]
        for templates in self.profiles.values():
            for t in templates:t['reference_equivalences']=blocks
        reads,_,_=self.extract([self.rec(),self.rec(flag=256,chrom='chrOther',mapq=0)])
        self.assertEqual(len(deduplicate_reads(reads)[0]),1)

    def test_reverse_equivalence_coordinates(self):
        block=dict(alt_contig='alt',alt_start=10,alt_end=110,primary_contig='chrT',
                   primary_start=1000,primary_end=1100,strand='-')
        self.assertEqual(a.normalized_position('alt',10,[block],31),('chrT',1069,True))
        self.assertEqual(a.normalized_position('alt',109,[block]),('chrT',1000,True))
        self.assertEqual(a.normalized_position('alt',100,[block],31),('alt',100,False))

    def test_equivalence_rejects_nonidentical_sequences(self):
        class Reference:
            def fetch(ref,chrom,start,end):return ('A' if chrom=='alt' else 'C')*(end-start)
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'blocks.tsv'
            path.write_text('alt_contig\talt_start\talt_end\tprimary_contig\tprimary_start\tprimary_end\tstrand\tevidence_url\n'
                            'alt\t0\t100\tchrT\t0\t100\t+\thttps://example.org/placement\n')
            with self.assertRaisesRegex(ValueError,'base-identical'):
                a.load_equivalences(path,Reference())


if __name__=='__main__':unittest.main()
