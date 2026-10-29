"""命名空间去重与结构备选 token 歧义的回归验证。"""
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))

from vntrgraph import Locus, load_source_namespaces
from vntrreads import read_bam, collect_reads, row_molecule_key
from vntrstructure import seed_dictionary, segment, summarize, learn_variants, emit
from vntrsummary import structural_summary


class Reference:
    def fetch(self, chrom, start, end):return 'A'*(end-start)
    def contig_length(self, chrom):return 1000


class SamProcess:
    def __init__(self, sam):self.stdout = io.StringIO(sam)
    def wait(self):return 0


def sam(libraries, groups):
    header = ['@SQ\tSN:chr1\tLN:1000']
    header += ['@RG\tID:'+group+('\tLB:'+library if library is not None else '') for group, library in libraries]
    records = ['\t'.join(['same','0','chr1','101','60','8M','*','0','0','ACGTACGT','IIIIIIII']+
                         (['RG:Z:'+group] if group else [])) for group in groups]
    return '\n'.join(header+records)+'\n'


def annotated_row(row_id, namespace=None, sequence='CTG'*4+'TTG'+'CTG'*4, platform='HiFi'):
    dictionary = seed_dictionary(['CTG','TTG'])
    fragment = dict(id=row_id+'_F1',sequence=sequence,qstart=10,qend=10+len(sequence),
                    coverage='spans_candidate_interval',anchor='both',anchor_scores=[1,1],
                    mapq=60,alignment_type='primary',flags=[],tokens=segment((sequence,dictionary,platform)))
    row = dict(id=row_id,read_id='same',sample='S1',platform=platform,bam=row_id+'.bam',
               fragments=[fragment],record_count=1,records=[],status='displayed')
    if namespace is not None:row['namespace']=namespace
    return row


class MoleculeIdentityTests(unittest.TestCase):
    def read(self, text, namespace=None):
        output = io.StringIO()
        with patch('vntrreads.subprocess.Popen',return_value=SamProcess(text)), patch('vntrreads.near_anchor',return_value=None):
            rows, audit = read_bam(Path('input.bam'),'S1','HiFi',Locus('X','chr1',100,108,()),
                                   Reference(),'samtools',100,output,namespace)
        return rows, audit, [json.loads(line) for line in output.getvalue().splitlines()]

    def test_same_qname_different_libraries_are_distinct(self):
        rows, audit, records = self.read(sam([('g1','lib1'),('g2','lib2')],['g1','g2']))
        self.assertEqual(audit['read_rows'],2)
        self.assertEqual({tuple(json.loads(row_molecule_key(r))) for r in rows},{('S1','lib1','same'),('S1','lib2','same')})
        self.assertEqual(len({r['id'] for r in rows}),2)
        self.assertEqual([r['library'] for r in records],['lib1','lib2'])
        self.assertTrue(all(r['molecule_key']==row_molecule_key(r) for r in rows))

    def test_different_read_groups_same_library_are_one_molecule(self):
        rows, audit, records = self.read(sam([('g1','lib'),('g2','lib')],['g1','g2']))
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['record_count'],2)
        self.assertEqual(len(rows[0]['fragments']),2)
        self.assertEqual([r['read_group'] for r in records],['g1','g2'])

    def test_explicit_namespace_including_empty_overrides_library(self):
        for override in ('merged',''):
            with self.subTest(namespace=override):
                rows, _, records = self.read(sam([('g1','lib1'),('g2','lib2')],['g1','g2']),override)
                self.assertEqual(len(rows),1)
                self.assertEqual(rows[0]['namespace'],override)
                self.assertEqual({r['namespace'] for r in records},{override})

    def test_missing_library_does_not_invent_bam_or_platform_namespace(self):
        rows, _, _ = self.read(sam([('g1',None)],['g1','']))
        self.assertEqual(len(rows),1)
        self.assertEqual(json.loads(row_molecule_key(rows[0])),['S1','','same'])

    def test_manifest_relative_paths_empty_values_and_duplicate_rejection(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as directory:
            path=Path(directory)/'sources.tsv'
            path.write_text('bam\tnamespace\na.bam\tchosen\nb.bam\t\n',encoding='utf-8')
            self.assertEqual(load_source_namespaces(path),{str((path.parent/'a.bam').resolve()):'chosen',str((path.parent/'b.bam').resolve()):''})
            path.write_text('bam\tnamespace\na.bam\tchosen\n./a.bam\tother\n',encoding='utf-8')
            with self.assertRaisesRegex(ValueError,'Duplicate BAM'):load_source_namespaces(path)

    def test_collect_reads_preserves_rows_but_uses_manifest_identity(self):
        bams=[('S1','HiFi',Path('first.bam')),('S1','ONT',Path('second.bam'))]
        mappings={str(bam.resolve()):'shared' for _,_,bam in bams}
        processes=[SamProcess(sam([('g','lib'+str(i))],['g'])) for i in range(2)]
        with patch('vntrreads.subprocess.Popen',side_effect=processes),patch('vntrreads.near_anchor',return_value=None):
            rows, _ = collect_reads(bams,Locus('X','chr1',100,108,()),Reference(),'samtools',100,io.StringIO(),mappings)
        self.assertEqual(len(rows),2)
        self.assertEqual(len({row_molecule_key(r) for r in rows}),1)
        self.assertEqual(len({r['id'] for r in rows}),2)

    def test_structure_nodes_edges_and_classes_share_namespace_identity(self):
        dictionary=seed_dictionary(['CTG','TTG'])
        rows=[annotated_row('a','lib1'),annotated_row('b','lib2'),annotated_row('copy','lib1',platform='ONT')]
        graph=summarize(rows,dictionary);summary=structural_summary(rows,minimum=2)
        self.assertTrue(all(n['support']==2 for n in graph['structure_nodes']))
        self.assertTrue(all(e['support']==2 for e in graph['edges']))
        self.assertEqual(summary['classes'][0]['support'],2)
        self.assertEqual(summary['classes'][0]['hifi_support'],2)
        self.assertEqual(summary['classes'][0]['status'],'supported')
        self.assertEqual(summary['audit']['eligible_hifi_by_sample'],{'S1':2})
        self.assertTrue(all(m['support']==2 for m in summary['modules']))
        self.assertTrue(all(json.loads(a['molecule'])==['S1',a['namespace'],'same'] for a in summary['assignments']))

    def test_conflicting_paths_only_conflict_within_same_namespace(self):
        rows=[annotated_row('a','lib1'),annotated_row('b','lib2',sequence='CTG'*7+'TTG'+'CTG'*4)]
        summary=structural_summary(rows)
        self.assertEqual(summary['audit']['conflicting_molecules'],0)
        self.assertEqual(summary['classes'][0]['support'],2)
        rows[1]['namespace']='lib1'
        summary=structural_summary(rows)
        self.assertEqual(summary['audit']['conflicting_molecules'],1)
        self.assertFalse(summary['classes'])

    def test_legacy_rows_and_explicit_empty_namespace_deduplicate(self):
        rows=[annotated_row('legacy'),annotated_row('new','')]
        self.assertEqual(row_molecule_key(rows[0]),row_molecule_key(rows[1]))
        summary=structural_summary(rows)
        self.assertEqual(summary['classes'][0]['support'],1)
        self.assertTrue(all(m['support']==1 for m in summary['modules']))

    def test_variant_learning_counts_namespaces_not_display_rows(self):
        motif='AACCGTACGATGAACC';variant='AACCGTTCGATGAACC';dictionary=seed_dictionary([motif])
        rows=[annotated_row(str(i),'lib'+str(i),sequence=variant*2) for i in range(3)]
        for r in rows:r['fragments'][0]['tokens']=segment((variant*2,dictionary,'HiFi'))
        learned, _ = learn_variants(rows,dictionary)
        self.assertEqual(len(learned),2)
        self.assertEqual(learned[1]['support'],3)
        for r in rows:r['namespace']='shared'
        learned, _ = learn_variants(rows,dictionary)
        self.assertEqual(len(learned),1)

    def test_structure_exports_include_molecular_identity(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as directory:
            destination=Path(directory);rows=[annotated_row('a','lib1'),annotated_row('b','lib2')]
            sequence=rows[0]['fragments'][0]['sequence']
            reference=dict(sequence=sequence,tokens=segment((sequence,seed_dictionary(['CTG','TTG']),'HiFi')))
            audit=emit(destination,Locus(destination.name,'chr1',100,127,('CTG','TTG')),rows,seed_dictionary(['CTG','TTG']),dict(read_rows=2),reference)
            self.assertEqual(audit['independent_molecules'],2)
            for filename, field in [('reads.tsv','molecule_key'),('read_nodes.tsv','molecule_key'),('read_assignments.tsv','molecule')]:
                with (destination/filename).open(encoding='utf-8',newline='') as handle:exported=list(csv.DictReader(handle,delimiter='\t'))
                self.assertEqual({r['namespace'] for r in exported},{'lib1','lib2'})
                self.assertTrue(all(json.loads(r[field])==['S1',r['namespace'],'same'] for r in exported))
            for verifier in ('verify_structure_results.py','verify_summary_results.py'):
                result=subprocess.run([sys.executable,str(ROOT/'tests'/verifier),str(destination)],capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)


class StructureAmbiguityTests(unittest.TestCase):
    def test_both_tied_paths_mark_only_different_tokens(self):
        sequence='TTGCACGT'
        primary, audit=segment((sequence,seed_dictionary(['TTGC','ACAT','ACCT']),'HiFi'),return_details=True)
        self.assertTrue(audit['annotation_ambiguity'])
        self.assertEqual(audit['score_gap'],0)
        for tokens in (primary,audit['alternative']):
            self.assertEqual([t['annotation_ambiguity'] for t in tokens],[False,True])
            self.assertEqual(''.join(t['exact_bases'] for t in tokens),sequence)

    def test_unequal_paths_do_not_mark_ambiguity(self):
        primary,audit=segment(('ACAT',seed_dictionary(['ACAT','ACCT']),'HiFi'),return_details=True)
        self.assertFalse(audit['annotation_ambiguity'])
        self.assertGreater(audit['score_gap'],0)
        self.assertTrue(audit['alternative'])
        self.assertFalse(any(t['annotation_ambiguity'] for t in primary+audit['alternative']))


if __name__=='__main__':unittest.main()
