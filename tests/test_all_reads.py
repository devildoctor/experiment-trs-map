"""Regression contracts for complete read inventory and partial read graphs."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrgraph import Locus
from vntrreads import cigar_interval, extract_record, annotate, write_outputs


class ReadViewTests(unittest.TestCase):
    def fields(self,sequence,cigar,pos=91,flag=0,mapq=60):
        return ['read1',str(flag),'chr1',str(pos),str(mapq),cigar,'*','0','0',sequence,'I'*len(sequence)]

    def test_large_internal_and_boundary_insertions_are_preserved(self):
        query, ref, flags, overlap=cigar_interval(90,'10M20I10M30I10M',100,110,80)
        self.assertEqual(query,(10,70))
        self.assertEqual(ref,(100,110))
        self.assertTrue(overlap)

    def test_single_left_anchor_is_displayed_and_partial(self):
        locus=Locus('X','chr1',100,120,())
        seq='AACCGGTTAA'+'ACGTACGTAC'
        with patch('vntrreads.near_anchor',side_effect=[(0,10,1),None]):
            r=extract_record(self.fields(seq,'20M'),locus,'AACCGGTTAA','TTGGAACCGG')
        self.assertEqual(r['status'],'displayed')
        self.assertEqual(r['anchor'],'left_only')
        self.assertEqual(r['coverage'],'left_partial')
        self.assertEqual(r['sequence'],'ACGTACGTAC')

    def test_single_right_anchor(self):
        locus=Locus('X','chr1',100,120,())
        with patch('vntrreads.near_anchor',side_effect=[None,(10,20,1)]):
            r=extract_record(self.fields('ACGTACGTAC'+'TTGGAACCGG','20M',111),locus,'AACCGGTTAA','TTGGAACCGG')
        self.assertEqual(r['anchor'],'right_only')
        self.assertEqual(r['coverage'],'right_partial')
        self.assertEqual(r['sequence'],'ACGTACGTAC')

    def test_no_anchor_low_mapq_duplicate_qc_fail_retained(self):
        locus=Locus('X','chr1',100,120,())
        with patch('vntrreads.near_anchor',return_value=None):
            r=extract_record(self.fields('ACGTACGT','8M',105,1024|512,0),locus,'AAA','TTT')
        self.assertEqual(r['status'],'displayed')
        self.assertEqual(r['anchor'],'cigar_only')
        self.assertEqual(r['coverage'],'internal_partial')
        self.assertIn('low_mapq',r['flags'])
        self.assertIn('duplicate_flag',r['flags'])

    def test_negative_strand_stored_sequence_not_reversed(self):
        locus=Locus('X','chr1',100,108,())
        with patch('vntrreads.near_anchor',return_value=None):
            a=extract_record(self.fields('AAGGCTTA','8M',101,0),locus,'','')
            b=extract_record(self.fields('AAGGCTTA','8M',101,16),locus,'','')
        self.assertEqual(a['sequence'],b['sequence'])

    def test_reference_skip_never_called_full(self):
        locus=Locus('X','chr1',100,120,())
        with patch('vntrreads.near_anchor',return_value=None):
            r=extract_record(self.fields('AACCGGTTAA','5M10N5M',101),locus,'','')
        self.assertEqual(r['coverage'],'split_reference_coverage')

    def test_soft_clip_not_invented_as_repeat(self):
        query,ref,flags,overlap=cigar_interval(100,'5S8M12S',100,120,25)
        self.assertEqual(query,(5,13))
        self.assertIn('soft_clipped_sequence_not_assigned_to_repeat',flags)

    def test_flank_only_and_unavailable_retained_as_status(self):
        locus=Locus('X','chr1',100,120,())
        r=extract_record(self.fields('AAAA','4M',11),locus,'','')
        self.assertEqual(r['status'],'flank_only')
        r=extract_record(self.fields('*','20M',101),locus,'','')
        self.assertEqual(r['status'],'sequence_unavailable')

    def test_observed_bases_survive_errors_and_unknown(self):
        for sequence in ('ACGT'*4,'ACGTACATACGT','NNNACGTACGTNN','AC','GATTACA'):
            tokens=annotate((sequence,('ACGT',),'ONT'))
            self.assertEqual(''.join(t['exact_bases'] for t in tokens),sequence)
            self.assertTrue(all(t['exact_bases']==sequence[t['start']:t['end']] for t in tokens))

    def test_outputs_include_empty_rows_and_each_fragment_with_same_node_ids(self):
        locus=Locus('X','chr1',100,112,())
        with patch('vntrreads.near_anchor',return_value=None):
            f=extract_record(self.fields('ACGT'*3,'12M',101),locus,'','')
        f.update(id='R1_F1',tokens=annotate((f['sequence'],('ACGT',),'HiFi')))
        rows=[dict(id='R1',read_id='read<1>',sample='S',platform='HiFi',bam='b.bam',record_count=1,records=[],fragments=[f],status='displayed'),
              dict(id='R2',read_id='read2',sample='S',platform='HiFi',bam='b.bam',record_count=1,records=[],fragments=[],status='flank_only')]
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)
            audit=dict(read_rows=2)
            write_outputs(dest,locus,rows,audit)
            gfa=(dest/'X.reads.gfa').read_text(encoding='utf-8')
            svg=(dest/'svg'/'X.reads.0001.svg').read_text(encoding='utf-8')
            ET.fromstring(svg)
            self.assertIn('read2',svg)
            self.assertIn('R1_F1_N1',gfa)
            self.assertIn('R1_F1_N1',svg)
            html=(dest/'X.reads.html').read_text(encoding='utf-8')
            self.assertNotIn('__DATA__',html)
            self.assertIn('read\\u003c1\\u003e',html)
            payload=json.loads((dest/'X.reads.json').read_text(encoding='utf-8'))
            self.assertEqual(len(payload['rows']),2)
            nodes={line.split('\t')[1]:line.split('\t')[2] for line in gfa.splitlines() if line.startswith('S\t')}
            paths=[line.split('\t')[2] for line in gfa.splitlines() if line.startswith('P\t')]
            self.assertEqual(''.join(nodes[n[:-1]] for n in paths[0].split(',')),'ACGT'*3)


if __name__=='__main__':
    unittest.main()
