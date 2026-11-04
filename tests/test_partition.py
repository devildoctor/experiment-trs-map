import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrgraph import read_catalog
from vntrstructure import segment, annotate_rows, summarize
from vntrpartition import load_model, unit_rows, check_partition, expression, export_core
from vntrsummary import is_primary
from test_structure import row

ROOT=Path(__file__).resolve().parents[1]


class PartitionTests(unittest.TestCase):
    def model(self, name='ATXN2'):
        locus=next(l for l in read_catalog(ROOT/'catalogs/grch38_partition_pilot_v1.bed') if l.name==name)
        model,dictionary=load_model(ROOT/f'catalogs/{name}.partition.v1.json',locus)
        return locus,model,dictionary

    def test_model_preserves_explicit_bases_and_lengths(self):
        _,_,atxn=self.model();_,_,abca=self.model('ABCA7')
        self.assertEqual([m['sequence'] for m in atxn],['CTG','TTG'])
        self.assertEqual([len(m['sequence']) for m in abca],[25,26])
        self.assertTrue(all(m['match_orientations']==['+'] for m in atxn+abca))

    def test_model_rejects_coordinate_and_length_drift(self):
        locus,model,_=self.model()
        for field,value in [('start',locus.start+1),('orientation','reverse')]:
            changed=dict(model,**{field:value})
            with tempfile.TemporaryDirectory() as tmp:
                p=Path(tmp)/'model.json';p.write_text(json.dumps(changed))
                with self.assertRaises(ValueError):load_model(p,locus)
        model['motifs'][0]['length']=4
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'model.json';p.write_text(json.dumps(model))
            with self.assertRaises(ValueError):load_model(p,locus)

    def test_added_samples_do_not_change_old_partition_or_model(self):
        _,_,dictionary=self.model();before=copy.deepcopy(dictionary)
        old=row('old','CTG'*4+'TTG'+'CTG'*3,dictionary)
        annotate_rows([old],dictionary,1);tokens=copy.deepcopy(old['fragments'][0]['tokens'])
        rows=[old,row('new','CTG'*3+'CCG'+'CTG'*2,dictionary)]
        annotate_rows(rows,dictionary,1)
        self.assertEqual(dictionary,before)
        self.assertEqual(old['fragments'][0]['tokens'],tokens)

    def test_interruption_position_is_retained_even_at_same_length(self):
        _,_,dictionary=self.model()
        a=unit_rows(segment(('CTG'*2+'TTG'+'CTG'*2,dictionary,'HiFi')))
        b=unit_rows(segment(('CTG'+'TTG'+'CTG'*3,dictionary,'HiFi')))
        self.assertNotEqual(expression(a),expression(b))
        self.assertEqual([u['start'] for u in a if u['motif']=='M2'],[6])
        self.assertEqual([u['start'] for u in b if u['motif']=='M2'],[3])
        self.assertEqual(expression(a,True),'(CAG)×2 – CAA – (CAG)×2')

    def test_short_substitutions_are_not_lost_before_dp(self):
        _,_,dictionary=self.model();dictionary=dictionary[:1]
        for i,base in enumerate('CTG'):
            for alternative in set('ACGT')-{base}:
                seq='CTG'[:i]+alternative+'CTG'[i+1:]
                tokens=segment((seq,dictionary,'HiFi'))
                self.assertEqual([(t['status'],t['edit_distance']) for t in tokens],[('unit',1)],seq)
                self.assertEqual(tokens[0]['exact_bases'],seq)

    def test_observed_indel_length_is_not_template_length(self):
        _,_,dictionary=self.model('ABCA7');motif=dictionary[0]['sequence']
        seq=motif+motif[:13]+motif[14:]+motif
        units=unit_rows(segment((seq,dictionary[:1],'HiFi')))
        check_partition(seq,units)
        self.assertTrue(any(u['template_length']==25 and u['observed_length']==24 for u in units))

    def test_mapq_unavailable_cannot_supply_primary_support(self):
        _,_,dictionary=self.model();r=row('r','CTG'*4,dictionary)
        r['fragments'][0]['mapq']=255
        self.assertFalse(is_primary(r['fragments'][0]))
        self.assertEqual(summarize([r],dictionary)['structure_nodes'][0]['primary_support'],0)

    def test_core_outputs_reconstruct_and_count_molecules_once(self):
        locus,model,dictionary=self.model();seq='CTG'*3+'TTG'+'CTG'*2
        rows=[row('same',seq,dictionary,copies=2),row('same',seq,dictionary,bam='copy.bam')]
        annotate_rows(rows,dictionary,1)
        audit=dict(dictionary_mode='frozen_apply',model_version=model['version'],bam_count=2,
                   mapping_context='supplied_BAM_only_not_verified_complete')
        reference=dict(sequence=seq,tokens=segment((seq,dictionary,'HiFi')))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);meta=export_core(p,locus,rows,dictionary,reference,audit)
            self.assertEqual(meta['molecules'],1)
            self.assertEqual(meta['reconstruction_checks'],4)
            text=(p/'index.html').read_text(encoding='utf-8')
            self.assertNotIn('__DATA__',text)
            payload=json.loads(text.split('<script id="data" type="application/json">')[1].split('</script>')[0])
            self.assertTrue(all(v['molecules']==1 for v in payload['observed_units']))
            self.assertEqual(meta['reference_expression'],'(CTG)×3 – TTG – (CTG)×2')

    def test_reconstruction_rejects_missing_or_modified_bases(self):
        _,_,dictionary=self.model();units=unit_rows(segment(('CTG'*3,dictionary,'HiFi')))
        with self.assertRaises(ValueError):check_partition('CTG'*3,units[:-1])
        units[0]['bases']='CAA'
        with self.assertRaises(ValueError):check_partition('CTG'*3,units)

    def test_reference_alternative_is_visible(self):
        locus,model,dictionary=self.model('ABCA7');seq='C'*10
        tokens,annotation=segment((seq,dictionary,'HiFi'),return_details=True)
        # Exercise the export contract independently of this sequence's score.
        annotation['annotation_ambiguity']=True
        annotation['alternative']=copy.deepcopy(tokens)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            audit=dict(dictionary_mode='frozen_apply',model_version=model['version'],bam_count=0,
                       mapping_context='supplied_BAM_only_not_verified_complete')
            reference=dict(sequence=seq,tokens=tokens,annotation=annotation)
            meta=export_core(p,locus,[],dictionary,reference,audit)
            self.assertTrue(meta['reference_ambiguity'])
            text=(p/'index.html').read_text(encoding='utf-8')
            payload=json.loads(text.split('<script id="data" type="application/json">')[1].split('</script>')[0])
            self.assertTrue(payload['observations'][0]['alternatives'])


if __name__=='__main__':unittest.main()
