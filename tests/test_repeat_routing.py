"""STR keeps its CTG/TTG decomposition when using the shared runner."""
import json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrgraph import Locus
from vntrpartition import repeat_type
from vntrpathstructure import annotate_paths
from vntrpathstructure_export import export_structure
from run_efficient_family import run
from run_hmm_parameters import rebuild
from test_pathstructure import record

ROOT=Path(__file__).resolve().parents[1]
MODEL_PATH=ROOT/'catalogs/ATXN2.structure.v1.json'
LOCUS=Locus('ATXN2','chr12',111598950,111599019,(),boundary_status='confirmed')


def records():
    two='CTG'*8+'TTG'+'CTG'*4+'TTG'+'CTG'*8
    one='CTG'*8+'TTG'+'CTG'*13
    rs=[record('P0001',two),record('P0002',two),record('P0003',one)]
    rs[0].update(role='reference',molecules=[],hifi_molecules=[],samples={},observations=[],
                 reference_interval=dict(chrom='chr12',start=111598950,end=111599019))
    for r in rs[1:]:r['role']='observation'
    return rs


class RepeatRoutingTests(unittest.TestCase):
    def test_existing_models_and_explicit_mixed_model_type(self):
        for name,expected in [('ATXN2','STR'),('D1S80','VNTR')]:
            self.assertEqual(repeat_type(json.loads((ROOT/'catalogs'/(name+'.structure.v1.json')).read_text())),expected)
        mixed=dict(motifs=[dict(sequence='CTG'),dict(sequence='ACGTACGT')])
        with self.assertRaisesRegex(ValueError,'Mixed'):repeat_type(mixed)
        self.assertEqual(repeat_type(dict(mixed,repeat_type='STR')),'STR')

    def test_shared_runner_restores_legacy_interruptions_without_vntr_modules(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);source=folder/'input.json';source.write_text(json.dumps(records()))
            with (patch('run_efficient_family.annotate',side_effect=AssertionError('VNTR annotation used')),
                 patch('run_efficient_family.select_motifs',side_effect=AssertionError('Dictionary compressed')),
                 patch('run_efficient_family.select_families',side_effect=AssertionError('Family merge used'))):
                result=run(source,MODEL_PATH,folder/'out',q=.1,fraction=.03,render=False)
            output=folder/'out';self.assertEqual(result[0]['motifs'],2)
            self.assertEqual(result[0]['total_edits_bp'],0)
            self.assertEqual(json.loads((output/'processing.json').read_text())['decoder'],'seed_dp')
            features={r['path']:r for r in map(json.loads,(output/'str_legacy/path_features.jsonl').read_text().splitlines())}
            self.assertEqual([r['unit'] for r in features['P0002']['interruptions']],[9,14])
            self.assertEqual([r['unit'] for r in features['P0003']['interruptions']],[9])
            self.assertNotEqual(features['P0002']['signature'],features['P0003']['signature'])
            self.assertEqual((output/'str_legacy/paths.original.fa').read_bytes(),(output/'str_legacy/paths.represented.fa').read_bytes())
            self.assertFalse((output/'locus_model.efficient.json').exists())

    def test_main_export_with_hmm_policy_routes_str_to_legacy(self):
        with tempfile.TemporaryDirectory() as temp:
            with patch('vntrhmm.segment_hmm',side_effect=AssertionError('HMM used for STR')):
                result=export_structure(Path(temp),LOCUS,records(),[],dict(backbone_path='P0001'),MODEL_PATH,
                    ROOT/'catalogs/D1S80.structure-policy.hmm.json',.01,0,None)
            self.assertEqual(result['loss_bp'],0)
            self.assertEqual(json.loads((Path(temp)/'structure_policy.json').read_text())['schema_version'],1)

    def test_explicit_hmm_and_skip_rebuild_do_not_process_str(self):
        with self.assertRaisesRegex(ValueError,'VNTR-only'):
            annotate_paths(records(),MODEL_PATH,LOCUS,decoder='hmm_single')
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);(folder/'source_members.json').write_text(json.dumps(records()))
            (folder/'locus_model.json').write_bytes(MODEL_PATH.read_bytes())
            with self.assertRaisesRegex(ValueError,'VNTR-only'):rebuild(folder,folder/'out',False)


if __name__=='__main__':unittest.main()
