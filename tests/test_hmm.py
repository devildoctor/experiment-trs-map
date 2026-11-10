"""Sequence-level checks for HMM segmentation and graph integration."""
import sys,unittest,json,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrhmm import segment_hmm,decode
from vntrpathstructure import annotate_paths
from vntrlocal import classify,assign_slots,plan_edits
from vntrpathstructure import restore_original

U1='CCTGGTGGTCCTCCTT'
U2='CCTGGTGGTTCTCCTT'


def annotate(sequence,motifs=(U1,U2)):
    dictionary=[dict(id='U%02d'%(i+1),base='U%02d'%(i+1),sequence=s,level='base') for i,s in enumerate(motifs)]
    return segment_hmm((sequence,dictionary,'HiFi'),return_details=True)


class HmmTests(unittest.TestCase):
    def test_exact_units_and_fixed_types(self):
        ts,a=annotate(U1*2+U2*3)
        self.assertEqual([t['motif'] for t in ts],['U01']*2+['U02']*3)
        self.assertTrue(all(t['edit_distance']==0 for t in ts))
        self.assertNotIn('optimal_dag',a)
        self.assertEqual(a['alternative'],[])

    def test_skip_between_motif_runs(self):
        sequence=U1*2+'A'*40+U2*2
        ts,_=annotate(sequence)
        self.assertEqual(''.join(t['exact_bases'] for t in ts),sequence)
        self.assertEqual([t['exact_bases'] for t in ts if t['status']=='unknown'],['A'*40])

    def test_indel_is_retained_inside_unit(self):
        variant='CCTGGTGGTCCCTCCTT'
        ts,_=annotate(U1+variant+U2)
        self.assertEqual([t['exact_bases'] for t in ts],[U1,variant,U2])
        self.assertEqual(ts[1]['edit_distance'],1)
        self.assertEqual(ts[1]['motif'],'U01')

    def test_deleted_bases_are_not_synthesized(self):
        variant=U1[:9]+U1[10:]
        ts,_=annotate(U1+variant+U2)
        self.assertEqual(''.join(t['exact_bases'] for t in ts),U1+variant+U2)
        self.assertTrue(any(t['status']=='unit' and len(t['exact_bases'])==15 for t in ts))

    def test_unknown_only_and_empty(self):
        ts,_=annotate('A'*20)
        self.assertEqual([(t['status'],t['exact_bases']) for t in ts],[('unknown','A'*20)])
        self.assertEqual(annotate('')[0],[])

    def test_tie_deterministic_one_result(self):
        ts,a=annotate(U1,(U1,U1))
        self.assertEqual([t['motif'] for t in ts],['U01'])
        self.assertEqual(a['ambiguity_scope'],'not_enumerated_single_trace')

    def test_real_short_unknowns_preserve_coordinates(self):
        for fragment in ['CTCAGGTCCTCCTT','CCTGGTGGTCCT']:
            sequence=U1+fragment+U2
            ts,_=annotate(sequence)
            self.assertEqual(''.join(t['exact_bases'] for t in ts),sequence)
            self.assertEqual([t['start'] for t in ts],[0]+[t['end'] for t in ts[:-1]])

    def test_frozen_terminal_and_loss_pipeline(self):
        root=Path(__file__).resolve().parents[1]
        model_path=root/'catalogs/D1S80.structure.v1.json'
        terminal=json.loads(model_path.read_text())['terminal_templates'][0]['sequence']
        sequences=[U1*8,U1*7+U2,U1*8+terminal]
        records=[dict(id='P%04d'%(i+1),sequence=s,role='reference' if i==0 else 'observation',
                      samples=['LCL5'],observations=[],molecules=[] if i==0 else ['m'+str(i)],
                      hifi_molecules=[] if i==0 else ['m'+str(i)]) for i,s in enumerate(sequences)]
        from vntrgraph import Locus
        locus=Locus('D1S80','chr1',2459279,2459743,(),boundary_status='confirmed')
        model,_,features=annotate_paths(records,model_path,locus,decoder='hmm_single')
        self.assertEqual(features['P0003']['tokens'][-1]['role'],'right_terminal')
        tokens={p:classify(f['tokens'],[]) for p,f in features.items()}
        slots={p:assign_slots(ts,tokens['P0001'])[0] for p,ts in tokens.items()}
        ts,views,represented,losses,_,_=plan_edits(records,tokens,slots,model,.01)
        self.assertEqual(represented['P0001'],sequences[0])
        by_id={r['id']:r for r in records}
        for loss in losses:
            original=by_id[loss['path']]['sequence']
            self.assertLessEqual(loss['loss_bp'],int(len(original)*.01))
            self.assertEqual(restore_original(represented[loss['path']],loss['events']),original)


if __name__=='__main__':unittest.main()
