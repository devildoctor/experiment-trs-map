from copy import deepcopy
from pathlib import Path
import sys
import unittest
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrdisplay_loss import (build_unit_graph,represented_graph,display_signature,
    compare_for_merge,canonical_edits,merge_paths,align_units)
from vntrdisplay import build_display,validate_display,add_support
from vntrpathstructure import restore_original
from vntrstructure import segment
from test_pathstructure import A,B,MODEL,feature,record

POLICY=dict(schema_version=2,version='toy.v2',model_version=MODEL['version'],unit_max_edits=2,terminal_max_edits=2)


class DisplayV2Tests(unittest.TestCase):
    def test_counts_parameterized_and_not_interpolated(self):
        records=[];fs={}
        for p,n in [('P0001',8),('P0002',12)]:
            s,f=feature('A'*n);r=record(p,s);r['role']='reference' if p=='P0001' else 'observation';records.append(r);fs[p]=f
        exact=build_unit_graph(records,fs,MODEL);reps={r['id']:r['id'] for r in records}
        d=build_display(exact,fs,reps,MODEL)
        self.assertEqual(len(d['nodes']),1)
        self.assertEqual(next(iter(d['nodes'].values()))['counts'],[8,12])
        self.assertEqual([d['paths'][p][0]['count'] for p in sorted(reps)],[8,12])
        validate_display(d,exact,{r['id']:r['sequence'] for r in records})

    def test_structure_changes_are_protected(self):
        for x,y in [('A'*8,'A'*9),('AABB','ABAB'),('AAA','ABA')]:
            a,fa=feature(x);b,fb=feature(y)
            self.assertEqual(compare_for_merge(a,b,fa,fb,POLICY,.05)[1],'protected_structure')

    def test_partial_completion_and_repeated_edits_are_allowed(self):
        target,tf=feature('A'*8)
        source,sf=feature('A'*8,{0:A[1:]})
        sf['tokens'][0]['status']='partial';sf['annotation']['annotation_ambiguity']=True
        es,_=compare_for_merge(source,target,sf,tf,POLICY,.01)
        self.assertEqual([e['op'] for e in es],['I']);self.assertEqual(es[0]['terminal_completion'],'left')
        self.assertEqual(restore_original(target,es),source)
        variant=A[:7]+'C'+A[8:]
        source,sf=feature('A'*8,{1:variant,3:variant},tied=True)
        es,_=compare_for_merge(source,target,sf,tf,POLICY,.03)
        self.assertEqual(len(es),2)

    def test_unit_budget_and_original_denominator(self):
        target,tf=feature('A'*8)
        source,sf=feature('A'*8,{1:'TTT'+A[3:]})
        self.assertEqual(compare_for_merge(source,target,sf,tf,POLICY,.05)[1],'unit_budget')
        source,sf=feature('A'*8,{0:A[1:]});sf['tokens'][0]['status']='partial'
        self.assertIsNone(compare_for_merge(source,target,sf,tf,POLICY,1/128)[0])

    def test_edit_placement_is_left_canonical_and_reversible(self):
        es=canonical_edits('CAAAAG','CAAAG')
        self.assertEqual(es[0]['ostart'],1)
        self.assertEqual(restore_original('CAAAG',es),'CAAAAG')

    def test_annotation_identity_does_not_block_literal_end_sharing(self):
        records=[];fs={}
        for p in ['P0001','P0002']:
            s,f=feature('ABA');records.append(record(p,s));fs[p]=f
        fs['P0002']['tokens'][-1].update(motif='RT',role='right_terminal')
        graph=build_unit_graph(records,fs,MODEL)
        self.assertEqual(graph['paths']['P0001'][-1],graph['paths']['P0002'][-1])

    def test_shared_detail_requires_a_batch_and_reference_stays_fixed(self):
        variant=A[:7]+'C'+A[8:];records=[];fs={}
        for p,labels,actual,count in [('P0001','A'*8,None,0),('P0002','A'*8,None,5),
            ('P0003','A'*8,{3:variant},3),('P0004','A'*7+'B',{3:variant},2),('P0005','A'*7+'B',None,5)]:
            s,f=feature(labels,actual);r=record(p,s,count);r['role']='reference' if p=='P0001' else 'observation';records.append(r);fs[p]=f
        graph=build_unit_graph(records,fs,MODEL)
        represented,reps,losses,audit=merge_paths(records,fs,MODEL,POLICY,.01,graph,build_display)
        reversed_graph=build_unit_graph(list(reversed(records)),fs,MODEL)
        reverse_result=merge_paths(list(reversed(records)),fs,MODEL,POLICY,.01,reversed_graph,build_display)
        self.assertEqual((represented,reps,losses,audit),reverse_result)
        self.assertEqual(reps['P0001'],'P0001');self.assertEqual(reps['P0003'],'P0002');self.assertEqual(reps['P0004'],'P0005')
        self.assertTrue(all(b['delta_C']>0 for b in audit['accepted_batches']))
        for loss in losses:self.assertEqual(restore_original(represented[loss['path']],loss['events']),next(r['sequence'] for r in records if r['id']==loss['path']))
        d=build_display(represented_graph(graph,reps),fs,reps,MODEL)
        validate_display(d,represented_graph(graph,reps),represented)

    def test_no_graph_benefit_does_not_spend_loss(self):
        s,fs=feature('A'*8,{3:A[:7]+'C'+A[8:]});t,ft=feature('A'*8)
        rs=[record('P0001',t),record('P0002',t,5),record('P0003',s)]
        rs[0]['role']='reference';features={'P0001':ft,'P0002':ft,'P0003':fs}
        graph=build_unit_graph(rs,features,MODEL)
        def constant(*args,**kwargs):return {'metrics':{'C':10}}
        _,reps,losses,_=merge_paths(rs,features,MODEL,POLICY,.05,graph,constant)
        self.assertFalse(losses);self.assertEqual(reps['P0003'],'P0003')

    def test_global_gap_objective_and_input_order(self):
        _,a=feature('A'*8);_,b=feature('A'*12)
        steps,audit=align_units(a['tokens'],b['tokens'])
        self.assertTrue(audit['had_primary_ties']);self.assertEqual(audit['gap_runs'],1)
        self.assertEqual(steps[:4],[(None,0),(None,1),(None,2),(None,3)])
        dictionary=[dict(id=k,base=k,sequence=v,level='base',match_orientations=['+']) for k,v in [('X','ACGA'),('Y','ACGC'),('Z','TTAA')]]
        x,au=segment(('TTAAACGTTTAA',dictionary,'HiFi'),True,select_representative=True)
        y,_=segment(('TTAAACGTTTAA',list(reversed(dictionary)),'HiFi'),True,select_representative=True)
        self.assertEqual([(t['motif'],t['start'],t['end']) for t in x],[(t['motif'],t['start'],t['end']) for t in y])
        self.assertEqual(au['ambiguity_intervals'],[[4,8]])

    def test_export_reconstruction_local_support_and_selection(self):
        from vntrdisplay_export import export_budget,select_budget
        s,sf=feature('A'*8,{0:A[1:]});sf['tokens'][0]['status']='partial'
        t,tf=feature('A'*8);rs=[record('P0001',t,0),record('P0002',t,5),record('P0003',s,3)]
        rs[0]['role']='reference';fs={'P0001':tf,'P0002':tf,'P0003':sf};graph=build_unit_graph(rs,fs,MODEL)
        with tempfile.TemporaryDirectory() as tmp:
            report=export_budget(Path(tmp),rs,fs,MODEL,POLICY,.01,graph,render=False)
            self.assertTrue(report['validation']['txt_roundtrip']);self.assertEqual(report['insertions'],1)
            self.assertIn('original_sequence\trepresented_sequence',(Path(tmp)/'D1S80.path_sequences.txt').read_text())
            import json
            display=json.loads((Path(tmp)/'display_graph.json').read_text())
            self.assertTrue(any(e['raw_exact_molecules']==e['route_molecules'] for e in display['edges']) or not display['edges'])
        reports=[dict(C=c,total_edits_bp=e,loss_budget_fraction=f) for c,e,f in [(100,0,0),(80,2,.01),(79,8,.03),(79,10,.05)]]
        self.assertEqual(select_budget(reports,dict(selection_gain_fraction=.9))['selected_fraction'],.01)


if __name__=='__main__':unittest.main()
