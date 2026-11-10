import sys,unittest
from pathlib import Path
from copy import deepcopy
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrlocal import classify,path_view,topology,plan_edits,replace_region,assign_slots
from vntrlocal_export import literal_graph,materialize,validate_parameters
from vntrstructure import segment

A='CCTGGTGGTCCTCCTT';B='CCTGGTGGTTCTCCTT';C='GCCGGTGGTCCTCCTT'


def tokens(parts):
    out=[];pos=0
    for motif,seq in parts:
        status='unknown' if motif=='U' else 'unit';template={'A':A,'B':B,'C':C,'U':''}[motif]
        out.append(dict(start=pos,end=pos+len(seq),exact_bases=seq,motif=motif,base=motif,
            template=template,status=status,orientation='+' if motif!='U' else '?',
            edit_distance=0 if seq==template else 1,annotation_ambiguity=False))
        pos+=len(seq)
    return classify(out,[])


def record(p,ts):
    s=''.join(t['exact_bases'] for t in ts)
    return dict(id=p,sequence=s,role='reference' if p=='P0001' else 'observation',
        molecules=[] if p=='P0001' else ['M'+p],hifi_molecules=[] if p=='P0001' else ['M'+p])


class LocalTests(unittest.TestCase):
    def test_A_ordered_variants_B_counts(self):
        ref=tokens([('A',A),('A',A),('B',B)])
        query=tokens([('A',A[:9]+'C'+A[9:]),('A',A),('A',A),('B',B)])
        ss,_=assign_slots(query,ref);vs=path_view(query,ss,ref)
        self.assertEqual(vs[0]['count'],3)
        self.assertEqual(''.join(vs[0]['variant_vector']),''.join(t['exact_bases'] for t in query[:3]))
        self.assertEqual(vs[0]['key'],path_view(ref,assign_slots(ref,ref)[0],ref)[0]['key'])

    def test_D_variable_lengths_and_no_rewrite(self):
        ref=tokens([('A',A),('A',A),('B',B)])
        q1=tokens([('A',A),('U','CTCAGGTCCTCCTT'),('B',B)])
        q2=tokens([('A',A),('U','CCTGGTGGTCCT'),('B',B)])
        ts={'P0001':ref,'P0002':q1,'P0003':q2};ss={p:assign_slots(v,ref)[0] for p,v in ts.items()}
        vs={p:path_view(v,ss[p],ref) for p,v in ts.items()}
        self.assertEqual(vs['P0002'][1]['key'],vs['P0003'][1]['key'])
        rs=[record(p,v) for p,v in ts.items()]
        rt,views,represented,losses,_,_=plan_edits(rs,ts,ss,{},.10)
        self.assertEqual(losses,[])
        g=literal_graph(rs,rt,ss);d,*_=materialize(rs,views,g,losses)
        validate_parameters(d,g,represented)
        self.assertEqual(represented['P0002'],rs[1]['sequence'])

    def test_C_cross_whole_structure_and_per_path_budget(self):
        ts={'P0001':tokens([('C',C),('A',A),('A',A),('C',C)]),
            'P0002':tokens([('C',C),('A',A),('A',A),('C',C)]),
            'P0003':tokens([('C',C),('B',B),('A',A),('C',C),('B',B)])}
        ref=ts['P0001'];ss={p:assign_slots(v,ref)[0] for p,v in ts.items()};rs=[record(p,v) for p,v in ts.items()]
        self.assertFalse(plan_edits(rs,ts,ss,{},.01)[3])
        result=plan_edits(rs,ts,ss,{},.03)
        self.assertTrue(result[3]);self.assertEqual(result[3][0]['loss_bp'],1)
        reverse=plan_edits(list(reversed(rs)),dict(reversed(list(ts.items()))),ss,{},.03)
        self.assertEqual(result[2:],reverse[2:])
        self.assertEqual(len(result[0]['P0003']),len(ts['P0003']))

    def test_E_split_unknown_preserves_all_flanking_bases(self):
        ts=tokens([('U','AACCGGTT')]);alt=tokens([('U','CCGG')]);alt[0].update(start=2,end=6)
        result=replace_region(ts,2,6,alt)
        self.assertEqual(result[0]['exact_bases'],'AACCGGTT')

    def test_all_optimal_DAG_export_keeps_legacy_selection(self):
        dictionary=[dict(id='A',base='A',sequence=A,level='base',match_orientations=['+']),
                    dict(id='B',base='B',sequence=B,level='base',match_orientations=['+'])]
        original=segment((A+B,dictionary,'HiFi'),True,True,True)
        new=segment((A+B,dictionary,'HiFi'),True,True,True,True)
        self.assertEqual(original[0],new[0]);self.assertIn('optimal_dag',new[1])
        self.assertEqual(new[1]['optimal_dag']['length'],32)

    def test_unknown_different_context_not_shared(self):
        ref=tokens([('A',A),('A',A),('B',B),('C',C)])
        a=tokens([('A',A),('U','AAA'),('B',B),('C',C)])
        b=tokens([('A',A),('U','AAA'),('C',C),('C',C)])
        va=path_view(a,assign_slots(a,ref)[0],ref);vb=path_view(b,assign_slots(b,ref)[0],ref)
        self.assertNotEqual(va[1]['key'],vb[1]['key'])

    def test_D_context_is_frozen_across_loss_budgets(self):
        ref=tokens([('A',A),('A',A),('B',B)])
        raw=tokens([('A',A),('U','AAA'),('B',B)]);ss=assign_slots(raw,ref)[0]
        represented=deepcopy(raw);represented[0].update(motif='B',exact_bases=B)
        self.assertEqual(path_view(raw,ss,ref)[1]['key'],path_view(represented,ss,ref,raw)[1]['key'])

    def test_ambiguous_and_weak_units_not_edited(self):
        ts={'P0001':tokens([('C',C),('A',A),('C',C)]),'P0002':tokens([('C',C),('A',A),('C',C)]),
            'P0003':tokens([('C',C),('B',B),('C',C)])}
        ts['P0003'][1]['ambiguity_kind']='type_only';ts['P0003'][1]['candidate_types']=['A','B']
        ss={p:assign_slots(v,ts['P0001'])[0] for p,v in ts.items()};rs=[record(p,v) for p,v in ts.items()]
        self.assertEqual(plan_edits(rs,ts,ss,{},.10)[3],[])


if __name__=='__main__':unittest.main()
