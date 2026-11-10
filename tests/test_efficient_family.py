import json,sys,tempfile,unittest
from pathlib import Path
from copy import deepcopy
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntreffmotifs import select_motifs
from vntrfamilies import select_families,apply_families
from vntrlocal import topology,path_view,plan_edits
from vntrlocal_export import export_budget
from vntrhmm import segment_hmm
from test_local_display import A,B,C,tokens,record


class EffectiveMotifTests(unittest.TestCase):
    def test_d1s80_q_point_one_matches_vamos_objective(self):
        model=json.loads((Path(__file__).resolve().parents[1]/'catalogs/D1S80.structure.v1.json').read_text())
        counts={'U01':2668,'U02':506,'U03':1314,'U04':988,'U05':264,'U06':393,'U07':523,'U08':193,'U09':590}
        ts={'P0001':[dict(status='unit',motif=m) for m,n in counts.items() for _ in range(n)]}
        r=select_motifs(model,ts,.1)
        self.assertEqual(r['delta'],743);self.assertEqual(r['weighted_template_cost'],657)
        self.assertEqual(r['selected_ids'],['U01','U02','U03','U04','U07','U08','U09'])
        self.assertEqual(r['mapping']['U05'],'U03');self.assertEqual(r['mapping']['U06'],'U04')

    def test_rare_distant_type_retained_and_terminal_excluded(self):
        model=dict(motifs=[dict(id='A',sequence='AAAA'),dict(id='B',sequence='AAAT'),dict(id='C',sequence='CCCC')],terminal_templates=[dict(id='RT',sequence='GGGG')])
        ts={'p':[dict(status='unit',motif=m) for m,n in [('A',100),('B',4),('C',5),('RT',100)] for _ in range(n)]}
        r=select_motifs(model,ts,.1)
        self.assertEqual(r['selected_ids'],['A','C']);self.assertNotIn('RT',r['counts'])
        self.assertEqual(select_motifs(model,ts,0)['selected_ids'],['A','B','C'])

    def test_effective_dictionary_actually_changes_hmm_choices_not_bases(self):
        dictionary=[dict(id='A',base='A',sequence=A,level='base')]
        result=segment_hmm((A+B+A,dictionary,'HiFi'))
        self.assertEqual(''.join(t['exact_bases'] for t in result),A+B+A)
        self.assertEqual({t['motif'] for t in result if t['status']=='unit'},{'A'})
        self.assertTrue(any(t['edit_distance'] for t in result if t['status']=='unit'))


def atom(kind,key,position):
    ts=tokens([(kind,{'A':A,'B':B,'C':C}[kind])])
    return dict(key=key,category=('unit',kind,'+','internal'),position=position,tokens=ts,
        count=1,variant_vector=[ts[0]['exact_bases']],occurrence_indices=[0],slots=[(0,1,0)],original_interval=[0,16])


class FamilyTests(unittest.TestCase):
    def test_same_position_exclusive_types_share_and_payload_is_unchanged(self):
        views={'p1':[atom('A','a',('slot',(0,1,0)))],'p2':[atom('B','b',('slot',(0,1,0)))]}
        original=deepcopy(views)
        family,rules,report=select_families(views,dict(motifs=[dict(id='A',sequence=A),dict(id='B',sequence=B)]))
        self.assertEqual(report['metrics']['N_d'],1);self.assertEqual(views,original)
        for p in views:
            self.assertEqual(family[p][0]['tokens'],views[p][0]['tokens'])
            self.assertEqual(family[p][0]['variant_vector'],views[p][0]['variant_vector'])
        self.assertEqual(len(report['families']),1)
        # A subsequent loss edit must not make two visits collapse to one.
        repeated=apply_families([views['p1'][0],views['p2'][0]],rules)
        self.assertEqual([v['key'] for v in repeated],['a','b'])

    def test_different_positions_and_cooccurrence_stay_separate(self):
        model=dict(motifs=[dict(id='A',sequence=A),dict(id='B',sequence=B)])
        different={'p1':[atom('A','a',('slot',(0,1,0)))],'p2':[atom('B','b',('slot',(1,1,0)))]}
        _,rules,_=select_families(different,model);self.assertEqual(rules,{})
        together={'p1':[atom('A','a',('slot',(0,1,0))),atom('B','b',('slot',(0,1,0)))],
                  'p2':[atom('A','a',('slot',(0,1,0)))],'p3':[atom('B','b',('slot',(0,1,0)))]}
        _,rules,_=select_families(together,model);self.assertEqual(rules,{})

    def test_similarity_chain_requires_one_direct_representative(self):
        motifs=[dict(id=m,sequence=s) for m,s in [('A','AAAA'),('B','AAAT'),('C','AATT'),('D','ATTT')]]
        views={m:[dict(atom('A',m,('slot',(0,1,0))),category=('unit',m,'+','internal'))] for m in 'ABCD'}
        _,_,r=select_families(views,dict(motifs=motifs))
        self.assertFalse(any(len(b)==4 for b in r['partition']))
        self.assertGreaterEqual(r['metrics']['N_d'],2)

    def test_three_percent_budget_and_atxn2_named_exports(self):
        ts={'P0001':tokens([('C',C),('A',A),('A',A),('C',C)]),
            'P0002':tokens([('C',C),('A',A),('A',A),('C',C)]),
            'P0003':tokens([('C',C),('B',B),('A',A),('C',C),('B',B)])}
        from vntrlocal import assign_slots
        ss={p:assign_slots(v,ts['P0001'])[0] for p,v in ts.items()};rs=[record(p,v) for p,v in ts.items()]
        builder=lambda p,v:path_view(v,ss[p],ts['P0001'],ts[p],parameterize=True)
        rep,vs,sequences,losses,_,decisions=plan_edits(rs,ts,ss,{},.03,view_builder=builder)
        self.assertTrue(all(x['fraction']<=.03 for x in losses))
        with tempfile.TemporaryDirectory() as directory:
            out=Path(directory)
            report=export_budget(out,rs,ts,rep,ss,vs,sequences,losses,.03,decisions,False,dict(locus='ATXN2',motifs=[dict(id='A',sequence=A)]))
            self.assertTrue((out/'ATXN2.represented.gfa').is_file())
            self.assertFalse((out/'D1S80.represented.gfa').exists())
            self.assertTrue(report['validation']['original_reconstruction'])


if __name__=='__main__':unittest.main()
