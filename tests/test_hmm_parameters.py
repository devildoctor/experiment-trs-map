import unittest
from copy import deepcopy
from test_local_display import A,B,C,tokens,record
from vntrlocal import classify,path_view,topology
from vntrlocal_export import literal_graph,materialize,validate_parameters


class ParameterDisplayTests(unittest.TestCase):
    def test_weak_run_keeps_each_string_and_annotation(self):
        ref=tokens([('A',A),('A',A)])
        query=tokens([('A','AAA'+A[3:]),('A',A)])
        query[0]['edit_distance']=3;query=classify(query,[])
        slots=[(0,1,0),(1,1,0)]
        self.assertEqual(len(path_view(query,slots,ref)),2)
        views={p:path_view(ts,slots,ref,parameterize=True)
               for p,ts in [('P0001',ref),('P0002',query)]}
        self.assertEqual(len(views['P0002']),1)
        self.assertEqual(views['P0001'][0]['key'],views['P0002'][0]['key'])
        original={'P0001':ref,'P0002':query};ss={p:slots for p in original}
        rs=[record(p,ts) for p,ts in original.items()];g=literal_graph(rs,original,ss)
        d,*_=materialize(rs,views,g,[],original)
        part=d['paths']['P0002'][0]['parts'][0]
        self.assertEqual(part['count'],2)
        self.assertEqual([a['weak_approximate'] for a in part['annotation']],[True,False])
        self.assertEqual([a['edit_distance'] for a in part['annotation']],[3,0])
        validate_parameters(d,g,{r['id']:r['sequence'] for r in rs})
        self.assertFalse(any(n['lossy_member_ids'] for n in d['nodes'].values()))

    def test_short_skip_context_is_preserved_as_visit_parameter(self):
        ref=tokens([('A',A),('A',A),('B',B),('C',C)])
        qa=tokens([('A',A),('U','AAA'),('B',B),('C',C)])
        qb=tokens([('C',C),('U','AA'),('C',C),('C',C)])
        slots=[(i,1,0) for i in range(4)]
        va=path_view(qa,slots,ref);vb=path_view(qb,slots,ref)
        self.assertNotEqual(va[1]['key'],vb[1]['key'])
        views={p:path_view(ts,slots,ref,parameterize=True)
               for p,ts in [('P0001',ref),('P0002',qa),('P0003',qb)]}
        self.assertEqual(views['P0002'][1]['key'],views['P0003'][1]['key'])
        self.assertNotEqual(views['P0002'][1]['skip_context'],views['P0003'][1]['skip_context'])
        original={'P0001':ref,'P0002':qa,'P0003':qb};ss={p:slots for p in original}
        rs=[record(p,ts) for p,ts in original.items()];g=literal_graph(rs,original,ss)
        d,*_=materialize(rs,views,g,[],original)
        for p in ['P0002','P0003']:
            part=next(part for visit in d['paths'][p] for part in visit['parts'] if 'skip_context' in part)
            self.assertEqual(part['skip_context'],views[p][1]['skip_context'])
        validate_parameters(d,g,{r['id']:r['sequence'] for r in rs})
        # Short skips in different comparison positions remain separate.
        shifted=list(slots);shifted[1]=(2,0,0)
        self.assertNotEqual(views['P0002'][1]['key'],path_view(qa,shifted,ref,parameterize=True)[1]['key'])

    def test_full_unit_length_skip_retains_context_and_literal_identity(self):
        ref=tokens([('A',A),('A',A),('B',B)])
        qa=tokens([('A',A),('U','A'*16),('B',B)])
        qb=tokens([('C',C),('U','A'*16),('B',B)])
        qc=tokens([('A',A),('U','C'*16),('B',B)])
        slots=[(i,1,0) for i in range(3)]
        ka=path_view(qa,slots,ref,parameterize=True)[1]['key']
        self.assertNotEqual(ka,path_view(qb,slots,ref,parameterize=True)[1]['key'])
        self.assertNotEqual(ka,path_view(qc,slots,ref,parameterize=True)[1]['key'])
        self.assertEqual(ka,path_view(qa,slots,ref)[1]['key'])

    def test_repeated_skip_visits_do_not_become_one_cycle(self):
        ref=tokens([('A',A),('B',B),('C',C)])
        query=tokens([('A',A),('U','T'),('B',B),('U','G'),('C',C)])
        slots=[(0,1,0),(1,0,0),(1,1,0),(1,0,0),(2,1,0)]
        views=path_view(query,slots,ref,parameterize=True)
        self.assertNotEqual(views[1]['key'],views[3]['key'])
        chains,_,paths,_=topology({'P0002':views})
        expanded=[key for node in paths['P0002'] for key in chains[node]]
        self.assertEqual(expanded,[v['key'] for v in views])

    def test_loss_marks_only_changed_node_including_boundary_indels(self):
        ref=tokens([('C',C),('A',A),('C',C)])
        for changed in [B[:-1],B+'T','G'+B[1:]]:
            with self.subTest(changed=changed):
                raw=tokens([('C',C),('B',B),('C',C)])
                original={'P0001':ref,'P0002':raw};represented=deepcopy(original)
                represented['P0002'][1]['exact_bases']=changed
                ss={p:[(i,1,0) for i in range(3)] for p in original}
                views={p:path_view(ts,ss[p],ref,original[p],parameterize=True) for p,ts in represented.items()}
                rs=[record(p,ts) for p,ts in original.items()];g=literal_graph(rs,represented,ss)
                d,*_=materialize(rs,views,g,[dict(path='P0002')],original)
                visits=d['paths']['P0002'];self.assertEqual(len(visits),3)
                self.assertEqual([d['nodes'][v['node']]['lossy_member_ids'] for v in visits],[[],['P0002'],[]])
                validate_parameters(d,g,{p:''.join(t['exact_bases'] for t in ts) for p,ts in represented.items()})


if __name__=='__main__':unittest.main()
