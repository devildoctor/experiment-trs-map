"""Reference identity, complete admission, support and full-canvas contracts."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrpathgraph import assign_reference_paths, graph_support, validate_graph
from vntrpathstructure import annotate_paths, build_structural_graph, merge_paths
from vntrpathstructure_export import export_structure
from test_pathstructure import record,feature,MODEL,POLICY,A
from test_pathgraph import row
from vntrpathgraph import select_paths


class AllFlankTests(unittest.TestCase):
    def test_reference_and_matching_observation_are_separate_zero_support_paths(self):
        selected,_=select_paths([row('good','ACGT'),row('low','ACGT',mapq=0),
                                  row('supp','ACGA',alignment_type='supplementary',flag=2048)])
        paths=assign_reference_paths(selected,'ACGT',dict(chrom='chr1',start=5,end=9))
        self.assertEqual([r['id'] for r in paths],['P0001','P0002','P0003'])
        self.assertEqual(paths[0]['molecules'],[])
        self.assertEqual(paths[1]['sequence'],paths[0]['sequence'])
        self.assertEqual(len(paths[1]['molecules']),2)
        again=assign_reference_paths(list(reversed(selected)),'ACGT',{})
        self.assertEqual([(r['id'],r['sequence']) for r in paths],[(r['id'],r['sequence']) for r in again])

    def test_reference_does_not_inherit_sample_specific_terminal(self):
        paths=assign_reference_paths([],A,{})
        model=dict(MODEL,terminal_templates=[dict(id='RT',sequence=A[-4:],samples=['S'])])
        token=dict(start=0,end=len(A),exact_bases=A,motif='unknown',status='unknown')
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'model.json';path.write_text('{}')
            with patch('vntrpathstructure.load_model',return_value=(model,[])),patch(
                    'vntrpathstructure.segment',return_value=([token],{})) as segment:
                _,_,features=annotate_paths(paths,path,'toy')
            self.assertEqual(segment.call_args.args[0][0],A)
            self.assertNotIn('terminal_scope',features['P0001']['annotation'])

    def test_reference_backbone_wins_over_observed_support_and_cannot_be_merged(self):
        a,fa=feature('AAA');b,fb=feature('AAAA')
        records=assign_reference_paths([record('old',b,20)],a,{})
        features={'P0001':fa,'P0002':fb};seqs={r['id']:r['sequence'] for r in records}
        graph=build_structural_graph(records,features,MODEL,seqs,{p:p for p in seqs})
        self.assertEqual(graph['backbone_path'],'P0001')
        validate_graph(graph,seqs)
        ns,_=graph_support(graph,records)
        self.assertTrue(all(v['molecules']==0 for n,v in ns.items() if v['paths']==['P0001']))
        represented,reps,losses,_=merge_paths(records,features,MODEL,POLICY,2)
        self.assertEqual(reps['P0001'],'P0001');self.assertEqual(represented['P0001'],a)

    def test_full_export_has_more_than_ten_paths_and_exact_txt_roundtrip(self):
        observed=[];by_sequence={}
        for i in range(1,13):
            sequence,f=feature('A'*i)
            observed.append(record('old'+str(i),sequence));by_sequence[sequence]=f
        reference,rf=feature('AB');by_sequence[reference]=rf
        records=assign_reference_paths(observed,reference,dict(chrom='chr1',start=0,end=len(reference)))
        features={r['id']:by_sequence[r['sequence']] for r in records}
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);policy=root/'policy.json';policy.write_text(json.dumps(POLICY))
            locus=type('Locus',(),dict(name='TOY'))()
            with patch('vntrpathstructure_export.annotate_paths',return_value=(MODEL,'hash',features)):
                result=export_structure(root,locus,records,[],dict(backbone_path='P0001'),
                                        root/'model.json',policy,0,0,'missing_dot_for_test')
            self.assertEqual(result['paths'],13);self.assertEqual(result['observed_paths'],12)
            self.assertEqual(result['reference_paths'],1);self.assertEqual(result['preview']['omitted_paths'],0)
            self.assertEqual(len(result['preview']['paths']),13)
            with (root/'TOY.path_sequences.txt').open() as handle:
                rows=list(csv.DictReader((s for s in handle if not s.startswith('#')),delimiter='\t'))
            self.assertEqual({r['name']:r['sequence'] for r in rows},{r['id']:r['sequence'] for r in records})
            self.assertEqual(rows[0]['molecules'],'0')


if __name__=='__main__':unittest.main()
