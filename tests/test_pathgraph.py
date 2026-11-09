"""Scientific invariants for literal-sequence, path-preserving GFA export."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import vntrpathgraph as g
from vntrpathgraph_cli import export_graphs


def record(name,seq,molecule=None):
    return dict(id=name,sequence=seq,sha256=g.digest(seq),molecules=[molecule or name],
                hifi_molecules=[molecule or name],samples={'S1':1},platforms={'HiFi':1},observations=[])


def row(read='q1',sequence='ACGT',**fragment_changes):
    f=dict(id='f'+read,sequence=sequence,alignment_type='primary',mapq=60,flag=0,
           coverage='spans_candidate_interval',anchor='both',anchor_scores=[1.,1.],flags=[],qstart=100,qend=100+len(sequence))
    from vntrflankalign import VERSION
    pair=dict(qstart=100,qend=100+len(sequence))
    f.update(extraction_version=VERSION,core_status='success',
             core_anchor_evidence=dict(complete=True,status='success',representative=pair,pairs=[pair]))
    f.update(fragment_changes)
    return dict(id='r'+read,read_id=read,sample='S1',namespace='lib',platform='HiFi',bam='x.bam',fragments=[f])


class SelectionTests(unittest.TestCase):
    def test_bam_copies_and_identical_sequences(self):
        a=row();b=copy.deepcopy(a);b['bam']='copy.bam';b['id']='copy'
        paths,excluded=g.select_paths([a,b,row('q2')])
        self.assertFalse(excluded);self.assertEqual(len(paths),1)
        self.assertEqual(len(paths[0]['molecules']),2)
        self.assertEqual(len(paths[0]['observations']),3)

    def test_same_qname_different_namespace_is_distinct(self):
        b=row();b['namespace']='another'
        paths,_=g.select_paths([row(),b])
        self.assertEqual(len(paths[0]['molecules']),2)

    def test_partial_low_quality_and_mapq255_never_become_paths(self):
        rows=[row('partial',core_status='unavailable'),row('unavailable',mapq=255),
              row('weak',core_anchor_evidence={}),row('supplementary',alignment_type='supplementary'),
              row('duplicate',flag=1024),row('skip',flags=['reference_skip']),row('good')]
        paths,excluded=g.select_paths(rows)
        self.assertEqual(len(paths),1);self.assertEqual(len(paths[0]['molecules']),1)
        self.assertEqual(len(excluded),6)

    def test_conflicting_complete_molecule_is_excluded(self):
        paths,excluded=g.select_paths([row(sequence='AAAA'),row(sequence='CCCC')])
        self.assertEqual(paths,[])
        self.assertEqual({e['reason'] for e in excluded},{'conflicting_complete_sequences'})


class GraphTests(unittest.TestCase):
    def test_substitution_bubble_keeps_bases(self):
        seqs={'P1':'AACGTT','P2':'AATGTT'}
        graph=g.build_graph(seqs,seqs['P1'],{'P1':'6=','P2':'2=1X3='})
        compact=g.compact_graph(graph)
        self.assertEqual(len(compact['nodes']),4)
        self.assertEqual(sorted(compact['nodes'].values()),sorted(['AA','C','T','GTT']))
        g.validate_graph(compact,seqs)

    def test_insertions_deletions_and_repeated_visits_are_not_counts(self):
        seqs={'P1':'ACGT','P2':'ACACGT','P3':'AGT'}
        graph=g.build_graph(seqs,'ACGT',{'P1':'4=','P2':'2I4=','P3':'1=1D2='})
        compact=g.compact_graph(graph)
        g.validate_graph(compact,seqs)
        self.assertEqual({p:sum(len(compact['nodes'][n]) for n in w) for p,w in compact['paths'].items()},
                         {'P1':4,'P2':6,'P3':3})

    def test_same_base_different_positions_not_globally_collapsed(self):
        graph=g.build_graph({'P':'AAAA'},'AAAA',{'P':'4='})
        self.assertEqual(len(graph['nodes']),4)
        self.assertEqual(len(set(graph['paths']['P'])),4)
        self.assertEqual(list(g.compact_graph(graph)['nodes'].values()),['AAAA'])

    def test_internal_path_end_prevents_compaction_past_it(self):
        seqs={'short':'A','long':'AC'}
        graph=g.build_graph(seqs,'AC',{'short':'1=1D','long':'2='})
        self.assertEqual(len(g.compact_graph(graph)['nodes']),2)
        g.validate_graph(g.compact_graph(graph),seqs)

    def test_gfa_roundtrip_and_corruption(self):
        seqs={'P1':'AACGTT','P2':'AATGTT'}
        graph=g.compact_graph(g.build_graph(seqs,seqs['P1'],{'P1':'6=','P2':'2=1X3='}))
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'x.gfa';records=[record(p,s) for p,s in seqs.items()]
            g.write_gfa(path,graph,records,'test',seqs)
            self.assertEqual(g.validate_gfa(path,seqs)['paths'],2)
            lines=path.read_text().splitlines()
            first=next(i for i,l in enumerate(lines) if l.startswith('L\t'))
            del lines[first];path.write_text('\n'.join(lines)+'\n')
            with self.assertRaises(ValueError):g.validate_gfa(path,seqs)

    def test_bad_cigar_rejected(self):
        for cigar in ['3=','2=1D','0=4=','2=1X1=']:
            with self.subTest(cigar=cigar),self.assertRaises(ValueError):g.alignment_columns('ACGT','ACGT',cigar)

    def test_node_support_deduplicates_molecule(self):
        seqs={'P1':'AC','P2':'AT'}
        graph=g.build_graph(seqs,'AC',{'P1':'2=','P2':'1=1X'})
        ns,_=g.graph_support(graph,[record('P1','AC','same'),record('P2','AT','same')])
        shared=next(n for n,v in ns.items() if len(v['paths'])==2)
        self.assertEqual(ns[shared]['molecules'],1)
        self.assertEqual(ns[shared]['visits'],2)

    def test_empty_graph_is_explicitly_valid(self):
        graph=g.compact_graph(g.build_graph({},''))
        self.assertEqual(g.validate_graph(graph,{}),dict(paths=0,nodes=0,edges=0,reconstructed_bases=0))


class LossTests(unittest.TestCase):
    def setUp(self):
        self.records=[record('P1','A'*100),record('P2','A'*50+'C'+'A'*49)]
        self.rule=dict(source_sha256=self.records[1]['sha256'],target_sha256=self.records[0]['sha256'],
                       allowed_positions=[50],reason='synthetic approved one-base substitution')
        self.policy=dict(schema_version=1,locus='D1S80',rules=[self.rule])

    def test_explicit_rule_preserves_ids_and_lengths(self):
        seqs,ledger=g.apply_merge_rules(self.records,self.policy,.01,'D1S80')
        self.assertEqual(set(seqs),{'P1','P2'});self.assertEqual(len(set(seqs.values())),1)
        self.assertEqual(ledger[0]['positions'],[50]);self.assertEqual(ledger[0]['substitutions'],1)

    def test_budget_unapproved_position_and_unknown_target_rejected(self):
        with self.assertRaises(ValueError):g.apply_merge_rules(self.records,self.policy,.005,'D1S80')
        policy=copy.deepcopy(self.policy);policy['rules'][0]['allowed_positions']=[]
        with self.assertRaises(ValueError):g.apply_merge_rules(self.records,policy,.01,'D1S80')
        policy=copy.deepcopy(self.policy);policy['rules'][0]['target_sha256']='absent'
        with self.assertRaises(ValueError):g.apply_merge_rules(self.records,policy,.01,'D1S80')

    def test_indel_and_chained_replacement_rejected(self):
        records=self.records+[record('P3','A'*101)]
        policy=copy.deepcopy(self.policy);policy['rules'][0]['target_sha256']=records[-1]['sha256']
        with self.assertRaises(ValueError):g.apply_merge_rules(records,policy,1,'D1S80')
        policy=copy.deepcopy(self.policy);policy['rules'].append(dict(source_sha256=self.records[0]['sha256'],
                         target_sha256=self.records[1]['sha256'],allowed_positions=[50],reason='cycle'))
        with self.assertRaises(ValueError):g.apply_merge_rules(self.records,policy,1,'D1S80')


@unittest.skipUnless(importlib.util.find_spec('edlib'),'edlib required; these tests run in Slurm')
class EdlibExportTests(unittest.TestCase):
    def test_actual_edlib_indel_semantics_and_sequence_identity(self):
        for a,b in [('AAC','AC'),('AC','AAC'),('CAG'*6,'CAG'*4),('ATCG','ATGG'),('NAC','AAC')]:
            graph=g.compact_graph(g.build_graph({'P':a},b))
            g.validate_graph(graph,{'P':a})

    def test_export_paths_with_multiple_variants_and_preview(self):
        A='ACGTGCTAGTACCGTA';B='ACGTGCTAGTATCGTA';C='ACGTGCTAGTAACCGTA'
        seqs=[A*4+B*2+A*3,A*6+B*2+A*3,A*7,A*4+B*2+C+A*3]
        records=[record(f'P{i:04d}',s) for i,s in enumerate(seqs,1)]
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            report=export_graphs(out,'TEST',records,[],{},preview_limit=2,dot='missing_dot_for_test')
            self.assertEqual(report['paths'],4);self.assertEqual(report['preview']['omitted_paths'],2)
            self.assertEqual(report['substituted_bases'],0)
            self.assertLess(report['graph_stats']['exact']['nodes'],report['graph_stats']['base']['nodes'])
            g.validate_gfa(out/'TEST.structure.gfa',{r['id']:r['sequence'] for r in records})

    def test_lossy_export_keeps_original_identity_and_reports_loss(self):
        a='ACGT'*25;b=a[:50]+'A'+a[51:]
        records=[record('P1',a),record('P2',b)]
        policy=dict(schema_version=1,locus='TEST',rules=[dict(source_sha256=g.digest(b),target_sha256=g.digest(a),
                    allowed_positions=[50],reason='synthetic')])
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            report=export_graphs(out,'TEST',records,[],{},merge_policy=policy,loss_budget=.01,dot='missing_dot_for_test')
            self.assertEqual(report['paths'],2);self.assertEqual(report['structure_routes'],1)
            self.assertEqual(report['substituted_bases'],1)
            self.assertEqual(report['graph_stats']['structure']['paths'],2)
            self.assertNotEqual((out/'paths.original.fa').read_text(),(out/'paths.represented.fa').read_text())


if __name__=='__main__':unittest.main()
