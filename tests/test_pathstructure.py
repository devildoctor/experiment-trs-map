"""Targeted contracts for structure protection, internal loss and graph paths."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrpathstructure import (summarize,compare_for_merge,merge_paths,restore_original,
                               build_structural_graph,select_subgraph,edit_script)
from vntrpathgraph import digest,validate_graph,write_gfa,validate_gfa
from vntrstructure import segment
from vntrmotifs import definition

A='ACGTCAGTACGATGCA'
B='G'+A[1:]
MODEL=dict(version='toy.v1',motifs=[dict(id='A',sequence=A),dict(id='B',sequence=B)],
           display_labels={},interruptions=['B'])
POLICY=dict(schema_version=1,version='toy.loss.v1',model_version='toy.v1',
            allowed_template_positions={'A':list(range(1,15)),'B':list(range(1,15))})


def feature(labels,actual=None,tied=False,templates=None):
    templates=templates or {'A':A,'B':B}; tokens=[]; start=0
    for i,label in enumerate(labels):
        sequence=(actual or {}).get(i,templates[label])
        edits=len(edit_script(sequence,templates[label])[0])
        tokens.append(dict(start=start,end=start+len(sequence),exact_bases=sequence,
            template=templates[label],motif=label,status='unit',orientation='+',phase=0,edit_distance=edits))
        start+=len(sequence)
    sequence=''.join(t['exact_bases'] for t in tokens)
    return sequence,summarize(sequence,tokens,dict(annotation_ambiguity=tied),MODEL)


def record(name,sequence,count=1):
    molecules=[name+'.'+str(i) for i in range(count)]
    return dict(id=name,sequence=sequence,sha256=digest(sequence),molecules=molecules,
                hifi_molecules=molecules,samples={'S':count},platforms={'HiFi':count},observations=[])


class StructureTests(unittest.TestCase):
    def test_all_optimal_partitions_localize_a_realistic_subtype_tie(self):
        dictionary=[dict(id=name,base=name,sequence=seq,match_orientations=['+'],**definition(seq))
                    for name,seq in [('X','ACGA'),('Y','ACGC'),('Z','TTAA')]]
        tokens,audit=segment(('TTAAACGTTTAA',dictionary,'HiFi'),True,collect_optimal=True)
        self.assertTrue(audit['annotation_ambiguity'])
        self.assertEqual(audit['ambiguity_intervals'],[[4,8]])
        self.assertEqual([t['annotation_ambiguity'] for t in tokens],[False,True,False])

    def test_real_atxn2_interruption_is_not_a_neutral_snp(self):
        a,fa=feature('A'*8+'B'+'A'*4+'B'+'A'*8,templates={'A':'CTG','B':'TTG'})
        b,fb=feature('A'*8+'B'+'A'*13,templates={'A':'CTG','B':'TTG'})
        self.assertEqual(sum(x!=y for x,y in zip(a,b)),1)
        self.assertEqual([x['unit'] for x in fa['interruptions']],[9,14])
        self.assertEqual(compare_for_merge(a,b,fa,fb,POLICY,2)[1],'protected_structure')

    def test_internal_substitution_and_indel_restore_original(self):
        target,tf=feature('AAABAA')
        variants=[A[:7]+'C'+A[8:],A[:8]+'C'+A[8:],A[:8]+A[9:]]
        for variant,expected in zip(variants,['S','D','I']):
            with self.subTest(expected=expected):
                original,sf=feature('AAABAA',{1:variant})
                events,reason=compare_for_merge(original,target,sf,tf,POLICY,2)
                self.assertEqual(reason,'allowed_internal_edit')
                self.assertEqual([e['op'] for e in events],[expected])
                self.assertEqual(restore_original(target,events),original)

    def test_count_and_order_remain_different(self):
        for left,right in [('A'*8,'A'*9),('AABB','ABAB')]:
            a,fa=feature(left); b,fb=feature(right)
            self.assertEqual(compare_for_merge(a,b,fa,fb,POLICY,2)[1],'protected_structure')

    def test_terminal_and_uncertain_partitions_stay_separate(self):
        target,tf=feature('AAA')
        for index,tied,reason in [(0,False,'protected_end'),(1,True,'unresolved_partition')]:
            original,sf=feature('AAA',{index:A[:7]+'C'+A[8:]},tied=tied)
            self.assertEqual(compare_for_merge(original,target,sf,tf,POLICY,2)[1],reason)

    def test_repeated_edited_subtype_preserves_pattern(self):
        variant=A[:7]+'C'+A[8:]
        source,sf=feature('AAAAA',{1:variant,3:variant});target,tf=feature('AAAAA')
        self.assertEqual(compare_for_merge(source,target,sf,tf,POLICY,2)[1],'repeated_edited_subtype')

    def test_any_optimal_forbidden_indel_placement_rejects_merge(self):
        events,unsafe=edit_script('CAAAAG','CAAAG',{2,3})
        self.assertEqual(len(events),1);self.assertTrue(unsafe)

    def test_no_transitive_budget_expansion(self):
        variants=[A,A[:4]+'G'+A[5:8]+'G'+A[9:],
                  ''.join('G' if i in (4,8,12,14) else b for i,b in enumerate(A))]
        records=[];features={}
        for i,v in enumerate(variants):
            s,f=feature('AAA',{1:v});p='P'+str(i);records.append(record(p,s,3-i));features[p]=f
        represented,reps,losses,_=merge_paths(records,features,MODEL,POLICY,2)
        self.assertEqual(reps,{'P0':'P0','P1':'P0','P2':'P2'})
        self.assertEqual(losses[0]['loss_bp'],2)

    def test_whole_run_nodes_preserve_count_without_single_base_bubbles(self):
        records=[];features={}
        for p,labels in [('P1','A'*8),('P2','A'*9)]:
            s,f=feature(labels);records.append(record(p,s));features[p]=f
        seqs={r['id']:r['sequence'] for r in records}
        graph=build_structural_graph(records,features,MODEL,seqs,{p:p for p in seqs})
        self.assertEqual(sorted(map(len,graph['nodes'].values())),[128,144])
        self.assertEqual(set(graph['labels'].values()),{'A^8','A^9'})
        validate_graph(graph,seqs)

    def test_shared_blocks_skip_and_unresolved_path_keep_all_sequences(self):
        records=[];features={}
        for i,labels in enumerate(['AAAABBAAA','AAAAAABBAAA','AAAAAAA','AAAABABAAA']):
            p='P'+str(i);s,f=feature(labels,tied=i==3);records.append(record(p,s,4-i));features[p]=f
        seqs={r['id']:r['sequence'] for r in records}
        graph=build_structural_graph(records,features,MODEL,seqs,{p:p for p in seqs})
        validate_graph(graph,seqs)
        self.assertTrue(any('PARTITION UNRESOLVED' in s for s in graph['labels'].values()))
        self.assertFalse(any(a==b for a,b in graph['edges']))
        sub=select_subgraph(graph,['P0','P1'])
        validate_graph(sub,{p:seqs[p] for p in sub['paths']})
        self.assertTrue(set(sub['nodes']).issubset(graph['nodes']))

    def test_gfa_records_real_represented_length_and_loss(self):
        s,sf=feature('AAA',{1:A[:8]+'C'+A[8:]});t,tf=feature('AAA')
        records=[record('P1',t,4),record('P2',s)];features={'P1':tf,'P2':sf}
        represented,reps,losses,_=merge_paths(records,features,MODEL,POLICY,2)
        graph=build_structural_graph(records,features,MODEL,represented,reps)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'test.gfa';write_gfa(path,graph,records,'structure',represented,losses)
            validate_gfa(path,represented)
            p=next(line for line in path.read_text().splitlines() if line.startswith('P\tP2\t'))
            self.assertIn('\tld:i:1',p);self.assertIn('\tol:i:49',p);self.assertIn('\trl:i:48',p)


if __name__=='__main__':unittest.main()
