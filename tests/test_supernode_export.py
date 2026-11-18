import json
import sys
import tempfile
import unittest
from unittest.mock import patch
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from vntrgraph import Locus
from vntrpathgraph_cli import structure_budget
from vntrsupernode_export import export_supernodes,validate_bundle
from vntrpathgraph import validate_gfa


def fixture():
    m=json.loads((ROOT/'catalogs/D1S80.structure.v1.json').read_bytes())
    a,b=[m['motifs'][i]['sequence'] for i in (0,3)]
    records=[]
    for i,s in enumerate([(a+b)*3,(a+b)*3,(a+b)*5]):
        mol=[] if i==0 else ['shared']
        obs=[dict(sample='sample',molecule='shared',core_status='success',core_ambiguous=False,
                  boundary_status='confirmed',platform='HiFi',fragment=str(i))] if i else []
        records.append(dict(id=f'P{i+1:04d}',role='reference' if i==0 else 'observation',sequence=s,
             sha256=sha256(s.encode()).hexdigest(),molecules=mol,observations=obs,samples={'sample':1} if i else {}))
    return m,records


def export(folder,records):
    m,_=fixture()
    locus=Locus(m['locus'],m['chrom'],m['start'],m['end'],(),boundary_status='confirmed')
    return export_supernodes(folder,locus,records,ROOT/'catalogs/D1S80.structure.v1.json',
                             ROOT/'catalogs/D1S80.supernodes.v1.json',render=False)


class ExportTests(unittest.TestCase):
    def test_bundle_and_hash_binding(self):
        _,records=fixture()
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)
            summary=export(out,records)
            graph=json.loads((out/'display_graph.json').read_bytes())
            n=next(n for n in graph['nodes'].values() if n['kind']=='tandem_module')
            self.assertEqual(n['observed_repeat_counts'],[3,5])
            self.assertEqual((n['path_count'],n['molecule_count'],n['visit_count']),(3,1,3))
            self.assertEqual(summary['molecules'],1)
            self.assertEqual(validate_bundle(out)['paths'],3)
            other=out/'shuffled';other.mkdir()
            export(other,records[::-1])
            for name in ('display_graph.json','path_visits.jsonl','D1S80.original.gfa','D1S80.represented.gfa'):
                self.assertEqual((other/name).read_bytes(),(out/name).read_bytes())
            with (out/'D1S80.original.gfa').open('a') as f:f.write('# mismatch\n')
            with self.assertRaisesRegex(ValueError,'hash mismatch'): validate_bundle(out)

    def test_reference_only_and_iupac(self):
        _,records=fixture()
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp);r=records[0]
            r['sequence']='NRY';r['sha256']=sha256(b'NRY').hexdigest()
            result=export(out,[r])
            self.assertEqual(result['molecules'],0)
            self.assertEqual(result['edited_bp'],0)
            self.assertEqual(validate_gfa(out/'D1S80.original.gfa',{'P0001':'NRY'})['paths'],1)

    def test_rehashed_semantic_corruption(self):
        m,records=fixture()
        a,b=[m['motifs'][i]['sequence'] for i in (0,3)]
        records[1]['molecules']=['m1','m2','shared']
        for r in records:
            r['sequence']=a*3+b+a+b
            if r['id']=='P0003':r['sequence']=r['sequence'][:3]+'A'+r['sequence'][4:]
            r['sha256']=sha256(r['sequence'].encode()).hexdigest()
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            self.assertEqual(export(out,records)['edited_bp'],1)
            def reject(name,mutate):
                path=out/name
                raw=path.read_bytes()
                jsonl=name.endswith('.jsonl')
                value=[json.loads(x) for x in raw.splitlines()] if jsonl else json.loads(raw)
                mutate(value)
                content=''.join(json.dumps(x)+'\n' for x in value) if jsonl else json.dumps(value)
                path.write_text(content)
                manifest=json.loads((out/'provenance.json').read_bytes())
                manifest['files'][name]=sha256(path.read_bytes()).hexdigest()
                (out/'provenance.json').write_text(json.dumps(manifest))
                try:
                    with self.assertRaises(ValueError): validate_bundle(out)
                finally:
                    path.write_bytes(raw)
                    manifest['files'][name]=sha256(raw).hexdigest()
                    (out/'provenance.json').write_text(json.dumps(manifest))
            def move_boundary(rows):
                first,second=rows[0]['visits'][:2]
                second['original_gfa_walk'].insert(0,first['original_gfa_walk'].pop())
            reject('path_visits.jsonl',move_boundary)
            reject('path_visits.jsonl',lambda rows:rows.append(rows[0]))
            reject('path_visits.jsonl',lambda rows:rows[0]['visits'][1].update(visit_id=rows[0]['visits'][0]['visit_id']))
            reject('display_graph.json',lambda g:next(iter(g['nodes'].values())).update(molecule_count=90))
            reject('display_graph.json',lambda g:g['edges'][0].update(visit_count=90))
            reject('display_graph.json',lambda g:next(iter(g['nodes'].values())).update(lossy_visit_ids=['extra']))
            reject('losses.jsonl',lambda es:es.append(es[0]))
            reject('losses.jsonl',lambda es:es[0].update(rstart=0))
            reject('losses.jsonl',lambda es:es[0].update(path='P0001'))
            reject('losses.jsonl',lambda es:es[0].update(representative_occurrence='P0002:u1'))
            reject('path_visits.jsonl',lambda rows:next(v for row in rows for v in row['visits'] if v['loss_event_ids']).update(position_assessment='path_local'))
            manifest=json.loads((out/'provenance.json').read_bytes())
            del manifest['files']['paths.original.fa']
            (out/'provenance.json').write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,'manifest'):validate_bundle(out)

    def test_repeated_iupac_failure_rejects_tool_run_after_report(self):
        import verify_supernodes_tools as verify
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'tools'
            (root/'summary.json').write_text('{"locus":"X"}')
            for mode in ('original','represented'):
                (root/('paths.'+mode+'.fa')).write_text('>P1\nA\n')
                (root/('X.'+mode+'.gfa')).touch()
            def fake(tool,verb,*args,**kwargs):
                if verb=='paths':
                    return '>Repeated\nAAA\n' if 'repeated_iupac' in str(args[1]) else '>P1\nA\n'
                if verb=='image':Path(args[1]).write_text('<svg/>')
                return 'ok'
            argv=['verify','--bundle',str(root),'--out',str(out),'--odgi','odgi','--bandage','Bandage']
            with patch.object(sys,'argv',argv),patch.object(verify,'command',fake):
                with self.assertRaisesRegex(ValueError,'repeated_iupac'):verify.main()
            checks=json.loads((out/'tool_validation.json').read_bytes())['checks']
            self.assertTrue(checks['original']['odgi_sequence_match'])
            self.assertFalse(checks['repeated_iupac']['odgi_sequence_match'])

    def test_schema_budget_routing(self):
        self.assertEqual(structure_budget(1,0,2),2)
        for schema in (2,3,4):
            self.assertEqual(structure_budget(schema,.01,0),.01)
