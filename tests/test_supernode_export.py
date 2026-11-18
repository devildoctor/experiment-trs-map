import json
import sys
import tempfile
import unittest
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

    def test_schema_budget_routing(self):
        self.assertEqual(structure_budget(1,0,2),2)
        for schema in (2,3,4):
            self.assertEqual(structure_budget(schema,.01,0),.01)
