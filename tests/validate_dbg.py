"""Independent GFA/SVG and provenance validation; run on a compute node."""
from collections import defaultdict
from hashlib import sha256
from pathlib import Path
import argparse
import json
import xml.etree.ElementTree as ET


def check(folder):
    original={r['id']:r for r in json.loads((folder/'input_paths.json').read_text())}
    results=[]
    for path in sorted(folder.glob('k*/*.gfa')):
        graph=json.loads(path.with_suffix('.json').read_text())
        k=graph['k']
        nodes,edges,paths,counts={},set(),{},{}
        for line in path.read_text(encoding='ascii').splitlines():
            f=line.split('\t')
            if f[0]=='S':
                assert f[1] not in nodes
                nodes[f[1]]=f[2]
                counts['node',f[1]]=int(next(x[5:] for x in f[3:] if x.startswith('ms:i:')))
            elif f[0]=='L':
                assert f[2]==f[4]=='+' and f[5]==f'{k-1}M'
                edge=(f[1],f[3])
                assert edge not in edges
                edges.add(edge)
                counts['edge',edge]=int(next(x[5:] for x in f[6:] if x.startswith('ms:i:')))
            elif f[0]=='P':
                assert f[1] not in paths and f[1] not in nodes
                assert all(n.endswith('+') for n in f[2].split(','))
                route=[n[:-1] for n in f[2].split(',')]
                assert f[3]==(','.join([f'{k-1}M']*(len(route)-1)) or '*')
                paths[f[1]]=route
        assert nodes=={n:v['sequence'] for n,v in graph['nodes'].items()}
        assert edges=={(e['source'],e['target']) for e in graph['edges']}
        assert set(paths)=={p['id'] for p in graph['paths']}
        for a,b in edges:
            assert nodes[a][-(k-1):]==nodes[b][:k-1]
        node_support,edge_support=defaultdict(set),defaultdict(set)
        used=set()
        raw_nodes=set()
        raw_edges=set()
        source_records=list(original.values())
        for p in source_records:
            w=[p['sequence'][i:i+k] for i in range(len(p['sequence'])-k+1)]
            raw_nodes.update(w)
            raw_edges.update(zip(w,w[1:]))
        for value in nodes.values():
            w=[value[i:i+k] for i in range(len(value)-k+1)]
            assert set(w)<=raw_nodes and set(zip(w,w[1:]))<=raw_edges
        for a,b in edges:
            assert (nodes[a][-k:],nodes[b][:k]) in raw_edges
        suffix='.supported' if graph['view']=='supported' else ''
        for p in graph['paths']:
            route=paths[p['id']]
            assert route==p['walk']
            adjacent=set(zip(route,route[1:]))
            assert adjacent<=edges
            used.update(adjacent)
            sequence=nodes[route[0]]+''.join(nodes[n][k-1:] for n in route[1:])
            assert sequence==p['sequence']
            parent=original[p.get('parent_id',p['id'])]
            expected=parent['sequence'][p.get('parent_start',0):p.get('parent_end',len(parent['sequence']))]
            assert sequence==expected
            assert sha256(sequence.encode()).hexdigest()==p['sequence_sha256']
            if p['kind']=='observed':
                for n in set(route):
                    node_support[n].add(p['molecule_key'])
                for edge in adjacent:
                    edge_support[edge].add(p['molecule_key'])
            svg=folder/path.parent.name/('paths-svg'+suffix)/(path.stem+'.'+p['id']+'.svg')
            root=ET.parse(svg).getroot()
            boxes=root.findall(".//*[@class='graph-node']")
            assert [x.attrib['data-node'] for x in boxes]==route
            assert sum(int(x.attrib['data-added']) for x in boxes)==len(sequence)
            links=root.findall(".//*[@class='graph-edge']")
            assert [(x.attrib['data-source'],x.attrib['data-target']) for x in links]==list(zip(route,route[1:]))
        assert used==edges
        for name,value in graph['nodes'].items():
            assert value['molecule_ids']==sorted(node_support[name])
            assert value['molecule_support']==counts['node',name]==len(node_support[name])
        for edge in graph['edges']:
            key=edge['source'],edge['target']
            assert edge['molecule_ids']==sorted(edge_support[key])
            assert edge['molecule_support']==counts['edge',key]==len(edge_support[key])
        svg=ET.parse(path.with_suffix('.svg')).getroot()
        ns={'s':'http://www.w3.org/2000/svg'}
        titles=lambda cls:{el.find('s:title',ns).text for el in svg.findall(f".//*[@class='{cls}']")}
        assert titles('node')==set(nodes)
        assert titles('edge')=={a+'->'+b for a,b in edges}
        allpaths=ET.parse(path.with_suffix('.paths.svg')).getroot()
        assert {p.attrib['data-path'] for p in allpaths.findall(".//*[@class='graph-path']")}==set(paths)
        results.append(dict(file=str(path.relative_to(folder)),status='PASS',nodes=len(nodes),edges=len(edges),exact_paths=len(paths),path_svgs=len(paths)))
    assert len(results)>0
    report=dict(status='PASS',graphs=results,checks=['GFA_1_S_L_P','exact_overlap','input_path_reconstruction','no_invented_raw_links','unique_molecule_support','all_topology_SVG_nodes_and_links','all_per_path_SVG_order_and_spelled_length','all_paths_SVG_inventory'])
    (folder/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('folder',type=Path)
    check(parser.parse_args().folder)
