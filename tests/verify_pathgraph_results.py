"""Independent file audit: GFA -> FASTA -> actual stored SAM SEQ slices.

This intentionally does not import the builder or its GFA parser.
"""
import csv
from collections import Counter, defaultdict
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET


def fasta(path):
    result = {}; name = None
    for line in path.read_text(encoding='ascii').splitlines():
        if line.startswith('>'):
            name = line[1:]; assert name not in result
            result[name] = ''
        else:
            assert name is not None
            result[name] += line
    return result


def tags(fields):
    return {x.split(':', 2)[0]: x.split(':', 2)[2] for x in fields}


def verify_gfa(path, expected, molecules):
    nodes = {}; edges = {}; walks = {}; path_tags = {}; node_tags = {}
    for line in path.read_text(encoding='ascii').splitlines():
        f = line.split('\t')
        if f[0] == 'S':
            assert f[1] not in nodes and f[2] != '*'
            nodes[f[1]] = f[2]; node_tags[f[1]] = tags(f[3:])
            assert int(node_tags[f[1]]['LN']) == len(f[2])
        elif f[0] == 'L':
            assert f[2] == f[4] == '+' and f[5] == '0M'
            edge = (f[1], f[3]); assert edge not in edges
            edges[edge] = tags(f[6:])
        elif f[0] == 'P':
            assert f[1] not in walks and f[3] == '*'
            steps = f[2].split(','); assert all(x.endswith('+') for x in steps)
            walks[f[1]] = [x[:-1] for x in steps]; path_tags[f[1]] = tags(f[4:])
    assert set(walks) == set(expected)
    node_members = defaultdict(set); edge_members = defaultdict(set); visits = Counter()
    for p, walk in walks.items():
        assert ''.join(nodes[n] for n in walk) == expected[p], (path, p)
        assert int(path_tags[p]['mc']) == len(molecules[p])
        for n in walk: node_members[n].add(p); visits[n] += 1
        for edge in zip(walk, walk[1:]):
            assert edge in edges
            edge_members[edge].add(p)
    assert set(nodes) == set(node_members) and set(edges) == set(edge_members)
    for n, ps in node_members.items():
        t = node_tags[n]
        assert int(t['mc']) == len(set().union(*(molecules[p] for p in ps)))
        assert int(t['pc']) == len(ps) and int(t['vc']) == visits[n]
    for edge, ps in edge_members.items():
        assert set(edges[edge]['ps'].split(',')) == ps
        assert int(edges[edge]['mc']) == len(set().union(*(molecules[p] for p in ps)))
    return dict(paths=len(walks), nodes=len(nodes), edges=len(edges),
                reconstructed_bases=sum(map(len, expected.values())))


def verify(root):
    summary = json.loads((root/'summary.json').read_text(encoding='utf-8'))
    assert summary['status'] == 'complete' and summary['paths'] > 0
    locus = summary['locus']
    original = fasta(root/'paths.original.fa'); represented = fasta(root/'paths.represented.fa')
    records = json.loads((root/'source_members.json').read_text(encoding='utf-8'))
    members = {r['id']: set(r['molecules']) for r in records}
    assert set(members) == set(original) == set(represented)
    assert len(members) == summary['paths']
    assert sum(map(len, members.values())) == len(set().union(*members.values())) == summary['molecules']
    stored = {}
    with (root/'records.jsonl').open(encoding='utf-8') as handle:
        for line in handle:
            r = json.loads(line); key = (r['row_id'], r['record_id'])
            assert key not in stored, 'Colliding original record identifiers'
            stored[key] = r
    observations = 0
    for r in records:
        p = r['id']; assert r['sequence'] == original[p]
        assert r['sha256'] == sha256(original[p].encode('ascii')).hexdigest()
        assert {o['molecule'] for o in r['observations']} == members[p]
        for o in r['observations']:
            raw = stored[(o['row_id'], o['record_id'])]
            assert raw['query_sequence'].upper()[o['qstart']:o['qend']] == original[p]
            assert raw['molecule_key'] == o['molecule']
            assert raw['alignment_type'] == 'primary' and not raw['flag'] & (4|256|512|1024|2048)
            assert 20 <= raw['mapq'] <= 254
            assert raw['coverage'] == 'spans_candidate_interval' and raw['anchor'] == 'both'
            cutoff = .90 if o['platform'] == 'HiFi' else .80
            assert all(x >= cutoff for x in raw['anchor_scores'])
            observations += 1
    results = {}
    for name in ['base', 'exact', 'structure', 'preview']:
        expected = represented if name in ['structure', 'preview'] else original
        if name == 'preview': expected = {p: represented[p] for p in summary['preview']['paths']}
        results[name] = verify_gfa(root/f'{locus}.{name}.gfa', expected, members)
        assert results[name] == summary['graph_stats'][name]
    with (root/'losses.tsv').open(encoding='utf-8', newline='') as handle:
        ledger = {r['path']: r for r in csv.DictReader(handle, delimiter='\t')}
    changed = {}
    for p in original:
        a, b = original[p], represented[p]; assert len(a) == len(b)
        positions = [i for i, (x, y) in enumerate(zip(a,b)) if x != y]
        if positions:
            loss = ledger[p]
            assert json.loads(loss['positions']) == positions
            assert int(loss['substitutions']) == len(positions)
            assert len(positions)/len(a) <= summary['loss_budget'] + 1e-12
            changed[p] = len(positions)
    assert set(changed) == set(ledger)
    assert sum(changed.values()) == summary['substituted_bases']
    assert len(changed) == summary['lossy_paths']
    assert summary['preview']['status'] == 'rendered'
    preview_nodes=set();preview_edges=set()
    for line in (root/f'{locus}.preview.gfa').read_text().splitlines():
        f=line.split('\t')
        if f[0]=='S':preview_nodes.add(f[1])
        if f[0]=='L':preview_edges.add((f[1],f[3]))
    svg=ET.parse(root/f'{locus}.preview.svg')
    ns={'s':'http://www.w3.org/2000/svg'}
    drawn_nodes=Counter();drawn_edges=Counter();continuations=0
    for element in svg.findall('.//s:g',ns):
        title=element.find('s:title',ns)
        if title is None:continue
        label=title.text
        if element.get('class')=='node' and re.fullmatch(r's\d+',label):drawn_nodes[label[1:]]+=1
        if element.get('class')=='edge':
            a,b=label.split('->')
            if a.startswith('port_in_'):
                source,target=a[len('port_in_'):].split('_');assert b=='s'+target
                continuations+=1
            elif b.startswith('port_out_'):
                source,target=b[len('port_out_'):].split('_');assert a=='s'+source
                continuations+=1
            else:
                assert a.startswith('s') and b.startswith('s');source,target=a[1:],b[1:]
            drawn_edges[(source,target)]+=1
    assert set(drawn_nodes)==preview_nodes and set(drawn_nodes.values())=={1}
    assert set(drawn_edges)==preview_edges and all(n in [1,2] for n in drawn_edges.values())
    assert continuations==2*sum(n==2 for n in drawn_edges.values())
    result = dict(locus=locus, independent_file_audit_passed=True, original_records=len(stored),
                  selected_observations=observations, molecules=summary['molecules'],
                  graphs=results, changed_bases=sum(changed.values()),
                  static_svg_topology_passed=True, static_edge_continuations=continuations//2,
                  scope='File reconstruction, selected SAM slices, molecule support and loss accounting; not biological truth')
    (root/'verification.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False), flush=True)


if __name__ == '__main__': verify(Path(sys.argv[1]))
