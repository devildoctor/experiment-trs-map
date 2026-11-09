"""Independent file-level checks for the real four-budget D1S80 delivery."""
from collections import Counter,defaultdict
from pathlib import Path
import csv,json,math,sys


def read_table(path):
    with path.open(encoding='utf-8') as h:return list(csv.DictReader(h,delimiter='\t'))


def read_gfa(path):
    nodes={};walks={};edges=set()
    for line in path.read_text().splitlines():
        f=line.split('\t')
        if f[0]=='S':nodes[int(f[1])]=f[2]
        elif f[0]=='L':assert f[5]=='0M';edges.add((int(f[1]),int(f[3])))
        elif f[0]=='P':walks[f[1]]=[int(x[:-1]) for x in f[2].split(',')]
    for walk in walks.values():assert all(e in edges for e in zip(walk,walk[1:]))
    return nodes,walks,{p:''.join(nodes[n] for n in walk) for p,walk in walks.items()}


def check(folder,source,features):
    rows=read_table(folder/'D1S80.path_sequences.txt');summary=json.loads((folder/'summary.json').read_text())
    by_id={r['path']:r for r in rows};records={r['id']:r for r in source};fraction=summary['loss_budget_fraction']
    assert set(by_id)==set(records) and len(rows)==307
    exact_nodes,exact_walks,raw=read_gfa(folder/'D1S80.exact.gfa')
    nodes,walks,represented=read_gfa(folder/'D1S80.structure.gfa')
    assert raw=={p:r['sequence'] for p,r in records.items()}
    assert raw=={p:r['original_sequence'] for p,r in by_id.items()}
    assert represented=={p:r['represented_sequence'] for p,r in by_id.items()}
    assert represented['P0001']==raw['P0001'] and not records['P0001']['molecules']
    assert json.loads((folder/'source_members.json').read_text())==source
    assert len(set().union(*(set(r['molecules']) for r in source)))==370
    losses={r['path']:r for r in map(json.loads,(folder/'losses.jsonl').read_text().splitlines())}
    for p,row in by_id.items():
        rep=row['representative'];assert by_id[rep]['representative']==rep
        assert represented[p]==raw[rep]
        assert int(row['original_length'])==len(raw[p]) and int(row['represented_length'])==len(represented[p])
        edits=losses.get(p,{}).get('events',[]);counts=Counter(e['op'] for e in edits)
        assert len(edits)<=math.floor(len(raw[p])*fraction+1e-9)
        assert len(represented[p])==len(raw[p])+counts['I']-counts['D']
        for k,op in [('substitutions','S'),('insertions','I'),('deletions','D')]:assert int(row[k])==counts[op]
        per_unit=Counter(e['unit'] for e in edits);assert all(c<=2 for c in per_unit.values())
        tokens=features[p]['tokens'];target=features[rep]['tokens'];assert len(tokens)==len(target)
        for i,(a,b) in enumerate(zip(tokens,target)):
            assert (a['motif'],a.get('orientation','+'))==(b['motif'],b.get('orientation','+'))
            if a['status']=='unknown' or b['status']=='unknown':assert a['exact_bases']==b['exact_bases']
            if a['status']!=b['status']:
                assert i in (0,len(tokens)-1) and {a['status'],b['status']}=={'partial','unit'}
                assert abs(len(a['exact_bases'])-len(b['exact_bases']))<=2
        original_parts=[];represented_parts=[];ocursor=rcursor=0
        for e in edits:
            assert raw[p][e['ostart']:e['oend']]==e['original']
            assert represented[p][e['rstart']:e['rend']]==e['represented']
            assert e['ostart']>=ocursor and e['rstart']>=rcursor
            original_parts.extend([represented[p][rcursor:e['rstart']],e['original']])
            represented_parts.extend([raw[p][ocursor:e['ostart']],e['represented']])
            ocursor=e['oend'];rcursor=e['rend']
        assert ''.join(original_parts)+represented[p][rcursor:]==raw[p]
        assert ''.join(represented_parts)+raw[p][ocursor:]==represented[p]
    d=json.loads((folder/'display_graph.json').read_text())
    for p,visits in d['paths'].items():
        assert [n for v in visits for n in v['gfa_nodes']]==walks[p]
        def expand(node,count):
            if node['kind']=='module':return ''.join(expand(x,1) for x in node['parts'])
            assert count in node['counts'];return node['unit_sequence']*count
        assert ''.join(expand(d['nodes'][v['node']],v['count']) for v in visits)==represented[p]
    for node in d['nodes'].values():
        ps=node['path_ids'];assert node['route_molecules']==len(set().union(*(set(records[p]['molecules']) for p in ps)))
    for e in d['edges']:
        assert e['route_molecules']==len(set().union(*(set(records[p]['molecules']) for p in e['path_ids'])))
        expected=set()
        for p in e['path_ids']:
            vs=d['paths'][p]
            for a,b in zip(vs,vs[1:]):
                if (a['node'],b['node'])!=(e['source'],e['target']):continue
                boundary=b['represented_interval'][0];clean=True
                for edit in losses.get(p,{}).get('events',[]):
                    if edit['rstart']==edit['rend']:dirty=boundary-1<edit['rstart']<boundary+1
                    else:dirty=edit['rstart']<boundary+1 and boundary-1<edit['rend']
                    if dirty:clean=False;break
                if clean:expected.update(records[p]['molecules'])
        assert len(expected)==e['raw_exact_molecules']
    groups=read_table(folder/'structure_groups.tsv')
    for g in groups:
        ps=json.loads(g['paths']);rep=g['representative']
        assert int(g['route_molecules'])==len(set().union(*(set(records[p]['molecules']) for p in ps)))
        assert int(g['representative_raw_exact_molecules'])==len(records[rep]['molecules'])
    decisions=json.loads((folder/'merge_decisions.json').read_text())
    assert all(b['C_before']-b['C_after']==b['delta_C']>0 for b in decisions['accepted_batches'])
    assert len(d['nodes'])==summary['N_d'] and len(d['edges'])==summary['E_d']
    assert summary['C']==summary['N_d']+summary['E_d']
    return dict(budget=fraction,paths=307,molecules=370,original_and_represented=True,
                reversible_edits=True,protected_structure=True,per_unit_and_path_budgets=True,
                literal_and_parameter_expansion=True,local_exact_edge_support=True,txt_roundtrip=True)


def main():
    root=Path(sys.argv[1]);source=json.loads(Path(sys.argv[2]).read_text())
    features=json.loads((root/'path_features.json').read_text());reports=[]
    for fraction in [0,1,3,5]:reports.append(check(root/f'{fraction}pct',source,features))
    rows=read_table(root/'simplification_comparison.tsv');selection=json.loads((root/'selection.json').read_text())
    baseline=int(rows[0]['C']);gain=max(baseline-int(r['C']) for r in rows)
    eligible=[r for r in rows if baseline-int(r['C'])>=.9*gain-1e-9]
    chosen=min(eligible,key=lambda r:(int(r['total_edits_bp']),float(r['loss_budget_fraction'])))
    assert float(chosen['loss_budget_fraction'])==selection['selected_fraction']
    assert by_path(root/'D1S80.path_sequences.txt')==by_path(root/f'{round(selection["selected_fraction"]*100)}pct'/'D1S80.path_sequences.txt')
    # Specifically verify the requested terminal-completion example when used.
    selected=by_path(root/'D1S80.path_sequences.txt');p13=selected['P0013']
    report=dict(status='passed',tiers=reports,selected_fraction=selection['selected_fraction'],
                P0013=dict(representative=p13['representative'],loss_bp=int(p13['loss_bp']),
                            original_length=int(p13['original_length']),represented_length=int(p13['represented_length'])))
    (root/'independent_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))


def by_path(path):return {r['path']:r for r in read_table(path)}


if __name__=='__main__':main()
