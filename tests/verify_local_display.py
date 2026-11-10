"""Independent all-P acceptance checks for the completed server v3 run."""
import argparse,csv,json,math,sys
from collections import Counter,defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrpathstructure import restore_original
from vntrpathgraph import validate_gfa
from vntrcore import edit_distance
from vntrdisplay import graph_metrics


def main():
    ap=argparse.ArgumentParser();ap.add_argument('result',type=Path);args=ap.parse_args();root=args.result
    records=json.loads((root/'source_members.json').read_text());originals={r['id']:r['sequence'] for r in records}
    by_id={r['id']:r for r in records};chosen=json.loads((root/'selected_annotations.json').read_text());tokens=chosen['tokens'];slots=chosen['slots']
    assert len(originals)==307 and sum(r.get('role')=='reference' for r in records)==1
    assert len({m for r in records for m in r['molecules']})==370
    reports=[]
    for pct in [0,1,3,5,10]:
        folder=root/('%dpct'%pct);display=json.loads((folder/'display_graph.json').read_text())
        with (folder/'D1S80.path_sequences.txt').open() as f:rows={r['path']:r for r in csv.DictReader(f,delimiter='\t')}
        assert set(rows)==set(originals)==set(display['paths'])
        losses={x['path']:x for x in [json.loads(l) for l in (folder/'losses.jsonl').read_text().splitlines()]}
        represented={p:r['represented_sequence'] for p,r in rows.items()}
        originals_check={p:r['original_sequence'] for p,r in rows.items()};assert originals_check==originals
        validate_gfa(folder/'D1S80.original.gfa',originals);validate_gfa(folder/'D1S80.represented.gfa',represented)
        touched=0;unknown=0
        for p,visits in display['paths'].items():
            expanded=''.join(display['sequence_variants'][v] for visit in visits for part in visit['parts'] for v in part['variant_ids'])
            assert expanded==represented[p]
            assert restore_original(expanded,losses.get(p,{}).get('events',[]))==originals[p]
            vector=[(i,vid,a) for visit in visits for part in visit['parts'] for i,vid,a in zip(part['occurrence_indices'],part['variant_ids'],part['annotation'])]
            assert [i for i,v,a in vector]==list(range(len(tokens[p])))
            edit_sum=0
            for i,vid,a in vector:
                old=tokens[p][i];new=display['sequence_variants'][vid];d=edit_distance(old['exact_bases'],new);edit_sum+=d
                if old['status']=='unknown':unknown+=1;assert old['exact_bases']==new and a['status']=='unknown'
                if d:
                    touched+=1
                    assert p!='P0001' and old['status']=='unit' and old['ambiguity_kind']=='none'
                    assert not old['weak_approximate'] and old.get('role','internal')=='internal' and d<=2
                    assert any(q!='P0001' and slots[q][j]==slots[p][i] and t['exact_bases']==new and t['motif']==a['motif']
                        for q,ts in tokens.items() for j,t in enumerate(ts))
            assert edit_sum==int(rows[p]['loss_bp'])==losses.get(p,{}).get('loss_bp',0)
            assert edit_sum<=math.floor(len(originals[p])*pct/100+1e-9)
        assert represented['P0001']==originals['P0001']
        expected_edges={(a['node'],b['node']) for vs in display['paths'].values() for a,b in zip(vs,vs[1:])}
        assert expected_edges=={(e['source'],e['target']) for e in display['edges']}
        metrics=graph_metrics(display['nodes'],expected_edges)
        summary=json.loads((folder/'summary.json').read_text())
        assert all(summary[k]==v for k,v in metrics.items())
        for n,v in display['nodes'].items():
            ps={p for p,vs in display['paths'].items() if any(w['node']==n for w in vs)}
            assert ps==set(v['path_ids'])
            assert len({m for p in ps for m in by_id[p]['molecules']})==v['route_molecules']
        for e in display['edges']:
            ps={p for p,vs in display['paths'].items() if any((a['node'],b['node'])==(e['source'],e['target']) for a,b in zip(vs,vs[1:]))}
            assert ps==set(e['path_ids'])
            assert len({m for p in ps for m in by_id[p]['molecules']})==e['route_molecules']
        reports.append(dict(budget_pct=pct,paths=len(rows),unknown_occurrences_preserved=unknown,
            edited_occurrences=touched,edits=sum(int(r['loss_bp']) for r in rows.values()),**metrics))
        print('BUDGET_VERIFIED',pct,flush=True)
    regions=json.loads((root/'annotation_choices.json').read_text());cache=json.loads((root/'annotation_cache.json').read_text())['features']
    for p,rs in regions.items():
        assert math.prod(r['number'] for r in rs)==cache[p]['annotation']['optimal_dag']['path_count']
        for r in rs:
            assert len(r['alternatives'])==r['number']
            for alt in r['alternatives']:
                assert ''.join(t['exact_bases'] for t in alt)==originals[p][r['start']:r['end']]
            if p!='P0001':
                selected=r['candidate_graph_metrics'][r['selected']]
                assert (selected['N_d'],selected['E_d'])==min((x['N_d'],x['E_d']) for x in r['candidate_graph_metrics'])
                assert selected['N_d']==reports[0]['N_d'] and selected['E_d']==reports[0]['E_d']
    result=dict(status='passed',all_original_P_and_370_molecular_sources_preserved=True,
        full_GFA_and_TXT_reconstruction=True,all_display_parameters_reconstruct=True,
        loss_ledgers_and_per_unit_per_P_budgets=True,unknown_sequences_and_ambiguity_preserved=True,
        E_all_regional_optimal_alternatives_and_final_minimum_verified=True,
        E_global_joint_optimum_proven=False,node_and_edge_support_deduplicated=True,reports=reports)
    (root/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
