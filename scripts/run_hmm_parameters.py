"""Rebuild display sharing while freezing an existing HMM 1% representation."""
import argparse,json,subprocess
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
from vntrlocal import path_view,topology
from vntrlocal_export import export_budget
from vntrdisplay_export import dump,table

VERSION='hmm-weak-skip-parameters-v1'


def read(path):return json.loads(path.read_text(encoding='utf-8'))


def rebuild(source,out,render=True):
    records=read(source/'source_members.json');model=read(source/'locus_model.json')
    selected=read(source/'selected_annotations.json');original=selected['tokens']
    slots={p:[tuple(x) for x in ss] for p,ss in selected['slots'].items()}
    previous=source/'1pct';old_display=read(previous/'display_graph.json');before=read(previous/'summary.json')
    assert before['loss_budget_fraction']==.01 and before['decoder']=='single_viterbi'
    losses=[json.loads(line) for line in (previous/'losses.jsonl').read_text(encoding='utf-8').splitlines()]
    assert all(x['fraction']<=.01 for x in losses)
    templates={m['id']:m['sequence'] for m in model['motifs']+model.get('terminal_templates',[])}
    represented_tokens=deepcopy(original);represented={}
    for p,visits in old_display['paths'].items():
        seen=[];strings=[]
        for visit in visits:
            for part in visit['parts']:
                for i,variant,a in zip(part['occurrence_indices'],part['variant_ids'],part['annotation']):
                    t=represented_tokens[p][i];t.update(a)
                    t['exact_bases']=old_display['sequence_variants'][variant]
                    if t['status']=='unit':t['template']=templates[t['motif']];t['base']=t['motif']
                    seen.append(i);strings.append(t['exact_bases'])
        assert seen==list(range(len(original[p])))
        represented[p]=''.join(strings)
    # First reproduce the predecessor topology, then change display grouping only.
    base={p:path_view(ts,slots[p],original['P0001'],original[p]) for p,ts in represented_tokens.items()}
    baseline_metrics=topology(base)[3]
    assert all(baseline_metrics[k]==before[k] for k in ['N_d','E_d','C'])
    views={p:path_view(ts,slots[p],original['P0001'],original[p],parameterize=True) for p,ts in represented_tokens.items()}
    for p,vs in views.items():assert ''.join(s for v in vs for s in v['variant_vector'])==represented[p]
    out.mkdir(parents=True,exist_ok=False)
    for name in ['source_members.json','locus_model.json','selected_annotations.json','original_annotations.json','alignment_choices.json']:
        (out/name).write_bytes((source/name).read_bytes())
    report=export_budget(out/'1pct',records,original,represented_tokens,slots,views,represented,losses,.01,
        dict(mode=VERSION,loss_optimization='not_rerun',frozen_source=str(source)),render,model)
    # Same strings, same ledger, same original/represented literal graphs.
    unchanged=['losses.jsonl','D1S80.path_sequences.txt','paths.original.fa','paths.represented.fa',
               'D1S80.original.gfa','D1S80.represented.gfa']
    for name in unchanged:assert (out/'1pct'/name).read_bytes()==(previous/name).read_bytes(),name
    assert report['total_edits_bp']==before['total_edits_bp']
    report.update(decoder='single_viterbi',parameterization=VERSION,E_selection=False,
                  loss_optimization='not_rerun',loss_marker_scope='node_local_string_change')
    dump(out/'1pct/summary.json',report)
    after_display=read(out/'1pct/display_graph.json')
    comparison=[]
    for label,r,g in [('previous',before,old_display),('parameterized',report,after_display)]:
        comparison.append(dict(version=label,**{k:r[k] for k in ['N_d','E_d','C','total_edits_bp','lossy_paths','max_path_loss_fraction','paths','molecules']},
            marked_nodes=sum(bool(n['lossy_member_ids']) for n in g['nodes'].values())))
    dump(out/'comparison.json',comparison);table(out/'comparison.tsv',comparison)
    dump(out/'provenance.json',dict(version=VERSION,source=str(source),unchanged_files=unchanged,
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[1]).decode().strip(),
        frozen_loss_sha256=sha256((previous/'losses.jsonl').read_bytes()).hexdigest(),
        changes=['weak_in_visit_parameters','short_skip_context_in_visit_parameters','node_local_loss_marker']))
    print(json.dumps(comparison,indent=2),flush=True)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--no-render',action='store_true');args=parser.parse_args()
    rebuild(args.source,args.out,not args.no_render)


if __name__=='__main__':main()
