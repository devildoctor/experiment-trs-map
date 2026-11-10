"""Run the approved local A-E display on all frozen B2 paths."""
import argparse,json,subprocess,time
from pathlib import Path
from hashlib import sha256
from vntrgraph import Locus
from vntrpathstructure import annotate_paths
from vntrdisplay_export import dump,table
from vntrlocal import select_interpretations,plan_edits,VERSION
from vntrlocal_export import export_budget


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--model',type=Path,required=True)
    p.add_argument('--policy',type=Path,default=Path(__file__).resolve().parents[1]/'catalogs/D1S80.structure-policy.v3.json')
    p.add_argument('--out',type=Path,required=True);p.add_argument('--prepare-only',action='store_true')
    p.add_argument('--no-render',action='store_true');args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    records=json.loads(args.input.read_text(encoding='utf-8'));model=json.loads(args.model.read_text(encoding='utf-8'))
    policy=json.loads(args.policy.read_text(encoding='utf-8'))
    assert policy['version']==VERSION and policy['model_version']==model['version']
    assert policy['requested_budgets']==[.01,.03,.05,.10] and policy['unit_max_edits']==2
    metadata=dict(input_sha256=sha256(args.input.read_bytes()).hexdigest(),model_sha256=sha256(args.model.read_bytes()).hexdigest(),version=VERSION)
    metadata['source_commit']=subprocess.check_output(['git','rev-parse','HEAD']).decode().strip()
    metadata['policy_sha256']=sha256(args.policy.read_bytes()).hexdigest()
    cache=args.out/'annotation_cache.json'
    if cache.exists():
        data=json.loads(cache.read_text());assert data['input_sha256']==metadata['input_sha256'] and data['model_sha256']==metadata['model_sha256']
        features=data['features']
    else:
        locus=Locus(model['locus'],model['chrom'],model['start'],model['end'],(),boundary_status=model['boundary_status'])
        features={};start=time.monotonic()
        for i,r in enumerate(records,1):
            _,_,f=annotate_paths([r],args.model,locus,select_representative=True,export_optimal=True);features.update(f)
            if i%20==0:print('ANNOTATED',i,'elapsed',round(time.monotonic()-start),flush=True)
        dump(cache,dict(metadata,features=features))
    dump(args.out/'source_members.json',records);dump(args.out/'locus_model.json',model);dump(args.out/'provenance.json',metadata)
    dump(args.out/'structure_policy.json',policy)
    if args.prepare_only:return
    choices=args.out/'selected_annotations.json'
    if choices.exists():
        data=json.loads(choices.read_text());tokens=data['tokens'];slots={p:[tuple(s) for s in ss] for p,ss in data['slots'].items()}
    else:
        print('CHOOSING_INTERPRETATIONS',flush=True)
        tokens,slots,alignments,regions,decisions=select_interpretations(features)
        dump(choices,dict(tokens=tokens,slots=slots))
        dump(args.out/'alignment_choices.json',alignments);dump(args.out/'annotation_choices.json',regions)
        dump(args.out/'interpretation_selection.json',decisions)
    reports=[]
    # Every budget starts from the same 0-loss state; no cross-budget state leakage.
    for fraction in [0,.01,.03,.05,.10]:
        print('BUDGET',fraction,flush=True)
        ts,views,represented,losses,mapping,decisions=plan_edits(records,tokens,slots,model,fraction)
        out=args.out/('%dpct'%round(fraction*100))
        report=export_budget(out,records,tokens,ts,slots,views,represented,losses,fraction,decisions,not args.no_render)
        reports.append(report);print(json.dumps({k:v for k,v in report.items() if k!='validation'}),flush=True)
    baseline=reports[0]['C'];requested=reports[1:];maximum=max(baseline-r['C'] for r in requested)
    feasible=[r for r in requested if baseline-r['C']>=.9*maximum-1e-9]
    chosen=min(feasible,key=lambda r:(r['total_edits_bp'],r['loss_budget_fraction']))
    dump(args.out/'comparison.json',reports);table(args.out/'comparison.tsv',reports,[k for k in reports[0] if k!='validation'])
    dump(args.out/'selection.json',dict(selected_fraction=chosen['loss_budget_fraction'],baseline_C=baseline,
        maximum_reduction=maximum,achieved_reduction=baseline-chosen['C'],
        criterion='lowest actual loss among requested budgets attaining >=90% of maximum incremental C reduction over v3 zero-loss'))
    print('LOCAL_V3_DONE',flush=True)


if __name__=='__main__':main()
