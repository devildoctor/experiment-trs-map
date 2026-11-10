"""Build all-P D1S80 graphs from one HMM/Viterbi annotation per sequence."""
import argparse,json,subprocess,time
from pathlib import Path
from hashlib import sha256
from vntrgraph import Locus
from vntrpathstructure import annotate_paths
from vntrlocal import classify,assign_slots,plan_edits
from vntrlocal_export import export_budget,render_graph
from vntrdisplay_export import dump,table
from vntrhmm import VERSION,PARAMETERS


def prepare(records,model_path,locus):
    model,model_hash,features=annotate_paths(records,model_path,locus,decoder='hmm_single')
    tokens={p:classify(f['tokens'],[]) for p,f in features.items()}
    slots={};alignments={}
    for p,ts in tokens.items():slots[p],alignments[p]=assign_slots(ts,tokens['P0001'])
    return model,model_hash,features,tokens,slots,alignments


def export_hmm(destination,locus,records,model_path,fraction=.01,render=True,dot='dot'):
    destination.mkdir(parents=True,exist_ok=True)
    model,model_hash,features,tokens,slots,alignments=prepare(records,model_path,locus)
    dump(destination/'source_members.json',records);dump(destination/'locus_model.json',model)
    dump(destination/'selected_annotations.json',dict(tokens=tokens,slots=slots))
    dump(destination/'original_annotations.json',features);dump(destination/'alignment_choices.json',alignments)
    ts,views,represented,losses,_,decisions=plan_edits(records,tokens,slots,model,fraction)
    report=export_budget(destination,records,tokens,ts,slots,views,represented,losses,fraction,decisions,False,model)
    if render:render_graph(destination,json.loads((destination/'display_graph.json').read_text()),dot)
    report.update(annotation_version=VERSION,decoder='single_viterbi',E_selection=False,model_sha256=model_hash)
    dump(destination/'summary.json',report)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--model',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--no-render',action='store_true')
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    records=json.loads(args.input.read_text());model=json.loads(args.model.read_text())
    locus=Locus(model['locus'],model['chrom'],model['start'],model['end'],(),boundary_status=model['boundary_status'])
    start=time.monotonic()
    model,model_hash,features,tokens,slots,alignments=prepare(records,args.model,locus)
    print('HMM_ANNOTATED',len(tokens),'seconds',round(time.monotonic()-start,2),flush=True)
    for name,data in [('source_members',records),('locus_model',model),('original_annotations',features),
                      ('selected_annotations',dict(tokens=tokens,slots=slots)),('alignment_choices',alignments)]:
        dump(args.out/(name+'.json'),data)
    dump(args.out/'provenance.json',dict(version=VERSION,decoder='single_viterbi',parameters=PARAMETERS,
        input_sha256=sha256(args.input.read_bytes()).hexdigest(),model_sha256=model_hash,
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[1]).decode().strip(),
        E_selection=False,requested_loss_fraction=.01))
    reports=[]
    for fraction in [0,.01]:
        ts,views,represented,losses,_,decisions=plan_edits(records,tokens,slots,model,fraction)
        report=export_budget(args.out/('%dpct'%round(fraction*100)),records,tokens,ts,slots,views,
            represented,losses,fraction,decisions,not args.no_render,model)
        report.update(annotation_version=VERSION,decoder='single_viterbi',E_selection=False)
        dump(args.out/('%dpct'%round(fraction*100))/'summary.json',report)
        reports.append(report);print(json.dumps({k:v for k,v in report.items() if k!='validation'}),flush=True)
    dump(args.out/'comparison.json',reports);table(args.out/'comparison.tsv',reports,[k for k in reports[0] if k!='validation'])
    print('HMM_DISPLAY_DONE',flush=True)


if __name__=='__main__':main()
