"""Leave-one-molecule discovery sensitivity with unchanged per-read anchor evidence."""
import argparse
import json
from pathlib import Path

from vntrflankbounds import load_snapshot,build_profile,discover,anchor_support,FETCH,SAMPLES
from vntrgraph import read_catalog,IndexedFastaReference
from vntrreads import write_table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',required=True,type=Path);p.add_argument('--project',required=True,type=Path)
    p.add_argument('--locus',required=True,choices=['ATXN2','D1S80'])
    args=p.parse_args()
    source=Path(__file__).resolve().parents[1]
    locus=next(l for l in read_catalog(source/'catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed') if l.name==args.locus)
    records=load_snapshot(args.run/'comparison/B0'/locus.name/'records.jsonl')
    fl,fr=FETCH[locus.name]
    with IndexedFastaReference(args.project/'reference/GRCh38_full_analysis_set_plus_decoy_hla.fa') as reference:
        dna=reference.fetch(locus.chrom,fl,fr).upper()
    full=build_profile(records,[fl,fr],dna);cache={};rows=[]
    full_models,_=discover(locus,full)
    for index,removed in enumerate(full['records']):
        key=removed['molecule_key']
        profile=build_profile([r for r in records if r['molecule_key']!=key],[fl,fr],dna)
        def lookup(boundary,side):
            k=(boundary,side)
            if k not in cache:cache[k]=anchor_support(full,boundary,side)
            base=cache[k]
            members={s:[m for m in base['members'][s] if m!=key] for s in SAMPLES}
            incomplete={s:[m for m in base.get('incomplete',{}).get(s,[]) if m!=key] for s in SAMPLES}
            counts={s:len(set(members[s])) for s in SAMPLES}
            status='pass' if min(counts.values())>=5 else 'anchor_support_insufficient' if any(counts[s]+len(incomplete[s])<5 for s in SAMPLES) else 'search_incomplete'
            return dict(base,status=status,members=members,counts=counts,incomplete=incomplete,
                        evidence=[e for e in base.get('evidence',[]) if e['molecule']!=key])
        models,audits=discover(locus,profile,support_lookup=lookup)
        for group,model in models.items():
            rows.append(dict(locus=locus.name,removed_molecule=key,sample=removed['sample'],group=group,
                full_interval=full_models[group]['analysis_interval'],leave_one_interval=model['analysis_interval'],
                changed=full_models[group]['analysis_interval']!=model['analysis_interval'],
                left_status=model['sides']['left']['status'],right_status=model['sides']['right']['status'],
                left_distance=model['sides']['left'].get('distance'),right_distance=model['sides']['right'].get('distance')))
        if (index+1)%20==0:print('LEAVE_ONE',locus.name,index+1,len(full['records']),flush=True)
    out=args.run/'audit'/locus.name
    write_table(out/'leave_one_discovery.tsv',rows,list(rows[0]))
    summary={g:dict(rounds=sum(r['group']==g for r in rows),changed=sum(r['group']==g and r['changed'] for r in rows),
         undetermined=sum(r['group']==g and r['leave_one_interval'] is None for r in rows)) for g in full_models}
    (out/'leave_one_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('LEAVE_ONE_COMPLETE',locus.name,json.dumps(summary),flush=True)


if __name__=='__main__':main()
