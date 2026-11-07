"""Fixed-denominator flank comparison and same-target search sensitivity audit."""
import argparse
from collections import Counter, defaultdict
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import statistics
from types import SimpleNamespace

import numpy as np

from vntrflankalign import locate_interval, digest
from vntrflankbounds import load_snapshot, quality, FETCH
from vntrreads import write_table
from vntrgraph import IndexedFastaReference,read_catalog,flank_hit

GROUPS=['B0','B1','B2']+[f'W{i}' for i in range(1,9)]
PERTURBATIONS=[dict(name='base',padding=[120,120],shifts=[0,0])]
for side in (0,1):
    for shift,pad in [(-10,120),(10,120),(-25,120),(25,120),(0,60),(0,240)]:
        pads=[120,120];shifts=[0,0];pads[side]=pad;shifts[side]=shift
        PERTURBATIONS.append(dict(name=f'{"left" if side==0 else "right"}_{shift}_{pad}',padding=pads,shifts=shifts))


def dump(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def load_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def compact_records(folder):
    result={}
    for line in (folder/'records.jsonl').open(encoding='utf-8'):
        r=json.loads(line);key=r['molecule_key']
        sequence=r.get('sequence','')
        if sequence:
            assert sequence==r['query_sequence'][r['qstart']:r['qend']].upper()
        if r.get('core_outer_relation')=='consistent':
            assert r['left_context']+sequence+r['right_context']==r['analysis_sequence']
        identity=sha256(json.dumps([r[k] for k in ['read_id','flag','chrom','pos','mapq','cigar','query_sequence','query_quality']],separators=(',',':')).encode()).hexdigest()
        item=dict(molecule=key,record_id=r['record_id'],qstart=r.get('qstart'),qend=r.get('qend'),
                  sequence_sha256=digest(sequence),sequence_length=len(sequence),raw_identity=identity,
                  status=r['status'],core_status=r.get('core_status'),extraction_status=r.get('extraction_status'),
                  core_ambiguous=r.get('core_ambiguous',False),analysis_status=r.get('analysis_status'),
                  analysis_ambiguous=r.get('analysis_ambiguous',False),core_outer_relation=r.get('core_outer_relation'),
                  threshold_path_dependent=r.get('threshold_path_dependent',False),
                  core_pairs=[[p['qstart'],p['qend']] for p in r.get('core_anchor_evidence',{}).get('pairs',[])],
                  core_representative=r.get('core_anchor_evidence',{}).get('representative'),
                  analysis_query_interval=r.get('analysis_query_interval'),mapq=r['mapq'],flag=r['flag'])
        result[r['record_id']]=item
    return result


def path_members(folder):
    paths=load_json(folder/'source_members.json');members={}
    for path in paths:
        for o in path['observations']:
            value=dict(path=path['id'],sequence_sha256=path['sha256'],length=len(path['sequence']),
                       sample=o['sample'],platform=o['platform'],qstart=o['qstart'],qend=o['qend'],record_id=o['record_id'])
            key=o['molecule']
            if key in members:
                assert members[key]['sequence_sha256']==value['sequence_sha256']
            members[key]=value
    return paths,members


def compare(run,locus):
    out=run/'audit'/locus;out.mkdir(parents=True,exist_ok=False)
    folders={g:run/'comparison'/g/locus for g in GROUPS}
    raw={g:compact_records(folder) for g,folder in folders.items()}
    path_data={g:path_members(folder) for g,folder in folders.items()}
    before=raw['B0'];base=path_data['B0'][1]
    source=load_snapshot(folders['B0']/'records.jsonl')
    identifiers={r['molecule_key']:dict(sample=r['sample'],platform=r['platform'],read_id=r['read_id'],mapq=r['mapq']) for r in source}
    inventory={r['record_id']:r for r in source}
    summaries=[];changes=[];transitions=[];strata=[]
    for group in GROUPS:
        now=raw[group];paths,members=path_data[group]
        assert set(now)==set(before),'Original record disappeared'
        assert all(now[k]['raw_identity']==before[k]['raw_identity'] for k in before),'Original input record changed'
        if group.startswith('W'):
            assert all((now[k]['core_pairs'],now[k]['core_representative'])==
                       (raw['B2'][k]['core_pairs'],raw['B2'][k]['core_representative']) for k in now),'Window changed independent core locator'
        summary=load_json(folders[group]/'summary.json')
        assert summary['molecules']==len(members) and summary['paths']==len(paths)
        common=set(base)&set(members)
        cut_changed={m for m in common if (base[m]['qstart'],base[m]['qend'])!=(members[m]['qstart'],members[m]['qend'])}
        seq_changed={m for m in common if base[m]['sequence_sha256']!=members[m]['sequence_sha256']}
        item=dict(locus=locus,group=group,input_records=len(now),input_molecules=len(identifiers),
            retained_molecules=len(members),hifi_molecules=sum(identifiers[m]['platform']=='HiFi' for m in members),
            ont_molecules=sum(identifiers[m]['platform']=='ONT' for m in members),
            ontul_molecules=sum(identifiers[m]['platform']=='ONTUL' for m in members),paths=len(paths),
            common_with_B0=len(common),newly_retained=len(set(members)-set(base)),lost=len(set(base)-set(members)),
            changed_cuts_in_common=len(cut_changed),changed_sequences_in_common=len(seq_changed),
            core_success=sum(r['core_status']=='success' for r in now.values()) if group!='B0' else None,
            core_ambiguous=sum(r['core_ambiguous'] for r in now.values()),
            retained_ambiguous=sum(now[members[m]['record_id']]['core_ambiguous'] for m in members),
            core_outer_conflicts=sum(r['core_outer_relation']=='conflict' for r in now.values()),
            analysis_success=sum(r['analysis_status']=='success' for r in now.values()),
            analysis_ambiguous=sum(r['analysis_ambiguous'] for r in now.values()),
            core_search_incomplete=sum(r['core_status']=='search_incomplete' for r in now.values()),
            analysis_search_incomplete=sum(r['analysis_status']=='search_incomplete' for r in now.values()),
            graph_reconstruction=all(v['paths']==len(paths) for k,v in summary['graph_stats'].items() if k!='preview'))
        if group!='B0':
            model=load_json(folders[group]/'boundary_model.json')
            item['analysis_interval']=model['analysis_interval']
        summaries.append(item)
        for sample in sorted({x['sample'] for x in identifiers.values()}):
            for platform in ['HiFi','ONT','ONTUL']:
                selected={m for m,x in identifiers.items() if x['sample']==sample and x['platform']==platform}
                if selected:
                    strata.append(dict(locus=locus,group=group,sample=sample,platform=platform,input_molecules=len(selected),
                         retained=len(selected&set(members)),common=len(selected&common),newly_retained=len(selected&(set(members)-set(base))),
                         lost=len(selected&(set(base)-set(members))),changed_cuts=len(selected&cut_changed),changed_sequences=len(selected&seq_changed)))
        transition_counts=defaultdict(list)
        for molecule,identity in sorted(identifiers.items()):
            old,new=base.get(molecule),members.get(molecule)
            record_id=(new or old or {}).get('record_id')
            if record_id is None:
                record_id=next(k for k,v in before.items() if v['molecule']==molecule)
            r=now[record_id];b=before[record_id]
            state=('retained_changed' if molecule in cut_changed or molecule in seq_changed else 'retained_unchanged') if old and new else 'recovered' if new else 'lost' if old else 'not_in_P'
            changes.append(dict(locus=locus,group=group,molecule=molecule,**identity,record_id=record_id,state=state,
                old_path=old['path'] if old else '',new_path=new['path'] if new else '',
                old_qstart=b['qstart'],old_qend=b['qend'],new_qstart=r['qstart'],new_qend=r['qend'],
                left_shift=r['qstart']-b['qstart'] if r['qstart'] is not None and b['qstart'] is not None else None,
                right_shift=r['qend']-b['qend'] if r['qend'] is not None and b['qend'] is not None else None,
                old_length=b['sequence_length'],new_length=r['sequence_length'],
                old_sequence_sha256=b['sequence_sha256'],new_sequence_sha256=r['sequence_sha256'],
                core_status=r['core_status'],extraction_status=r['extraction_status'],core_ambiguous=r['core_ambiguous'],
                analysis_status=r['analysis_status'],analysis_ambiguous=r['analysis_ambiguous'],
                core_outer_relation=r['core_outer_relation'],core_candidates=r['core_pairs']))
            transition_counts[(old['path'] if old else '',new['path'] if new else '',
                               old['sequence_sha256'] if old else '',new['sequence_sha256'] if new else '')].append(molecule)
        for key,mm in transition_counts.items():
            if not key[0] and not key[1]:continue
            transitions.append(dict(locus=locus,group=group,old_path=key[0],new_path=key[1],old_sha256=key[2],new_sha256=key[3],
                                    molecules=len(mm),members=mm))
    dump(out/'comparison_summary.json',summaries)
    write_table(out/'window_comparison.tsv',summaries,list(summaries[-1]))
    write_table(out/'sample_platform_retention.tsv',strata,list(strata[0]))
    write_table(out/'molecule_changes.tsv',changes,list(changes[0]))
    write_table(out/'path_transitions.tsv',transitions,list(transitions[0]))
    dump(out/'input_checks.json',dict(all_raw_records_identical=True,all_core_slices_exact=True,
         all_context_concatenations_exact=True,all_window_core_candidates_equal_B2=True,all_graph_counts_match=True,
         input_records_sha256=sha256((folders['B0']/'records.jsonl').read_bytes()).hexdigest()))
    # Every changed cut (including display-only), gain/loss, and ambiguous/conflict
    # case is retained. Evidence labels are deliberately not inferred from the
    # new locator's score or CIGAR agreement.
    review_records=set()
    reasons=defaultdict(set)
    for c in changes:
        if c['group']=='B0':continue
        if (c['old_qstart'],c['old_qend'])!=(c['new_qstart'],c['new_qend']):
            reasons[c['record_id']].add('cut_changed_including_display')
        if c['state'] in ('recovered','lost'):reasons[c['record_id']].add(c['state'])
        if c['core_ambiguous']:reasons[c['record_id']].add('core_ambiguity')
        if c['core_outer_relation']=='conflict':reasons[c['record_id']].add('core_outer_conflict')
        if c['core_status'] in ('search_incomplete','unavailable'):reasons[c['record_id']].add('core_not_determined')
    review_records.update(reasons)
    controls=defaultdict(list)
    for r in source:
        if r['record_id'] not in review_records and quality(r,20):
            controls[(r['sample'],r['platform'],'20-49' if r['mapq']<50 else '50-254')].append(r)
    for group,rr in controls.items():
        chosen=sorted(rr,key=lambda r:sha256(('flank-review-v1-20261108'+r['molecule_key']).encode()).hexdigest())[:10]
        for r in chosen:reasons[r['record_id']].add('stratified_unchanged_control')
    review=[]
    for record_id in sorted(reasons):
        r=inventory[record_id];candidates={g:raw[g][record_id] for g in GROUPS}
        any_ambiguous=any(c['core_ambiguous'] for c in candidates.values())
        item=dict(record_id=record_id,molecule=r['molecule_key'],sample=r['sample'],platform=r['platform'],mapq=r['mapq'],
                  triggers=sorted(reasons[record_id]),coordinate_reconstruction='verified',
                  boundary_label='multiple_model_supported_cuts' if any_ambiguous else 'insufficient_independent_boundary_evidence',
                  biological_repair_confirmed=False,biological_regression_confirmed=False,
                  evidence_basis='Raw slices and complete affine paths were replayed; neither the new score nor original CIGAR is an independent truth label.',
                  input_flag=r['flag'],input_cigar=r['cigar'],input_reference_start=r['pos'],configurations=candidates)
        cuts={p for c in candidates.values() for p in [c['qstart'],c['qend']] if p is not None}
        item['raw_cut_contexts']={str(c):dict(start=max(0,c-60),sequence=r['query_sequence'][max(0,c-60):c+60]) for c in sorted(cuts)}
        review.append(item)
    (out/'case_review.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=True)+'\n' for r in review),encoding='utf-8')
    write_table(out/'case_review.tsv',[{k:v for k,v in r.items() if k not in ('configurations','raw_cut_contexts')} for r in review],
                [k for k in review[0] if k not in ('configurations','raw_cut_contexts')] if review else ['record_id'])
    return out,source


def legacy_locator(module,record,interval,fetch,perturbation):
    calls=[0]
    def near(template,sequence,expected,identity):
        side=calls[0];calls[0]+=1
        exact=sequence.find(template)
        if exact>=0 and sequence.find(template,exact+1)<0:return exact,exact+len(template),1.0
        padding=perturbation['padding'][side];shift=perturbation['shifts'][side]
        a=max(0,expected+shift-padding);b=min(len(sequence),expected+shift+len(template)+padding)
        hit=flank_hit(template,sequence[a:b],expected-a,identity)
        return (hit[0]+a,hit[1]+a,hit[2]) if hit else None
    module.near_anchor=near
    s,e=interval
    locus=SimpleNamespace(chrom=record['chrom'],start=s,end=e,boundary_status='catalog')
    fields=[record['read_id'],str(record['flag']),record['chrom'],str(record['pos']+1),str(record['mapq']),record['cigar'],
            '*','0','0',record['query_sequence'],record['query_quality']]
    result=module.extract_record(fields,locus,fetch(s-100,s),fetch(e,e+100))
    threshold=.9 if record['platform']=='HiFi' else .8
    valid=(result.get('coverage')=='spans_candidate_interval' and result.get('anchor')=='both' and
           all(x is not None and x>=threshold for x in result.get('anchor_scores',[])) and
           not set(result.get('flags',[]))&{'reference_skip','anchor_order_conflict_cigar_fallback'})
    return dict(status='success' if valid else 'unavailable',cuts=[result.get('qstart'),result.get('qend')] if valid else None,
                candidates=[[result['qstart'],result['qend']]] if valid else [])


def sensitivity(run,project,locus_name,source,output):
    catalog=Path(__file__).resolve().parents[1]/'catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed'
    locus=next(l for l in read_catalog(catalog) if l.name==locus_name)
    fl,fr=FETCH[locus.name]
    with IndexedFastaReference(project/'reference/GRCh38_full_analysis_set_plus_decoy_hla.fa') as reference:
        dna=reference.fetch(locus.chrom,fl,fr).upper()
    fetch=lambda a,b:dna[a-fl:b-fl] if fl<=a<=b<=fr else ''
    spec=importlib.util.spec_from_file_location('legacy_vntrreads_sensitivity',run/'baseline_source/scripts/vntrreads.py')
    legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy)
    targets={tuple([locus.start,locus.end]):['core']}
    for group in GROUPS[3:]:
        model=load_json(run/'models'/locus.name/f'{group}.json')
        if model['analysis_interval']:
            targets.setdefault(tuple(model['analysis_interval']),[]).append(group)
    details=[];metrics=[];status_rows=[]
    for interval,labels in targets.items():
        print('SENSITIVITY_TARGET',locus.name,interval,labels,flush=True)
        per_method={method:{} for method in ('legacy100','affine100','affine250')}
        for record in source:
            if not quality(record,20):continue
            for method in per_method:
                results=[]
                for config in PERTURBATIONS:
                    if method=='legacy100':
                        result=legacy_locator(legacy,record,interval,fetch,config)
                    else:
                        value=locate_interval(record,interval,fetch,100 if method=='affine100' else 250,
                             padding=config['padding'],shifts=config['shifts'])
                        rep=value.get('representative')
                        result=dict(status=value['status'],cuts=[rep['qstart'],rep['qend']] if rep else None,
                                    candidates=[[p['qstart'],p['qend']] for p in value['pairs']],
                                    threshold_path_dependent=value.get('threshold_path_dependent',False))
                    results.append(dict(condition=config['name'],**result))
                per_method[method][record['record_id']]=dict(record=record,results=results)
                details.append(dict(record_id=record['record_id'],molecule=record['molecule_key'],sample=record['sample'],
                    platform=record['platform'],target=list(interval),target_labels=labels,method=method,conditions=results))
        for platform in ['HiFi','ONT','ONTUL']:
            keys={r['record_id'] for r in source if r['platform']==platform and quality(r,20)}
            for method,observations in per_method.items():
                subset=[observations[k] for k in keys]
                status_rows.append(dict(target=list(interval),labels=labels,platform=platform,method=method,
                    quality_records=len(keys),base_success=sum(x['results'][0]['status']=='success' for x in subset),
                    all_conditions_success=sum(all(r['status']=='success' for r in x['results']) for x in subset),
                    conditions_without_result=sum(r['status']!='success' for x in subset for r in x['results'])))
            common={k for k in keys if all(all(r['status']=='success' for r in per_method[m][k]['results']) for m in per_method)}
            for method in per_method:
                left_moves=[];right_moves=[];changes=0;maxima=[]
                for k in sorted(common):
                    rr=per_method[method][k]['results'];base=rr[0]['cuts']
                    lm=max(abs(r['cuts'][0]-base[0]) for r in rr);rm=max(abs(r['cuts'][1]-base[1]) for r in rr)
                    left_moves.append(lm);right_moves.append(rm);maxima.append(max(lm,rm));changes+=int(max(lm,rm)>0)
                metrics.append(dict(locus=locus.name,target=list(interval),labels=labels,platform=platform,method=method,
                    common_denominator=len(common),changed_molecules=changes,p95_max_move=float(np.percentile(maxima,95)) if maxima else None,
                    maximum_move=max(maxima,default=None),maximum_left_move=max(left_moves,default=None),maximum_right_move=max(right_moves,default=None)))
    (output/'sensitivity_records.jsonl').write_text(''.join(json.dumps(x,ensure_ascii=True)+'\n' for x in details),encoding='utf-8')
    write_table(output/'sensitivity_summary.tsv',metrics,list(metrics[0]))
    write_table(output/'sensitivity_status.tsv',status_rows,list(status_rows[0]))
    dump(output/'sensitivity_summary.json',metrics)
    dump(output/'sensitivity_settings.json',dict(perturbations=PERTURBATIONS,original_cigar_hints_unchanged=True,
        comparison_denominator='All three methods successful in all 13 conditions; unavailable cases separately enumerated',
        outer_controls='The same three locators are evaluated on each fixed outer target, not compared across shifted targets.'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',required=True,type=Path);parser.add_argument('--project',required=True,type=Path)
    parser.add_argument('--locus',required=True,choices=['ATXN2','D1S80'])
    parser.add_argument('--sensitivity',action='store_true')
    args=parser.parse_args()
    out,source=compare(args.run,args.locus)
    if args.sensitivity:sensitivity(args.run,args.project,args.locus,source,out)
    print('AUDIT_COMPLETE',args.locus,flush=True)


if __name__=='__main__':main()
