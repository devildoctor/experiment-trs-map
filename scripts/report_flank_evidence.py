"""Derive review tables from completed flank runs without relabeling biological truth."""
import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
from pathlib import Path

from audit_flank_results import GROUPS, PERTURBATIONS, load_json, path_members
from vntrflankbounds import load_snapshot, discovery_records, FETCH, SAMPLES
from vntrflankalign import flank_search
from vntrreads import write_table
from vntrgraph import IndexedFastaReference


def dump(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def table(path,rows):
    write_table(path,rows,list(rows[0]) if rows else ['no_rows'])


def json_lines(path):
    with path.open(encoding='utf-8') as stream:
        return [json.loads(line) for line in stream]


def read_table(path):
    with path.open(encoding='utf-8',newline='') as stream:
        return list(csv.DictReader(stream,delimiter='\t'))


def mapq_bin(value):
    return 'unknown_255' if value==255 else 'below_20' if value<20 else '20-49' if value<50 else '50-254'


def compact_evidence(record):
    result={k:record.get(k) for k in ['qstart','qend','sequence','core_status','analysis_status','extraction_status',
             'core_outer_relation','core_ambiguous','analysis_ambiguous','normalization_reason','anchor','anchor_scores']}
    for layer in ['core_anchor_evidence','analysis_anchor_evidence']:
        e=record.get(layer,{})
        result[layer]={k:e.get(k) for k in ['reference_interval','status','complete','pairs','representative',
                'ambiguous','distinct_sequence_count','threshold_path_dependent','normalization_reason']}
        for side in ['left','right']:
            s=e.get(side,{})
            result[layer][side]={k:s.get(k) for k in ['status','reference_interval','original_hint',
                'requested_search_interval','search_interval','physical_read_edge_clipped','optimal_traces']}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',required=True,type=Path)
    parser.add_argument('--project',required=True,type=Path)
    args=parser.parse_args();run=args.run
    out=run/'report';out.mkdir(exist_ok=False)
    all_stats=[];key_cases=[];paths_summary=[];discovery=[];strata=[];status_counts=[];control_sampling=[]
    sensitivity_cases=[];sensitivity_stats=[];lost_cases=[];blind_key=[];replay=[];complete_review=[]
    for locus in ['ATXN2','D1S80']:
        source=load_snapshot(run/'comparison/B0'/locus/'records.jsonl')
        byid={r['record_id']:r for r in source};bymolecule={r['molecule_key']:r for r in source}
        raw={g:{r['record_id']:r for r in json_lines(run/'comparison'/g/locus/'records.jsonl')} for g in GROUPS}
        paths={g:path_members(run/'comparison'/g/locus) for g in GROUPS}
        baseline_hashes={p['sha256'] for p in paths['B0'][0]}
        reviews=json_lines(run/'audit'/locus/'case_review.jsonl')
        reasons={r['record_id']:set(r['triggers'])-{'stratified_unchanged_control'} for r in reviews}
        for rid in byid:
            extra=set()
            if any(raw[g][rid].get('analysis_ambiguous') for g in GROUPS[1:]):extra.add('analysis_ambiguity')
            if any(raw[g][rid].get(layer,{}).get('threshold_path_dependent') for g in GROUPS[1:]
                   for layer in ['core_anchor_evidence','analysis_anchor_evidence']):extra.add('threshold_path_dependent')
            if extra:reasons.setdefault(rid,set()).update(extra)
        reasons={rid:triggers for rid,triggers in reasons.items() if triggers}
        pool=defaultdict(list)
        for r in source:
            if r['record_id'] not in reasons and 20<=r['mapq']<255 and not r['flag']&(4|256|512|1024|2048):
                pool[(r['sample'],r['platform'],mapq_bin(r['mapq']))].append(r)
        for rr in pool.values():
            for r in sorted(rr,key=lambda r:sha256(('flank-review-v1-20261108'+r['molecule_key']).encode()).hexdigest())[:10]:
                reasons[r['record_id']]={'stratified_unchanged_control'}
        reasons={rid:sorted(triggers) for rid,triggers in reasons.items()}
        reviews=[dict(locus=locus,record_id=rid,triggers=triggers,
            boundary_label='multiple_model_supported_cuts' if any('ambiguity' in t for t in triggers) else 'insufficient_independent_boundary_evidence',
            biological_repair_confirmed=False,biological_regression_confirmed=False,
            coordinate_reconstruction='verified',independent_biological_review='not_available') for rid,triggers in sorted(reasons.items())]
        complete_review.extend(reviews)
        original=load_json(run/'audit'/locus/'comparison_summary.json')
        final=load_json(run/'final_replay/audit'/locus/'comparison_summary.json')
        assert original==final,'Final code changed primary comparison outcomes'
        before=read_table(run/'audit'/locus/'molecule_changes.tsv')
        after=read_table(run/'final_replay/audit'/locus/'molecule_changes.tsv')
        assert before==after,'Final code changed individual molecule outcomes'
        for group in GROUPS:
            rows=raw[group];pp,mm=paths[group]
            summary=load_json(run/'comparison'/group/locus/'summary.json')
            latest=load_json(run/'final_replay/comparison'/group/locus/'summary.json')
            assert summary['graph_stats']==latest['graph_stats']
            latest_paths=path_members(run/'final_replay/comparison'/group/locus)
            assert {p['sha256'] for p in pp}=={p['sha256'] for p in latest_paths[0]}
            for rid,record in rows.items():
                assert record.get('sequence','')==record['query_sequence'][record.get('qstart',0):record.get('qend',0)].upper() if record.get('sequence') else True
            if group!='B0':
                final_records={r['record_id']:r for r in json_lines(run/'final_replay/comparison'/group/locus/'records.jsonl')}
                assert set(final_records)==set(rows)
                assert all(compact_evidence(rows[k])==compact_evidence(final_records[k]) for k in rows)
            replay.append(dict(locus=locus,group=group,identical_molecule_table=True,identical_evidence=True,identical_graph_stats=True))
            hashes={p['sha256'] for p in pp}
            paths_summary.append(dict(locus=locus,group=group,paths=len(pp),same_sequence_as_B0=len(hashes&baseline_hashes),
                added_sequence_paths=len(hashes-baseline_hashes),absent_old_sequence_paths=len(baseline_hashes-hashes),
                singleton_paths=sum(len(p['molecules'])==1 for p in pp),no_hifi_paths=sum(not p['hifi_molecules'] for p in pp),
                hifi_supported_paths=sum(bool(p['hifi_molecules']) for p in pp),
                changed_or_new_sequences_are_not_confirmed_alleles=True))
            counts=Counter()
            for rid,r in rows.items():
                identity=byid[rid]
                key=(identity['sample'],identity['platform'],mapq_bin(r['mapq']))
                counts[key+('input',)]+=1
                if r['molecule_key'] in mm:counts[key+('retained',)]+=1
                core=r.get('core_anchor_evidence',{})
                side_state='legacy'
                if core:
                    l=bool(core.get('left',{}).get('endpoints'));rr=bool(core.get('right',{}).get('endpoints'))
                    side_state='both_ordered' if core.get('status')=='success' else 'both_without_ordered_pair' if l and rr else 'left_only' if l else 'right_only' if rr else 'neither'
                status_counts.append(dict(locus=locus,group=group,record_id=rid,molecule=r['molecule_key'],
                    sample=identity['sample'],platform=identity['platform'],mapq=r['mapq'],mapq_bin=mapq_bin(r['mapq']),
                    input_status=r['status'],core_status=r.get('core_status','legacy'),side_state=side_state,
                    analysis_status=r.get('analysis_status','legacy'),extraction_status=r.get('extraction_status','legacy'),
                    in_P=r['molecule_key'] in mm,core_ambiguous=r.get('core_ambiguous',False)))
            for sample,platform,mq in sorted({k[:3] for k in counts}):
                strata.append(dict(locus=locus,group=group,sample=sample,platform=platform,mapq_bin=mq,
                    input_molecules=counts[(sample,platform,mq,'input')],retained_molecules=counts[(sample,platform,mq,'retained')]))
        picked,_=discovery_records(source)
        for group in GROUPS[3:]:
            model=load_json(run/'models'/locus/f'{group}.json')
            for sample in SAMPLES:
                left=set(model['sides']['left']['support']['members'][sample]);right=set(model['sides']['right']['support']['members'][sample])
                paired={r['molecule_key'] for r in picked if r['sample']==sample and raw[group][r['record_id']].get('analysis_status')=='success'}
                discovery.append(dict(locus=locus,group=group,sample=sample,discovery_molecules=sum(r['sample']==sample for r in picked),
                    left_support=len(left),right_support=len(right),intersection=len(left&right),ordered_pair_support=len(paired),
                    left_boundary=model['sides']['left']['boundary'],right_boundary=model['sides']['right']['boundary'],
                    fixed_window_mean_depth=model['depths'][sample]))
        for group in ['B1','B2','W1']:
            old,new=paths['B0'][1],paths[group][1]
            excluded=defaultdict(list)
            for r in read_table(run/'comparison'/group/locus/'excluded.tsv'):excluded[r['molecule']].append(r['reason'])
            for molecule in sorted(set(old)|set(new)):
                a,b=old.get(molecule),new.get(molecule);identity=bymolecule[molecule];rid=identity['record_id']
                if a and b and (a['qstart'],a['qend'])!=(b['qstart'],b['qend']):
                    key_cases.append(dict(locus=locus,group=group,record_id=rid,sample=identity['sample'],platform=identity['platform'],
                        mapq=identity['mapq'],old_cuts=[a['qstart'],a['qend']],new_cuts=[b['qstart'],b['qend']],
                        left_shift=b['qstart']-a['qstart'],right_shift=b['qend']-a['qend'],old_length=a['length'],new_length=b['length'],
                        old_path=a['path'],new_path=b['path'],core_candidates=[[p['qstart'],p['qend']] for p in raw[group][rid]['core_anchor_evidence']['pairs']],
                        old_sequence=raw['B0'][rid].get('sequence'),new_sequence=raw[group][rid].get('sequence'),
                        biological_repair_confirmed=False))
                if a and not b:
                    lost_cases.append(dict(locus=locus,group=group,record_id=rid,molecule=molecule,sample=identity['sample'],
                        platform=identity['platform'],mapq=identity['mapq'],old_path=a['path'],old_cuts=[a['qstart'],a['qend']],
                        new_cuts=[raw[group][rid].get('qstart'),raw[group][rid].get('qend')],
                        new_core_status=raw[group][rid].get('core_status'),new_state=raw[group][rid].get('extraction_status'),
                        original_cigar=identity['cigar'],query_length=len(identity['query_sequence']),
                        left_reason=raw[group][rid].get('core_anchor_evidence',{}).get('left',{}).get('status'),
                        right_reason=raw[group][rid].get('core_anchor_evidence',{}).get('right',{}).get('status'),
                        exclusion_reasons=sorted(set(excluded[molecule])),independent_correctness_label='unresolved'))
        control_groups=defaultdict(list)
        for r in source:
            if (r['record_id'] not in reasons or reasons[r['record_id']]==['stratified_unchanged_control']) and 20<=r['mapq']<255 and not r['flag']&(4|256|512|1024|2048):
                control_groups[(r['sample'],r['platform'],mapq_bin(r['mapq']))].append(r)
        for key,rr in control_groups.items():
            selected=[r['record_id'] for r in rr if 'stratified_unchanged_control' in reasons.get(r['record_id'],[])]
            control_sampling.append(dict(locus=locus,sample=key[0],platform=key[1],mapq_bin=key[2],eligible_controls=len(rr),
                selected=len(selected),sampling_fraction=len(selected)/len(rr),record_ids=selected))
        fl,fr=FETCH[locus]
        with IndexedFastaReference(args.project/'reference/GRCh38_full_analysis_set_plus_decoy_hla.fa') as reference:
            dna=reference.fetch(source[0]['chrom'],fl,fr).upper()
        blind=[]
        for review in reviews:
            rid=review['record_id'];record=byid[rid]
            case=sha256(('flank-blind-v1-'+locus+rid).encode()).hexdigest()[:16]
            order=sorted(['B0','B1','B2','W1'],key=lambda g:sha256((case+g).encode()).hexdigest())
            candidates={chr(65+i):compact_evidence(raw[g][rid]) for i,g in enumerate(order)}
            blind_key.extend(dict(locus=locus,case=case,display=chr(65+i),method=g,record_id=rid) for i,g in enumerate(order))
            # Method names and automatic labels are omitted from the review face.
            item=dict(case=case,sample=record['sample'],platform=record['platform'],mapq=record['mapq'],
                query_length=len(record['query_sequence']),query_sequence=record['query_sequence'],
                original_cigar=record['cigar'],reference_start=record['pos'],flag=record['flag'],
                reference_window_start=fl,reference_window_sequence=dna,alternatives=candidates)
            blind.append(item)
        blind.sort(key=lambda r:r['case'])
        (out/f'{locus}.blind_review.jsonl').write_text(''.join(json.dumps(x,ensure_ascii=True)+'\n' for x in blind),encoding='utf-8')
        payload=json.dumps(blind,ensure_ascii=True).replace('</','<\\/')
        page='''<!doctype html><meta charset="utf-8"><title>切点复核</title><style>
body{font-family:Arial,"Microsoft YaHei",sans-serif;margin:24px;color:#24364b}pre{white-space:pre-wrap;word-break:break-all;background:#f3f6fa;padding:12px}button,input,select{font:inherit;padding:7px}summary{padding:8px;cursor:pointer}td,th{padding:6px;border-bottom:1px solid #ddd}section{margin-top:16px}</style>
<h1>切点复核材料</h1><p>方法身份已隐藏，四种不同配置的展示顺序由固定哈希决定。相同结果的 W2–W8 可在完整原始审计中核查。材料中有逐字 SEQ、原 CIGAR、参考上下文、候选和可回放路径。方法分数与 CIGAR 一致性不构成独立真值。尚未完成独立生物学裁定。</p>
<input id="filter" placeholder="病例编号 / 样本 / 平台"><select id="cases"></select><div id="view"></div><script>
const rows=PAYLOAD; const sel=document.getElementById('cases'),filter=document.getElementById('filter'),view=document.getElementById('view');
const esc=x=>String(x).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function render(){const r=rows.find(x=>x.case===sel.value);if(!r){view.textContent='无匹配病例';return;}
view.innerHTML='<h2>'+esc(r.case)+'</h2><p>'+esc([r.sample,r.platform,'MAPQ '+r.mapq,'read 长 '+r.query_length].join(' · '))+'</p><p>原 CIGAR：'+esc(r.original_cigar)+'</p>'+
Object.entries(r.alternatives).map(([k,v])=>'<section><h3>候选 '+k+'</h3><p>核心切点 '+esc(JSON.stringify([v.qstart,v.qend]))+' · 状态 '+esc(v.extraction_status||'未提供分层状态')+'</p><pre>'+esc(v.sequence||'无核心片段')+'</pre><details><summary>候选与路径证据</summary><pre>'+esc(JSON.stringify(v,null,2))+'</pre></details></section>').join('')+
'<details><summary>完整原始 read 与参考窗口</summary><pre>'+esc(JSON.stringify({query:r.query_sequence,reference_start:r.reference_window_start,reference:r.reference_window_sequence},null,2))+'</pre></details>';}
function list(){const q=filter.value.toLowerCase();sel.innerHTML=rows.filter(r=>[r.case,r.sample,r.platform].join(' ').toLowerCase().includes(q)).map(r=>'<option value="'+r.case+'">'+esc(r.case+' '+r.sample+' '+r.platform)+'</option>').join('');render();}filter.oninput=list;sel.onchange=render;list();</script>'''.replace('PAYLOAD',payload)
        (out/f'{locus}.blind_review.html').write_text(page,encoding='utf-8')
        details=json_lines(run/'audit'/locus/'sensitivity_records.jsonl')
        aggregates=defaultdict(list)
        for row in details:
            record=byid[row['record_id']];base=row['conditions'][0];moves=[];candidate_changed=False
            for config,r in zip(PERTURBATIONS,row['conditions']):
                candidates=r['candidates'];strings={record['query_sequence'][a:b] for a,b in candidates}
                spans=[max((p[i] for p in candidates),default=0)-min((p[i] for p in candidates),default=0) for i in (0,1)]
                domains={}
                if row['method']!='legacy100':
                    length=100 if row['method']=='affine100' else 250;a,b=row['target']
                    for i,(side,interval) in enumerate([('left',[a-length,a]),('right',[b,b+length])]):
                        domains[side]=flank_search(record,interval,side,config['padding'][i],config['shifts'][i])
                changed=(r['cuts']!=base['cuts']) if base['cuts'] is not None and r['cuts'] is not None else None
                move=max(abs(r['cuts'][i]-base['cuts'][i]) for i in (0,1)) if changed is not None else None
                if move is not None:moves.append(move)
                candidate_changed |= candidates!=base['candidates']
                sensitivity_cases.append(dict(locus=locus,record_id=row['record_id'],sample=row['sample'],platform=row['platform'],
                    target=row['target'],method=row['method'],condition=r['condition'],status=r['status'],cuts=r['cuts'],
                    changed_vs_base=changed,maximum_move=move,candidate_count=len(candidates),candidate_span=spans,
                    distinct_candidate_strings=len(strings),candidate_set_changed=candidates!=base['candidates'],search_domains=domains))
            aggregates[(tuple(row['target']),row['platform'],row['method'])].append(dict(
                base_success=base['status']=='success',all_success=all(r['status']=='success' for r in row['conditions']),
                any_output_lost=base['status']=='success' and any(r['status']!='success' for r in row['conditions']),
                any_output_recovered=base['status']!='success' and any(r['status']=='success' for r in row['conditions']),
                changed_any_comparable=any(v>0 for v in moves),max_move=max(moves,default=None),candidate_set_changed=candidate_changed))
        for (target,platform,method),items in aggregates.items():
            sensitivity_stats.append(dict(locus=locus,target=target,platform=platform,method=method,quality_records=len(items),
                **{k:sum(x[k] for x in items) for k in ['base_success','all_success','any_output_lost','any_output_recovered','changed_any_comparable','candidate_set_changed']},
                maximum_observed_move=max((x['max_move'] for x in items if x['max_move'] is not None),default=None)))
        all_stats.append(dict(locus=locus,review_cases=len(reviews),changed_case_counts=dict(Counter(t for r in reviews for t in r['triggers'])),
            labels=dict(Counter(r['boundary_label'] for r in reviews)),known_independent_repair=0,independent_biological_review='not_available'))
    for name,rows in [('stratified_retention',strata),('record_states',status_counts),('discovery_support',discovery),
          ('path_summary',paths_summary),('changed_common_P',key_cases),('lost_B0_molecules',lost_cases),('review_control_sampling',control_sampling),
          ('case_review_complete',complete_review),('blind_review_key',blind_key),('sensitivity_extended',sensitivity_cases),
          ('sensitivity_full_denominator',sensitivity_stats),('final_replay_checks',replay)]:
        table(out/f'{name}.tsv',rows)
    for name,rows in [('path_summary',paths_summary),('changed_common_P',key_cases),('lost_B0_molecules',lost_cases),
          ('review_summary',all_stats),('sensitivity_full_denominator',sensitivity_stats),('discovery_support',discovery)]:dump(out/f'{name}.json',rows)
    print('REPORT_EVIDENCE_COMPLETE',json.dumps(all_stats),flush=True)


if __name__=='__main__':main()
