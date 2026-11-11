"""Files and static view for the structure-first mode of vntrpathgraph_cli."""
from collections import Counter, defaultdict
from hashlib import sha256
from html import escape
import json
import csv

from vntrreads import write_table
from vntrpathgraph import write_gfa, validate_gfa, graph_support, digest
from vntrpathstructure import (VERSION, annotate_paths, merge_paths, build_structural_graph,
                              select_subgraph, restore_original)
from vntrpartition import repeat_type,str_policy_path

def write_path_sequences_txt(path,records,metadata):
    """Write complete original strings, using exactly the GFA path names."""
    old=metadata.get('previous_name_by_sha256',{})
    with path.open('w',encoding='utf-8',newline='') as handle:
        handle.write('# P0001 is the reference core; observation paths start at P0002.\n')
        handle.write('# Admission: complete B2 flank evidence; MAPQ/FLAG retained as provenance.\n')
        handle.write('# Coordinates are 0-based half-open; sequences exclude flank templates.\n')
        reference=next((r for r in records if r.get('role')=='reference'),None)
        if reference: handle.write('# Reference interval: '+json.dumps(reference['reference_interval'])+'\n')
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(['name','role','length_bp','molecules','hifi_molecules','old_name','sequence'])
        for r in records:
            writer.writerow([r['id'],r.get('role','observation'),len(r['sequence']),len(r['molecules']),
                             len(r['hifi_molecules']),old.get(r['sha256'],'') if r.get('role')!='reference' else '',r['sequence']])


def export_structure(destination,locus,records,excluded,metadata,model_path,policy_path,budget,
                     preview_limit,dot,imported_rows=None,imported_audit=None):
    from vntrpathgraph_cli import dump, render_preview
    policy=json.loads(policy_path.read_text(encoding='utf-8'))
    if policy.get('schema_version')==4 and repeat_type(json.loads(model_path.read_text(encoding='utf-8')))=='STR':
        return export_structure(destination,locus,records,excluded,metadata,model_path,
            str_policy_path(model_path,locus.name),0,preview_limit,dot,imported_rows,imported_audit)
    if policy.get('schema_version')==4:
        from run_hmm_display import export_hmm
        result=export_hmm(destination,locus,records,model_path,budget,render=bool(dot),dot=dot)
        result.update(locus=locus.name,status='complete',graph_stats=result['validation'],
            substituted_bases=result['substitutions'],preview={'paths':sorted(r['id'] for r in records)},metadata=metadata)
        dump(destination/'structure_policy.json',policy);dump(destination/'summary.json',result)
        return result
    if policy.get('schema_version')==3:
        from vntrlocal import select_interpretations,plan_edits,VERSION as LOCAL_VERSION
        from vntrlocal_export import export_budget,render_graph
        model,model_hash,features=annotate_paths(records,model_path,locus,select_representative=True,export_optimal=True)
        if policy['version']!=LOCAL_VERSION or policy['model_version']!=model['version']:raise ValueError('v3 policy/model mismatch')
        tokens,slots,alignments,regions,selection=select_interpretations(features)
        ts,views,represented,losses,_,decisions=plan_edits(records,tokens,slots,model,budget)
        result=export_budget(destination,records,tokens,ts,slots,views,represented,losses,budget,decisions,render=False,model=model)
        display=json.loads((destination/'display_graph.json').read_text(encoding='utf-8'));render_graph(destination,display,dot)
        for name,data in [('source_members',records),('locus_model',model),('structure_policy',policy),
            ('selected_annotations',dict(tokens=tokens,slots=slots)),('alignment_choices',alignments),
            ('annotation_choices',regions),('interpretation_selection',selection),('original_annotations',features)]:dump(destination/(name+'.json'),data)
        result.update(locus=locus.name,status='complete',graph_stats=result['validation'],
            substituted_bases=result['substitutions'],preview={'paths':sorted(represented)},metadata=metadata,model_sha256=model_hash)
        dump(destination/'summary.json',result)
        return result
    if policy.get('schema_version')==2:
        from vntrdisplay_loss import build_unit_graph
        from vntrdisplay_export import export_budget
        model,_,features=annotate_paths(records,model_path,locus,select_representative=True)
        exact=build_unit_graph(records,features,model)
        result=export_budget(destination,records,features,model,policy,budget,exact,metadata,dot)
        result.update(locus=locus.name,status='complete',graph_stats=result['validation'],
                      substituted_bases=result['substitutions'],preview={'paths':list(exact['paths'])})
        dump(destination/'summary.json',result)
        return result
    model,model_hash,features=annotate_paths(records,model_path,locus,imported_rows,imported_audit)
    policy=json.loads(policy_path.read_text(encoding='utf-8'))
    represented,reps,losses,decisions=merge_paths(records,features,model,policy,budget)
    originals={r['id']:r['sequence'] for r in records}
    backbone=metadata.get('backbone_path')
    exact=build_structural_graph(records,features,model,originals,{p:p for p in originals},backbone)
    structure=build_structural_graph(records,features,model,represented,reps,backbone) if losses else exact
    selected=records if preview_limit==0 else records[:preview_limit]
    preview=select_subgraph(structure,[r['id'] for r in selected])
    validation={}
    for name,graph,seqs in [('exact',exact,originals),('structure',structure,represented),
                           ('preview',preview,{r['id']:represented[r['id']] for r in selected})]:
        path=destination/f'{locus.name}.{name}.gfa'
        write_gfa(path,graph,records,'structure_first_'+name,seqs,[] if name=='exact' else losses)
        validation[name]=validate_gfa(path,seqs)
    loss_by_id={x['path']:x for x in losses}
    for p,sequence in represented.items():
        if restore_original(sequence,loss_by_id.get(p,{}).get('events',[]))!=originals[p]:
            raise ValueError('Original sequence not recovered: '+p)
        if features[p]['signature']!=features[reps[p]]['signature']: raise ValueError('Protected structure changed')
    for name,data in [('original',originals),('represented',represented)]:
        (destination/f'paths.{name}.fa').write_text(''.join(f'>{p}\n{s}\n' for p,s in data.items()),encoding='ascii')
    groups=defaultdict(list)
    for p in represented: groups[(reps[p],digest(represented[p]))].append(p)
    group_id={p:f'G{i:04d}' for i,(_,members) in enumerate(sorted(groups.items()),1) for p in members}
    path_rows=[]
    for r in records:
        p=r['id']; f=features[p]; loss=loss_by_id.get(p,{})
        path_rows.append(dict(path=p,role=r.get('role','observation'),sequence_sha256=r['sha256'],length=len(originals[p]),
            represented_length=len(represented[p]),molecules=len(r['molecules']),hifi_molecules=len(r['hifi_molecules']),
            structure_group=group_id[p],representative=reps[p],expression=f['expression'],
            count_status=f['count_status'],complete_units=f['complete_units'],interruptions=f['interruptions'],
            loss_bp=loss.get('loss_bp',0),substitutions=loss.get('substitutions',0),insertions=loss.get('insertions',0),
            deletions=loss.get('deletions',0),loss_fraction=loss.get('fraction',0),
            exact_walk=exact['paths'][p],structure_walk=structure['paths'][p]))
    path_fields=['path','role','sequence_sha256','length','represented_length','molecules','hifi_molecules','structure_group',
                 'representative','expression','count_status','complete_units','interruptions','loss_bp','substitutions',
                 'insertions','deletions','loss_fraction','exact_walk','structure_walk']
    write_table(destination/'paths.tsv',path_rows,path_fields)
    write_table(destination/'excluded.tsv',excluded,['molecule','fragment','reason'])
    write_table(destination/'members.tsv',[dict(path=r['id'],**{k:v for k,v in o.items() if k!='sequence'})
        for r in records for o in r['observations']],['path','molecule','sample','platform','namespace','read_id','row_id',
        'record_id','fragment','bam','qstart','qend','mapq','flag','alignment_type','flags','anchor_scores','strand','reference_interval','boundary_status']+
        (['extraction_version','boundary_model_sha256','extraction_status','core_status',
          'core_ambiguous','threshold_path_dependent',
          'normalization_reason','core_candidate_pairs'] if metadata.get('extraction_version') else []))
    write_table(destination/'merge_decisions.tsv',decisions,['path','target_path','accepted','reason'])
    write_table(destination/'losses.tsv',losses,['path','target_path','loss_bp','substitutions','insertions','deletions',
                'original_length','represented_length','fraction','positions','events','policy_version','reason'])
    (destination/'losses.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in losses),encoding='utf-8')
    (destination/'path_features.jsonl').write_text(''.join(json.dumps(dict(path=p,**f),ensure_ascii=False)+'\n'
        for p,f in features.items()),encoding='utf-8')
    nstats,estats=graph_support(structure,records)
    write_table(destination/'nodes.tsv',[dict(id=n,length=len(s),sequence=s,**structure['node_features'][n],**nstats[n])
        for n,s in sorted(structure['nodes'].items())],['id','length','sequence','expression','unit_count','uncertain',
        'unit_keys','paths','molecules','visits'])
    write_table(destination/'edges.tsv',[dict(source=a,target=b,**v) for (a,b),v in sorted(estats.items())],
                ['source','target','paths','molecules'])
    dump(destination/'source_members.json',records)
    write_path_sequences_txt(destination/f'{locus.name}.path_sequences.txt',records,metadata)
    dump(destination/'locus_model.json',model);dump(destination/'structure_policy.json',policy)
    dump(destination/'alignment_audit.json',structure['alignments'])
    group_rows=[]; by_id={r['id']:r for r in records}
    for (rep,_),ps in sorted(groups.items()):
        group_rows.append(dict(group=group_id[rep],representative=rep,paths=ps,
            molecules=len(set().union(*(set(by_id[p]['molecules']) for p in ps))),
            representative_exact_molecules=len(by_id[rep]['molecules'])))
    write_table(destination/'structure_groups.tsv',group_rows,
                ['group','representative','paths','molecules','representative_exact_molecules'])
    preview_report=render_preview(destination,locus.name,preview,selected,records,losses,dot) if records else dict(status='empty',paths=[])
    molecules=set().union(*(set(r['molecules']) for r in records)) if records else set()
    summary=dict(version=VERSION,status='complete' if records else 'no_eligible_complete_paths',locus=locus.name,
        paths=len(records),reference_paths=sum(r.get('role')=='reference' for r in records),
        observed_paths=sum(r.get('role')!='reference' for r in records),
        backbone=structure.get('backbone_path'),molecules=len(molecules),structure_routes=len(groups),graph_stats=validation,
        lossy_paths=len(losses),loss_bp=sum(x['loss_bp'] for x in losses),substituted_bases=sum(x['substitutions'] for x in losses),
        inserted_bases=sum(x['insertions'] for x in losses),deleted_bases=sum(x['deletions'] for x in losses),max_loss_bp=budget,
        annotation_status=dict(Counter(f['count_status'] for f in features.values())),model_sha256=model_hash,
        policy_sha256=sha256(policy_path.read_bytes()).hexdigest(),preview=dict(preview_report,node_ids='same_as_full_structure_graph',
        selected_molecules=len(set().union(*(set(r['molecules']) for r in selected))) if selected else 0,
        total_molecules=len(molecules),omitted_paths=len(records)-len(selected)),metadata=metadata,
        validation=dict(original_reconstruction=True,protected_structure=True,molecule_identity=True),
        limitations=['Fixed-model annotations and comparison positions are not biological truth or unique homology.',
            'P0001 is reference when role=reference; other P records are observations; other traversable combinations are not observations.',
            'Display loss is not sequencing error correction; original and represented lengths are distinct.',
            'Source BAM selection is unchanged; no genome-wide competitive mapping or genotype claim.'])
    dump(destination/'summary.json',summary);dump(destination/'validation.json',validation)
    table=''.join('<tr>'+''.join('<td>'+escape(str(r[k]))+'</td>' for k in
        ['path','structure_group','expression','length','represented_length','molecules','loss_bp','count_status'])+'</tr>' for r in path_rows)
    html=f'''<!doctype html><meta charset="utf-8"><title>{locus.name} 结构图</title>
<style>body{{font-family:Arial,"Microsoft YaHei",sans-serif;margin:30px;color:#253748}}img{{max-width:100%}}td,th{{padding:8px;border-bottom:1px solid #ddd;text-align:left}}p{{line-height:1.7}}</style>
<h1>{locus.name}：重复段与模块结构图</h1>
<p>{len(records)} 条 path（{summary['reference_paths']} 条参考、{summary['observed_paths']} 条观测），{len(molecules)} 个独立分子。累计隐藏 {summary['loss_bp']} bp。</p>
<p>静态图展示 {len(selected)} 条 path；完整 GFA 保留全部路径。~ 表示实际单元与模板不同，未必被有损替换；* 表示路线有明确记录的表示损失。UNRESOLVED 是划分未定，ALIGNMENT ? 是区块对应未定。</p>
<a href="{locus.name}.preview.svg"><img src="{locus.name}.preview.svg"></a>
<p><a href="{locus.name}.structure.gfa">结构 GFA</a> · <a href="{locus.name}.path_sequences.txt">名称与完整序列 TXT</a> · <a href="paths.tsv">路径表</a> · <a href="structure_groups.tsv">路线成员与精确支持</a></p>
<table><tr><th>Path</th><th>路线</th><th>原始结构表达</th><th>原始 bp</th><th>表示 bp</th><th>分子</th><th>损失 bp</th><th>注释状态</th></tr>{table}</table>
<p>图中节点存固定实际碱基。P0001 为参考，不增加分子支持；其余完整 P 为观测，其他可走组合不代表额外序列。侧翼通过不代表全基因组唯一定位；模型单元数不是基因型或已验证拷贝数。</p>'''
    if metadata.get('extraction_version'):
        html+='<p><a href="boundary_model.json">冻结边界模型</a> · <a href="extraction_audit.jsonl">切点及歧义证据</a></p>'
    (destination/'index.html').write_text(html,encoding='utf-8')
    return summary
