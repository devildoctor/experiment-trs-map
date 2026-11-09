#!/usr/bin/env python3
"""Build literal-sequence GFA graphs from complete observed locus sequences."""
import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
from hashlib import sha256
from html import escape
import json
from pathlib import Path
import shutil
import subprocess
import textwrap

from vntrgraph import IndexedFastaReference, read_catalog, discover_bams, load_source_namespaces, locus_output_dir
from vntrreads import collect_reads, write_table
from vntrpathgraph import (VERSION, select_paths, build_graph, compact_graph, apply_merge_rules,
                           graph_support, write_gfa, validate_gfa, digest, assign_reference_paths)


def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')


def sequence_label(sequence):
    for size in range(1, min(32, len(sequence)//2)+1):
        if len(sequence)%size == 0 and sequence == sequence[:size]*(len(sequence)//size):
            return f'{sequence[:size]} x{len(sequence)//size}'
    return sequence if len(sequence) <= 26 else sequence[:12]+'...'+sequence[-10:]


def number_ranges(values):
    groups=[]; start=last=values[0]
    for n in values[1:]+[None]:
        if n==last+1: last=n; continue
        groups.append(str(start) if start==last else f'{start}-{last}')
        start=last=n
    return ','.join(groups)


def render_preview(destination, locus, graph, records, all_records, losses, dot):
    """Draw every supplied node and edge on one continuous canvas."""
    names=sorted(graph['paths'])
    nstats,estats=graph_support(graph,records)
    reference=next((r['id'] for r in records if r.get('role')=='reference'),None)
    reference_nodes=set(graph['paths'].get(reference,[]))
    reference_edges=set(zip(graph['paths'].get(reference,[]),graph['paths'].get(reference,[])[1:]))
    observed=sum(r.get('role')!='reference' for r in records)
    molecules=set().union(*(set(r['molecules']) for r in records)) if records else set()
    title=f'{locus} | {len(records)}/{len(all_records)} paths | {observed} observed strings | {len(molecules)} molecules'
    if reference: title+='\nP0001 = reference (blue); reference adds zero molecule support.'
    title+='\nAll supplied paths are drawn. Hover nodes/edges for path names; see sequence TXT.'
    if graph.get('structure_mode'):
        title+='\n~ = differs from template; UNRESOLVED = partition uncertain; ALIGNMENT ? = correspondence uncertain.'
        legend=graph.get('dictionary_legend',[])
        title+='\n'+'\n'.join(' ; '.join(legend[i:i+3]) for i in range(0,len(legend),3))
    lines=['digraph G {',
           'graph [rankdir=LR, bgcolor="white", nodesep=0.18, ranksep=0.3, fontname="DejaVu Sans", fontsize=13];',
           'node [shape=box, style="rounded,filled", fillcolor="#EEF3F7", color="#526978", fontname="DejaVu Sans", fontsize=10];',
           'edge [fontname="DejaVu Sans", fontsize=8, arrowsize=0.5];',
           'graph [labelloc=t, label='+json.dumps(title)+'];']
    for n,bases in sorted(graph['nodes'].items()):
        expression=graph.get('labels',{}).get(n,sequence_label(bases))
        label=f'{n} | {len(bases)} bp\n'+'\n'.join(textwrap.wrap(expression,42))
        uncertain=graph.get('node_features',{}).get(n,{}).get('uncertain')
        fill='#FFF4D6' if uncertain else '#EEF3F7'
        stroke='#0072B2' if n in reference_nodes else '#526978'
        tooltip=f'{len(nstats[n]["paths"])} paths; {nstats[n]["molecules"]} molecules: '+','.join(nstats[n]['paths'])
        lines.append(f's{n} [label='+json.dumps(label)+', fillcolor='+json.dumps(fill)+
                     ', color='+json.dumps(stroke)+', penwidth='+('2.5' if n in reference_nodes else '1')+
                     ', tooltip='+json.dumps(tooltip)+'];')
    for (a,b),v in sorted(estats.items()):
        ids=v['paths'];is_reference=(a,b) in reference_edges
        label=','.join(ids) if len(names)<=10 else f'{len(ids)} P'
        tooltip=f'{v["molecules"]} molecules: '+','.join(ids)
        lines.append(f's{a} -> s{b} [color='+json.dumps('#0072B2' if is_reference else '#8A949F')+
                     ', penwidth='+('2.5' if is_reference else '1')+', label='+json.dumps(label)+
                     ', tooltip='+json.dumps(tooltip)+'];')
    lines.append('}')
    path=destination/f'{locus}.preview.dot'
    path.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    if not dot or not shutil.which(dot): return dict(status='dot_unavailable',paths=names)
    subprocess.run([dot,'-Tsvg',str(path),'-o',str(destination/f'{locus}.preview.svg')],
                   check=True,timeout=600,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    shutil.copyfile(destination/f'{locus}.preview.svg',destination/f'{locus}.preview.linear.svg')
    subprocess.run([dot,'-Gsize=30,30','-Gdpi=96','-Tpng',str(path),'-o',str(destination/f'{locus}.preview.png')],
                   check=True,timeout=600,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    return dict(status='rendered',paths=names,layout='single_canvas',panels=1)


def export_graphs(destination, locus, records, excluded, metadata, preview_limit=0,
                  merge_policy=None, loss_budget=0, dot='dot'):
    originals={r['id']:r['sequence'] for r in records}
    represented, losses=apply_merge_rules(records,merge_policy,loss_budget,locus)
    # 首条为按 HiFi 分子数、总分子数、长度及哈希确定的实际观测骨架。
    backbone=next((r['sequence'] for r in records if r.get('role')=='reference'),records[0]['sequence'] if records else '')
    raw=build_graph(originals,backbone)
    exact=compact_graph(raw)
    structure=compact_graph(build_graph(represented,backbone)) if losses else exact
    validation={}
    for name,graph,expected in [('base',raw,originals),('exact',exact,originals),('structure',structure,represented)]:
        path=destination/f'{locus}.{name}.gfa'
        mode='structure_lossy' if name=='structure' and losses else name+'_lossless'
        write_gfa(path,graph,records,mode,expected)
        validation[name]=validate_gfa(path,expected)
    from vntrpathstructure_export import write_path_sequences_txt
    write_path_sequences_txt(destination/f'{locus}.path_sequences.txt',records,metadata)
    (destination/'paths.original.fa').write_text(''.join(f'>{p}\n{s}\n' for p,s in originals.items()),encoding='ascii')
    (destination/'paths.represented.fa').write_text(''.join(f'>{p}\n{s}\n' for p,s in represented.items()),encoding='ascii')
    loss_by_id={x['path']:x for x in losses}
    groups=defaultdict(list)
    for p,s in represented.items(): groups[digest(s)].append(p)
    group_id={p:f'G{i:04d}' for i,(_,members) in enumerate(sorted(groups.items()),1) for p in members}
    path_rows=[]
    for r in records:
        path_rows.append(dict(path=r['id'],sequence_sha256=r['sha256'],length=len(r['sequence']),
            molecules=len(r['molecules']),hifi_molecules=len(r['hifi_molecules']),samples=r['samples'],platforms=r['platforms'],
            structure_group=group_id[r['id']],loss_bp=loss_by_id.get(r['id'],{}).get('substitutions',0),
            loss_fraction=loss_by_id.get(r['id'],{}).get('fraction',0),
            exact_walk=','.join(map(str,exact['paths'][r['id']])),structure_walk=','.join(map(str,structure['paths'][r['id']]))))
    columns=['path','sequence_sha256','length','molecules','hifi_molecules','samples','platforms','structure_group','loss_bp','loss_fraction','exact_walk','structure_walk']
    write_table(destination/'paths.tsv',path_rows,columns)
    write_table(destination/'members.tsv',[dict(path=r['id'],**{k:v for k,v in o.items() if k!='sequence'})
                for r in records for o in r['observations']],
                ['path','molecule','sample','platform','read_id','namespace','bam','row_id','record_id','fragment',
                 'qstart','qend','mapq','anchor_scores','strand','reference_interval','boundary_status']+
                (['extraction_version','boundary_model_sha256','extraction_status','core_status',
                  'core_ambiguous','threshold_path_dependent',
                  'normalization_reason','core_candidate_pairs'] if metadata.get('extraction_version') else []))
    write_table(destination/'excluded.tsv',excluded,['molecule','fragment','reason'])
    write_table(destination/'losses.tsv',losses,['path','target_path','substitutions','fraction','positions','source_bases','target_bases','reason'])
    nstats,estats=graph_support(structure,records)
    write_table(destination/'nodes.tsv',[dict(id=n,length=len(s),sequence=s,**nstats[n],
                alignment_columns=structure['coordinates'][n]) for n,s in sorted(structure['nodes'].items())],
                ['id','length','sequence','paths','molecules','visits','alignment_columns'])
    write_table(destination/'edges.tsv',[dict(source=a,target=b,**v) for (a,b),v in sorted(estats.items())],
                ['source','target','paths','molecules'])
    selected=records if preview_limit==0 else records[:preview_limit]
    preview_seq={r['id']:represented[r['id']] for r in selected}
    preview_graph=compact_graph(build_graph(preview_seq,backbone))
    write_gfa(destination/f'{locus}.preview.gfa',preview_graph,selected,'preview_selected_paths',preview_seq)
    validation['preview']=validate_gfa(destination/f'{locus}.preview.gfa',preview_seq)
    preview=render_preview(destination,locus,preview_graph,selected,records,losses,dot) if selected else dict(status='empty',paths=[])
    selected_molecules=set().union(*(set(r['molecules']) for r in selected)) if selected else set()
    all_molecules=set().union(*(set(r['molecules']) for r in records)) if records else set()
    summary=dict(version=VERSION,locus=locus,status='complete' if records else 'no_eligible_complete_paths',
        paths=len(records),structure_routes=len(groups),molecules=len(all_molecules),
        hifi_molecules=len(set().union(*(set(r['hifi_molecules']) for r in records))) if records else 0,
        excluded_fragments=len(excluded),excluded_reasons=dict(Counter(reason for x in excluded for reason in x['reason'].split(','))),
        backbone=next((r['id'] for r in records if r.get('role')=='reference'),records[0]['id'] if records else None),backbone_sha256=digest(backbone) if records else None,
        graph_stats=validation,lossy_paths=len(losses),substituted_bases=sum(x['substitutions'] for x in losses),
        maximum_loss_fraction=max((x['fraction'] for x in losses),default=0),loss_budget=loss_budget,
        preview=dict(preview,selected_molecules=len(selected_molecules),total_molecules=len(all_molecules),
                     omitted_paths=len(records)-len(selected),
                     node_ids='local_to_preview_gfa_not_full_graph'),
        mapping_context='supplied_BAM_only_not_verified_complete',
        alignment_model='observed_backbone_center_star_global_edit_distance; insertion nodes keyed by boundary and inserted length',
        alignment_ambiguity='one optimal alignment selected; equivalent alignments not enumerated',
        observation_definition='all independently B2-complete valid core strings; MAPQ and flags are provenance; reference P0001 adds zero molecule support',
        graph_walk_semantics='P role distinguishes reference from observation; recombining edges does not establish a new observation',
        metadata=metadata)
    dump(destination/'validation.json',validation);dump(destination/'summary.json',summary)
    dump(destination/'source_members.json',records)
    dump(destination/'alignment_cigars.json',raw['alignments'])
    if merge_policy: dump(destination/'merge_rules.json',merge_policy)
    preview_names=set(preview['paths'])
    table=''.join('<tr>'+''.join('<td>'+escape(str(x))+'</td>' for x in
                 [r['path'],'是' if r['path'] in preview_names else '',r['length'],r['molecules'],r['hifi_molecules'],r['structure_group'],r['loss_bp']])+'</tr>' for r in path_rows)
    html=f'''<!doctype html><meta charset="utf-8"><title>{escape(locus)} 结构图</title>
<style>body{{font-family:Arial,"Microsoft YaHei",sans-serif;max-width:1280px;margin:28px auto;padding:0 18px;color:#1f2937}}img{{width:100%;height:auto}}table{{border-collapse:collapse}}td,th{{padding:7px 14px;border-bottom:1px solid #ddd;text-align:left}}a{{color:#075d91}}p{{line-height:1.7}}</style>
<h1>{escape(locus)} 序列结构图</h1>
<p>{len(records)} 条完整区间序列 path，{len(all_molecules)} 个去重分子，{len(groups)} 条结构路线。人工压缩替换 {summary['substituted_bases']} bp。</p>
<p>静态图展示 {len(selected)} / {len(records)} 条 path，覆盖 {len(selected_molecules)} / {len(all_molecules)} 个合格分子；完整 GFA 保留全部合格路径。预览节点编号仅对应 preview.gfa。</p>
<p><a href="{locus}.structure.gfa">结构 GFA</a> · <a href="{locus}.exact.gfa">无损 GFA</a> · <a href="{locus}.preview.gfa">图中子图 GFA</a> · <a href="paths.tsv">全部路径表</a> · <a href="nodes.tsv">节点碱基</a> · <a href="excluded.tsv">排除记录</a></p>
'''
    if preview['status']=='rendered':html+=f'<a href="{locus}.preview.svg"><img src="{locus}.preview.svg" alt="{locus} 静态序列图，节点为实际碱基片段，边标出经过的路径"></a>'
    html+='<table><tr><th>Path</th><th>图中</th><th>长度 bp</th><th>分子数</th><th>HiFi 分子数</th><th>结构组</th><th>替换 bp</th></tr>'+table+'</table>'
    html+='<p>一个最优比对用于组织图的对应位置，并列对应未枚举；序列重建通过不证明划分或同源位置唯一。输入来自已有位点 BAM，未完成全基因组竞争定位核查。图中其他可走通组合不等于实际观测，不输出基因型。</p>'
    if metadata.get('extraction_version'):
        html+='<p><a href="boundary_model.json">冻结边界模型</a> · <a href="extraction_audit.jsonl">切点及歧义证据</a></p>'
    (destination/'index.html').write_text(html,encoding='utf-8')
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    inputs=p.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--bam-dir',type=Path);inputs.add_argument('--structures-json',type=Path)
    p.add_argument('--reference',type=Path);p.add_argument('--catalog',type=Path)
    p.add_argument('--locus',required=True);p.add_argument('--out',required=True,type=Path)
    p.add_argument('--samtools',default='samtools');p.add_argument('--sources',type=Path)
    p.add_argument('--platform',action='append',choices=['HiFi','ONT','ONTUL'])
    p.add_argument('--preview-paths',type=int,default=0,help='0 draws all paths; positive values explicitly limit display')
    p.add_argument('--merge-rules',type=Path);p.add_argument('--loss-budget',type=float,default=0)
    p.add_argument('--graph-mode',choices=['sequence','structure'],default='sequence')
    p.add_argument('--locus-model',type=Path);p.add_argument('--structure-policy',type=Path)
    p.add_argument('--max-loss-bp',type=int,default=0,choices=[0,1,2])
    p.add_argument('--dot',default='dot')
    p.add_argument('--previous-members',type=Path,help='Optional previous source_members.json for old-name mapping')
    args=p.parse_args()
    if args.preview_paths<0 or not 0<=args.loss_budget<=1:p.error('Invalid limits')
    if not args.reference or not args.catalog:p.error('B2 needs --reference and --catalog, including imports')
    from vntrflankalign import build_context
    selected=[l for l in read_catalog(args.catalog) if l.name==args.locus]
    if len(selected)!=1:p.error('Unknown B2 locus')
    with IndexedFastaReference(args.reference) as reference:
        flank_context=build_context(selected[0],reference)
        reference_sequence=reference.fetch(selected[0].chrom,selected[0].start,selected[0].end).upper()
    if args.graph_mode=='structure':
        if not args.locus_model or not args.structure_policy:p.error('Structure mode needs --locus-model and --structure-policy')
        structure_policy=json.loads(args.structure_policy.read_text(encoding='utf-8'))
        if structure_policy.get('schema_version')==2:
            if args.merge_rules or args.max_loss_bp or args.loss_budget>.05:p.error('v2 uses --loss-budget in [0, .05]')
        elif args.merge_rules or args.loss_budget:p.error('Use --max-loss-bp for v1 structure mode')
    elif args.locus_model or args.structure_policy or args.max_loss_bp:
        p.error('Structure options require --graph-mode structure')
    destination=locus_output_dir(args.out,args.locus);destination.mkdir(parents=True,exist_ok=False)
    dump(destination/'summary.json',dict(status='running',locus=args.locus))
    if args.structures_json:
        payload=json.loads(args.structures_json.read_text(encoding='utf-8'))
        if payload['locus']['name']!=args.locus:raise ValueError('Input locus mismatch')
        rows=payload['rows'];metadata=dict(locus=payload['locus'],source=str(args.structures_json),
                    source_sha256=sha256(args.structures_json.read_bytes()).hexdigest())
        if flank_context:
            if payload.get('audit',{}).get('boundary_model_sha256')!=flank_context.model_hash:
                raise ValueError('Imported flank boundary model mismatch')
            if any(f.get('extraction_version')!=flank_context.model['extraction_version'] or
                   f.get('boundary_model_sha256')!=flank_context.model_hash for r in rows for f in r['fragments']):
                raise ValueError('Imported fragment extraction version mismatch')
        if args.graph_mode=='structure':
            from vntrgraph import Locus
            locus=Locus(**payload['locus'])
    else:
        selected=[l for l in read_catalog(args.catalog) if l.name==args.locus]
        if len(selected)!=1:raise ValueError('Unknown locus')
        locus=selected[0];bams=discover_bams(args.bam_dir,args.locus)
        if not bams:raise ValueError('No locus BAMs')
        namespaces=load_source_namespaces(args.sources)
        with IndexedFastaReference(args.reference) as reference:
            with (destination/'records.jsonl').open('w',encoding='utf-8') as records_file:
                rows,audits=collect_reads(bams,locus,reference,args.samtools,records_file,namespaces,flank_context)
        metadata=dict(locus=asdict(locus),bam_audits=audits,bam_count=len(bams),
            reference=str(args.reference),catalog_sha256=sha256(args.catalog.read_bytes()).hexdigest(),
            input_files=[dict(path=str(b),bytes=b.stat().st_size,mtime_ns=b.stat().st_mtime_ns) for _,_,b in bams])
    records,excluded=select_paths(rows,set(args.platform) if args.platform else None)
    records=assign_reference_paths(records,reference_sequence,dict(chrom=selected[0].chrom,start=selected[0].start,end=selected[0].end))
    metadata.update(backbone_path='P0001',reference_core_sha256=digest(reference_sequence),
                    reference=str(args.reference),admission='all_independently_B2_complete')
    if args.previous_members:
        previous=json.loads(args.previous_members.read_text(encoding='utf-8'))
        metadata['previous_name_by_sha256']={r['sha256']:r['id'] for r in previous if r.get('role')!='reference'}
    if flank_context:
        from vntrflankalign import write_flank_outputs
        metadata.update(extraction_version=flank_context.model['extraction_version'],boundary_model_sha256=flank_context.model_hash)
        dump(destination/'boundary_model.json',flank_context.model)
        write_flank_outputs(destination,rows,records)
    policy=json.loads(args.merge_rules.read_text(encoding='utf-8')) if args.merge_rules else None
    metadata['parameters']={k:str(v) for k,v in vars(args).items()}
    metadata['source_hashes']={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                              for name in ['vntrpathgraph.py','vntrpathgraph_cli.py','vntrreads.py','vntrgraph.py']}
    if flank_context:
        for name in ['vntrflankalign.py']:
            metadata['source_hashes'][name]=sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
    if args.graph_mode=='structure':
        from vntrpathstructure_export import export_structure
        for name in ['vntrpathstructure.py','vntrpathstructure_export.py','vntrstructure.py','vntrpartition.py','vntrmotifs.py']:
            metadata['source_hashes'][name]=sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
        summary=export_structure(destination,locus,records,excluded,metadata,args.locus_model,args.structure_policy,
            args.loss_budget if structure_policy.get('schema_version')==2 else args.max_loss_bp,args.preview_paths,args.dot,rows if args.structures_json else None,
            payload.get('audit') if args.structures_json else None)
    else:
        summary=export_graphs(destination,args.locus,records,excluded,metadata,args.preview_paths,policy,args.loss_budget,args.dot)
    print(json.dumps({k:summary[k] for k in ['locus','status','paths','molecules','graph_stats','substituted_bases','preview']},ensure_ascii=False),flush=True)


if __name__=='__main__':main()
