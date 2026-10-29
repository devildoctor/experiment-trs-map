#!/usr/bin/env python3
"""All-read VNTR graphs. Single anchors and CIGAR overlap are display evidence.

No read is removed by MAPQ, duplicate flags, support thresholds or failed anchors.
One row per BAM/namespace/QNAME; supplementary and secondary observations remain separate
fragments, never concatenated into an invented complete allele.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import asdict
from hashlib import sha256
from html import escape
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from vntrgraph import IndexedFastaReference, read_catalog, discover_bams, flank_hit, write_tsv, source_hashes, load_source_namespaces, locus_output_dir
from vntrcore import canonical_motif, discover_motifs, edit_distance, compact

from vntrmotifs import primitive_unit, definition, base_metrics
from vntrcore import motif_orientation

VERSION = '3.1-phase1-all-reads'
CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')
COLORS = ('#3977ad', '#e18435', '#3a9460', '#a34e93', '#389baa', '#93812b')


# 把参考区间投影到 BAM 存储序列；边界插入保留，软剪切不自动归入重复区。
def cigar_interval(pos, cigar, start, end, query_length):
    """Map observed interval to stored SEQ. Keep insertions at either boundary.

    Return query bounds, aligned reference bounds, and evidence flags. Never
    assume that a query interval from a soft clip belongs to the VNTR.
    """
    ops = [(int(n), op) for n, op in CIGAR.findall(cigar)]
    if ''.join(f'{n}{op}' for n, op in ops) != cigar or not ops:
        raise ValueError('invalid_cigar')
    if sum(n for n, op in ops if op in 'MIS=X') != query_length:
        raise ValueError('cigar_sequence_length_mismatch')
    r, q = pos, 0
    spans, reference_spans, flags = [], [], []
    for n, op in ops:
        if op in 'M=X':
            a, b = max(r, start), min(r+n, end)
            if a < b:
                spans.append((q+a-r, q+b-r))
                reference_spans.append((a, b))
            r += n
            q += n
        elif op == 'I':
            if start <= r <= end:
                spans.append((q, q+n))
            q += n
        elif op in 'DN':
            a, b = max(r, start), min(r+n, end)
            if a < b:
                reference_spans.append((a, b))
                flags.append('reference_skip' if op == 'N' else 'deletion_in_interval')
            r += n
        elif op == 'S':
            flags.append('soft_clipped_sequence_not_assigned_to_repeat')
            q += n
        elif op == 'H':
            flags.append('hard_clipped_bases_unavailable')
    overlap = pos < end and r > start
    if not spans:
        return None, (max(pos, start), min(r, end)), sorted(set(flags)), overlap
    return (min(a for a, b in spans), max(b for a, b in spans)), (min(a for a,b in reference_spans), max(b for a,b in reference_spans)) if reference_spans else (start, start), sorted(set(flags)), overlap


# 这里只做局部显示定位；锚点命中不等于全基因组唯一，也不是分型证据。
def near_anchor(flank, sequence, expected, identity):
    """A local anchor helps delimit sequence; it is not a uniqueness gate."""
    if not flank:
        return None
    exact = sequence.find(flank)
    if exact >= 0 and sequence.find(flank, exact+1) < 0:
        return exact, exact+len(flank), 1.0
    # Restrict expensive approximate matching to a CIGAR-guided neighborhood.
    a = max(0, expected-120)
    b = min(len(sequence), expected+len(flank)+120)
    hit = flank_hit(flank, sequence[a:b], expected-a, identity)
    return (hit[0]+a, hit[1]+a, hit[2]) if hit else None


# 逐条比对记录提取实际片段；锚点不可用时回退到 CIGAR，并保留限制状态。
def extract_record(fields, locus, left, right):
    name, flag, chrom, pos, mapq, cigar = fields[:6]
    flag, pos, mapq = int(flag), int(pos)-1, int(mapq)
    sequence = fields[9].upper()
    meta = dict(flag=flag, chrom=chrom, pos=pos, mapq=mapq, cigar=cigar,
                strand='-' if flag & 16 else '+',
                alignment_type='secondary' if flag & 256 else 'supplementary' if flag & 2048 else 'primary')
    meta['flags'] = ([x for bit, x in ((1024,'duplicate_flag'),(512,'qc_fail'),(256,'secondary_alignment'),(2048,'supplementary_alignment')) if flag & bit]
                     + (['low_mapq'] if mapq < 20 else []))
    if flag & 4 or chrom != locus.chrom:
        return dict(meta, status='unmapped' if flag & 4 else 'other_contig')
    if sequence == '*':
        return dict(meta, status='sequence_unavailable')
    try:
        query, ref, flags, overlap = cigar_interval(pos, cigar, locus.start, locus.end, len(sequence))
    except ValueError as exc:
        return dict(meta, status=str(exc))
    meta['flags'] += flags
    if not overlap:
        return dict(meta, status='flank_only')
    if query is None:
        return dict(meta, status='no_observed_repeat_bases', reference_interval=ref)
    q0, q1 = query
    left_hit = near_anchor(left, sequence, q0-len(left), .72)
    right_hit = near_anchor(right, sequence, q1, .72)
    # A match must agree with the locus-projected neighborhood. Reject implausible
    # coordinates as anchors, not as reasons to discard the record.
    if left_hit and abs(left_hit[1]-q0) > 120:
        left_hit = None
    if right_hit and abs(right_hit[0]-q1) > 120:
        right_hit = None
    a = left_hit[1] if left_hit else q0
    b = right_hit[0] if right_hit else q1
    if a >= b:
        a, b = q0, q1
        left_hit = right_hit = None
        meta['flags'].append('anchor_order_conflict_cigar_fallback')
    coverage = ('spans_candidate_interval' if ref[0] <= locus.start and ref[1] >= locus.end
                else 'left_partial' if ref[0] <= locus.start else 'right_partial' if ref[1] >= locus.end
                else 'internal_partial')
    if 'reference_skip' in meta['flags']:
        coverage = 'split_reference_coverage'
    anchor = 'both' if left_hit and right_hit else 'left_only' if left_hit else 'right_only' if right_hit else 'cigar_only'
    return dict(meta, status='displayed', sequence=sequence[a:b], qstart=a, qend=b,
                reference_interval=ref, anchor=anchor, coverage=coverage,
                anchor_scores=[left_hit[2] if left_hit else None, right_hit[2] if right_hit else None],
                boundary_status=locus.boundary_status)


# 缺少 namespace 的旧结果按空字符串兼容；BAM 和平台不是分子身份。
def row_molecule_key(row):
    return compact([row['sample'], row.get('namespace', ''), row['read_id']])


# 读取整个 BAM，按 namespace/QNAME 建行；补充/次级比对独立保留，不拼接。
def read_bam(bam, sample, platform, locus, reference, samtools, anchor_length, records_file, namespace=None):
    left = reference.fetch(locus.chrom, max(0,locus.start-anchor_length), locus.start)
    right = reference.fetch(locus.chrom, locus.end, min(reference.contig_length(locus.chrom),locus.end+anchor_length))
    grouped, audit, read_groups = {}, Counter(), {}
    with tempfile.TemporaryFile(mode='w+t', encoding='utf-8') as errors:
        process = subprocess.Popen([samtools,'view','-h',str(bam)], stdout=subprocess.PIPE, stderr=errors, text=True)
        header_ok = False
        try:
            for record_index, line in enumerate(process.stdout, 1):
                if line.startswith('@'):
                    if line.startswith('@RG\t'):
                        tags = dict(v.split(':',1) for v in line.strip().split('\t')[1:])
                        read_groups[tags['ID']] = tags
                    if line.startswith('@SQ\t'):
                        tags = dict(v.split(':',1) for v in line.strip().split('\t')[1:])
                        if tags.get('SN') == locus.chrom:
                            if int(tags['LN']) != reference.contig_length(locus.chrom):
                                raise ValueError(f'BAM/reference length mismatch: {bam}')
                            header_ok = True
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 11:
                    raise ValueError(f'Malformed SAM record in {bam}: {record_index}')
                audit['alignment_records'] += 1
                tags = {p[0]:p[2] for value in f[11:] if len(p:=value.split(':',2)) == 3}
                read_group = tags.get('RG', '')
                library = read_groups.get(read_group, {}).get('LB', '')
                # 显式来源覆盖优先，空字符串也有效；否则沿用主流程的 RG/LB。
                ns = namespace if namespace is not None else library
                group_key = (ns, f[0])
                if group_key not in grouped:
                    key = compact([sample, platform, bam.name, ns, f[0]])
                    grouped[group_key] = dict(id='R'+sha256(key.encode()).hexdigest()[:16], read_id=f[0], sample=sample,
                                             namespace=ns, molecule_key=compact([sample, ns, f[0]]),
                                             platform=platform, bam=bam.name, fragments=[], records=[], record_count=0)
                row = grouped[group_key]
                observation = extract_record(f,locus,left,right)
                observation.update(namespace=ns, molecule_key=row['molecule_key'], read_group=read_group, library=library)
                observation['record_index'] = record_index
                record_id = row['id']+f'_rec{row["record_count"]+1}'
                observation['record_id'] = record_id
                # Complete stored SEQ remains available for clipping/interval audits.
                records_file.write(compact(dict(observation, bam=str(bam), read_id=f[0], row_id=row['id'],
                                                query_sequence=f[9], query_quality=f[10]))+'\n')
                row['record_count'] += 1
                row['records'].append({k:v for k,v in observation.items() if k != 'sequence'})
                audit[observation['status']] += 1
                if observation['status'] == 'displayed':
                    # Keep all alignment alternatives without stitching their bases.
                    observation['id'] = row['id']+f'_F{len(row["fragments"])+1}'
                    row['fragments'].append(observation)
        finally:
            process.stdout.close()
            code = process.wait()
        if code:
            errors.seek(0)
            raise RuntimeError(f'samtools failed for {bam}: {errors.read()}')
        if not header_ok:
            raise ValueError(f'Target contig missing in {bam}')
    rows = list(grouped.values())
    for row in rows:
        row['status'] = 'displayed' if row['fragments'] else ','.join(sorted({r['status'] for r in row['records']}))
    audit['read_rows'] = len(rows)
    audit['rows_with_repeat_sequence'] = sum(bool(r['fragments']) for r in rows)
    return rows, dict(audit)


# 依次收集各 BAM 的全部行和数量审计；跨 BAM 的重复行在此不删除。
def collect_reads(bams, locus, reference, samtools, anchor_length, records_file, namespaces=None):
    rows, audits = [], []
    namespaces = namespaces or {}
    for sample, platform, bam in bams:
        incoming, audit = read_bam(bam, sample, platform, locus, reference, samtools, anchor_length, records_file,
                                   namespaces.get(str(bam.resolve())))
        rows.extend(incoming)
        audits.append(dict(bam=str(bam), sample=sample, platform=platform, **audit))
        print(f'{locus.name} {sample}/{platform}: {len(incoming)} read rows, {audit["rows_with_repeat_sequence"]} with locus sequence', flush=True)
    return rows, audits


# 串行与进程池都保持输入顺序，避免注释结果配到另一条片段。
def annotation_results(annotator, tasks, workers):
    """Keep input order and the pool lifetime identical for both read views."""
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            yield from pool.map(annotator, tasks, chunksize=1)
    else:
        yield from map(annotator, tasks)


# 合并可用种子与目录模板；只有缺少种子时才尝试从参考发现模板。
def motif_templates(locus, reference_sequence, seed_root):
    records = []
    if seed_root:
        path = seed_root/locus.name/'motifs.tsv'
        if path.exists():
            with path.open(encoding='utf-8',newline='') as handle:
                for r in csv.DictReader(handle,delimiter='\t'):
                    if r.get('validated') == 'True' and r.get('level') == 'base':
                        records.append(dict(motif=r['motif'], source=str(path)))
    for motif in locus.motifs:
        records.append(dict(motif=motif, source='catalog'))
    if not records:
        records.extend(dict(motif=m,source='reference_discovery') for m in discover_motifs([reference_sequence]))
    unique = {}
    for r in list(records):
        r.update(definition(r['motif']))
        if r['primitive_sequence'] != r['motif']:
            records.append(dict(motif=r['primitive_sequence'],source='derived_primitive',**definition(r['primitive_sequence'])))
    for r in records:
        if r['motif'] not in unique:
            unique[r['motif']] = r
    return list(unique.values())


# 用精确短种子限定候选位置，再做近似单元分解；未解释碱基仍完整保存。
def annotate(task):
    """Seed-gated DP over true read slices; unknown intervals retain every base.

    This is exploratory read annotation, not a genotype or calibrated score.
    Runtime is bounded by candidate starts identified by exact short seeds;
    missing a noisy seed only yields UNKNOWN, never omission of a read/base.
    """
    sequence, templates, platform = task
    templates = tuple(dict.fromkeys(primitive_unit(m) for m in templates))
    n = len(sequence)
    candidates = defaultdict(set)
    for motif in templates:
        length = len(motif)
        allowance = 0 if length == 1 else max(1,min(10,math.ceil(length*(.12 if platform == 'HiFi' else .22))))
        k = min(length,5, max(2,length//4))
        for offset in sorted({0,max(0,length//2-k//2),length-k}):
            seed = motif[offset:offset+k]
            hit = sequence.find(seed)
            while hit >= 0:
                for delta in range(-allowance,allowance+1):
                    start = hit-offset+delta
                    if 0 <= start < n:
                        candidates[start].add(motif)
                hit = sequence.find(seed,hit+1)
    costs = [float('inf')]*(n+1)
    back = [None]*(n+1)
    ties = [False]*(n+1)
    costs[0] = 0
    def update(end, value, step):
        if value < costs[end]-1e-8:
            costs[end],back[end],ties[end] = value,step,ties[step[0]]
        elif abs(value-costs[end]) < 1e-8 and back[end] != step:
            ties[end] = True
    for i in range(n):
        weak_terminal = any(len(m)>1 and ((i==0 and sequence[i]==m[-1]) or (i==n-1 and sequence[i]==m[0])) for m in templates)
        update(i+1,costs[i]+(.8 if weak_terminal else 2),(i,'','unknown',0))
        for motif in sorted(candidates.get(i,())):
            length = len(motif)
            allowance = 0 if length == 1 else max(1,min(10,math.ceil(length*(.12 if platform == 'HiFi' else .22))))
            for width in range(max(1 if length==1 else 2,length-allowance),min(n-i,length+allowance)+1):
                exact = sequence[i:i+width]
                distance = edit_distance(motif,exact)
                terminal_partial = width < length and ((i+width == n and motif.startswith(exact)) or (i == 0 and motif.endswith(exact)))
                if distance <= allowance and not terminal_partial:
                    update(i+width,costs[i]+.5+distance*1.5,(i,motif,'unit',distance))
        # Exact terminal partials are not called complete copies.
        if i == 0:
            for motif in templates:
                for width in range(2,min(len(motif),n+1)):
                    if sequence[:width] == motif[-width:]:
                        update(width,.8,(0,motif,'partial',0))
        for motif in templates:
            if 2 <= n-i < len(motif) and sequence[i:] == motif[:n-i]:
                update(n,costs[i]+.8,(i,motif,'partial',0))
    reverse, end = [], n
    while end:
        start,motif,status,distance = back[end]
        reverse.append(dict(start=start,end=end,exact_bases=sequence[start:end],status=status,
                            template=motif,family=canonical_motif(motif) if motif else 'UNKNOWN',
                            edit_distance=distance,annotation_ambiguity=ties[end]))
        end = start
    tokens = []
    for t in reversed(reverse):
        if t['status'] == 'unknown' and tokens and tokens[-1]['status'] == 'unknown':
            tokens[-1]['end'] = t['end']
            tokens[-1]['exact_bases'] += t['exact_bases']
            tokens[-1]['annotation_ambiguity'] |= t['annotation_ambiguity']
        else:
            tokens.append(t)
    assert ''.join(t['exact_bases'] for t in tokens) == sequence
    for t in tokens:
        t['phase'],t['orientation']=motif_orientation(t['template'],t['family']) if t['template'] else (0,'+')
        t['phase_policy']='fixed_template_phase'
    return tokens


# 只写入声明的字段，复杂值的序列化沿用主程序的 TSV 写出规则。
def write_table(path, rows, fields):
    write_tsv(path, fields, ({k:v for k,v in row.items() if k in fields} for row in rows))


# 逐片段绘图；节点坐标来自实际切片，不把不同片段连成完整等位序列。
def svg_read(row, units, width=1560):
    """Nodes share IDs with GFA; a read can have several unjoined fragments."""
    chunks = []
    y = 28
    def text(x, yy, label, size=13, color='#243449'):
        chunks.append(f'<text x="{x}" y="{yy}" font-size="{size}" fill="{color}">{escape(str(label))}</text>')
    text(20,y,f'{row["sample"]} / {row["platform"]} / {row["read_id"]}',15)
    y += 23
    text(20,y,f'{row["id"]} | {row["bam"]} | namespace={row.get("namespace", "")} | records={row["record_count"]}')
    y += 25
    if not row['fragments']:
        text(20,y,'No extractable locus sequence: '+row['status'],13,'#955422')
        y += 30
    for fragment in row['fragments']:
        text(20,y,f'{fragment["id"]} | {fragment["alignment_type"]} | {fragment["coverage"]} | anchor={fragment["anchor"]} | MAPQ={fragment["mapq"]} | {len(fragment["sequence"])} observed bp')
        y += 22
        text(20,y,f'Stored-read coordinates [{fragment["qstart"]},{fragment["qend"]}); boundary={fragment["boundary_status"]}; observed units only')
        y += 18
        if 'str_annotation' in fragment:
            import vntrstr
            diagram, height = vntrstr.svg_blocks(fragment['str_annotation'],y=y,width=width-40)
            chunks.append(diagram)
            y += height
            a = fragment['str_annotation']
            text(20,y,f'Exact STR fragment; literal bases={a["literal_bases"]}; model match={a["model_match"]}; ambiguous={a["annotation_ambiguity"]}',12)
            y += 25
            for start in range(0,len(', '.join(fragment['flags'])),170):
                text(20,y,', '.join(fragment['flags'])[start:start+170],11,'#955422')
                y += 17
            continue
        # Bound visual width; never omit nodes or truncate the path.
        columns = 8
        for i, token in enumerate(fragment['tokens']):
            col, level = i%columns,i//columns
            x, top = 20+col*190,y+level*100
            color = '#707780' if token['status']=='unknown' else '#9b7537' if token['status']=='partial' else COLORS[int(sha256(token['family'].encode()).hexdigest()[:4],16)%len(COLORS)]
            node = token['node_id']
            title = f'{node}\n{token["unit_id"]}\nread [{token["qstart"]},{token["qend"]})\n{token["status"]}\n{token["exact_bases"]}'
            chunks.append(f'<g data-node="{node}"><title>{escape(title)}</title><rect x="{x}" y="{top}" width="176" height="80" rx="5" fill="{color}"/></g>')
            text(x+7,top+18,node,11,'white')
            text(x+7,top+35,f'{token["unit_id"]} | {len(token["exact_bases"])} bp | {token["status"]}',11,'white')
            text(x+7,top+52,f'[{token["qstart"]},{token["qend"]})',11,'white')
            bases = token['exact_bases']
            text(x+7,top+69,bases if len(bases)<=20 else bases[:8]+'…'+bases[-8:],11,'white')
            if col < columns-1 and i+1 < len(fragment['tokens']):
                text(x+178,top+44,'→',13)
        y += math.ceil(len(fragment['tokens'])/columns)*100+8
        if fragment['flags']:
            # Wrap all flags rather than silently slicing a warning string.
            words = ', '.join(fragment['flags'])
            for start in range(0,len(words),170):
                text(20,y,words[start:start+170],11,'#955422')
                y += 17
        y += 15
    return f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{y+10}" viewBox="0 0 {width} {y+10}" font-family="monospace"><rect width="100%" height="100%" fill="white"/>'+''.join(chunks)+'</svg>'


HTML = '''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>__LOCUS__ 全部读段节点图</title>
<style>body{font:15px system-ui;margin:24px;color:#253449}header{position:sticky;top:0;background:white;padding:12px 0;z-index:2;border-bottom:1px solid #ddd}input,select,button{padding:7px;margin:4px}article{border:1px solid #ccd5dd;margin:18px 0;padding:12px;border-radius:8px}svg{max-width:100%;height:auto}summary{cursor:pointer;padding:8px}pre{white-space:pre-wrap;word-break:break-all;background:#f2f5f8;padding:12px}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ddd;padding:8px;text-align:left}.seq{font:13px monospace;word-break:break-all}small{color:#627184}.warn{color:#955422}</style>
<header><h2>__LOCUS__：全部 BAM 读段及节点组成</h2><p>每行是一组 BAM/namespace/read ID，显示行数不等于独立分子数。单侧锚点、CIGAR 定位和部分覆盖均保留；不抽样，不按支持数或 MAPQ 删除。补充/次级比对分别展示，不能相加为完整等位。</p>
<input id="search" placeholder="搜索 read ID / BAM / namespace"><select id="sample"><option value="">所有样本</option></select><select id="platform"><option value="">所有平台</option></select><select id="coverage"><option value="">所有覆盖状态</option><option>displayed</option><option>left_only</option><option>right_only</option><option>cigar_only</option><option>both</option><option>flank_only</option></select>
<button id="prev">上一页</button><span id="page"></span><button id="next">下一页</button><input id="goto" type="number" min="1" style="width:65px" aria-label="页码"><button id="jump">跳页</button><div id="count"></div></header>
<details><summary>处理逻辑、总量核对和输出文件</summary><pre id="audit"></pre><p><a href="__LOCUS__.reads.gfa">全部读段 GFA</a> · <a href="reads.tsv">读段清单</a> · <a href="read_nodes.tsv">节点组成表</a> · <a href="records.jsonl">原始比对与完整存储序列</a> · <a href="svg/index.html">全部静态 SVG 分页</a></p></details>
<main id="rows"></main><script id="data" type="application/json">__DATA__</script>
<script>
const data=JSON.parse(document.getElementById('data').textContent);let page=0;const size=20;
if(data.rows.some(r=>r.fragments.some(f=>f.str_annotation))){const p=document.createElement('p');p.innerHTML='<a href="str_paths.tsv">STR 表达表</a> · <a href="str_paths.json">STR 完整切分</a> · <a href="str-svg/index.html">全部 STR 图</a>';document.querySelector('header').append(p);}
const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
for(const field of ['sample','platform'])for(const value of [...new Set(data.rows.map(r=>r[field]))].sort()){let op=document.createElement('option');op.value=value;op.textContent=value;document.getElementById(field).append(op)}
document.getElementById('audit').textContent=JSON.stringify(data.audit,null,2)+'\\n流程：读取全部记录 → 按 BAM/namespace/read ID 分组 → 单侧锚点或 CIGAR 截取 → 逐 read 分解真实碱基 → 保存全部节点 → 重构核验。UNKNOWN 表示未完成 motif 注释，碱基仍保留。坐标为 BAM 存储 SEQ 的 0-based 半开坐标。';
function render(){const q=document.getElementById('search').value.toLowerCase(),s=document.getElementById('sample').value,p=document.getElementById('platform').value,c=document.getElementById('coverage').value;
const selected=data.rows.filter(r=>(!q||(r.read_id+' '+r.bam+' '+(r.namespace||'')).toLowerCase().includes(q))&&(!s||r.sample===s)&&(!p||r.platform===p)&&(!c||r.status===c||r.fragments.some(f=>f.anchor===c)));
page=Math.max(0,Math.min(page,Math.ceil(selected.length/size)-1));document.getElementById('count').textContent=`共 ${data.rows.length} 行 read；筛选后 ${selected.length} 行。每页 ${size} 行，所有页均来自完整数据。`;document.getElementById('page').textContent=`${page+1} / ${Math.max(1,Math.ceil(selected.length/size))}`;
document.getElementById('rows').innerHTML=selected.slice(page*size,(page+1)*size).map(r=>`<article>${r.svg}<details><summary>展开该 read 的全部节点与完整碱基</summary>${r.fragments.map(f=>`<h4>${esc(f.id)} · ${esc(f.coverage)} · ${esc(f.anchor)}</h4>${f.str_annotation?`<p><b>精确 STR 表达：</b><span class="seq">${esc(f.str_annotation.expression)}</span></p><p class="warn">仅描述观察到的片段；${f.str_annotation.annotation_ambiguity?'存在同分切分歧义':'未发现同分切分'}。模型 ${esc(f.str_annotation.structure)}；字面碱基 ${f.str_annotation.literal_bases} bp。</p><details><summary>STR 切分与备选</summary><pre>${esc(JSON.stringify(f.str_annotation,null,2))}</pre></details>`:''}<p>有序路径：${f.tokens.map(t=>esc(t.node_id)).join(' → ')}</p><table><tr><th>节点 / 单元</th><th>read 坐标</th><th>状态 / 模板编辑距离</th><th>实际碱基</th></tr>${f.tokens.map(t=>`<tr><td>${esc(t.node_id)} / ${esc(t.unit_id)}</td><td>[${t.qstart},${t.qend})</td><td>${esc(t.status)} / ${t.edit_distance}${t.annotation_ambiguity?' / 注释歧义':''}</td><td class="seq">${esc(t.exact_bases)}</td></tr>`).join('')}</table><details><summary>这段实际序列 ${f.sequence.length} bp</summary><pre>${esc(f.sequence)}</pre></details>`).join('')}<details><summary>全部比对记录状态</summary><pre>${esc(JSON.stringify(r.records,null,2))}</pre></details></details></article>`).join('');}
for(const id of ['search','sample','platform','coverage'])document.getElementById(id).addEventListener('input',()=>{page=0;render()});document.getElementById('prev').onclick=()=>{page--;render()};document.getElementById('next').onclick=()=>{page++;render()};document.getElementById('jump').onclick=()=>{page=(Number(document.getElementById('goto').value)||1)-1;render()};render();
</script></html>'''


# 输出独立片段路径、逐节点真实碱基和全部读段分页，并核对路径重建。
def write_outputs(destination,locus,rows,audit):
    if locus.str_structure:
        import vntrstr
        str_records = vntrstr.annotate_rows(rows,locus.str_structure)
        vntrstr.write_outputs(destination,locus.name,locus.str_structure,str_records)
    units = {}
    node_rows = []
    for row in rows:
        for fragment in row['fragments']:
            for index,t in enumerate(fragment['tokens'],1):
                key = (t['exact_bases'],t['status'],t['family'])
                if key not in units:
                    units[key] = f'U{len(units)+1}'
                t.update(node_id=f'{fragment["id"]}_N{index}',unit_id=units[key],
                         qstart=fragment['qstart']+t['start'],qend=fragment['qstart']+t['end'])
                node_rows.append(dict(t,row_id=row['id'],fragment_id=fragment['id'],read_id=row['read_id'],
                                      sample=row['sample'],namespace=row.get('namespace',''),molecule_key=row_molecule_key(row),order=index))
            assert ''.join(t['exact_bases'] for t in fragment['tokens']) == fragment['sequence']
            assert all(t['exact_bases']==fragment['sequence'][t['start']:t['end']] for t in fragment['tokens'])
    with (destination/f'{locus.name}.reads.gfa').open('w',encoding='utf-8') as handle:
        handle.write('H\tVN:Z:1.0\tPG:Z:vntrgraph-3.0-all-reads\n')
        for row in rows:
            for f in row['fragments']:
                for t in f['tokens']:
                    handle.write(f'S\t{t["node_id"]}\t{t["exact_bases"]}\tRC:i:1\tUT:Z:{t["unit_id"]}\tST:Z:{t["status"]}\n')
                for a,b in zip(f['tokens'],f['tokens'][1:]):
                    handle.write(f'L\t{a["node_id"]}\t+\t{b["node_id"]}\t+\t0M\n')
                ids=','.join(t['node_id']+'+' for t in f['tokens'])
                if ids:
                    handle.write(f'P\t{f["id"]}\t{ids}\t*\tRN:Z:{row["read_id"]}\tCV:Z:{f["coverage"]}\tAN:Z:{f["anchor"]}\n')
    write_table(destination/'read_nodes.tsv',node_rows,['row_id','read_id','fragment_id','order','node_id','unit_id','start','end','qstart','qend','status','family','template','edit_distance','annotation_ambiguity','exact_bases','sample','namespace','molecule_key'])
    write_table(destination/'reads.tsv',[dict(r,namespace=r.get('namespace',''),molecule_key=row_molecule_key(r),fragment_count=len(r['fragments']),anchors=sorted({f['anchor'] for f in r['fragments']})) for r in rows],['id','read_id','sample','platform','bam','record_count','fragment_count','status','anchors','namespace','molecule_key'])
    svg_dir=destination/'svg'
    svg_dir.mkdir()
    links=[]
    for start in range(0,len(rows),20):
        group=rows[start:start+20]
        offset=60
        body=[]
        for row in group:
            row['svg']=svg_read(row,units)
            height=int(re.search(r'height="(\d+)"',row['svg'])[1])
            body.append(f'<g transform="translate(0,{offset})">{row["svg"]}</g>')
            offset+=height+15
        number=start//20+1
        name=f'{locus.name}.reads.{number:04d}.svg'
        (svg_dir/name).write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="1560" height="{offset}" viewBox="0 0 1560 {offset}"><rect width="100%" height="100%" fill="white"/><text x="20" y="30" font-size="20">{locus.name} all reads | page {number} | rows {start+1}-{start+len(group)} / {len(rows)}</text>'+''.join(body)+'</svg>',encoding='utf-8')
        links.append(f'<li><a href="{name}">第 {number} 页，read {start+1}–{start+len(group)}</a></li>')
    (svg_dir/'index.html').write_text('<!doctype html><meta charset="utf-8"><h2>全部 read SVG 分页，无抽样</h2><ol>'+''.join(links)+'</ol>',encoding='utf-8')
    audit['svg_pages']=len(links)
    audit['graph_fragments']=sum(len(r['fragments']) for r in rows)
    audit['graph_nodes']=len(node_rows)
    audit['sequence_reconstruction_checks_passed']=audit['graph_fragments']
    (destination/'status.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    # Escape script delimiters in arbitrary BAM names/read IDs.
    payload=compact(dict(rows=rows,audit=audit)).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
    (destination/f'{locus.name}.reads.html').write_text(HTML.replace('__LOCUS__',locus.name).replace('__DATA__',payload),encoding='utf-8')
    for row in rows:
        row.pop('svg',None)
    (destination/f'{locus.name}.reads.json').write_text(compact(dict(rows=rows,audit=audit)),encoding='utf-8')


# 单个位点入口：读取全部记录、准备模板、注释片段、输出并记录参数。
def run_locus(locus,args,reference):
    destination=locus_output_dir(args.out,locus.name)
    destination.mkdir(parents=True,exist_ok=False)
    bams=discover_bams(args.bam_dir,locus.name)
    if not bams:
        raise ValueError(f'No BAM files: {locus.name}')
    with (destination/'records.jsonl').open('w',encoding='utf-8') as record_file:
        rows,bam_audits=collect_reads(bams,locus,reference,args.samtools,args.anchor_length,record_file,args.source_namespaces)
    records=motif_templates(locus,reference.fetch(locus.chrom,locus.start,locus.end),args.motif_seed_dir)
    write_table(destination/'motifs.tsv',records,['motif','source','level','primitive_sequence','repeat_factor','phase_policy'])
    templates=tuple(r['motif'] for r in records)
    fragments=[(f,row['platform']) for row in rows for f in row['fragments']]
    tasks=[(f['sequence'],templates,platform) for f,platform in fragments]
    print(f'{locus.name}: annotating {len(tasks)} fragments in {len(rows)} reads; {len(templates)} templates; no coverage gate',flush=True)
    results=annotation_results(annotate,tasks,args.workers)
    for i,((fragment,_),tokens) in enumerate(zip(fragments,results),1):
        fragment['tokens']=tokens
        fragment['base_metrics']=base_metrics(tokens)
        if args.workers>1 and i%25==0:
            print(f'{locus.name}: annotated {i}/{len(tasks)}',flush=True)
    audit=dict(version=VERSION,locus=asdict(locus),status='complete',bam_count=len(bams),read_rows=len(rows),
               rows_with_repeat_sequence=sum(bool(r['fragments']) for r in rows),
               rows_without_extractable_repeat=sum(not r['fragments'] for r in rows),
               alignment_records=sum(r['record_count'] for r in rows),
               anchors=dict(Counter(f['anchor'] for f,_ in fragments)),
               coverage=dict(Counter(f['coverage'] for f,_ in fragments)),
               bam_audits=bam_audits,read_row_key='BAM + namespace + QNAME; display rows are not independent molecular support',
               independent_molecules=len({row_molecule_key(r) for r in rows}),
               molecule_key='sample + namespace + QNAME; --sources overrides RG/LB, otherwise empty when LB is absent',
               filtering='no MAPQ, flag, minimum-support or anchor-count read exclusion',
               annotation='seed-gated DP; unrecognized bases preserved; read units are not definitive allele copies')
    assert audit['read_rows']==sum(b['read_rows'] for b in bam_audits)
    write_outputs(destination,locus,rows,audit)
    (destination/'parameters.json').write_text(json.dumps(dict(version=VERSION,args={k:str(v) for k,v in vars(args).items()},
                 source_hashes=source_hashes(),
                 catalog_sha256=sha256(args.catalog.read_bytes()).hexdigest()),ensure_ascii=False,indent=2),encoding='utf-8')
    print(compact({k:v for k,v in audit.items() if k!='bam_audits'}),flush=True)


# 命令行校验位点和参数；输出目录必须是新的，避免覆盖已有结果。
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam-dir',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--catalog',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--locus',action='append',required=True)
    parser.add_argument('--samtools',default='samtools')
    parser.add_argument('--anchor-length',type=int,default=100)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--motif-seed-dir',type=Path)
    parser.add_argument('--sources',type=Path,help='Optional bam/namespace TSV; overrides RG/LB, including an empty namespace')
    args=parser.parse_args()
    try:args.source_namespaces=load_source_namespaces(args.sources)
    except ValueError as exc:parser.error(str(exc))
    loci=[l for l in read_catalog(args.catalog) if l.name in args.locus]
    if len(loci)!=len(set(args.locus)):
        parser.error('Unknown locus')
    if args.workers<1 or args.anchor_length<20:
        parser.error('workers >= 1 and anchor-length >= 20 required')
    for locus in loci:
        if locus_output_dir(args.out,locus.name).exists():
            parser.error('Use a new output directory; existing results are preserved')
    with IndexedFastaReference(args.reference) as reference:
        for locus in loci:
            run_locus(locus,args,reference)


if __name__=='__main__':
    main()
