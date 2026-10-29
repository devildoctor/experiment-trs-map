#!/usr/bin/env python3
"""Exact, oriented compact de Bruijn graphs; no allele inference or STR changes."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from html import escape
from pathlib import Path
import csv
import json
import math
import re
import subprocess

VERSION = 'vntr_compact_dbg_v1'
DNA = set('ACGT')


def digest(text):
    return sha256(text.encode('ascii')).hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')


def tsv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or (list(rows[0]) if rows else ['id'])
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False, sort_keys=True)
                             if isinstance(value, (list, dict, tuple)) else value
                             for key, value in row.items() if key in fields})


# 按原序列方向依次取 k-mer，重复出现的位置仍保留在路径中。
def words(sequence, k):
    return [sequence[i:i+k] for i in range(len(sequence)-k+1)]


# 原始图只连接输入序列中真实相邻的 k-mer；REF/CAND 参与拓扑但不增加分子票。
def raw_graph(records, k):
    """Edges must occur in a real input k+1-mer; overlap alone is insufficient."""
    nodes, edges, paths = set(), set(), {}
    ns, es = defaultdict(set), defaultdict(set)
    molecules = {}
    for record in sorted(records, key=lambda r:r['id']):
        sequence = record['sequence']
        if set(sequence)-DNA:
            raise ValueError('Split ambiguous bases before graph construction')
        if len(sequence)<k:
            continue
        path = words(sequence, k)
        paths[record['id']] = path
        nodes.update(path)
        links = set(zip(path,path[1:]))
        edges.update(links)
        # 用集合去重：同一分子多次经过节点或边，独立支持仍只记一次。
        if record['kind']=='observed':
            key=record['molecule_key']
            meta=molecules.setdefault(key, dict(sample=record['sample'],platforms=set()))
            if meta['sample']!=record['sample']:
                raise ValueError('Molecule key crosses samples')
            meta['platforms'].update(record['platforms'])
            for node in set(path):
                ns[node].add(key)
            for edge in links:
                es[edge].add(key)
    return dict(nodes=nodes,edges=edges,paths=paths,node_support=ns,edge_support=es,molecules=molecules)


# 支持门槛按单个样本计算；过滤断点必须拆成片段，不能把两端重新连起来。
def supported_records(records, raw, k, minimum):
    """Keep edges supported by >=minimum unique molecules in at least one sample.

    Reference is retained as explicitly zero-evidence context. Removed edges
    split paths; no false full-length paths or bridges across removed regions.
    """
    def enough(keys):
        counts=Counter(raw['molecules'][key]['sample'] for key in keys)
        return max(counts.values(), default=0)>=minimum
    keep_nodes={n for n in raw['nodes'] if enough(raw['node_support'][n])}
    keep_edges={e for e in raw['edges'] if enough(raw['edge_support'][e])}
    retained=[]
    for record in records:
        if record['kind']=='reference':
            retained.append(dict(record))
            continue
        path=raw['paths'].get(record['id'],[])
        intervals=[]
        start=None
        for i,node in enumerate(path):
            if node not in keep_nodes:
                if start is not None:
                    intervals.append((start,i))
                start=None
            elif start is None:
                start=i
            elif (path[i-1],node) not in keep_edges:
                intervals.append((start,i))
                start=i
        if start is not None:
            intervals.append((start,len(path)))
        for index,(a,b) in enumerate(intervals,1):
            full=a==0 and b==len(path)
            # k-mer 索引区间 [a,b) 对应碱基区间 [a,b+k-1)。
            seq=record['sequence'][a:b+k-1]
            item=dict(record, id=record['id'] if full else f"{record['id']}_F{index:03}",
                      parent_id=record['id'],parent_start=a,parent_end=b+k-1,
                      sequence=seq,sequence_sha256=digest(seq),
                      coverage=record['coverage'] if full else 'filtered_fragment',
                      filter_complete=full)
            retained.append(item)
    return retained, keep_nodes, keep_edges


# 将无分支链压缩为精确 unitig，同时保留所有输入路径的起止位置。
def compact_graph(records, raw, k, locus, view, minimum=3):
    incoming,outgoing=defaultdict(set),defaultdict(set)
    for a,b in raw['edges']:
        outgoing[a].add(b)
        incoming[b].add(a)
    starts={p[0] for p in raw['paths'].values()}
    ends={p[-1] for p in raw['paths'].values()}
    # Node-disjoint chains. Preserve every path start/end, and stop at branches.
    successor={a:next(iter(outgoing[a])) for a in raw['nodes']
               if len(outgoing[a])==1 and a not in ends}
    successor={a:b for a,b in successor.items() if len(incoming[b])==1 and b not in starts}
    predecessor={b:a for a,b in successor.items()}
    chains,seen=[],set()
    for seed in sorted(raw['nodes']-set(predecessor))+sorted(raw['nodes']):
        if seed in seen:
            continue
        chain=[]
        node=seed
        while node not in seen:
            chain.append(node)
            seen.add(node)
            if node not in successor:
                break
            node=successor[node]
        chains.append(chain)
    chains.sort(key=lambda chain:tuple(chain))
    nodes,owner,internal={}, {}, set()
    for index,chain in enumerate(chains,1):
        name=f'U{index:05}'
        sequence=chain[0]+''.join(word[-1] for word in chain[1:])
        for offset,word in enumerate(chain):
            owner[word]=(name,offset)
        internal.update(zip(chain,chain[1:]))
        nodes[name]=dict(sequence=sequence,length=len(sequence),kmers=chain,raw_kmers=len(chain))
    edges={(owner[a][0],owner[b][0]) for a,b in raw['edges'] if (a,b) not in internal}
    paths=[]
    for record in sorted(records,key=lambda r:(r['kind']!='reference',r['kind']=='observed',r['id'])):
        path=raw['paths'].get(record['id'])
        if not path:
            continue
        route=[]
        i=0
        while i<len(path):
            name,offset=owner[path[i]]
            if offset:
                raise AssertionError('Path begins inside compacted segment')
            chain=nodes[name]['kmers']
            if path[i:i+len(chain)]!=chain:
                raise AssertionError('Compaction crossed a path endpoint or branch')
            route.append(name)
            i+=len(chain)
        # 相邻 unitig 重叠 k-1 bp，重建时只追加后续节点的新碱基。
        result=nodes[route[0]]['sequence']+''.join(nodes[n]['sequence'][k-1:] for n in route[1:])
        assert result==record['sequence']
        paths.append(dict(record,walk=route,sequence_sha256=digest(result)))
    ns,es=defaultdict(set),defaultdict(set)
    nplatform,eplatform=defaultdict(lambda:defaultdict(set)),defaultdict(lambda:defaultdict(set))
    nt,et=defaultdict(set),defaultdict(set)
    for path in paths:
        for n in set(path['walk']):
            nt[n].add(path['kind'])
            if path['kind']=='observed':
                ns[n].add(path['molecule_key'])
                for platform in path['platforms']:
                    nplatform[n][platform].add(path['molecule_key'])
        for edge in set(zip(path['walk'],path['walk'][1:])):
            et[edge].add(path['kind'])
            if path['kind']=='observed':
                es[edge].add(path['molecule_key'])
                for platform in path['platforms']:
                    eplatform[edge][platform].add(path['molecule_key'])
    # Platform support uses actual paths across the object, not every platform
    # ever associated with that molecule elsewhere in the locus.
    # 平台小计可能重叠；分样本和总支持均来自独立分子集合。
    def evidence(keys, kind, obj):
        platforms=nplatform[obj] if kind=='node' else eplatform[obj]
        samples=dict(sorted(Counter(raw['molecules'][key]['sample'] for key in keys).items()))
        return dict(molecule_ids=sorted(keys),molecule_support=len(keys),
                    sample_support=samples,max_sample_support=max(samples.values(),default=0),
                    platform_support={p:len(v) for p,v in sorted(platforms.items())})
    for n,value in nodes.items():
        value.update(evidence(ns[n],'node',n),origins=sorted(nt[n]))
    links=[dict(source=a,target=b,overlap=k-1,origins=sorted(et[a,b]),**evidence(es[a,b],'edge',(a,b)))
           for a,b in sorted(edges)]
    return dict(version=VERSION,locus=locus,k=k,view=view,overlap=k-1,
                min_support=minimum,support_scope='max_unique_molecules_within_one_sample',
                encoding='reference_oriented_exact_string',nodes=nodes,edges=links,paths=paths)


# GFA 的 S 为精确序列，L 保存重叠，P 保存实际输入路径；不能枚举拓扑来宣称等位基因。
def write_gfa(path, graph):
    k=graph['k']
    lines=[f"H\tVN:Z:1.0\tgm:Z:{VERSION}\tkk:i:{k}\tvw:Z:{graph['view']}"]
    for name,node in graph['nodes'].items():
        lines.append(f"S\t{name}\t{node['sequence']}\tLN:i:{node['length']}\tms:i:{node['molecule_support']}\tor:Z:{','.join(node['origins'])}")
    for e in graph['edges']:
        lines.append(f"L\t{e['source']}\t+\t{e['target']}\t+\t{k-1}M\tms:i:{e['molecule_support']}\tor:Z:{','.join(e['origins'])}")
    for p in graph['paths']:
        overlaps=','.join([f'{k-1}M']*(len(p['walk'])-1)) or '*'
        lines.append(f"P\t{p['id']}\t{','.join(n+'+' for n in p['walk'])}\t{overlaps}\tpt:Z:{p['kind']}\tcv:Z:{p['coverage']}")
    Path(path).write_text('\n'.join(lines)+'\n',encoding='ascii')


def color(node):
    if node['molecule_support']==0:
        return '#e4e7ec' if 'reference' in node['origins'] else '#f7dfba'
    return '#cee9df' if node['platform_support'].get('HiFi',0)>=2 else '#dce7f6'


# Graphviz 仅负责布局，节点尺寸不代表 bp 长度，边粗细仅用于展示支持量。
def topology(path, graph):
    count=len(graph['nodes'])
    minimum=graph.get('min_support',3)
    engine='dot' if count<=450 else 'neato'
    # VPSC works with the server Graphviz build and keeps all graph objects.
    layout='overlap=false, splines=true' if engine=='dot' else 'overlap=vpsc, splines=line, maxiter=100'
    lines=['digraph DBG {',
           f'graph [bgcolor="white", {layout}, outputorder=edgesfirst, pad=0.4, start=42];',
           'node [shape=box, style="rounded,filled", fontname="monospace", fontsize=9];',
           'edge [arrowsize=0.55, fontname="sans", fontsize=8];',
           f'label="{graph["locus"]} k={graph["k"]} {graph["view"]} | S=exact unitig; L overlap={graph["overlap"]} bp; ms=unique molecules; edge color=max sample support >= {minimum}; edge width=total ms; node size not bp scale";',
           'labelloc=t;']
    if engine=='dot':
        lines.append('rankdir=LR;')
    for name,n in graph['nodes'].items():
        seq=n['sequence']
        short=seq if len(seq)<=24 else seq[:11]+'...'+seq[-10:]
        label=f'{name} | {len(seq)} bp\n{short}\nms={n["molecule_support"]}'
        lines.append(f'{name} [label={json.dumps(label)}, fillcolor="{color(n)}", tooltip={json.dumps(seq)}];')
    for e in graph['edges']:
        support=e['molecule_support']
        # 颜色使用同一样本最大独立支持和本次门槛，不能跨样本凑票。
        sample_max=max(e['sample_support'].values(),default=0)
        line_color='#b8c0ca' if sample_max<minimum else '#47566e'
        width=0.6+min(3,math.log2(1+support)*0.45)
        lines.append(f'{e["source"]} -> {e["target"]} [color="{line_color}", penwidth={width:.2f}, tooltip="ms={support}; max_sample_support={sample_max}; min_support={minimum}; overlap={graph["overlap"]}"];')
    lines.append('}')
    dotpath=Path(path).with_suffix('.dot')
    dotpath.write_text('\n'.join(lines)+'\n',encoding='ascii')
    subprocess.run([engine,'-Tsvg',str(dotpath),'-o',str(path)],check=True,timeout=240)
    graph['layout_engine']=engine
    graph['layout_options']=layout


# 展示明确保存的路径；重复访问同一图节点不等于重复单元拷贝数。
def paths_svg(path, graph, paths=None, annotations=None):
    paths=graph['paths'] if paths is None else paths
    annotations=annotations or {}
    body=[]
    y=68
    width=1300
    for p in paths:
        walk=p['walk']
        body.append(f'<g class="graph-path" data-path="{p["id"]}" data-walk="{",".join(walk)}">')
        title=f"{p['id']} | {p['kind']} | {p.get('sample','reference')} | {p['coverage']} | {len(p['sequence'])} bp"
        body.append(f'<text x="24" y="{y}" font-size="15" font-weight="bold">{escape(title)}</text>')
        y+=18
        tokens=annotations.get(p.get('parent_id',p['id']),[])
        # 过滤后的片段不能直接套用完整路径的单元坐标。
        if tokens and p.get('filter_complete',True):
            palette=['#dce7f6','#cce5df','#f4ddaf','#e8d4e8','#d5e4bb']
            base=y+8
            for i,t in enumerate(tokens):
                x=24+t['start']/len(p['sequence'])*1200
                w=(t['end']-t['start'])/len(p['sequence'])*1200
                fill='#dddddd' if t['status']=='unknown' else palette[int(digest(t['exact_bases'])[:6],16)%len(palette)]
                body.append(f'<rect x="{x:.2f}" y="{base}" width="{w:.2f}" height="16" fill="{fill}" stroke="white"><title>{escape(t["family"])} [{t["start"]},{t["end"]}) {t["status"]}: {t["exact_bases"]}</title></rect>')
            body.append(f'<text x="24" y="{base+31}" font-size="10">Motif annotation along sequence coordinates; colors identify exact slices, not graph nodes. Hover for bases and boundaries.</text>')
            y+=50
        rows=(len(walk)+7)//8
        for i,n in enumerate(walk):
            x=24+(i%8)*155
            yy=y+(i//8)*94
            node=graph['nodes'][n]
            seq=node['sequence']
            short=seq if len(seq)<=18 else seq[:8]+'...'+seq[-7:]
            added=len(seq) if i==0 else len(seq)-graph['overlap']
            body.append(f'<g class="graph-node" data-node="{n}" data-index="{i}" data-added="{added}"><title>{seq}</title>'
                        f'<rect x="{x}" y="{yy}" width="140" height="68" rx="5" fill="{color(node)}" stroke="#627087"/>'
                        f'<text x="{x+7}" y="{yy+17}" font-size="12">{n} | {len(seq)} bp</text>'
                        f'<text x="{x+7}" y="{yy+36}" font-size="10" font-family="monospace">{short}</text>'
                        f'<text x="{x+7}" y="{yy+55}" font-size="10">+{added} bp | ms={node["molecule_support"]}</text></g>')
            if i+1<len(walk):
                if i%8==7:
                    shape=f'<text x="{x+72}" y="{yy+84}" font-size="10">continue next row</text>'
                else:
                    shape=f'<path d="M{x+141},{yy+34} h12" stroke="#526175" marker-end="url(#a)"/>'
                body.append(f'<g class="graph-edge" data-source="{n}" data-target="{walk[i+1]}" data-index="{i}">{shape}</g>')
        y+=rows*94+30
        body.append('</g>')
    header=(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{y+12}" viewBox="0 0 {width} {y+12}" font-family="sans-serif">'
            '<rect width="100%" height="100%" fill="white"/><defs><marker id="a" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto"><path d="M0 0 L10 5 L0 10z" fill="#526175"/></marker></defs>'
            f'<text x="24" y="27" font-size="20">{graph["locus"]} | k={graph["k"]} | {graph["view"]} | explicit paths</text>'
            f'<text x="24" y="46" font-size="11">Each arrow overlaps {graph["overlap"]} bp. +bp is newly spelled sequence. Repeated visits are not motif copy counts.</text>')
    Path(path).write_text(header+''.join(body)+'</svg>\n',encoding='utf-8')


# 同步导出拓扑、逐路径图、GFA、JSON 和独立支持表。
def write_outputs(folder, graph, annotations):
    suffix='' if graph['view']=='all_evidence' else '.supported'
    stem=f"{graph['locus']}.k{graph['k']}{suffix}"
    write_gfa(folder/(stem+'.gfa'),graph)
    topology(folder/(stem+'.svg'),graph)
    paths_svg(folder/(stem+'.paths.svg'),graph,annotations=annotations)
    paths_svg(folder/(stem+'.summary.paths.svg'),graph,
              [p for p in graph['paths'] if p['kind']!='observed'],annotations)
    detail=folder/('paths-svg'+suffix)
    detail.mkdir()
    for p in graph['paths']:
        paths_svg(detail/f"{stem}.{p['id']}.svg",graph,[p],annotations)
    save_json(folder/(stem+'.json'),graph)
    tsv(folder/(stem+'.paths.tsv'),graph['paths'],['id','kind','sample','coverage','molecule_key','platforms','parent_id','parent_start','parent_end','filter_complete','sequence_sha256','walk','source_records'])
    tsv(folder/(stem+'.support.tsv'),
        [dict(object='node',id=n,**v) for n,v in graph['nodes'].items()]+
        [dict(object='edge',id=e['source']+'>'+e['target'],**e) for e in graph['edges']],
        ['object','id','molecule_support','sample_support','platform_support','molecule_ids','origins','max_sample_support'])


# 这些指标比较图压缩和证据保留程度，不是分型准确率。
def metrics(raw, graph, kept_nodes, kept_edges, records, k):
    incoming,outgoing=defaultdict(set),defaultdict(set)
    for a,b in kept_edges:
        incoming[b].add(a)
        outgoing[a].add(b)
    candidates=[r for r in records if r['kind']=='candidate']
    totals=covered=repeated=windows=0
    for p in candidates:
        route=words(p['sequence'],k)
        counts=Counter(route)
        repeated+=sum(counts[n]>1 for n in route)
        windows+=len(route)
        edges=list(zip(route,route[1:]))
        totals+=len(edges)
        covered+=sum(e in kept_edges for e in edges)
    result=dict(k=k,raw_kmer_nodes=len(raw['nodes']),raw_kmer_edges=len(raw['edges']),
                compact_nodes=len(graph['nodes']),compact_edges=len(graph['edges']),
                raw_single_molecule_edge_fraction=sum(len(raw['edge_support'][e])==1 for e in raw['edges'])/max(1,len(raw['edges'])),
                supported_kmer_nodes=len(kept_nodes),supported_kmer_edges=len(kept_edges),
                supported_branch_nodes=sum(len(incoming[n])>1 or len(outgoing[n])>1 for n in kept_nodes),
                supported_merge_and_split_nodes=sum(len(incoming[n])>1 and len(outgoing[n])>1 for n in kept_nodes),
                candidate_edge_retention=covered/max(1,totals),
                candidate_repeated_kmer_fraction=repeated/max(1,windows))
    for platform in ['HiFi','ONT','ONTUL']:
        pp=[p for p in records if p['kind']=='observed' and platform in p['platforms'] and len(p['sequence'])>=k]
        full=sum(all(e in kept_edges for e in zip(words(p['sequence'],k),words(p['sequence'],k)[1:])) and all(n in kept_nodes for n in words(p['sequence'],k)) for p in pp)
        result[platform+'_observed_paths']=len(pp)
        result[platform+'_fully_retained']=full
    return result


# 每个 k 独立生成全证据图和支持图，二者都保留各自的明确路径。
def build_all(records, annotations, locus, ks, out, minimum=3):
    comparison=[]
    for k in ks:
        folder=out/f'k{k}'
        folder.mkdir()
        raw=raw_graph(records,k)
        kept,kn,ke=supported_records(records,raw,k,minimum)
        graph=compact_graph(records,raw,k,locus,'all_evidence',minimum)
        print(f'k={k}: {len(graph["nodes"])} compact nodes, {len(graph["edges"])} links, {len(graph["paths"])} paths',flush=True)
        write_outputs(folder,graph,annotations)
        filtered_raw=raw_graph(kept,k)
        filtered=compact_graph(kept,filtered_raw,k,locus,'supported',minimum)
        write_outputs(folder,filtered,annotations)
        row=metrics(raw,graph,kn,ke,records,k)
        row.update(min_support=minimum,supported_compact_nodes=len(filtered['nodes']),supported_compact_edges=len(filtered['edges']),
                   exported_paths=len(graph['paths']),supported_exported_fragments=len(filtered['paths']),
                   supported_candidate_full_paths=sum(p['kind']=='candidate' and p.get('filter_complete',True) for p in filtered['paths']))
        comparison.append(row)
        # 短于 k 的输入无法形成该图节点，单独列出以避免静默遗漏。
        tsv(folder/'excluded_short_paths.tsv',[dict(id=r['id'],length=len(r['sequence']),reason='shorter_than_k') for r in records if len(r['sequence'])<k],['id','length','reason'])
        save_json(folder/'graph_semantics.json',dict(version=VERSION,k=k,min_support=minimum,support_scope='max_unique_molecules_within_one_sample',node='exact oriented unitig, not a repeat unit',edge=f'exact {k-1}-bp overlap observed in an input path',support='unique observed molecules; reference/candidate do not add molecule support',supported_view=f'at least {minimum} molecules within one sample per raw node/edge; reference retained explicitly; all interrupted paths split',path='exact reference, observed fragment or uncertain candidate; no topology-enumerated allele',reference_context='catalog repeat interval only; no synthetic flanks attached',orientation='stored BAM sequence on target reference orientation; no canonical reverse-complement collapse',compaction='node-disjoint nonbranching chains; all path endpoints preserved',scope='locus union; sample and platform evidence retained; local support does not phase a full haplotype'))
        save_json(folder/'metrics.json',row)
    tsv(out/'k_comparison.tsv',comparison)
    save_json(out/'k_comparison.json',comparison)
    return comparison


# 命令行入口：检查参数及输出目录，再冻结输入并逐 k 构图。
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-dir',type=Path,required=True)
    parser.add_argument('--bam-dir',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--catalog',type=Path,required=True)
    parser.add_argument('--models',type=Path,required=True)
    parser.add_argument('--locus',required=True)
    parser.add_argument('--k',type=int,nargs='+',default=[17,21,31,33,49,65])
    parser.add_argument('--samtools',required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--min-support',type=int,default=3)
    args=parser.parse_args()
    if any(k<2 for k in args.k) or len(set(args.k))!=len(args.k) or args.min_support<1:
        parser.error('Unique k>=2 and positive support are required')
    if args.out.exists() and any(args.out.iterdir()):
        parser.error('Refusing to overwrite a nonempty output directory')
    args.out.mkdir(parents=True,exist_ok=True)
    from vntrdbg_inputs import prepare
    records,annotations=prepare(args)
    build_all(records,annotations,args.locus,args.k,args.out,args.min_support)


if __name__=='__main__':
    main()
