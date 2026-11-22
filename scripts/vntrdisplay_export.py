"""Four-budget v2 exports, verifiable expansion, and a single connected canvas."""
from collections import Counter,defaultdict
from hashlib import sha256
from html import escape
from pathlib import Path
import csv
import json
import math
import re
import shutil
import subprocess

from vntrpathstructure import annotate_paths,restore_original
from vntrpathgraph import write_gfa,validate_gfa,validate_graph
from vntrdisplay_loss import build_unit_graph,represented_graph,merge_paths,display_signature
from vntrdisplay import build_display,validate_display,add_support


def dump(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def table(path,rows,fields=None):
    fields=fields or (list(rows[0]) if rows else [])
    with path.open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for row in rows:
            w.writerow({k:json.dumps(row.get(k),ensure_ascii=False) if isinstance(row.get(k),(dict,list,tuple)) else row.get(k,'') for k in fields})


def node_expression(node):
    if node['kind']=='module':return ' / '.join(node_expression(p) for p in node['parts'])
    types='|'.join(node['protected_types'])
    if node['kind']=='segment':return f'{types} [{len(node["unit_sequence"])}bp]'
    counts=node['counts'];count=str(counts[0]) if len(counts)==1 else '{'+','.join(map(str,counts))+'}'
    return f'{types} x {count}'


def render_display(destination,display,losses,records,originals,represented,groups,dot='dot'):
    name='D1S80';loss_by_id={x['path']:x for x in losses};ref_nodes={v['node'] for v in display['paths']['P0001']}
    lines=['digraph D1S80 {','graph [rankdir=LR, bgcolor="white", pad=0.3, nodesep=0.14, ranksep=0.6, splines=polyline, outputorder=edgesfirst, nslimit=2, nslimit1=2, mclimit=0.2, searchsize=5];',
           'node [shape=box, style="rounded,filled", fontname="DejaVu Sans", fontsize=10, margin="0.10,0.07"];',
           'edge [color="#94a3b8", arrowsize=0.4, penwidth=0.7];']
    def quoted(s):return json.dumps(str(s),ensure_ascii=True)
    by_id={r['id']:r for r in records}
    for n,node in display['nodes'].items():
        expression=node_expression(node);parts=expression.split(' / ')
        label='\n'.join(' / '.join(parts[i:i+3]) for i in range(0,len(parts),3))
        label+='\nbp={'+','.join(map(str,node['represented_lengths']))+'}'
        label+=f'\n{n} | M={node["route_molecules"]}; exact={node["raw_exact_molecules"]}'
        if node['lossy_member_ids']:label+=' | LOSS*'
        if 'annotation_ambiguity' in node['uncertainty_flags']:label+=' | ?'
        ps=node['path_ids'];fractions=[loss_by_id.get(p,{}).get('fraction',0) for p in ps]
        tooltip=dict(expression=expression,groups=node['structure_group_ids'],paths=ps,
            original_lengths=sorted({len(originals[p]) for p in ps}),represented_lengths=sorted({len(represented[p]) for p in ps}),
            max_loss_fraction=max(fractions),mean_loss_fraction=sum(fractions)/len(fractions),
            route_molecules=node['route_molecules'],raw_exact_molecules=node['raw_exact_molecules'],
            lossy_members=node['lossy_member_ids'],uncertainty=node['uncertainty_flags'])
        fill='#fff1d6' if node['lossy_member_ids'] else '#eff6ff' if n in ref_nodes else '#f5f7fa'
        color='#2563eb' if n in ref_nodes else '#e59a28' if node['lossy_member_ids'] else '#94a3b8'
        lines.append(f'{n} [label={quoted(label)}, tooltip={quoted(json.dumps(tooltip,ensure_ascii=False))}, fillcolor="{fill}", color="{color}", penwidth={2 if n in ref_nodes else 1}];')
    ref_edges={(a['node'],b['node']) for a,b in zip(display['paths']['P0001'],display['paths']['P0001'][1:])}
    for e in display['edges']:
        a,b=e['source'],e['target'];gs=e['structure_group_ids']
        label=(','.join(gs) if len(gs)<=3 else f'{len(gs)}G')+f' M{e["route_molecules"]}'
        extra=', color="#2563eb", penwidth=1.7' if (a,b) in ref_edges else ''
        lines.append(f'{a} -> {b} [label={quoted(label)}, fontsize=7, tooltip={quoted(json.dumps(e,ensure_ascii=False))}{extra}];')
    lines.append('}')
    dotfile=destination/f'{name}.display.dot';dotfile.write_text('\n'.join(lines)+'\n',encoding='ascii')
    subprocess.run([dot,'-Tsvg',str(dotfile),'-o',str(destination/f'{name}.display.svg')],check=True,timeout=240)


def write_local_gfa_support(path,graph,records,losses):
    # Reuse the same local projection definition on each literal GFA segment.
    d=dict(nodes={str(n):{} for n in graph['nodes']},
        edges=[dict(source=str(a),target=str(b)) for a,b in sorted(graph['edges'])],
        paths={p:[dict(node=str(n),gfa_nodes=[n]) for n in walk] for p,walk in graph['paths'].items()})
    add_support(d,graph,records,losses,{r['id']:r['id'] for r in records})
    result=[]
    for line in path.read_text(encoding='ascii').splitlines():
        f=line.split('\t')
        if f[0]=='S':line+=f'\txc:i:{d["nodes"][f[1]]["raw_exact_molecules"]}'
        elif f[0]=='L':
            e=next(x for x in d['edges'] if x['source']==f[1] and x['target']==f[3])
            line+=f'\txc:i:{e["raw_exact_molecules"]}'
        result.append(line)
    path.write_text('\n'.join(result)+'\n',encoding='ascii')
    return d


def export_budget(destination,records,features,model,policy,fraction,exact,metadata=None,dot='dot',render=True):
    destination.mkdir(parents=True,exist_ok=True);metadata=metadata or {}
    represented,reps,losses,decisions=merge_paths(records,features,model,policy,fraction,exact,build_display)
    originals={r['id']:r['sequence'] for r in records};by_id={r['id']:r for r in records};loss_by_id={x['path']:x for x in losses}
    graph=represented_graph(exact,reps);validate_graph(graph,represented)
    group_names={rep:f'G{i:04d}' for i,rep in enumerate(sorted(set(reps.values())),1)}
    groups={p:group_names[q] for p,q in reps.items()}
    display=build_display(graph,features,reps,model);add_support(display,graph,records,losses,groups)
    validation=dict(display=validate_display(display,graph,represented),protected_structure=True,
        reference_unchanged=represented['P0001']==originals['P0001'],original_reconstruction=True,
        direct_representatives=True,budgets=True)
    for p,q in reps.items():
        loss=loss_by_id.get(p,{})
        assert display_signature(features[p])==display_signature(features[q]),p
        assert restore_original(represented[p],loss.get('events',[]))==originals[p],p
        assert loss.get('loss_bp',0)<=math.floor(len(originals[p])*fraction+1e-9),p
        assert reps[q]==q and (q==p or by_id[q].get('role')!='reference'),p
    for label,g,seqs,ledger in [('exact',exact,originals,[]),('structure',graph,represented,losses)]:
        f=destination/f'D1S80.{label}.gfa';write_gfa(f,g,records,'structure_v2_'+label,seqs,ledger)
        write_local_gfa_support(f,g,records,ledger);validation[label]=validate_gfa(f,seqs)
    for label,data in [('original',originals),('represented',represented)]:
        (destination/f'paths.{label}.fa').write_text(''.join(f'>{p}\n{s}\n' for p,s in sorted(data.items())),encoding='ascii')
    rows=[]
    for r in sorted(records,key=lambda r:r['id']):
        p=r['id'];loss=loss_by_id.get(p,{})
        rows.append(dict(path=p,role=r.get('role','observation'),structure_group=groups[p],representative=reps[p],
            molecules=len(r['molecules']),hifi_molecules=len(r['hifi_molecules']),original_length=len(originals[p]),
            represented_length=len(represented[p]),substitutions=loss.get('substitutions',0),insertions=loss.get('insertions',0),
            deletions=loss.get('deletions',0),loss_bp=loss.get('loss_bp',0),loss_fraction=loss.get('fraction',0),
            original_expression=features[p]['expression'],represented_expression=features[reps[p]]['expression'],
            original_sequence=originals[p],represented_sequence=represented[p]))
    table(destination/'D1S80.path_sequences.txt',rows);table(destination/'paths.tsv',rows)
    # Re-read the requested TXT, including full untruncated strings.
    with (destination/'D1S80.path_sequences.txt').open(encoding='utf-8') as h:
        readback=list(csv.DictReader(h,delimiter='\t'))
    assert {r['path']:r['original_sequence'] for r in readback}==originals
    assert {r['path']:r['represented_sequence'] for r in readback}==represented
    validation['txt_roundtrip']=True
    group_rows=[]
    for rep,g in sorted(group_names.items()):
        ps=sorted(p for p,q in reps.items() if q==rep);mols=set().union(*(set(by_id[p]['molecules']) for p in ps))
        group_rows.append(dict(group=g,representative=rep,paths=ps,route_molecules=len(mols),
            representative_raw_exact_molecules=len(by_id[rep]['molecules']),expression=features[rep]['expression']))
    table(destination/'structure_groups.tsv',group_rows)
    table(destination/'losses.tsv',losses,['path','target_path','loss_bp','substitutions','insertions','deletions','original_length','represented_length','fraction','events','policy_version','protected_structure'])
    (destination/'losses.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in losses),encoding='utf-8')
    dump(destination/'display_graph.json',display);dump(destination/'merge_decisions.json',decisions)
    dump(destination/'alignment_audit.json',exact['alignments'])
    dump(destination/'annotation_audit.json',{p:f['annotation'] for p,f in features.items()})
    table(destination/'nodes.tsv',list(display['nodes'].values()),['id','kind','protected_types','counts','represented_lengths',
        'uncertainty_flags','path_ids','structure_group_ids','route_molecules','raw_exact_molecules','lossy_member_ids'])
    table(destination/'edges.tsv',display['edges'])
    dump(destination/'source_members.json',records);dump(destination/'structure_policy.json',policy)
    dump(destination/'locus_model.json',model)
    molecules=set().union(*(set(r['molecules']) for r in records))
    weighted=sum(x['loss_bp']*len(by_id[x['path']]['molecules']) for x in losses)
    summary=dict(loss_budget_fraction=fraction,**display['metrics'],paths=len(records),molecules=len(molecules),
        structure_groups=len(group_names),lossy_paths=len(losses),total_edits_bp=sum(x['loss_bp'] for x in losses),
        substitutions=sum(x['substitutions'] for x in losses),insertions=sum(x['insertions'] for x in losses),
        deletions=sum(x['deletions'] for x in losses),max_path_loss_fraction=max((x['fraction'] for x in losses),default=0),
        max_path_loss_bp=max((x['loss_bp'] for x in losses),default=0),molecule_weighted_edit_bp=weighted,
        molecule_weighted_original_bp=sum(len(r['sequence'])*len(r['molecules']) for r in records),
        accepted_batches=len(decisions['accepted_batches']),literal_gfa_nodes=len(graph['nodes']),
        literal_gfa_edges=len(graph['edges']),metadata=metadata,validation=validation)
    dump(destination/'summary.json',summary);dump(destination/'validation.json',validation)
    if render:render_display(destination,display,losses,records,originals,represented,groups,dot)
    return summary


def select_budget(reports,policy):
    baseline=reports[0]['C'];maximum=max(baseline-r['C'] for r in reports)
    nondominated=[r for r in reports if not any(s['total_edits_bp']<=r['total_edits_bp'] and s['C']<=r['C'] and
        (s['total_edits_bp']<r['total_edits_bp'] or s['C']<r['C']) for s in reports)]
    qualifying=[r for r in nondominated if baseline-r['C']>=policy['selection_gain_fraction']*maximum-1e-9]
    chosen=min(qualifying,key=lambda r:(r['total_edits_bp'],r['loss_budget_fraction']))
    return dict(selected_fraction=chosen['loss_budget_fraction'],baseline_C=baseline,maximum_reduction=maximum,
                achieved_reduction=baseline-chosen['C'],target_fraction=policy['selection_gain_fraction'],
                nondominated_fractions=[r['loss_budget_fraction'] for r in nondominated],
                criterion='minimum_actual_total_loss_achieving_90_percent_of_maximum_C_reduction')


def comparison_html(destination,reports,selection):
    tr=''.join('<tr>'+''.join('<td>'+escape(str(r[k]))+'</td>' for k in ['loss_budget_fraction','N_d','E_d','C','branch_nodes','independent_blocks','structure_groups','lossy_paths','total_edits_bp','max_path_loss_fraction'])+'</tr>' for r in reports)
    data=json.loads((destination/'display_graph.json').read_text(encoding='utf-8'))
    source=json.loads((destination/'source_members.json').read_text(encoding='utf-8'))
    opts='<option value="">全部路径</option>'+''.join(f'<option value="{r["id"]}">{r["id"]} · {r.get("role","observation")} · {len(r["sequence"])} bp</option>' for r in source)
    payload=json.dumps(data,ensure_ascii=False).replace('</',r'<\/')
    svg=(destination/'D1S80.display.svg').read_text(encoding='utf-8')
    svg=re.sub(r'<\?xml[^>]*\?>|<!DOCTYPE[\s\S]*?>','',svg)
    with (destination/'paths.tsv').open(encoding='utf-8') as h:path_rows={r['path']:r for r in csv.DictReader(h,delimiter='\t')}
    row_payload=json.dumps(path_rows,ensure_ascii=False).replace('</',r'<\/')
    html=r'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>D1S80 损失与图复杂度</title>
<style>body{font-family:Arial,"Microsoft YaHei",sans-serif;color:#22354a;margin:28px}p{line-height:1.7}table{border-collapse:collapse}th,td{padding:9px;border-bottom:1px solid #ddd;text-align:right}button,select{padding:8px;margin:4px}#viewport{height:70vh;overflow:hidden;border:1px solid #cbd5e1;background:#fafcff;touch-action:none;cursor:grab}#graph,#graph>svg{display:block;width:100%;height:100%}.dim{opacity:.08}.active path,.active polygon{stroke:#ea770c!important;stroke-width:2!important}pre{white-space:pre-wrap;word-break:break-all;background:#f1f5f9;padding:14px}a{color:#2563eb}</style>
<h1>D1S80：损失与图复杂度</h1>
<p>全部307条P：1条参考、306种观测核心；370个分子。参考以蓝色突出，橙色*标明存在有损成员，?保留划分歧义。单元数、U01–U09类型及顺序受保护。次数集合仅列已观测值。</p>
<p>选定预算：SELECTED。损失按原始长度计算；原始序列始终保留。这是模型下的展示压缩，不是基因型或测序错误校正。主图为单张SVG，所有分支和独立块均保留。</p>
<table><tr><th>预算</th><th>节点</th><th>边</th><th>C</th><th>分叉</th><th>独立块</th><th>组</th><th>有损P</th><th>编辑bp</th><th>最大P损失</th></tr>ROWS</table>
<p><a href="D1S80.display.svg">新图SVG</a> · <a href="D1S80.path_sequences.txt">原始／表示序列TXT</a> · <a href="D1S80.structure.gfa">表示GFA</a> · <a href="D1S80.exact.gfa">原始GFA</a> · <a href="simplification_comparison.tsv">四档比较</a></p>
<div><button onclick="zoom(1.25)">放大</button><button onclick="zoom(.8)">缩小</button><button onclick="fit()">显示全图</button><button onclick="focusPath()">定位所选P</button><span>滚轮缩放、拖动平移；悬停查看成员、精确支持和损失。选择P可突出对应路线。</span></div>
<div id="viewport"><div id="graph">INLINE_SVG</div></div>
<h2>逐P追溯</h2><select id="path" onchange="showPath()">OPTIONS</select><pre id="detail"></pre>
<p>只有完整登记的P路线及每次访问参数属于本次表示。不能任意组合图中分支生成新的观测。小端部补齐计入损失，不改写原始B2切点。布局坐标不是唯一生物学事件位置。</p>
<script>const data=PAYLOAD,rows=ROW_DATA,svg=document.querySelector('#graph svg'),viewport=document.getElementById('viewport');
const original=svg.getAttribute('viewBox').split(/\s+/).map(Number);let box=original.slice(),drag=null;
function apply(){svg.setAttribute('viewBox',box.join(' '))}function fit(){box=original.slice();apply()}
function zoom(f,px=.5,py=.5){const w=box[2]/f,h=box[3]/f;box=[box[0]+(box[2]-w)*px,box[1]+(box[3]-h)*py,w,h];apply()}
function showPath(){const p=document.getElementById('path').value,vs=data.paths[p]||[],ns=new Set(vs.map(v=>v.node)),es=new Set(vs.slice(1).map((v,i)=>vs[i].node+'->'+v.node));document.querySelectorAll('#graph g.node,#graph g.edge').forEach(g=>{const t=g.querySelector('title').textContent,on=g.classList.contains('node')?ns.has(t):es.has(t);g.classList.toggle('dim',!!p&&!on);g.classList.toggle('active',!!p&&on)});document.getElementById('detail').textContent=p?JSON.stringify({path:p,...rows[p],visits:vs},null,2):'全部307条P均在图中；选择P查看原始／表示序列、损失和精确次数。'}
function focusPath(){const p=document.getElementById('path').value;if(!p){fit();return}const ns=new Set(data.paths[p].map(v=>v.node));let x=Infinity,y=Infinity,x2=-Infinity,y2=-Infinity;document.querySelectorAll('#graph g.node').forEach(g=>{if(!ns.has(g.querySelector('title').textContent))return;const b=g.getBBox(),m=svg.getScreenCTM().inverse().multiply(g.getScreenCTM());[[b.x,b.y],[b.x+b.width,b.y+b.height]].forEach(([a,c])=>{const pt=svg.createSVGPoint();pt.x=a;pt.y=c;const q=pt.matrixTransform(m);x=Math.min(x,q.x);y=Math.min(y,q.y);x2=Math.max(x2,q.x);y2=Math.max(y2,q.y)})});box=[x-30,y-30,Math.max(x2-x+60,100),Math.max(y2-y+60,100)];apply()}
viewport.addEventListener('wheel',e=>{e.preventDefault();const pt=svg.createSVGPoint();pt.x=e.clientX;pt.y=e.clientY;const q=pt.matrixTransform(svg.getScreenCTM().inverse());zoom(Math.exp(-e.deltaY*.0015),(q.x-box[0])/box[2],(q.y-box[1])/box[3])},{passive:false});
viewport.addEventListener('pointerdown',e=>{const m=svg.getScreenCTM();drag={x:e.clientX,y:e.clientY,box:box.slice(),sx:m.a,sy:m.d};viewport.setPointerCapture(e.pointerId)});viewport.addEventListener('pointermove',e=>{if(!drag)return;box=[drag.box[0]-(e.clientX-drag.x)/drag.sx,drag.box[1]-(e.clientY-drag.y)/drag.sy,drag.box[2],drag.box[3]];apply()});viewport.addEventListener('pointerup',()=>drag=null);viewport.addEventListener('pointercancel',()=>drag=null);showPath();</script></html>'''
    html=html.replace('SELECTED',f'{selection["selected_fraction"]:.0%}').replace('ROWS',tr).replace('OPTIONS',opts).replace('PAYLOAD',payload).replace('ROW_DATA',row_payload).replace('INLINE_SVG',svg)
    (destination/'index.html').write_text(html,encoding='utf-8')


def run_comparison(source,model_path,policy_path,locus,destination,dot='dot'):
    destination.mkdir(parents=True,exist_ok=False)
    records=json.loads(source.read_text(encoding='utf-8'));policy=json.loads(policy_path.read_text(encoding='utf-8'))
    assert len({r['id'] for r in records})==len(records)
    metadata=dict(source=str(source),source_sha256=sha256(source.read_bytes()).hexdigest(),
                  model_sha256=sha256(model_path.read_bytes()).hexdigest(),policy_sha256=sha256(policy_path.read_bytes()).hexdigest(),
                  source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=str(Path(__file__).resolve().parents[1])).decode().strip())
    print('ANNOTATE',len(records),flush=True)
    model,model_hash,features=annotate_paths(records,model_path,locus,select_representative=True)
    dump(destination/'path_features.json',features)
    print('ALIGN',flush=True);exact=build_unit_graph(records,features,model)
    reports=[]
    for fraction in policy['loss_fractions']:
        folder=destination/f'{round(fraction*100)}pct';print('BUDGET',fraction,flush=True)
        report=export_budget(folder,records,features,model,policy,fraction,exact,metadata,dot,render=False)
        reports.append(report);print(json.dumps({k:v for k,v in report.items() if k not in ('metadata','validation')}),flush=True)
    selection=select_budget(reports,policy);chosen=destination/f'{round(selection["selected_fraction"]*100)}pct'
    for f in chosen.iterdir():
        if f.is_file():shutil.copy2(f,destination/f.name)
    table(destination/'simplification_comparison.tsv',[{k:v for k,v in r.items() if k not in ('metadata','validation')} for r in reports])
    dump(destination/'selection.json',selection);dump(destination/'all_budget_summaries.json',reports)
    # Render every budget after candidate evaluation; never call a layout engine
    # from the greedy objective. The chosen one is also delivered at root.
    for report in reports:
        folder=destination/f'{round(report["loss_budget_fraction"]*100)}pct'
        d=json.loads((folder/'display_graph.json').read_text(encoding='utf-8'))
        with (folder/'paths.tsv').open(encoding='utf-8') as h:rows=list(csv.DictReader(h,delimiter='\t'))
        losses=[json.loads(s) for s in (folder/'losses.jsonl').read_text().splitlines()]
        render_display(folder,d,losses,records,{r['path']:r['original_sequence'] for r in rows},
            {r['path']:r['represented_sequence'] for r in rows},{r['path']:r['structure_group'] for r in rows},dot)
    for suffix in ('svg','dot'):shutil.copy2(chosen/f'D1S80.display.{suffix}',destination/f'D1S80.display.{suffix}')
    comparison_html(destination,reports,selection)
    print('COMPARISON_COMPLETE',json.dumps(selection),flush=True)
    return reports,selection
