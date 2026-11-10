"""Self-contained family/expanded graph viewer with per-path exact payloads."""
import json
from pathlib import Path


def write_viewer(out,locus):
    data={}
    for mode,directory in [('family','family_3pct'),('expanded','expanded_3pct')]:
        folder=out/directory;graph=json.loads((folder/'display_graph.json').read_text())
        svg=(folder/(locus+'.display.svg')).read_text(encoding='utf-8')
        data[mode]=dict(graph=graph,svg=svg[svg.index('<svg'):],summary=json.loads((folder/'summary.json').read_text()))
    payload=json.dumps(data,ensure_ascii=True).replace('</','<\\/')
    page=r'''<!doctype html>
<html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>__LOCUS__ 家族图</title>
<style>
*{box-sizing:border-box}body{margin:0;background:#f5f7fb;color:#142334;font:14px system-ui,"Microsoft YaHei",sans-serif}
header{padding:20px 26px;background:white;border-bottom:1px solid #dbe3ec}h1{font-size:23px;margin:0 0 8px}
.note{color:#53677c;margin:6px 0 12px}.bar{display:flex;gap:12px;align-items:center;flex-wrap:wrap}
button,select{padding:8px 12px;border:1px solid #aebed0;border-radius:6px;background:white;color:#142334}button{cursor:pointer}
button.active{background:#234f73;color:white}#metrics{font-weight:600;margin-left:auto}
main{display:grid;grid-template-columns:minmax(0,1fr) 390px;height:calc(100vh - 155px)}
#canvas{overflow:auto;padding:18px;background:white}#canvas svg{display:block;max-width:none;height:auto}
aside{overflow:auto;padding:18px;border-left:1px solid #dbe3ec}h2{font-size:18px;margin:0 0 10px}
.card{background:white;border:1px solid #dbe3ec;border-radius:7px;padding:12px;margin:10px 0}
pre{font:12px/1.6 Consolas,monospace;white-space:pre-wrap;overflow-wrap:anywhere;margin:8px 0 0}
.node{cursor:pointer}.node.selected polygon,.node.selected path{stroke:#d81b60;stroke-width:3px}
.muted{opacity:.16}@media(max-width:900px){main{grid-template-columns:1fr;height:auto}#canvas{height:65vh}aside{border-top:1px solid #dbe3ec}#metrics{margin-left:0}}
</style>
<header><h1>__LOCUS__ · 单元家族与路径</h1>
<div class="note">家族汇总相近单元；点击节点查看各条路径的实际类型、次数和字符串。橙色 LOSS* 表示该节点至少一次访问有局部改写。</div>
<div class="bar"><button id="family" class="active">家族视图</button><button id="expanded">展开全部家族</button>
<label>路径 <select id="path"><option value="">全部路径</option></select></label>
<label>缩放 <input id="zoom" type="range" min="0.05" max="1.5" step="0.05" value="0.25"></label><span id="metrics"></span></div></header>
<main><section id="canvas"></section><aside id="detail"><h2>选择一个节点</h2><p>可先选定一条路径，再查看它在家族节点中的具体单元。</p></aside></main>
<script id="data" type="application/json">__DATA__</script>
<script>
const data=JSON.parse(document.getElementById('data').textContent);
let mode='family',selected=null;
const canvas=document.getElementById('canvas'),detail=document.getElementById('detail'),pathSelect=document.getElementById('path');
for(const p of Object.keys(data.family.graph.paths).sort()){const o=document.createElement('option');o.value=p;o.textContent=p;pathSelect.append(o)}
function text(tag,value,parent){const e=document.createElement(tag);e.textContent=value;parent.append(e);return e}
function zoom(){const s=canvas.querySelector('svg');if(s)s.style.width=(s.viewBox.baseVal.width*Number(document.getElementById('zoom').value))+'px'}
function highlight(){const p=pathSelect.value,g=data[mode].graph,walk=p?g.paths[p]:[],nodes=new Set(walk.map(v=>v.node));
 const edges=new Set(walk.slice(1).map((v,i)=>walk[i].node+'->'+v.node));
 for(const n of canvas.querySelectorAll('g.node'))n.classList.toggle('muted',!!p&&!nodes.has(n.querySelector('title').textContent));
 for(const e of canvas.querySelectorAll('g.edge'))e.classList.toggle('muted',!!p&&!edges.has(e.querySelector('title').textContent));
 if(selected)show(selected)}
function show(id){selected=id;detail.replaceChildren();const g=data[mode].graph,n=g.nodes[id];if(!n)return;
 text('h2',id+' · '+n.path_ids.length+' 条路径',detail);
 for(const part of n.parts){const f=part.display_family;
  text('p',f?f.id+'：'+f.members.join(' / ')+'；代表 '+f.representative:'类型：'+part.types.join(' / '),detail)}
 const chosen=pathSelect.value;
 for(const [p,visits] of Object.entries(g.paths)){if(chosen&&chosen!==p)continue;
  for(const v of visits){if(v.node!==id)continue;const c=document.createElement('div');c.className='card';detail.append(c);text('strong',p,c);
   for(const part of v.parts){text('p','次数 '+part.count+'；类型顺序 '+part.types.join(' → '),c);
    text('pre',part.variant_ids.map((x,i)=>part.types[i]+'  '+g.sequence_variants[x]+(part.annotation[i].weak_approximate?'  [weak]':'')).join('\n'),c);
    if(part.skip_context)text('p','此片段保留逐次 skip 上下文。',c)}}}
 for(const el of canvas.querySelectorAll('g.node'))el.classList.toggle('selected',el.querySelector('title').textContent===id)}
function render(){selected=null;canvas.innerHTML=data[mode].svg;const s=data[mode].summary;
 document.getElementById('metrics').textContent=s.N_d+' 节点 · '+s.E_d+' 边 · 改写 '+s.total_edits_bp+' bp · 上限 3%';
 for(const m of ['family','expanded'])document.getElementById(m).classList.toggle('active',m===mode);
 detail.replaceChildren();text('h2','选择一个节点',detail);zoom();highlight();
 for(const n of canvas.querySelectorAll('g.node'))n.addEventListener('click',()=>show(n.querySelector('title').textContent))}
for(const m of ['family','expanded'])document.getElementById(m).onclick=()=>{mode=m;render()};
document.getElementById('zoom').oninput=zoom;pathSelect.onchange=highlight;render();
</script></html>'''
    (out/(locus+'.families.html')).write_text(page.replace('__LOCUS__',locus).replace('__DATA__',payload),encoding='utf-8')
