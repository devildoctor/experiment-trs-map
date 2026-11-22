const assert=require('node:assert/strict');
const {projectStateGraph,stateId,summaryId,copyId,containerId,stateContainerId}=require('../scripts/supernode_projection');
const paths={},variants={a:'A',c:'C',g:'G',ac:'AC',empty:'',ag:'AG'};
function visit(p,sn,n,original=['a'],represented=original){
  const vs=paths[p]??=[],id=p+':v'+vs.length;
  vs.push({visit_id:id,supernode_id:sn,repeat_count:n,original_variant_ids:original,
    represented_variant_ids:represented,ordered_occurrence_ids:original.map((_,i)=>id+':'+i),
    original_gfa_walk:original.map((_,i)=>i+1),represented_gfa_walk:original.map((_,i)=>i+11),
    loss_event_ids:original.join()===represented.join()?[]:['loss']});
}
for(const p of ['P1','P2','P3','P4','REF'])visit(p,'L',1);
visit('P1','X',2,['a','c']);visit('P1','L',1);visit('P1','X',2,['a','g'],['a','c']);
visit('P2','X',2,['a','c']);visit('P3','X',1);visit('REF','X',6,Array(6).fill('a'));
for(const p of Object.keys(paths))visit(p,'R',1);
visit('P2','AB',3,['a','g','a','g','a','g']);
const samples={P1:{S:['m1','m2']},P2:{S:['m1']},P3:{T:['m3']},P4:{S:['m4']},REF:{}};
const info=Object.fromEntries(Object.entries(samples).map(([p,s])=>[p,{role:p==='REF'?'reference':'observation',molecule_keys:Object.values(s).flat()}]));
function freeze(x){Object.values(x).forEach(v=>{if(v && typeof v==='object' && !Object.isFrozen(v))freeze(v)});return Object.freeze(x)}
const bundle=freeze({paths,variants,samples,path_info:info});
const sid=stateId('X',2),expandedSupernodes=new Set(['X','AB']);
const project=(extra={})=>projectStateGraph(bundle,{expandedSupernodes,...extra});
const r=project(),s=r.states.get(sid);
assert.deepEqual(r.blocks.get('X').state_ids.map(id=>r.states.get(id).repeat_count),[1,2,6]);
assert.deepEqual([s.visit_count,s.path_count,s.molecule_count,s.original_sequence_count,s.represented_sequence_count,s.changed_visit_ids.length],[3,2,2,2,1,1]);
assert.deepEqual(r.routes.get('P1').map(v=>v.display_node_ids[0]),[summaryId('L'),sid,summaryId('L'),sid,summaryId('R')]);
assert.deepEqual(r.routes.get('P4').map(v=>v.display_node_ids[0]),[summaryId('L'),summaryId('R')]);
assert.equal(r.nodes.get(sid).parent,containerId('X'));
assert.equal(new Set(r.edges.map(e=>e.id)).size,r.edges.length);
const incoming=r.edges.filter(e=>e.target===sid);
assert.equal(new Set(incoming.map(e=>e.inlet)).size,3);
const filtered=project({sample:'S',pathIds:['P1','P2','REF','P3']});
assert.deepEqual(filtered.support,{path_count:3,molecule_count:2});
assert.deepEqual(filtered.blocks.get('X').state_ids,[sid,stateId('X',6)]);
assert.equal(filtered.edges.find(e=>e.id===incoming[0].id).inlet,incoming[0].inlet);
assert.equal(project({pathIds:[]}).nodes.size,0);
const opts={expandedStates:new Set([sid,stateId('AB',3)])},opened=project(opts),rep=project({...opts,mode:'represented'});
assert.deepEqual(opened.routes,rep.routes);assert.deepEqual(opened.states,rep.states);assert.deepEqual(opened.edges,rep.edges);
assert.deepEqual([...opened.nodes.keys()],[...rep.nodes.keys()]);
assert.equal(opened.routes.get('P2').at(-1).display_node_ids.length,6);
const id=copyId('P1:v3','P1:v3:1');
assert.equal(opened.nodes.get(id).parent,stateContainerId('X',2));
assert.deepEqual([opened.nodes.get(id).dna,rep.nodes.get(id).dna,rep.nodes.get(id).source_node_id],['G','C',12]);
for(const [p,route] of opened.routes){
  const ids=route.flatMap(v=>v.display_node_ids);
  assert.deepEqual(opened.edges.filter(e=>e.path_id===p).map(e=>[e.source,e.target]),ids.slice(1).map((n,i)=>[ids[i],n]));
  assert.ok(ids.every(n=>!opened.nodes.get(n).view_kind.endsWith('container')));
}
// Alternative token boundaries can spell the same local string.
const twin=structuredClone(bundle);twin.variants.aa='AA';twin.variants.aac='AAC';
twin.paths.P1[1].original_variant_ids=['aa','c'];twin.paths.P2[1].original_variant_ids=['a','ac'];
assert.equal(projectStateGraph(twin,{pathIds:['P1','P2']}).states.get(sid).original_sequence_count,2);
console.log('PASS states, repeated visits and edges, bypass, stable ports, scoped support, exact DNA types, ordered copies and both modes');
