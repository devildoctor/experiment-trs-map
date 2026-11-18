(() => {
  function check(ok,msg){if(!ok)throw Error(msg)}
  const set=(id,value)=>{el(id).value=value;el(id).dispatchEvent(new Event('change'))};
  const recorded=new Map(Object.values(paths).flat().map(v=>[v.visit_id,v]));
  const sample=el('sample').options[1]?.value;
  if(sample){set('sample',sample);check(visiblePaths().every(p=>data.samples[p][sample]?.length),'sample filter');set('sample','')}
  const module=Array.from(recorded.values()).find(v=>nodes[v.supernode_id].kind==='tandem_module' && v.repeat_count>=2);
  check(module,'expected a registered tandem module');
  set('pathA',module.path_id);
  const other=Object.keys(paths).find(p=>p!==module.path_id);
  set('pathB',other);el('only').checked=true;el('only').dispatchEvent(new Event('change'));
  function checkHighlight(p,cls){
    const vs=paths[p],ids=vs.map(v=>v.supernode_id),edges=ids.slice(1).map((n,i)=>'e:'+ids[i]+'>'+n);
    check(cy.nodes('.'+cls).length===new Set(ids).size && ids.every(id=>cy.getElementById(id).hasClass(cls)),'node walk '+p);
    check(cy.edges('.'+cls).length===new Set(edges).size && edges.every(id=>cy.getElementById(id).hasClass(cls)),'edge walk '+p);
  }
  checkHighlight(module.path_id,'a');checkHighlight(other,'b');
  cy.getElementById(module.supernode_id).emit('tap');set('visit',module.visit_id);
  const before=cy.zoom(),support=el('info').textContent,meta=el('meta').textContent;
  const unrelated=cy.nodes().filter(n=>n.id()!==module.supernode_id)[0];
  const pos=unrelated?{...unrelated.position()}:null;
  el('toggle').click();
  const units=cy.nodes().filter(n=>n.data('visit')===module.visit_id);
  check(units.length===module.ordered_occurrence_ids.length,'expanded occurrences');
  check(units.every((n,i)=>n.data('source_node_id')===module.original_gfa_walk[i]),'ordered original walk');
  check(el('details').querySelectorAll('li').length===units.length,'unit details');
  check(el('info').textContent===support && el('meta').textContent===meta && cy.zoom()===before,'stable support and zoom');
  if(pos)check(unrelated.position().x===pos.x && unrelated.position().y===pos.y,'stable unrelated position');
  set('mode','represented');
  check(units.every((n,i)=>n.data('source_node_id')===module.represented_gfa_walk[i]),'represented walk');
  el('toggle').click();checkHighlight(module.path_id,'a');
  const loss=Array.from(recorded.values()).find(v=>v.loss_event_ids.length);
  if(loss){
    set('pathA',loss.path_id);detail(loss.supernode_id,loss.visit_id);el('toggle').click();
    check(el('details').textContent.includes(data.losses[loss.loss_event_ids[0]].representative_occurrence),'donor details');
    const changed=cy.nodes().filter(n=>n.data('visit')===loss.visit_id && n.data('flag')==='loss');
    check(changed.length>0 && changed.every(n=>n.data('label')===data.variants[loss.represented_variant_ids[loss.ordered_occurrence_ids.indexOf(n.id().slice(2))]]),'local substitutions');
  }
  const countSet=Object.values(nodes).find(n=>JSON.stringify(n.observed_repeat_counts)==='[3,5]');
  if(countSet)check(label(countSet).includes('{3,5}'),'literal observed count set');
  check(!performance.getEntriesByType('resource').some(r=>/^https?:/.test(r.name)),'offline dependencies');
  return {cytoscape:cytoscape.version,paths:Object.keys(paths).length,sample,compared:[module.path_id,other],
          expanded_visit:module.visit_id,pattern:pattern(nodes[module.supernode_id]),count:module.repeat_count,
          expanded_units:units.length,loss_visit:loss?.visit_id||null,checked_counts_3_5:!!countSet,
          offline:true,path_highlight:true,stable_support_and_viewport:true};
})()
