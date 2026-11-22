(() => {
  function check(value,message){if(!value)throw Error(message)}
  const change=(id,value)=>{el(id).value=value;el(id).dispatchEvent(new Event('change'))};
  const options=id=>Array.from(el(id).options,o=>o.value);
  const real=!!paths.P0188, sample=Object.keys(data.samples.P0002)[0];
  check(el('view').value==='five' && !el('fivePanel').hidden && el('fiveRoutes').querySelector('.structRow'), 'initial observed five-path view');
  change('view','overview');
  check(el('view').value==='overview' && el('tracks').hidden && !el('globalSample').hidden &&
    cy.nodes().length===Object.keys(nodes).length,'initial global overview');
  change('view','compare');
  check(el('legend').textContent.includes('请选择'),'missing path guidance');
  change('sampleA',real?'LCL5':sample);change('pathA','P0002');
  change('sampleB',real?'LCL6':sample);change('pathB','P0002');
  check(el('meta').textContent.startsWith('1 条登记路径'),'same registered path');
  const sameSummary=el('meta').textContent;
  if(real)check(sameSummary.includes('13 个独立分子') && el('statsA').textContent.startsWith('6 个') &&
    el('statsB').textContent.startsWith('7 个'),'sample-specific molecule counts');
  function ordered(track,p){
    const row=cy.nodes('.'+track.toLowerCase()).sort((a,b)=>a.position('x')-b.position('x'));
    check(row.length===paths[p].length && row.every((n,i)=>n.data('visit')===paths[p][i].visit_id),'ordered visits '+track);
    check(row.every((n,i)=>!i || n.position('x')>row[i-1].position('x')),'monotone layout');
    const edges=cy.edges('.'+track.toLowerCase());
    check(edges.length===row.length-1 && row.slice(1).every((n,i)=>edges.some(e=>e.source().id()===row[i].id() && e.target().id()===n.id())),'registered edges');
    return row;
  }
  ordered('A','P0002');ordered('B','P0002');
  check(el('statsA').parentElement.style.borderColor===el('statsB').parentElement.style.borderColor &&
    cy.edges('.a,.b').every(e=>e.data('color')===pathColors.get('P0002')),'one path color on both tracks');
  check(cy.edges('.correspondence').length===paths.P0002.length,'same-source connectors');
  const oldB=options('pathB').join(), chosenB=el('pathB').value;
  if(real){
    change('sampleA','LCL8');
    check(!el('pathA').value && options('pathB').join()===oldB && el('pathB').value===chosenB,'independent controls');
  }
  change('pathA','P0001');check(el('statsA').textContent.startsWith('0 个') && el('statsA').textContent.includes('预算 0 bp'),'reference scope');
  if(real){
    check(!cy.edges('.correspondence').length && el('legend').textContent.includes('尚未确定'),'uncertain position');
    change('sampleB','LCL5');change('pathB','P0188');
    const links=cy.edges('.correspondence');
    check(links.length===4 && links.every(e=>e.source().position('x')===e.target().position('x')),'real anchored columns');
    for(const [side,count] of [['A',2],['B',4]]){
      const row=ordered(side,side==='A'?'P0001':'P0188');
      check(row[2].data('label').includes('×'+count),'exact visit count');
    }
  }
  const moduleVisit=paths.P0001.find(v=>nodes[v.supernode_id].kind==='tandem_module');
  const visitNode=(t,v)=>cy.getElementById(JSON.stringify([t,v.path_id,v.visit_id]));
  visitNode('A',moduleVisit).emit('tap');
  check(el('name').textContent.includes('A 行') && el('name').textContent.includes(moduleVisit.visit_id),'visit detail target');
  check(el('details').querySelectorAll('li').length===moduleVisit.ordered_occurrence_ids.length,'ordered units');
  if(!real)check(moduleVisit.repeat_count===3 && moduleVisit.ordered_occurrence_ids.length===6,'AB cubed');
  const auditUnits=(v,mode)=>Array.from(el('details').querySelectorAll('li'),(li,i)=>{
    const text=li.textContent;
    return text.includes(v.ordered_occurrence_ids[i]+' ') && text.includes(data.variants[v[mode+'_variant_ids'][i]]+' · S'+v[mode+'_gfa_walk'][i]);
  }).every(Boolean);
  check(auditUnits(moduleVisit,'original'),'original walk detail');
  const snapshot=()=>JSON.stringify([cy.nodes().map(n=>[n.id(),n.position().x,n.position().y]),cy.pan(),cy.zoom(),selections(),
    el('meta').textContent,el('statsA').textContent,el('statsB').textContent]);
  const before=snapshot();change('mode','represented');
  check(before===snapshot() && auditUnits(moduleVisit,'represented'),'stable viewport, counts and dual walks');
  let lossVisit=null;
  if(real){
    change('sampleB','LCL7');change('pathB','P0007');
    const v=paths.P0007.find(v=>v.visit_id==='P0007:v3');lossVisit=v.visit_id;visitNode('B',v).emit('tap');
    check(el('details').textContent.includes('G → A') && el('details').textContent.includes('P0097:u6') &&
      el('statsB').textContent.includes('展示整串未见于原始 OBS'),'donor provenance and whole string status');
    check(el('details').querySelector('details pre').textContent.includes('LCL7'),'molecule origin');
  }
  change('view','overview');
  const sn=moduleVisit.supernode_id,sid=stateId(sn,moduleVisit.repeat_count);
  const peers=projection.states.get(sid).path_ids.slice(0,3);
  const scope=new Set(['P0001','P0002',...peers,real?'P0188':'P0003']);
  for(const opt of el('range').options)opt.selected=scope.has(opt.value);
  el('range').dispatchEvent(new Event('change'));
  const supportBefore=el('meta').textContent, viewport=()=>JSON.stringify([cy.zoom(),cy.pan()]);
  const camera=viewport(),tap=id=>{const n=cy.getElementById(id);check(n.length===1,'present '+id);n.emit('tap')};
  tap(summaryId(sn));el('toggle').click();
  const state=projection.states.get(sid);
  check(cy.getElementById(sid).parent().id()===containerId(sn),'state compound');
  tap(sid);
  check(JSON.stringify(options('visit'))===JSON.stringify(state.member_visit_ids),'scoped state members');
  change('visit',moduleVisit.visit_id);
  check(el('details').querySelectorAll('li').length===moduleVisit.ordered_occurrence_ids.length,'per-visit sequences');
  const trace=projection.edges.find(e=>e.target===sid && e.path_id==='P0001');
  if(trace){
    const stateNode=cy.getElementById(sid),incoming=stateNode.incomers('edge.trace').map(e=>e.targetEndpoint()).sort((a,b)=>a.y-b.y);
    check(incoming.every((end,i)=>Math.abs(end.x-(stateNode.position('x')-stateNode.outerWidth()/2))<1 &&
      (!i || end.y-incoming[i-1].y>=4)),'separate physical state entry ports');
    const outgoing=stateNode.outgoers('edge.trace').map(e=>e.sourceEndpoint()).sort((a,b)=>a.y-b.y);
    check(outgoing.every((end,i)=>Math.abs(end.x-(stateNode.position('x')+stateNode.outerWidth()/2))<1 &&
      (!i || end.y-outgoing[i-1].y>=4)),'separate physical state exit ports');
    const edge=cy.getElementById(trace.id);edge.data('color','#123456');overview();
    check(edge.data('color')===pathColors.get('P0001'),'existing edge data refreshed');
    tap(edge.id());check(el('info').textContent.includes(trace.from_visit_id),'edge sources');
  }
  change('tracePath','P0001');
  check(cy.getElementById(sid).hasClass('onRoute') && cy.edges('.onRoute').every(e=>e.data('path_id')==='P0001') &&
    cy.edges('.onRoute').length===projection.edges.filter(e=>e.path_id==='P0001').length,'registered highlight');
  tap(sid);el('toggle').click();
  const parentId=stateContainerId(sn,moduleVisit.repeat_count),parent=cy.getElementById(parentId);
  check(parent.parent().id()===containerId(sn),'nested compound');
  for(const vid of state.member_visit_ids){
    const v=byId.get(vid),children=v.ordered_occurrence_ids.map(oid=>cy.getElementById(copyId(vid,oid)));
    check(children.every((n,i)=>n.parent().id()===parentId && n.data('source_node_id')===v.represented_gfa_walk[i]),'copy mapping');
    check(children.every((n,i)=>!i || (n.position('x')>children[i-1].position('x') && n.position('y')===children[0].position('y'))),'copy order');
  }
  const first=vid=>cy.getElementById(copyId(vid,byId.get(vid).ordered_occurrence_ids[0]));
  check(new Set(state.member_visit_ids.map(vid=>first(vid).position('y'))).size===state.visit_count,'one row per visit');
  const copy=first(moduleVisit.visit_id);tap(copy.id());
  check(el('details').querySelectorAll('li').length===1 && el('details').textContent.includes(copy.data('occurrence_id')),'one occurrence detail');
  let adjacent=false;
  if(real){
    const other=paths.P0001[paths.P0001.indexOf(moduleVisit)+1].supernode_id;
    tap(summaryId(other));el('toggle').click();adjacent=true;
    const containers=[cy.getElementById(containerId(sn)),cy.getElementById(containerId(other))];
    function intersects(a,b){return a.x1<b.x2 && b.x1<a.x2 && a.y1<b.y2 && b.y1<a.y2}
    for(const box of containers){
      const rect=box.boundingBox();
      for(const sibling of cy.nodes().filter(n=>n.parent().empty() && n.id()!==box.id())){
        check(!intersects(rect,sibling.boundingBox()),'top-level bounding boxes');
      }
    }
  }
  const coordinates=()=>JSON.stringify(cy.nodes().map(n=>[n.id(),n.position('x'),n.position('y')]));
  const beforeMode=coordinates();change('mode','original');
  check(coordinates()===beforeMode && viewport()===camera && supportBefore===el('meta').textContent,'stable expanded layout and support');
  check(copy.data('source_node_id')===moduleVisit.original_gfa_walk[0],'updated GFA mapping');
  tap(containerId(sn));el('toggle').click();check(cy.getElementById(parentId).empty(),'collapse descendants');
  tap(summaryId(sn));el('toggle').click();check(cy.getElementById(parentId).length===1 && viewport()===camera,'restore nested expansion');
  check(!performance.getEntriesByType('resource').some(r=>/^https?:/.test(r.name)),'offline');
  return {initial_view:'overview',paths:Object.keys(paths).length,cytoscape:cytoscape.version,same_path:sameSummary,
    module_visit:moduleVisit.visit_id,ordered_occurrences:moduleVisit.ordered_occurrence_ids.length,
    loss_visit:lossVisit,dual_walks:true,stable_viewport:true,overview_expand_collapse:true,expanded_visits:state.visit_count,adjacent_blocks_checked:adjacent,registered_trace_highlight:true,offline:true};
})()
