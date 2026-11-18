(() => {
  function check(value,message){if(!value)throw Error(message)}
  const change=(id,value)=>{el(id).value=value;el(id).dispatchEvent(new Event('change'))};
  const options=id=>Array.from(el(id).options,o=>o.value);
  const real=!!paths.P0188, sample=Object.keys(data.samples.P0002)[0];
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
  const sn=moduleVisit.supernode_id;cy.getElementById(sn).emit('tap');change('visit',moduleVisit.visit_id);
  el('toggle').click();check(cy.nodes().some(n=>n.data('visit')===moduleVisit.visit_id),'overview inline expansion');
  el('toggle').click();
  check(!performance.getEntriesByType('resource').some(r=>/^https?:/.test(r.name)),'offline');
  return {initial_view:'overview',paths:Object.keys(paths).length,cytoscape:cytoscape.version,same_path:sameSummary,
    module_visit:moduleVisit.visit_id,ordered_occurrences:moduleVisit.ordered_occurrence_ids.length,
    loss_visit:lossVisit,dual_walks:true,stable_viewport:true,overview_expand_collapse:true,offline:true};
})()
