/* Display identities and edges are independent of the selected sequence representation. */
const displayId = (kind, ...parts) => JSON.stringify([kind, ...parts]);
const stateId = (sn, count) => displayId('state', sn, count);
const summaryId = sn => displayId('summary', sn);
const containerId = sn => displayId('container', sn);
const stateContainerId = (sn,count) => displayId('state-container',sn,count);
const copyId = (visit,occurrence) => displayId('copy',visit,occurrence);

function projectStateGraph(bundle, {pathIds, sample='', expandedSupernodes=new Set(),
                                    expandedStates=new Set(), mode='original'}={}) {
  const routes=new Map(), states=new Map(), blocks=new Map(), nodes=new Map(), edges=[];
  const members=new Map(), molecules=new Set(), ports=new Map(), total=new Map();
  const requested=pathIds==null ? null : new Set(pathIds);
  for(const p of Object.keys(bundle.paths).sort()){
    const info=bundle.path_info[p];
    const keys=info.role==='reference'?[]:sample?(bundle.samples[p][sample]||[]):info.molecule_keys;
    const visible=(!requested || requested.has(p)) && (!sample || info.role==='reference' || keys.length>0);
    if(visible){members.set(p,keys);keys.forEach(k=>molecules.add(k))}
    for(const v of bundle.paths[p]){
      const sn=v.supernode_id, sid=stateId(sn,v.repeat_count);
      // Rank in the full registered state, so filtering never renumbers the ports.
      ports.set(v.visit_id,total.get(sid)||0);
      total.set(sid,(total.get(sid)||0)+1);
      if(!visible)continue;
      if(!blocks.has(sn))blocks.set(sn,{sn,state_ids:new Set(),member_visit_ids:[],path_ids:new Set(),molecule_keys:new Set()});
      if(!states.has(sid))states.set(sid,{state_id:sid,supernode_id:sn,repeat_count:v.repeat_count,
        member_visit_ids:[],path_ids:new Set(),molecule_keys:new Set(),original:new Set(),represented:new Set(),changed_visit_ids:[]});
      const block=blocks.get(sn),state=states.get(sid);
      block.state_ids.add(sid);
      for(const group of [block,state]){
        group.member_visit_ids.push(v.visit_id);
        if(!group.path_ids.has(p)){group.path_ids.add(p);keys.forEach(k=>group.molecule_keys.add(k))}
      }
      for(const m of ['original','represented']){
        state[m].add(v[m+'_variant_ids'].map(id=>bundle.variants[id]).join(''));
      }
      if(v.loss_event_ids.length)state.changed_visit_ids.push(v.visit_id);
    }
  }
  for(const group of [...blocks.values(),...states.values()]){
    group.visit_count=group.member_visit_ids.length;
    group.path_count=group.path_ids.size;
    group.molecule_count=group.molecule_keys.size;
    group.path_ids=[...group.path_ids];group.molecule_keys=[...group.molecule_keys];
  }
  function node(id,view_kind,sn,extra={}){
    nodes.set(id,{id,view_kind,sn,...extra});
  }
  for(const [sn,block] of blocks){
    block.state_ids=[...block.state_ids].sort((a,b)=>states.get(a).repeat_count-states.get(b).repeat_count);
    if(expandedSupernodes.has(sn))node(containerId(sn),'supernode-container',sn);
    else node(summaryId(sn),'supernode-summary',sn);
    for(const sid of block.state_ids){
      const state=states.get(sid);
      state.original_sequence_count=state.original.size;
      state.represented_sequence_count=state.represented.size;
      delete state.original;delete state.represented;
      if(!expandedSupernodes.has(sn))continue;
      const opened=expandedStates.has(sid);
      node(opened?stateContainerId(sn,state.repeat_count):sid,opened?'state-container':'state',sn,
           {parent:containerId(sn),state_id:sid});
    }
  }
  for(const p of members.keys()){
    const route=[];routes.set(p,route);
    let previous=null,step=0;
    for(const v of bundle.paths[p]){
      const sn=v.supernode_id,sid=stateId(sn,v.repeat_count),ids=[];
      let offset=0;
      if(!expandedSupernodes.has(sn))ids.push(summaryId(sn));
      else if(!expandedStates.has(sid)){
        ids.push(sid);
        offset=80*((ports.get(v.visit_id)+0.5)/total.get(sid)-0.5);
      }else{
        v.ordered_occurrence_ids.forEach((oid,i)=>{
          const id=copyId(v.visit_id,oid);
          ids.push(id);
          node(id,'copy',sn,{parent:stateContainerId(sn,v.repeat_count),state_id:sid,
            visit_id:v.visit_id,path_id:p,occurrence_id:oid,index:i,
            dna:bundle.variants[v[mode+'_variant_ids'][i]],source_node_id:v[mode+'_gfa_walk'][i],
            changed:v.original_variant_ids[i]!==v.represented_variant_ids[i]});
        });
      }
      route.push({visit_id:v.visit_id,display_node_ids:ids});
      for(const id of ids){
        if(previous){
          edges.push({id:displayId('trace',p,step,previous.id,id),path_id:p,step_index:step++,
            source:previous.id,target:id,from_visit_id:previous.visit,to_visit_id:v.visit_id,
            outlet:'50% '+previous.offset+'%',inlet:'-50% '+offset+'%'});
        }
        previous={id,visit:v.visit_id,offset};
      }
    }
  }
  return {nodes,edges,routes,states,blocks,support:{path_count:members.size,molecule_count:molecules.size}};
}

if(typeof module!=='undefined' && module.exports){
  module.exports={projectStateGraph,displayId,stateId,summaryId,containerId,stateContainerId,copyId};
}
