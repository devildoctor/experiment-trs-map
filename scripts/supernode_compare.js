/* Pair projection uses registered visit order and the backend's position IDs. */
function correspondingVisits(left, right, samePath) {
  function positions(visits, key) {
    const result = new Map();
    visits.forEach((visit, i) => {
      const id = visit[key];
      result.set(id, result.has(id) ? -1 : i);
    });
    return result;
  }
  const key = samePath ? 'visit_id' : 'supernode_id';
  const a = positions(left, key), b = positions(right, key);
  const assessed = v => v.structure_assessment === 'invariant_within_model' &&
    v.boundary_assessment === 'invariant_within_model' && v.position_assessment === 'aligned';
  const candidates = [];
  left.forEach((v, i) => {
    const j = b.get(v[key]);
    if (a.get(v[key]) === i && j >= 0 && (samePath || (assessed(v) && assessed(right[j])))) {
      candidates.push([i,j]);
    }
  });
  // An inversion invalidates both endpoints, even when one is in a longer increasing subsequence.
  let minAfter = Infinity, maxBefore = -1;
  const safe = Array(candidates.length);
  for (let k = candidates.length - 1; k >= 0; k--) {
    const j = candidates[k][1];
    safe[k] = j < minAfter;
    minAfter = Math.min(minAfter, j);
  }
  return candidates.filter(([,j], k) => {
    const keep = safe[k] && j > maxBefore;
    maxBefore = Math.max(maxBefore, j);
    return keep;
  });
}

function layoutVisits(lengthA, lengthB, pairs, spacing) {
  const columns = [[],[]];
  let previous = [-1,-1], x = 0;
  for (const next of [...pairs,[lengthA,lengthB]]) {
    const width = Math.max(next[0]-previous[0], next[1]-previous[1])*spacing;
    for (const side of [0,1]) {
      const steps = next[side]-previous[side];
      for (let i = previous[side]+1; i < next[side]; i++) {
        columns[side][i] = x + width*(i-previous[side])/steps;
      }
      if (next[side] < [lengthA,lengthB][side]) columns[side][next[side]] = x + width;
    }
    x += width;
    previous = next;
  }
  return columns;
}

function selectionSupport(selections, info, samples) {
  const sets = selections.map(({path,sample}) => {
    if (!path || info[path].role === 'reference') return [];
    return sample ? (samples[path][sample] || []) : info[path].molecule_keys;
  });
  return {perTrack:sets, molecules:new Set(sets.flat()).size,
          paths:new Set(selections.map(s => s.path).filter(Boolean)).size};
}

function rawDifference(a,b) {
  if(a.length===b.length) {
    let mismatches=0;
    for(let i=0;i<a.length;i++) if(a[i]!==b[i]) mismatches++;
    return {lengthDifference:0,mismatches};
  }
  const [longer,shorter]=a.length>b.length?[a,b]:[b,a];
  let prefix=0,suffix=0;
  while(prefix<shorter.length && longer[prefix]===shorter[prefix]) prefix++;
  while(suffix<shorter.length && longer[longer.length-1-suffix]===shorter[shorter.length-1-suffix]) suffix++;
  const first=shorter.length-suffix;
  return {lengthDifference:longer.length-shorter.length,first,last:prefix,starts:Math.max(0,prefix-first+1)};
}

function untraversedJunctions(edges,paths) {
  const triples=new Set(),outgoing=new Map(),results=[];
  for(const visits of Object.values(paths))for(let i=2;i<visits.length;i++)
    triples.add(JSON.stringify([visits[i-2].supernode_id,visits[i-1].supernode_id,visits[i].supernode_id]));
  for(const edge of edges)if(edge.molecule_count>0){
    if(!outgoing.has(edge.source))outgoing.set(edge.source,[]);
    outgoing.get(edge.source).push(edge);
  }
  for(const first of edges)if(first.molecule_count>0)for(const second of outgoing.get(first.target)||[]){
    const route=[first.source,first.target,second.target];
    if(!triples.has(JSON.stringify(route)))results.push({route,left:first.molecule_count,right:second.molecule_count});
  }
  return results.sort((a,b)=>Math.min(b.left,b.right)-Math.min(a.left,a.right));
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = {correspondingVisits,layoutVisits,selectionSupport,rawDifference,untraversedJunctions};
}
