const assert = require('node:assert/strict');
const {correspondingVisits, layoutVisits, selectionSupport} = require('../scripts/supernode_compare');
const visits = (names,p) => names.map((name,i)=>({
  visit_id:p+':v'+i, supernode_id:name, structure_assessment:'invariant_within_model',
  boundary_assessment:'invariant_within_model', position_assessment:'aligned'
}));
const left=visits(['L','R'],'A'), right=visits(['L','A','B','C','D','R'],'B');
const pairs=correspondingVisits(left,right,false);
assert.deepEqual(pairs,[[0,0],[1,5]]);
const cols=layoutVisits(2,6,pairs,220);
for(const row of cols) assert.ok(row.every((x,i)=>i===0 || x-row[i-1]>=220));
assert.equal(cols[0][0],cols[1][0]);assert.equal(cols[0][1],cols[1][5]);
assert.deepEqual(correspondingVisits(visits([1,2,3,4],'A'),visits([1,3,2,4],'B'),false),[[0,0],[3,3]]);
assert.deepEqual(correspondingVisits(visits(['X','X'],'A'),visits(['X'],'B'),false),[]);
assert.deepEqual(correspondingVisits(visits(['alpha'],'A'),visits(['beta'],'B'),false),[]);
for(const [field,value] of [['boundary_assessment','ambiguous'],['structure_assessment','not_assessed'],['position_assessment','path_local']]){
  const x=visits(['X'],'A'),y=visits(['X'],'B');x[0][field]=value;
  assert.deepEqual(correspondingVisits(x,y,false),[]);
  assert.deepEqual(correspondingVisits(x,x,true),[[0,0]]);
}
assert.deepEqual(layoutVisits(0,2,[],200),[[],[200,400]]);
const samples={p:{s:['m1','m2']},q:{t:['m1','m3']},ref:{}};
const info={p:{role:'observation',molecule_keys:['m1','m2']},q:{role:'observation',molecule_keys:['m1','m3']},ref:{role:'reference',molecule_keys:[]}};
const sel=(a,b)=>selectionSupport([a,b],info,samples);
assert.equal(sel({path:'p',sample:'s'},{path:'q',sample:'t'}).molecules,3);
assert.equal(sel({path:'p',sample:'s'},{path:'p',sample:'s'}).paths,1);
assert.deepEqual(sel({path:'ref',sample:'s'},{path:'p',sample:'s'}).perTrack,[[],['m1','m2']]);
console.log('PASS correspondence, crossed/repeated/uncertain IDs, ordered layout and scoped support');
