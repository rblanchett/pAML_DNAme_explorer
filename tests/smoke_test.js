// Minimal offline runtime smoke test: no browser, Node 18+ required.
// Usage: node tests/smoke_test.js
const fs=require('fs'),vm=require('vm'),path=require('path');
const root=path.join(__dirname,'..');
class Element{
 constructor(tag='div'){this.tag=tag;this.children=[];this.attributes={};this.events={};this.hidden=false;this.value='';this.style={};this.dataset={};this.textContent='';this.innerHTML='';}
 setAttribute(k,v){this.attributes[k]=v;}append(...x){this.children.push(...x);}replaceChildren(...x){this.children=x;}addEventListener(type,fn){this.events[type]=fn;}querySelectorAll(){return [];}querySelector(){return this.children.find(c=>c.tag==='button')||null;}click(){this.events.click?.();}getBoundingClientRect(){return {width:100,height:100};}contains(other){return this.children.includes(other);}
}
const elems=new Map(),doc={getElementById:id=>{if(!elems.has(id))elems.set(id,new Element());return elems.get(id);},createElement:tag=>new Element(tag),createElementNS:(_,tag)=>new Element(tag),querySelectorAll:()=>[],addEventListener:()=>{}};
const sandbox={window:{innerWidth:1600,innerHeight:900},document:doc,Blob,URL,Math,Number,Map,Set,Object,String,Array,console};
vm.runInNewContext(fs.readFileSync(path.join(root,'data','atlas.js'),'utf8'),sandbox,{filename:'atlas.js'});
vm.runInNewContext(fs.readFileSync(path.join(root,'assets','app.js'),'utf8'),sandbox,{filename:'app.js'});
function assert(ok,msg){if(!ok)throw new Error(msg);console.log('PASS '+msg)}
assert(elems.get('graphTitle').textContent==='GATA3','Initial gene is GATA3');
assert(elems.get('nodeLayer').children.length===23,'All 23 original NMF nodes are drawn');
assert(elems.get('edgeLayer').children.length===3,'Expected three qualifying GATA3 arrows');
assert(elems.get('evidenceBody').children.length===13,'Underlying GATA3 evidence rows appear');
assert(elems.get('details').innerHTML.includes('GATA3'),'Gene-level details populate');
elems.get('contrast').events.change({target:{value:'IP'}});
assert(elems.get('nodeLayer').children.length===23,'Switching contrast updates network without exception');
elems.get('viewMode').events.change({target:{value:'cpg'}});
assert(elems.get('graphTitle').textContent.startsWith('GATA3 · cg'),'CpG mode updates selected probe');
elems.get('viewMode').events.change({target:{value:'gene'}});
elems.get('signature').events.change({target:{value:'conserved'}});
assert(elems.get('graphTitle').textContent==='GATA3','Signature filtering retains gene selection');
elems.get('signature').events.change({target:{value:'all'}});
console.log('SMOKE TEST COMPLETE');
