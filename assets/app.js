/* Integrated Pediatric AML NMF CpG molecular explorer v2. Static, no external calls. */
(()=>{'use strict';
const D=window.PAML_METHYLATION_ATLAS, get=id=>document.getElementById(id), NS='http://www.w3.org/2000/svg';
if(!D||D.schema!=='pAML-methylation-v2'){
  get('graphTitle').textContent='The atlas data could not be loaded';
  get('graphSubtitle').textContent='Check that data/atlas.js is uploaded beside index.html.';return;
}
const R={C:0,D:1,P:2,I:3,N:4,NP:5,NI:6,NN:7,QIN:8,QPI:9,QPN:10,QT:11,MASK:12,W:13,WSPEC:14,WSHARE:15,WSENS:16,SDP:17,SDI:18,SDN:19,PIN:20,PPI:21,PPN:22,PT:23,DIN:24,DPI:25,DPN:26,CONS:27};
const P={ID:0,CHR:1,START:2,END:3,STRAND:4,REGION:5,ISLAND:6,PROM:7,BODY:8,NEAREST:9,DIST:10};
const C={SOURCE:0,TARGET:1,POP:2,STATUS:3,NS:4,NB:5,NN:6,FEAT:7,LEVEL:8,TYPE:9,EFFECT:10,UNIT:11,PV:12,QW:13,QO:14,TEST:15,ASSOC:16,UNIQUE:17,RESTRICT:18,CONS:19,BOTH:20,ONLYR:21,ONLYC:22,NEITHER:23,HAS:24};
const groups=D.groups.map(g=>g[0]),counts=Object.fromEntries(D.groups),geneList=Object.keys(D.geneRows),geneIndex=new Map(geneList.map(g=>[g.toUpperCase(),g]));
const probeLookup=new Map(D.probes.map((p,i)=>[p[P.ID].toUpperCase(),i]));
const rowsByProbe=new Map();D.rows.forEach((r,i)=>{const p=r[R.C];if(!rowsByProbe.has(p))rowsByProbe.set(p,[]);rowsByProbe.get(p).push(i);});
const groupMeans=new Map(Object.entries(D.groupMeans).map(([group,vals])=>[group,new Map(vals.map(v=>[v[0],{beta:v[1],n:v[2]}]))]));
const clinical=new Map(D.clinical.map(r=>[`${r[C.SOURCE]}|${r[C.TARGET]}|${r[C.POP]}`,r]));
const flags=Object.fromEntries(D.flags.map((f,i)=>[f,i]));
const flagGroups={restricted:['Primary_Restricted'],conserved:['Neighborhood_Conserved'],restricted010:['Primary_Restricted_010'],conserved010:['Neighborhood_Conserved_010'],
 core:['Core_Target_06','Core_Source_06'],bridge:['Prior_Bridge_Source_to_Target_07a','Prior_Bridge_Target_to_Source_07a'],
 distinguishing:['Prior_06b_A_Higher_015','Prior_06b_B_Higher_015','Prior_06b_A_Defining_Concordant_Distinct_015','Prior_06b_B_Defining_Concordant_Distinct_015','Prior_06b_Target_Higher_015','Prior_06b_Source_Higher_015','Prior_06b_Target_Defining_Distinct_015','Prior_06b_Source_Defining_Distinct_015'],
 shared:['Prior_06b_Shared_Higher_015','Prior_06b_Shared_Lower_015','Prior_06b_Shared_Higher_Formal_015','Prior_06b_Shared_Lower_Formal_015'],
 '11c':['Prior_11C_Target_Restricted','Prior_11C_Target_Conserved']};
const colors={red:'#c84d5d',blue:'#3478bf',mixed:'#888999'};
const state={gene:'',directProbe:null,mode:'gene',cpg:null,region:'Promoter',island:'All',signature:'all',contrast:'IN',qMax:.05,effectMin:.15,wMin:0,source:'All',target:'All',selection:{type:'none',id:null}};
let view={ids:[],probes:[],edges:[],nodes:[]};
const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const clamp=(x,a,b)=>Math.min(b,Math.max(a,x));
const fmt=(x,d=3)=>x===null||x===undefined||!Number.isFinite(x)?'—':Number(x).toFixed(d);
const fq=x=>x===null||x===undefined||!Number.isFinite(x)?'—':x<.0001?'<0.0001':x.toFixed(4);
const pct=x=>x===null||!Number.isFinite(x)?'—':(x*100).toFixed(1)+'%';
const flag=(mask,key)=>Number.isInteger(flags[key])&&Boolean(mask & (1<<flags[key]));
const allFlags=mask=>D.flags.filter((name,i)=>(mask & (1<<i))!==0);
const shortFlags=mask=>{let list=[];if(flag(mask,'Primary_Restricted'))list.push('11D restricted');if(flag(mask,'Neighborhood_Conserved'))list.push('11D conserved');
 for(const [key,label] of [['core','06 core'],['bridge','07a bridge'],['distinguishing','06b distinguishing'],['shared','06b shared'],['11c','11C']])if(flagGroups[key].some(k=>flag(mask,k)))list.push(label);return list;};
const contrastText=()=>state.contrast==='IN'?'incoming bridge I − source nonbridge N':state.contrast==='IP'?'incoming bridge I − target primary P':'target primary P − source nonbridge N';
const effect=r=>state.contrast==='IN'?r[R.DIN]:state.contrast==='IP'?(r[R.I]===null||r[R.P]===null?null:r[R.I]-r[R.P]):r[R.DPN];
const adjusted=r=>r[state.contrast==='IN'?R.QIN:state.contrast==='IP'?R.QPI:R.QPN];
const qualifies=r=>Number.isFinite(effect(r))&&Number.isFinite(adjusted(r))&&adjusted(r)<=state.qMax+1e-11&&Math.abs(effect(r))+1e-11>=state.effectMin;
function colorBeta(b){if(!Number.isFinite(b))return '#e9edf2';const v=clamp(b,0,1),a=v<=.5?[165,203,236]:[248,249,251],z=v<=.5?[248,249,251]:[243,186,196],t=v<=.5?v*2:(v-.5)*2;
 return `rgb(${a.map((x,i)=>Math.round(x+(z[i]-x)*t)).join(',')})`;}
function svg(tag,attributes={},text=null){const a=document.createElementNS(NS,tag);for(const [k,v] of Object.entries(attributes))a.setAttribute(k,String(v));if(text!==null)a.textContent=text;return a;}
function select(type,id){state.selection={type,id};draw();}
function getRowsForGene(){return state.directProbe!==null?(rowsByProbe.get(state.directProbe)||[]):(D.geneRows[state.gene]||[]);}
function contextRows(){const ids=getRowsForGene();return ids.filter(i=>state.region==='All'||D.probes[D.rows[i][R.C]][P.REGION]===state.region);}
function activeRows(){return contextRows().filter(i=>{const r=D.rows[i],pr=D.probes[r[R.C]],d=D.directions[r[R.D]];
 return (state.mode!=='cpg'||r[R.C]===state.cpg)&&(state.island==='All'||pr[P.ISLAND]===state.island)&&
 (state.signature==='all'||flagGroups[state.signature].some(k=>flag(r[R.MASK],k)))&&
 (state.wMin===0||Number.isFinite(r[R.W])&&r[R.W]>=state.wMin/100-1e-10)&&
 (state.source==='All'||d[0]===state.source)&&(state.target==='All'||d[1]===state.target);});}
function computeView(){const ids=activeRows(),sProbes=new Set(),edges=new Map();for(const i of ids){const r=D.rows[i];sProbes.add(r[R.C]);if(!qualifies(r))continue;let di=r[R.D];if(!edges.has(di))edges.set(di,[]);edges.get(di).push(i);}
 const arrows=[];for(const [di,indices] of edges){const vals=indices.map(i=>effect(D.rows[i]));const pos=vals.filter(x=>x>0).length,neg=vals.filter(x=>x<0).length;
 arrows.push({di,indices,number:new Set(indices.map(i=>D.rows[i][R.C])).size,pos,neg,sign:pos&&neg?'mixed':pos?'red':neg?'blue':'mixed',avg:vals.reduce((a,b)=>a+b,0)/vals.length});}
 arrows.sort((a,b)=>b.number-a.number||a.di-b.di);
 const nodeData=groups.map(g=>{const gm=groupMeans.get(g)||new Map(),found=[...sProbes].map(id=>({id,entry:gm.get(id)})).filter(z=>z.entry);
 return {g,mean:found.length?found.reduce((t,z)=>t+z.entry.beta,0)/found.length:null,found,available:sProbes.size,count:counts[g]};});
 return {ids,probes:[...sProbes],edges:arrows,nodes:nodeData};}
function updateCpGSelect(){const ids=[...new Set(contextRows().map(i=>D.rows[i][R.C]))].sort((a,b)=>D.probes[a][P.ID].localeCompare(D.probes[b][P.ID]));
 const box=get('cpgSelect');box.replaceChildren();for(const id of ids){const pr=D.probes[id];const opt=document.createElement('option');opt.value=String(id);opt.textContent=`${pr[P.ID]} · chr${pr[P.CHR]}:${pr[P.START]}`;box.append(opt);}
 if(!ids.includes(state.cpg))state.cpg=ids.length?ids[0]:null;
 box.value=state.cpg===null?'':String(state.cpg);get('cpgSelectHolder').hidden=state.mode!=='cpg';}
function setGene(g){if(!Object.prototype.hasOwnProperty.call(D.geneRows,g))return;state.gene=g;state.directProbe=null;state.mode='gene';get('viewMode').value='gene';state.cpg=null;state.selection={type:'none',id:null};get('geneSearch').value=g;get('probeSearch').value='';get('probeStatus').textContent='';get('suggestions').hidden=true;updateCpGSelect();draw();}
function jumpProbe(){const raw=get('probeSearch').value.trim().toUpperCase(),id=probeLookup.get(raw);if(id===undefined){get('probeStatus').textContent='Exact CpG ID not found in the 11D selected suite.';return;}
 const pr=D.probes[id];state.directProbe=id;state.gene='';state.cpg=id;state.mode='cpg';get('viewMode').value='cpg';state.region=pr[P.REGION];get('region').value=state.region;
 get('geneSearch').value='';get('probeStatus').textContent=`${pr[P.ID]} · ${pr[P.REGION]} · chr${pr[P.CHR]}:${pr[P.START]}`;state.selection={type:'none',id:null};updateCpGSelect();draw();}
function groupPosition(g){const idx=groups.indexOf(g),a=-Math.PI/2+idx*2*Math.PI/groups.length;return{x:500+365*Math.cos(a),y:393+304*Math.sin(a)};}
function pathFor(s,t){let a=groupPosition(s),b=groupPosition(t),dx=b.x-a.x,dy=b.y-a.y,L=Math.hypot(dx,dy)||1,ux=dx/L,uy=dy/L,nx=-uy,ny=ux,
 start={x:a.x+ux*34,y:a.y+uy*34},end={x:b.x-ux*42,y:b.y-uy*42};let bend=Math.min(90,20+L*.085),cx=(start.x+end.x)/2+nx*bend,cy=(start.y+end.y)/2+ny*bend;
 return `M${start.x.toFixed(2)} ${start.y.toFixed(2)} Q${cx.toFixed(2)} ${cy.toFixed(2)} ${end.x.toFixed(2)} ${end.y.toFixed(2)}`;}
function tip(html,e){const b=get('tooltip');b.innerHTML=html;b.hidden=false;const rect=b.getBoundingClientRect();b.style.left=clamp(e.clientX+14,8,window.innerWidth-rect.width-12)+'px';b.style.top=clamp(e.clientY+10,8,window.innerHeight-rect.height-12)+'px';}
function drawNetwork(){const edgeLayer=get('edgeLayer'),nodeLayer=get('nodeLayer');edgeLayer.replaceChildren();nodeLayer.replaceChildren();
 for(const ed of view.edges){const [s,t,ns,nb]=D.directions[ed.di],active=state.selection.type==='edge'&&state.selection.id===ed.di;
 const g=svg('g',{class:'edge'+(active?' active':'')});const d=pathFor(s,t),width=2+Math.min(4,Math.log2(1+ed.number)*.9);
 g.append(svg('path',{d,stroke:colors[ed.sign],'stroke-width':width,'stroke-dasharray':ed.sign==='mixed'?'6 5':'none','marker-end':`url(#arrow-${ed.sign})`,opacity:active?1:.77}));
 const hit=svg('path',{d,stroke:'transparent','stroke-width':18,fill:'none',role:'button',tabindex:0,'aria-label':`${s} to ${t}: ${ed.number} qualifying CpGs`});
 hit.addEventListener('pointermove',e=>tip(`<strong>${esc(s)} → ${esc(t)}</strong><div>${ed.number} qualifying CpG(s); ${nb} pooled source bridge patients</div><div>Mean qualifying-site Δβ: ${fmt(ed.avg)}</div><div class="dim">${esc(contrastText())}</div>`,e));
 hit.addEventListener('pointerleave',()=>get('tooltip').hidden=true);hit.addEventListener('click',()=>select('edge',ed.di));hit.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();select('edge',ed.di);}});
 g.append(hit);edgeLayer.append(g);}
 for(const nd of view.nodes){let {x,y}=groupPosition(nd.g),active=state.selection.type==='node'&&state.selection.id===nd.g;
 const g=svg('g',{class:'node'+(nd.mean===null?' no-data':'')+(active?' active':''),role:'button',tabindex:0,'aria-label':`${nd.g}: ${nd.mean===null?'no mean available':'mean beta '+fmt(nd.mean)}`});
 g.append(svg('circle',{cx:x,cy:y,r:32,fill:colorBeta(nd.mean),stroke:active?'#1f4465':nd.mean===null?'#c9d4e0':'#879eb5','stroke-width':active?4:1.5}));
 g.append(svg('text',{x,y:y+1,'text-anchor':'middle'},nd.g));g.append(svg('text',{x,y:y+14,'text-anchor':'middle',class:'n-label'},nd.mean===null?'no data':'β '+fmt(nd.mean,2)));
 g.addEventListener('click',()=>select('node',nd.g));g.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();select('node',nd.g);}});
 g.addEventListener('pointermove',e=>tip(`<strong>${esc(nd.g)}</strong><div>Original primary mean β: ${fmt(nd.mean)}</div><div class="dim">${nd.found.length}/${nd.available} selected CpGs with means</div>`,e));g.addEventListener('pointerleave',()=>get('tooltip').hidden=true);nodeLayer.append(g);
 }
 if(view.edges.length===0){const text=svg('text',{x:500,y:392,'text-anchor':'middle',class:'selection-tip'},'No CpG arrows pass the active filters');edgeLayer.append(text);}
}
function metric(title,value,note=''){return `<div class="metric"><strong>${esc(title)}</strong><span class="value">${esc(value)}</span>${note?`<span class="micro">${esc(note)}</span>`:''}</div>`;}
function smallMetric(title,value){return `<div class="metric small"><strong>${esc(title)}</strong><span class="value">${esc(value)}</span></div>`;}
function displayGene(){return state.directProbe!==null?D.probes[state.directProbe][P.ID]:state.gene;}
function prettyLocus(pr){return `chr${pr[P.CHR]}:${pr[P.START]}–${pr[P.END]} (${pr[P.STRAND]??'?'} strand)`;}
function clinicalHtml(source,target){let html='<div class="clinical-list"><h3>Clinical bridge context (13C)</h3><p class="help">All / Female / Male summaries refer to clinical tests. CpG β and q-values above are pooled only.</p>';
 const all=clinical.get(`${source}|${target}|All`);
 if(all){html+=`<div class="note"><b>11D suite counts for this direction:</b> ${all[C.ASSOC]??0} association(s); ${all[C.UNIQUE]??0} unique CpGs; ${all[C.RESTRICT]??0} restricted; ${all[C.CONS]??0} conserved; ${all[C.BOTH]??0} both. These are pooled 11D counts, not clinical significance tests or sex-stratified CpG counts.</div>`;}
 for(const pop of ['All','Female','Male']){const v=clinical.get(`${source}|${target}|${pop}`);
 if(!v){html+=`<div class="clinical-card"><h4>${pop}</h4><div class="dimmed">No 13C pair record available.</div></div>`;continue;}
 const eligible=v[C.STATUS]==='eligible'&&v[C.NB]!==null&&v[C.NN]!==null;
 html+=`<div class="clinical-card"><h4>${esc(pop)} · ${eligible?'eligible':'not testable / redacted'}</h4>`;
 if(eligible){html+=`<div class="clinical-line">Source n=${v[C.NS]}, bridge n=${v[C.NB]}, nonbridge n=${v[C.NN]}</div>`;
 if(v[C.FEAT]!==null){html+=`<div class="clinical-line"><b>Top reported characteristic:</b> ${esc(v[C.FEAT])}${v[C.LEVEL]?' · '+esc(v[C.LEVEL]):''}</div>`;
 html+=`<div class="clinical-line">Effect (bridge − nonbridge): ${fmt(v[C.EFFECT])} ${esc(v[C.UNIT]??'')} · within-population BH q=${fq(v[C.QW])} · overall BH q=${fq(v[C.QO])}</div>`;
 }else html+='<div class="dimmed">13C has no ranked feature for this eligible direction.</div>';
 }else html+='<div class="dimmed">13C did not support subgroup testing; potentially small counts and top clinical results are withheld in the public export.</div>';
 html+='</div>';}
 html+='<div class="note">The “top characteristic” is the one retained by the 13C summary table, <b>not the complete set</b> of tested clinical features. q-values are from clinical tests, not the CpG contrast, and the selected feature may not pass q≤0.05.</div></div>';
 return html;}
function chosenEdgeHtml(di, pairMode=false){const d=D.directions[di],source=d?.[0],target=d?.[1];if(!d)return '<p class="help">No CpG direction record.</p>';const ed=view.edges.find(x=>x.di===di);
 let h=`<h3>${esc(source)} → ${esc(target)}</h3><p class="detail-subtitle">SOURCE original-primary patients evaluated for TARGET H-neighborhood inclusion</p>`;
 h+=metric('Qualifying CpGs for selected gene',ed?.number??0,contrastText());
 if(ed){h+=metric('Mean Δβ of qualifying sites',fmt(ed.avg),`${ed.pos} higher, ${ed.neg} lower`);
 h+=`<div class="note">${ed.sign==='mixed'?'Mixed signs: no single gene-level hyper/hypomethylation classification.':`All qualifying CpGs have ${ed.sign==='red'?'positive':'negative'} effects.`} Gene-level mean is descriptive and has no combined p-value.</div>`;
 h+='<h3>Highest |Δβ| sites</h3>'+ed.indices.slice().sort((a,b)=>Math.abs(effect(D.rows[b]))-Math.abs(effect(D.rows[a]))).slice(0,12).map(i=>{
 const r=D.rows[i],p=D.probes[r[R.C]];return `<div class="metric"><a href="#" class="clickable" data-site="${i}">${esc(p[P.ID])}</a> · ${esc(p[P.REGION])}<div>Δβ ${fmt(effect(r))} · BH q ${fq(adjusted(r))}</div><div class="micro">P/I/N β: ${fmt(r[R.P])} / ${fmt(r[R.I])} / ${fmt(r[R.N])}</div></div>`;}).join('');
 }else h+='<div class="note">No CpGs pass the current gene, region, and arrow filters for this direction. The clinical summary below remains available.</div>';
 return h+clinicalHtml(source,target);}
function rowDetail(i){const r=D.rows[i];if(!r)return '<p>Unknown CpG association.</p>';const pr=D.probes[r[R.C]],d=D.directions[r[R.D]];
 let h=`<h3>${esc(pr[P.ID])}</h3><p class="detail-subtitle">${esc(d[0])} → ${esc(d[1])} · ${esc(pr[P.REGION])} · ${esc(pr[P.ISLAND]??'unclassified')}</p>`;
 h+=`<div class="note"><b>${esc(prettyLocus(pr))}</b><div>Promoter overlap: ${esc(pr[P.PROM]||'—')}</div><div>Gene-body overlap: ${esc(pr[P.BODY]||'—')}</div><div>Nearest gene: ${esc(pr[P.NEAREST]||'—')} (${pr[P.DIST]===null?'unknown':pr[P.DIST]} bp)</div><div>Nearest does not imply regulation.</div></div>`;
 h+='<div class="detail-grid">'+[
 ['Mean β P',fmt(r[R.P],4)],['Mean β I',fmt(r[R.I],4)],['Mean β N',fmt(r[R.N],4)],
 ['Observed n P / I / N',`${r[R.NP]} / ${r[R.NI]} / ${r[R.NN]}`],
 ['SD P / I / N',`${fmt(r[R.SDP])} / ${fmt(r[R.SDI])} / ${fmt(r[R.SDN])}`],
 ['Δβ I − N',fmt(r[R.DIN])],['Δβ I − P',fmt(-r[R.DPI])],['Δβ P − N',fmt(r[R.DPN])],
 ['BH q I vs N',fq(r[R.QIN])],['BH q P vs I',fq(r[R.QPI])],['BH q P vs N',fq(r[R.QPN])],['BH q equivalence P~I',fq(r[R.QT])],
 ['Raw Welch p I vs N',fq(r[R.PIN])],['Raw Welch p P vs I',fq(r[R.PPI])],['Raw Welch p P vs N',fq(r[R.PPN])],['Raw TOST p P~I',fq(r[R.PT])],
 ['W target percentile',pct(r[R.W])],['W target specificity',fmt(r[R.WSPEC],4)],['W column share',fmt(r[R.WSHARE],6)],['W sensitivity weight',fmt(r[R.WSENS],4)]
 ].map(k=>smallMetric(...k)).join('')+'</div>';
 h+=`<h3>All preserved suite labels</h3><div>${allFlags(r[R.MASK]).map(k=>`<span class="bubble">${esc(k)}</span>`).join('')||'<span class="help">No prior flags retained</span>'}</div>`;
 h+=`<p class="help">Conserved-direction annotation: ${esc(r[R.CONS]||'not specified')}. BH tests were corrected separately by direction/test family; TOST q is not a differential-methylation q.</p>`;
 return h+clinicalHtml(d[0],d[1]);}
function showDetails(){const root=get('details');let h='';const s=state.selection;
 if(s.type==='edge'){h=chosenEdgeHtml(s.id);}else if(s.type==='pair'){
 const [source,target]=s.id.split('|');let di=D.directions.findIndex(d=>d[0]===source&&d[1]===target);
 h=di>=0?chosenEdgeHtml(di,true):`<h3>${esc(source)} → ${esc(target)}</h3><p class="help">No selected CpG association for this direction in 11D.</p>`+clinicalHtml(source,target);
 }else if(s.type==='site'){h=rowDetail(s.id);}else if(s.type==='node'){
 const node=view.nodes.find(v=>v.g===s.id);h=`<h3>${esc(node.g)}</h3><p class="detail-subtitle">Original NMF primary-group methylation for selected sites</p>`+
 metric('Mean β across represented unique CpGs',fmt(node.mean),'Not relative to controls or other groups')+
 metric('CpG coverage',`${node.found.length} / ${node.available}`,'CpGs represented in this original group / selected CpGs')+
 metric('Original group size',node.count===null?'Not available':String(node.count),'No data fabricated for ineligible groups');
 h+='<h3>Available CpG means</h3>'+node.found.slice(0,35).map(v=>`<div class="metric"><strong>${esc(D.probes[v.id][P.ID])}</strong><span class="value">β ${fmt(v.entry.beta)}</span><div class="micro">n=${v.entry.n} CpG-wise observed primary patients</div></div>`).join('');
 h+='<div class="note">Each unique CpG contributes once per original group. Means are reported for target primary patients or reconstructed from source bridge and nonbridge counts. Unequal CpG coverage can affect comparisons.</div>';
 }else{
 const grouped={Promoter:0,'Gene body':0,Intergenic:0};for(const p of new Set(getRowsForGene().map(i=>D.rows[i][R.C])))grouped[D.probes[p][P.REGION]]++;
 h=`<h3>${esc(displayGene()||'No gene selected')}</h3><p class="detail-subtitle">${esc(state.region)} · ${state.mode==='gene'?'gene-level descriptive view':'single CpG'}</p>`+
 metric('Unique CpGs after filters',String(view.probes.length))+
 metric('Arrows passing per-CpG tests',String(view.edges.length),contrastText())+
 metric('Original groups with mean β',`${view.nodes.filter(x=>x.mean!==null).length}/23`)+
 `<div class="note"><b>Coverage before filters:</b> ${grouped.Promoter} promoter CpGs, ${grouped['Gene body']} gene-body CpGs, ${grouped.Intergenic} intergenic CpGs (assigned by nearest gene).</div>`+
 '<div class="note">A CpG can carry several 06/06b/07a/11C/11D flags. Red/blue circles show absolute primary-group β; red/blue arrows show the selected signed group contrast. Neither demonstrates a time transition.</div>';
 if(view.edges.length){h+='<h3>Qualifying directions</h3><div class="dir-list">'+view.edges.slice(0,30).map(ed=>{const d=D.directions[ed.di];return `<button data-direction="${ed.di}"><span style="color:${colors[ed.sign]}">●</span><span>${esc(d[0])} → ${esc(d[1])}</span><b>${ed.number} CpGs</b></button>`;}).join('')+'</div>';}
 }
 root.innerHTML=h;for(const b of root.querySelectorAll('[data-direction]'))b.addEventListener('click',()=>select('edge',Number(b.dataset.direction)));
 for(const b of root.querySelectorAll('[data-site]'))b.addEventListener('click',e=>{e.preventDefault();select('site',Number(b.dataset.site));});
}
function renderEvidence(){let ids=view.ids;
 if(state.selection.type==='edge')ids=ids.filter(i=>D.rows[i][R.D]===state.selection.id);
 const ranked=ids.slice().sort((a,b)=>{let aQ=adjusted(D.rows[a]),bQ=adjusted(D.rows[b]);return (aQ??Infinity)-(bQ??Infinity)||Math.abs(effect(D.rows[b])??0)-Math.abs(effect(D.rows[a])??0);});
 get('tableCount').textContent=`${ranked.length.toLocaleString()} association row(s); showing ${Math.min(250,ranked.length)}`;
 const body=get('evidenceBody');body.replaceChildren();if(!ranked.length){let tr=document.createElement('tr');tr.innerHTML='<td colspan="12" class="empty-row">No selected CpG associations match these filters. Try All regions or a different signature label.</td>';body.append(tr);return;}
 for(const i of ranked.slice(0,250)){const r=D.rows[i],p=D.probes[r[R.C]],d=D.directions[r[R.D]],v=effect(r);
 const tr=document.createElement('tr');tr.className='evidence-row'+(state.selection.type==='site'&&state.selection.id===i?' selected':'');tr.tabIndex=0;
 tr.innerHTML=`<td>${esc(d[0])} → ${esc(d[1])}</td><td class="clickable">${esc(p[P.ID])}</td><td><span class="ctx-label">${esc(p[P.REGION])}</span></td><td>chr${esc(p[P.CHR])}:${p[P.START]}</td><td>${esc(p[P.ISLAND])}</td><td>${fmt(r[R.P])}</td><td>${fmt(r[R.I])}</td><td>${fmt(r[R.N])}</td><td class="${v>0?'effect-up':v<0?'effect-down':'effect-flat'}">${fmt(v)}</td><td>${fq(adjusted(r))}</td><td>${pct(r[R.W])}</td><td class="labels-cell">${shortFlags(r[R.MASK]).map(z=>`<span class="small-label">${esc(z)}</span>`).join('')||'—'}</td>`;
 tr.addEventListener('click',()=>select('site',i));tr.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();select('site',i);}});body.append(tr);}
}
function redrawText(){const raw=getRowsForGene(),allCpGs=new Set(raw.map(i=>D.rows[i][R.C]));const name=displayGene();
 get('geneStatus').textContent=raw.length?`${allCpGs.size.toLocaleString()} CpGs linked by genomic position; ${view.probes.length.toLocaleString()} selected after filters`:'No annotated CpGs in this suite for this gene.';
 get('graphTitle').textContent=name+(state.mode==='cpg'&&state.cpg!==null&&state.directProbe===null?' · '+D.probes[state.cpg][P.ID]:'');
 get('graphSubtitle').textContent=`${view.probes.length} unique CpGs · ${view.edges.length} qualifying arrow(s) · ${contrastText()} · ${state.region}`;
 get('edgeSummary').textContent=`${view.edges.filter(x=>x.sign==='red').length} higher · ${view.edges.filter(x=>x.sign==='blue').length} lower · ${view.edges.filter(x=>x.sign==='mixed').length} mixed`;
}
function draw(){view=computeView();if(state.selection.type==='edge'&&!view.edges.some(e=>e.di===state.selection.id))state.selection={type:'none',id:null};
 if(state.selection.type==='site'&&!view.ids.includes(state.selection.id))state.selection={type:'none',id:null};
 redrawText();drawNetwork();renderEvidence();showDetails();}
function bind(id,event,fn){get(id).addEventListener(event,e=>{fn(e.target.value);state.selection={type:'none',id:null};draw();});}
function setupSearch(){const input=get('geneSearch'),suggest=get('suggestions');
 const refresh=()=>{const value=input.value.trim().toUpperCase();suggest.replaceChildren();if(!value){suggest.hidden=true;return;}
 const starts=geneList.filter(g=>g.toUpperCase().startsWith(value)).slice(0,16),rest=geneList.filter(g=>!g.toUpperCase().startsWith(value)&&g.toUpperCase().includes(value)).slice(0,8);
 for(const g of starts.concat(rest)){const b=document.createElement('button');b.type='button';b.textContent=g;b.addEventListener('click',()=>setGene(g));suggest.append(b);}suggest.hidden=!suggest.children.length;};
 input.addEventListener('input',refresh);input.addEventListener('focus',()=>{if(input.value!==state.gene)refresh();});
 input.addEventListener('keydown',e=>{if(e.key==='Escape'){suggest.hidden=true;input.blur();}if(e.key==='Enter'){e.preventDefault();const option=suggest.querySelector('button');if(option)option.click();else{const actual=geneIndex.get(input.value.trim().toUpperCase());if(actual)setGene(actual);}}});
 document.addEventListener('pointerdown',e=>{if(!suggest.contains(e.target)&&e.target!==input)suggest.hidden=true;});
 document.querySelectorAll('[data-try]').forEach(b=>b.addEventListener('click',()=>{const gene=geneIndex.get(b.dataset.try.toUpperCase());if(gene)setGene(gene);}));
 get('goProbe').addEventListener('click',jumpProbe);get('probeSearch').addEventListener('keydown',e=>{if(e.key==='Enter')jumpProbe();});}
function setupFilters(){setupSearch();
 bind('region','change',v=>{state.region=v;updateCpGSelect();});
 bind('island','change',v=>state.island=v);
 bind('viewMode','change',v=>{state.mode=v;updateCpGSelect();});
 bind('cpgSelect','change',v=>state.cpg=v===''?null:Number(v));
 bind('signature','change',v=>state.signature=v);bind('contrast','change',v=>state.contrast=v);
 bind('qMax','change',v=>state.qMax=Number(v));
 bind('effectMin','input',v=>{state.effectMin=Number(v);get('effectMinLabel').textContent=fmt(state.effectMin);});
 bind('wMin','input',v=>{state.wMin=Number(v);get('wMinLabel').textContent=v+'%';});
 for(const id of ['sourceFilter','targetFilter']){let select=get(id);for(const g of groups){const o=document.createElement('option');o.value=g;o.textContent=g;select.append(o);}}
 bind('sourceFilter','change',v=>state.source=v);bind('targetFilter','change',v=>state.target=v);
 const allPairs=[...new Set(D.clinical.filter(v=>v[C.POP]==='All').map(v=>`${v[C.SOURCE]}|${v[C.TARGET]}`))].sort((a,b)=>{
 let x=a.split('|'),y=b.split('|');return Number(x[0].slice(3))-Number(y[0].slice(3))||Number(x[1].slice(3))-Number(y[1].slice(3));});
 for(const pair of allPairs){const o=document.createElement('option');o.value=pair;o.textContent=pair.replace('|',' → ');get('pairSelect').append(o);}
 get('inspectPair').addEventListener('click',()=>select('pair',get('pairSelect').value));
 get('exportCsv').addEventListener('click',exportCpGs);get('exportClinicalCsv').addEventListener('click',exportClinical);
}
function saveCsv(name,headers,lines){const x=[headers,...lines].map(row=>row.map(v=>'"'+String(v??'').replaceAll('"','""')+'"').join(',')).join('\r\n');
 const blob=new Blob(['\ufeff',x],{type:'text/csv;charset=utf-8;'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function exportCpGs(){let ids=view.ids;if(state.selection.type==='edge')ids=ids.filter(i=>D.rows[i][R.D]===state.selection.id);
 const head=['Gene_selected','CpG','Source','Target','Region','Island_context','Chromosome','Start_hg38','End_hg38','Strand','Promoter_overlaps','Gene_body_overlaps','Nearest_gene','Nearest_distance_bp',
 'Mean_beta_target_P','Mean_beta_bridge_I','Mean_beta_nonbridge_N','SD_P','SD_I','SD_N','Nobs_P','Nobs_I','Nobs_N',
 'Delta_I_minus_N','Delta_I_minus_P','Delta_P_minus_N','BH_q_I_vs_N','BH_q_P_vs_I','BH_q_P_vs_N','BH_q_TOST_P_equiv_I','P_Welch_I_vs_N','P_Welch_P_vs_I','P_Welch_P_vs_N','P_TOST',
 'W_target_percentile','W_target_specificity','W_target_column_share','W_sensitivity_weight','Conserved_direction','All_22_Provenance_Labels'];
 const lines=ids.map(i=>{const r=D.rows[i],p=D.probes[r[R.C]],d=D.directions[r[R.D]];return [displayGene(),p[P.ID],d[0],d[1],p[P.REGION],p[P.ISLAND],p[P.CHR],p[P.START],p[P.END],p[P.STRAND],p[P.PROM],p[P.BODY],p[P.NEAREST],p[P.DIST],
 r[R.P],r[R.I],r[R.N],r[R.SDP],r[R.SDI],r[R.SDN],r[R.NP],r[R.NI],r[R.NN],r[R.DIN],r[R.I]===null||r[R.P]===null?null:r[R.I]-r[R.P],r[R.DPN],r[R.QIN],r[R.QPI],r[R.QPN],r[R.QT],r[R.PIN],r[R.PPI],r[R.PPN],r[R.PT],r[R.W],r[R.WSPEC],r[R.WSHARE],r[R.WSENS],r[R.CONS],allFlags(r[R.MASK]).join(';')];});
 saveCsv(`pAML_CpG_${displayGene().replace(/[^\w.-]/g,'_')}_${state.region.replace(/\s+/g,'_')}.csv`,head,lines);}
function exportClinical(){
 const headers=['Source','Target','Population','13C_status','N_Source_if_eligible','N_Bridge_if_eligible','N_Nonbridge_if_eligible',
 'Top_Feature_if_eligible','Top_Level','Top_Feature_Type','Top_Effect_Bridge_Minus_NonBridge','Effect_Unit','Top_P','Top_FDR_Within_Population','Top_FDR_Overall','Top_Primary_Test',
 'N_11D_Directional_Associations','N_11D_Unique_CpGs','N_PrimaryRestricted_CpGs','N_NeighborhoodConserved_CpGs','N_Both','N_PrimaryOnly','N_ConservedOnly','N_Neither','Has_11D_CpGs'];
 saveCsv('pAML_13C_Clinical_Directional_Abbreviated_Public.csv',headers,D.clinical);
}
function boot(){get('datasetStats').innerHTML=`<span class="stat"><b>${D.meta.cohortPatients}</b> original patients</span><span class="stat"><b>${D.meta.uniqueCpGs.toLocaleString()}</b> unique CpGs</span><span class="stat"><b>${D.meta.associations.toLocaleString()}</b> CpG × direction rows</span><span class="stat"><b>${D.meta.strictPromoterCpGs.toLocaleString()}</b> promoter CpGs</span><span class="stat"><b>${D.meta.geneLabels.toLocaleString()}</b> genomic gene labels</span><span class="stat"><b>${D.meta.eligibleDirectionsWithCpGs}</b> directions with selected CpGs</span>`;
 setupFilters();const g=geneIndex.get('GATA3')||geneList[0];setGene(g);}
boot();
})();
