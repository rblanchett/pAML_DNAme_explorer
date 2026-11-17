/* Pediatric AML Promoter CpG Network | standalone/static, no requests or dependencies */
(()=>{'use strict';
const D=window.PAML_PROMOTER_ATLAS;
const el=id=>document.getElementById(id), NS='http://www.w3.org/2000/svg';
if(!D||D.schema!=='pAML-promoter-v1') {el('graphTitle').textContent='Could not load the promoter data';el('graphSubtitle').textContent='Ensure data/atlas.js is present in the same folder as this page.';return;}
const R={PROBE:0,DIR:1,P:2,I:3,N:4,NP:5,NI:6,NN:7,QIN:8,QPI:9,QPN:10,QTOST:11,MASK:12,W:13,WSPEC:14,WSHARE:15,DIN:16};
const flagIndex=Object.fromEntries(D.flags.map((x,i)=>[x,i]));
const groupList=D.groups.map(x=>x[0]);
const groupCounts=Object.fromEntries(D.groups);
const groupMeanMap=new Map(Object.entries(D.groupMeans).map(([key,rows])=>[key,new Map(rows.map(a=>[a[0],{beta:a[1],n:a[2]}]))]));
const geneList=Object.keys(D.geneRows);
const state={gene:'',mode:'gene',cpg:null,signature:'all',contrast:'IN',qMax:0.05,effectMin:.15,wMin:0,select:{type:'none',id:null}};
let current={rows:[],probes:[],edges:[],nodes:[]};
const colors={red:'#c84d5d',blue:'#3478bf',mixed:'#888999'};
const compareName=()=>state.contrast==='IN'?'incoming bridges − source nonbridges':'incoming bridges − target primary';
const clamp=(x,a,b)=>Math.min(b,Math.max(a,x));
const colorLerp=(a,b,t)=>`rgb(${a.map((v,i)=>Math.round(v+(b[i]-v)*t)).join(',')})`;
const betaColor=b=>{if(b===null||!Number.isFinite(b))return '#e9edf2';let v=clamp(b,0,1);return v<=.5?colorLerp([166,205,239],[247,249,251],v*2):colorLerp([247,249,251],[245,188,197],(v-.5)*2);};
const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt=(x,d=3)=>x===null||x===undefined||!Number.isFinite(x)?'—':Number(x).toFixed(d);
const fmtP=q=>q===null||q===undefined||!Number.isFinite(q)?'—':q<.0001?'<0.0001':q.toFixed(4);
const labelDir=di=>{let d=D.directions[di];return `${d[0]} → ${d[1]}`;};
const flag=(mask,key)=>(mask & (1<<flagIndex[key]))!==0;
const flagGroup={
 'restricted':['Primary_Restricted'],
 'conserved':['Neighborhood_Conserved'],
 'core':['Core_Target_06','Core_Source_06'],
 'bridge':['Prior_Bridge_Source_to_Target_07a','Prior_Bridge_Target_to_Source_07a'],
 'distinguishing':['Prior_06b_A_Higher_015','Prior_06b_B_Higher_015','Prior_06b_A_Defining_Concordant_Distinct_015','Prior_06b_B_Defining_Concordant_Distinct_015','Prior_06b_Target_Higher_015','Prior_06b_Source_Higher_015','Prior_06b_Target_Defining_Distinct_015','Prior_06b_Source_Defining_Distinct_015'],
 'shared':['Prior_06b_Shared_Higher_015','Prior_06b_Shared_Lower_015','Prior_06b_Shared_Higher_Formal_015','Prior_06b_Shared_Lower_Formal_015'],
 '11c':['Prior_11C_Target_Restricted','Prior_11C_Target_Conserved']
};
const fullFlags=mask=>D.flags.filter((k,i)=>(mask&(1<<i))!==0);
const shortFlags=mask=>{
  let out=[];if(flag(mask,'Primary_Restricted'))out.push('11D restricted');if(flag(mask,'Neighborhood_Conserved'))out.push('11D conserved');
  if(flagGroup.core.some(x=>flag(mask,x)))out.push('Core 06');if(flagGroup.bridge.some(x=>flag(mask,x)))out.push('Bridge 07a');
  if(flagGroup.distinguishing.some(x=>flag(mask,x)))out.push('Pair higher/distinct');
  if(flagGroup.shared.some(x=>flag(mask,x)))out.push('Pair shared');
  if(flagGroup['11c'].some(x=>flag(mask,x)))out.push('11C');return out;
};
function qualifyCategory(row){if(state.signature==='all')return true;return (flagGroup[state.signature]||[]).some(x=>flag(row[R.MASK],x));}
function diff(row){return state.contrast==='IN'?row[R.DIN]:row[R.I]===null||row[R.P]===null?null:row[R.I]-row[R.P];}
function qValue(row){return state.contrast==='IN'?row[R.QIN]:row[R.QPI];}
function qualifiesArrow(row){let v=diff(row),q=qValue(row);return Number.isFinite(v)&&Math.abs(v)+1e-10>=state.effectMin&&(state.qMax===1||q!==null&&q<=state.qMax+1e-10);}
function updateCpGs(){let records=D.geneRows[state.gene]||[];let ids=[...new Set(records.map(i=>D.rows[i][R.PROBE]))].sort((a,b)=>D.probes[a][0].localeCompare(D.probes[b][0]));
  let box=el('cpgSelect');box.replaceChildren();for(let p of ids){const o=document.createElement('option');o.value=p;o.textContent=D.probes[p][0]+' · chr'+D.probes[p][1]+':'+D.probes[p][2];box.append(o);}
  if(!ids.includes(state.cpg))state.cpg=ids[0]??null;box.value=state.cpg===null?'':String(state.cpg);el('cpgSelectHolder').hidden=state.mode!=='cpg';
}
function setGene(g){if(!D.geneRows[g])return;state.gene=g;state.cpg=null;state.select={type:'none',id:null};el('geneSearch').value=g;el('suggestions').hidden=true;updateCpGs();draw();}
function activeRows(){let ids=D.geneRows[state.gene]||[];return ids.filter(i=>{let r=D.rows[i];return (state.mode!=='cpg'||r[R.PROBE]===state.cpg)&&qualifyCategory(r)&&(r[R.W]!==null&&r[R.W]>=state.wMin/100-1e-10||state.wMin===0);});}
function summarize(){let rowIndices=activeRows(),byDir=new Map();let probes=new Set();for(const ri of rowIndices){let r=D.rows[ri];probes.add(r[R.PROBE]);if(!qualifiesArrow(r))continue;let di=r[R.DIR];if(!byDir.has(di))byDir.set(di,[]);byDir.get(di).push(ri);}
 let edges=[];for(let [id,ids] of byDir){let vals=ids.map(i=>diff(D.rows[i]));let pos=vals.filter(v=>v>0).length,neg=vals.filter(v=>v<0).length;
  let sign=pos&&neg?'mixed':pos?'red':neg?'blue':'mixed';let mean=vals.reduce((a,b)=>a+b,0)/vals.length;
  edges.push({id,ids,nCpG:new Set(ids.map(i=>D.rows[i][R.PROBE])).size,sign,mean,nPositive:pos,nNegative:neg});}
 const nodeSummary=[];for(const group of groupList){let gm=groupMeanMap.get(group)||new Map();let observations=[...probes].map(p=>({p,value:gm.get(p)})).filter(x=>x.value);
  let mean=observations.length?observations.reduce((s,v)=>s+v.value.beta,0)/observations.length:null;
  nodeSummary.push({group,mean,nProbes:observations.length,available:probes.size,count:groupCounts[group],obs:observations});}
 return {rows:rowIndices,probes:[...probes],edges:edges.sort((a,b)=>b.nCpG-a.nCpG),nodes:nodeSummary};
}
function svg(name,attrs={},text=null){let x=document.createElementNS(NS,name);for(let [k,v] of Object.entries(attrs))x.setAttribute(k,String(v));if(text!==null)x.textContent=text;return x;}
function locate(group){let i=groupList.indexOf(group), angle=-Math.PI/2+i*2*Math.PI/groupList.length;
 return {x:500+365*Math.cos(angle),y:393+304*Math.sin(angle)};}
function edgePath(g1,g2){let a=locate(g1),b=locate(g2),dx=b.x-a.x,dy=b.y-a.y,len=Math.hypot(dx,dy)||1,ux=dx/len,uy=dy/len,nx=-uy,ny=ux;
 let r=34,start={x:a.x+ux*r,y:a.y+uy*r},end={x:b.x-ux*(r+8),y:b.y-uy*(r+8)};
 let mid={x:(start.x+end.x)/2,y:(start.y+end.y)/2};let bend=Math.min(90,20+len*.085);
 let cx=mid.x+nx*bend+(.06*(500-mid.x)),cy=mid.y+ny*bend+(.06*(393-mid.y));return `M${start.x.toFixed(2)} ${start.y.toFixed(2)} Q${cx.toFixed(2)} ${cy.toFixed(2)} ${end.x.toFixed(2)} ${end.y.toFixed(2)}`;}
function tooltip(content,e){const t=el('tooltip');t.innerHTML=content;t.hidden=false;const rect=t.getBoundingClientRect();t.style.left=clamp(e.clientX+14,8,window.innerWidth-rect.width-12)+'px';t.style.top=clamp(e.clientY+12,8,window.innerHeight-rect.height-12)+'px';}
function hideTip(){el('tooltip').hidden=true;}
function selectEdge(id){state.select={type:'edge',id};draw();}
function selectNode(id){state.select={type:'node',id};draw();}
function populateNetwork(){const edgeLayer=el('edgeLayer'),nodeLayer=el('nodeLayer');edgeLayer.replaceChildren();nodeLayer.replaceChildren();
 for(let ed of current.edges){let [src,tgt,nSource,nBridge]=D.directions[ed.id],s=ed.sign,stroke=colors[s];let active=state.select.type==='edge'&&state.select.id===ed.id;
  let g=svg('g',{class:'edge'+(active?' active':''),'aria-label':`${src} to ${tgt}: ${ed.nCpG} qualifying CpGs`});
  let path=edgePath(src,tgt);let p=svg('path',{d:path,stroke,'stroke-width':2.1+Math.min(4,Math.log2(1+ed.nCpG)*.9),'stroke-dasharray':s==='mixed'?'7 5':'none','marker-end':`url(#arrow-${s})`,opacity:active?1:.75});
  let hit=svg('path',{d:path,stroke:'transparent','stroke-width':18,fill:'none',tabindex:0,role:'button'});
  let tip=()=>`<strong>${esc(src)} → ${esc(tgt)}</strong><div class="dim">${ed.nCpG} qualifying CpG(s); ${nBridge} bridge patients (before CpG-specific missingness)</div>${s==='mixed'?'Mixed directions':`Mean qualifying-site Δβ ${fmt(ed.mean)}`}<div class="dim">${esc(compareName())}</div>`;
  hit.addEventListener('pointermove',e=>tooltip(tip(),e));hit.addEventListener('pointerleave',hideTip);
  hit.addEventListener('click',()=>selectEdge(ed.id));hit.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();selectEdge(ed.id);}});g.append(p,hit);edgeLayer.append(g);
 }
 for(let nd of current.nodes){let {x,y}=locate(nd.group);let active=state.select.type==='node'&&state.select.id===nd.group;
  let g=svg('g',{class:'node'+(nd.mean===null?' no-data':'')+(active?' active':''),tabindex:0,role:'button','aria-label':nd.group+' mean beta '+(nd.mean===null?'unavailable':fmt(nd.mean))});
  g.append(svg('circle',{cx:x,cy:y,r:32,fill:betaColor(nd.mean),stroke:active?'#264d70':nd.mean===null?'#ced8e3':'#879eb5','stroke-width':active?4:1.5}));
  g.append(svg('text',{x,y:y+1,'text-anchor':'middle'},nd.group));g.append(svg('text',{x,y:y+14,'text-anchor':'middle',class:'n-label'},nd.mean!==null?'β '+fmt(nd.mean,2):'no data'));
  let tip=()=>`<strong>${esc(nd.group)}</strong><div>Original primary group mean β: ${nd.mean===null?'not represented':fmt(nd.mean)}</div><div class="dim">${nd.nProbes}/${nd.available} selected CpGs have group-level values<br>Primary group size: ${nd.count??'not supplied (not eligible for 11D)'}</div>`;
  g.addEventListener('pointermove',e=>tooltip(tip(),e));g.addEventListener('pointerleave',hideTip);g.addEventListener('click',()=>selectNode(nd.group));g.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();selectNode(nd.group);}});nodeLayer.append(g);
 }
}
function updateRows(){let subset=current.rows;if(state.select.type==='edge')subset=subset.filter(i=>D.rows[i][R.DIR]===state.select.id);
 let tb=el('evidenceBody');tb.replaceChildren();let sorted=subset.slice().sort((a,b)=>{
  const qa=qValue(D.rows[a]),qb=qValue(D.rows[b]);return (qa??999)-(qb??999)||Math.abs(diff(D.rows[b])??0)-Math.abs(diff(D.rows[a])??0);
 });
 for(const ri of sorted){let r=D.rows[ri],p=D.probes[r[R.PROBE]],q=qValue(r),del=diff(r),hit=qualifiesArrow(r);let row=document.createElement('tr');
  const cells=[labelDir(r[R.DIR]),p[0],p[1]+':'+p[2],fmt(r[R.P]),fmt(r[R.I]),fmt(r[R.N]),fmt(del),fmtP(q),r[R.W]===null?'—':Math.round(100*r[R.W])+'%',shortFlags(r[R.MASK]).join(' · ')];
  cells.forEach((v,j)=>{let cell=document.createElement('td');cell.textContent=v;if(j===6)cell.className=hit&&del!==null?(del>0?'effect-up':'effect-down'):'effect-flat';if(j===9){cell.className='labels-cell';cell.title=fullFlags(r[R.MASK]).join('; ');}row.append(cell);});tb.append(row);
 }
 if(!sorted.length){let tr=document.createElement('tr');let td=document.createElement('td');td.colSpan=10;td.className='empty-row';td.textContent='No CpGs match the current filters.';tr.append(td);tb.append(tr);}
 el('tableCount').textContent=sorted.length+' CpG–direction associations';el('tableExplainer').textContent='P = target primary; I = incoming from source; N = source nonbridge. Values are observed mean β. Bold Δβ means the row passes the active arrow cutoffs. q is CpG-wise BH FDR for the selected comparison.';
}
function metric(title,value,sub=''){return `<div class="metric"><span class="micro">${esc(title)}</span><strong class="value">${esc(value)}</strong>${sub?`<div class="help">${esc(sub)}</div>`:''}</div>`;}
function detailEdges(){if(!current.edges.length)return '<div class="note">No arrow meets the current CpG-level contrast and cutoff. Try a different gene, a smaller |Δβ| threshold, or broader FDR. Gray circles may indicate unrepresented CpGs, not unmethylated genes.</div>';
 return '<div class="dir-list"><h3>Qualifying directions</h3>'+current.edges.slice(0,26).map(ed=>{
  let d=D.directions[ed.id],active=state.select.type==='edge'&&ed.id===state.select.id;return `<button data-dir="${ed.id}" class="${active?'active':''}"><span style="color:${colors[ed.sign]}">●</span><span>${esc(d[0])} → ${esc(d[1])}</span><b>${ed.nCpG} CpG${ed.nCpG===1?'':'s'}</b></button>`;}).join('')+'</div>';
}
function showDetails(){const root=el('details');let html='';
 if(!state.gene){root.innerHTML='<p class="help">Select a gene from the search box.</p>';return;}
 let totalProbes=current.probes.length, totalEdges=current.edges.length;
 if(state.select.type==='edge'){
  let ed=current.edges.find(x=>x.id===state.select.id);
  if(!ed){state.select={type:'none',id:null};return showDetails();}
  let [src,tgt,nS,nB,nN]=D.directions[ed.id];html=`<h3>${esc(src)} → ${esc(tgt)}</h3><p class="detail-subtitle">${esc(state.gene)} | incoming H neighborhood of ${esc(tgt)}</p>`+
   metric('Qualifying promoter CpGs',ed.nCpG+' sites','CpG-level tests, not independent gene-level inference')+
   metric('Original source patients',nS+' total',nB+' bridge / '+nN+' nonbridge (13C pooled counts)')+
   metric('Mean Δβ among qualifying CpGs',fmt(ed.mean),compareName())+
   `<div class="note">${ed.sign==='mixed'?`<strong>Mixed methylation directions.</strong> ${ed.nPositive} CpGs higher and ${ed.nNegative} lower. No single red/blue gene-level direction is imposed.`:ed.sign==='red'?'Higher beta in the incoming population at every qualifying promoter CpG for this arrow.':'Lower beta in the incoming population at every qualifying promoter CpG for this arrow.'} CpG-level observed counts vary with missing values.</div>`;
  let hits=ed.ids.slice().sort((a,b)=>Math.abs(diff(D.rows[b]))-Math.abs(diff(D.rows[a]))).slice(0,16);
  html+='<h3>Strongest CpG evidence</h3>'+hits.map(i=>{let r=D.rows[i],p=D.probes[r[R.PROBE]];return `<div class="metric"><strong>${esc(p[0])}</strong> <span class="micro">chr${esc(p[1])}:${p[2]}</span><div>Δβ <b class="${diff(r)>0?'effect-up':'effect-down'}">${fmt(diff(r))}</b> · q ${fmtP(qValue(r))}</div><div class="micro">β P/I/N: ${fmt(r[R.P])} / ${fmt(r[R.I])} / ${fmt(r[R.N])}</div><div class="micro">Observed n P/I/N: ${r[R.NP]} / ${r[R.NI]} / ${r[R.NN]}</div><div class="help">${esc(shortFlags(r[R.MASK]).join(' · '))}</div></div>`;}).join('');
  html+='<div class="note">Only selected, promoter-annotated CpGs in the 12C/11D discovery suite appear. Passing a threshold here does not establish gene expression regulation or causality.</div>';
 }else if(state.select.type==='node'){
  let nd=current.nodes.find(x=>x.group===state.select.id);html=`<h3>${esc(nd.group)}: original primary patients</h3><p class="detail-subtitle">${esc(state.gene)}; group-level beta across represented promoter CpGs</p>`+
   metric('Mean β across available CpGs',fmt(nd.mean),'Descriptive arithmetic mean across unique CpGs; not a differential test')+
   metric('CpG coverage',nd.nProbes+' / '+nd.available,'Unrepresented sites are missing, not beta=0')+
   metric('Original primary group size',nd.count??'Not available','13C tested-group size; smaller/empty groups were not assessed');
  html+='<h3>Group-level CpG means</h3>'+nd.obs.slice(0,35).map(v=>`<div class="metric"><strong>${esc(D.probes[v.p][0])}</strong><span class="value">β ${fmt(v.value.beta)}</span><div class="micro">n=${v.value.n} observed original-primary patients</div></div>`).join('');
  html+='<div class="note">Values were either reported directly as target-primary means (11D) or reconstructed exactly from source bridge/nonbridge mean beta and CpG-specific observed counts. Different groups may have different represented CpGs, so compare coverage before interpreting differences.</div>';
 }else{
   html=`<h3>${esc(state.gene)}</h3><p class="detail-subtitle">GRCh38 promoter overlap; ${state.mode==='cpg'?'single probe':'gene view'}</p>`+
     metric('Selected distinct promoter CpGs',totalProbes)+metric('Directions passing CpG-level filters',totalEdges)+metric('NMF group node means available',current.nodes.filter(n=>n.mean!==null).length+' / 23','Source and target primary groups; four original groups were not eligible for 11D')+
     '<div class="note"><b>Reading the network:</b> Red/blue circles reflect absolute mean β (low ↔ high) for original primary members. Red/blue arrows reflect a signed bridge versus reference contrast (higher ↔ lower). Gene-level arrows summarize individual CpG tests; they are not gene-level significance tests.</div>';
 }
 if(state.select.type!=='edge')html+=detailEdges();
 root.innerHTML=html;root.querySelectorAll('[data-dir]').forEach(btn=>btn.addEventListener('click',()=>selectEdge(Number(btn.dataset.dir))));
}
function draw(){current=summarize();if(state.select.type==='edge'&&!current.edges.some(x=>x.id===state.select.id))state.select={type:'none',id:null};
 el('geneStatus').textContent=`${new Set((D.geneRows[state.gene]||[]).map(i=>D.rows[i][R.PROBE])).size} unique promoter CpGs for ${state.gene} before filters`;
 el('graphTitle').textContent=state.gene+(state.mode==='cpg'&&state.cpg!==null?` · ${D.probes[state.cpg][0]}`:'');
 el('graphSubtitle').textContent=`${current.probes.length} distinct selected promoter CpGs · ${current.edges.length} directional arrows with CpG-level evidence · ${compareName()}`;
 el('edgeSummary').textContent=current.edges.length?`${current.edges.filter(e=>e.sign==='red').length} higher · ${current.edges.filter(e=>e.sign==='blue').length} lower · ${current.edges.filter(e=>e.sign==='mixed').length} mixed`: 'No arrows meet the active thresholds';
 populateNetwork();updateRows();showDetails();
}
function setupSearch(){let inp=el('geneSearch'),list=el('suggestions');
 const sugg=()=>{let q=inp.value.trim().toUpperCase();list.replaceChildren();if(!q){list.hidden=true;return;}
  let starts=geneList.filter(g=>g.toUpperCase().startsWith(q)).slice(0,18),others=geneList.filter(g=>!g.toUpperCase().startsWith(q)&&g.toUpperCase().includes(q)).slice(0,12);
  for(let gene of starts.concat(others)){let b=document.createElement('button');b.textContent=gene;b.type='button';b.addEventListener('click',()=>setGene(gene));list.append(b);}list.hidden=!list.children.length;};
 inp.addEventListener('input',sugg);inp.addEventListener('focus',()=>{if(inp.value!==state.gene)sugg();});
 inp.addEventListener('keydown',e=>{if(e.key==='Escape'){list.hidden=true;inp.blur();}if(e.key==='Enter'){e.preventDefault();let first=list.querySelector('button');if(first)first.click();else{let found=geneList.find(g=>g.toUpperCase()===inp.value.trim().toUpperCase());if(found)setGene(found);}}});
 document.addEventListener('pointerdown',e=>{if(!list.contains(e.target)&&e.target!==inp)list.hidden=true;});
 document.querySelectorAll('[data-try]').forEach(b=>b.addEventListener('click',()=>setGene(b.dataset.try)));
}
function setupControls(){setupSearch();
 const bind=(id,event,cb)=>el(id).addEventListener(event,e=>{cb(e.target.value);state.select={type:'none',id:null};draw();});
 bind('viewMode','change',v=>{state.mode=v;updateCpGs();});
 bind('cpgSelect','change',v=>{state.cpg=Number(v);});
 bind('signature','change',v=>state.signature=v);
 bind('contrast','change',v=>state.contrast=v);
 bind('qMax','change',v=>state.qMax=Number(v));
 bind('effectMin','input',v=>{state.effectMin=Number(v);el('effectMinLabel').textContent=fmt(state.effectMin,3);});
 bind('wMin','input',v=>{state.wMin=Number(v);el('wMinLabel').textContent=v+'%';});
 el('exportCsv').addEventListener('click',exportCSV);
}
function exportCSV(){let rs=current.rows;if(state.select.type==='edge')rs=rs.filter(i=>D.rows[i][R.DIR]===state.select.id);
 const heads=['Gene','CpG','Chr','Start_GRCh38','Source','Target','N_Target_Primary_Observed','N_Source_Bridge_Observed','N_Source_Nonbridge_Observed','Mean_Beta_Target_Primary','Mean_Beta_Source_Bridge','Mean_Beta_Source_Nonbridge','Delta_Bridge_Minus_Source_Nonbridge','BH_Bridge_vs_SourceNonbridge','BH_TargetPrimary_vs_Bridge','BH_TargetPrimary_vs_SourceNonbridge','BH_TOST_Primary_Equiv_Bridge','W_Target_Percentile','W_Target_Specificity','All_Provenance_Flags'];
 const lines=[heads];for(let i of rs){let r=D.rows[i],p=D.probes[r[R.PROBE]],d=D.directions[r[R.DIR]];
   lines.push([state.gene,p[0],p[1],p[2],d[0],d[1],r[R.NP],r[R.NI],r[R.NN],r[R.P],r[R.I],r[R.N],r[R.DIN],r[R.QIN],r[R.QPI],r[R.QPN],r[R.QTOST],r[R.W],r[R.WSPEC],fullFlags(r[R.MASK]).join(';')]);}
 const contents=lines.map(cols=>cols.map(v=>`"${String(v??'').replaceAll('"','""')}"`).join(',')).join('\r\n');
 let blob=new Blob(['\ufeff',contents],{type:'text/csv;charset=utf-8;'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=`Promoter_${state.gene.replace(/[^\w.-]+/g,'_')}_${state.contrast}.csv`;a.click();URL.revokeObjectURL(url);
}
function boot(){el('datasetStats').innerHTML=`<span class="stat"><b>${D.meta.cohortPatients}</b> original patients</span><span class="stat"><b>${D.meta.uniqueCpGs.toLocaleString()}</b> promoter CpGs</span><span class="stat"><b>${D.meta.promoterGenes.toLocaleString()}</b> promoter gene labels</span><span class="stat"><b>${D.meta.eligibleDirections}</b> analyzed bridge directions with promoter CpGs</span><span class="stat"><b>19</b> eligible original primary groups</span>`;
 setupControls();setGene(D.geneRows.GATA3?'GATA3':geneList[0]);
}
boot();
})();
