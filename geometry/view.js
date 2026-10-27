/* Selection and page layout; the campaign SurfaceView renderer is unchanged. */
'use strict';
const ringEndings={ring1:[10],ring2:[20],ring3:[30],ring4:[40],ring5:[51,52,53],ring6:[60]};
const mainEndings=Object.values(ringEndings).flat();
function detectorRegion(det){
 if(det>=73001&&det<=73028)return det%4===3?'open':det%4===1?'closed':'transition';
 if(det>=110000&&det<112000)return 'transition';
 return det%200<100?'open':'closed';
}
function detectorMatches(c,det){
 const sm=det>=73001&&det<=73028;
 const main=((det>=110000&&det<112000)||(det>=150000&&det<152000))&&mainEndings.includes(det%100);
 if(c.system==='showermax'&&!sm)return false;
 if(c.system!=='all'&&c.system!=='showermax'&&!main)return false;
 if(ringEndings[c.system]&&!ringEndings[c.system].includes(det%100))return false;
 return c.region==='all'||((main||sm)&&detectorRegion(det)===c.region);
}
function buildCut(c){
 const parts=[];
 const mainSelection=()=>{
  const family=c.region==='transition'?'hit.det >= 110000 && hit.det < 112000':
   c.region==='all'?'((hit.det >= 110000 && hit.det < 112000) || (hit.det >= 150000 && hit.det < 152000))':
   `hit.det >= 150000 && hit.det < 152000 && hit.det%200 ${c.region==='open'?'<':'>='} 100`;
  const endings=ringEndings[c.system]||mainEndings;
  const ring=endings.map(n=>`hit.det%100 == ${n}`).join(' || ');
  return family+' && '+(endings.length===1?ring:`(${ring})`);
 };
 const smSelection=()=> 'hit.det >= 73001 && hit.det <= 73028'+(c.region==='all'?'':
  c.region==='transition'?' && hit.det%2 == 0':` && hit.det%4 == ${c.region==='open'?3:1}`);
 // Choose one detector expression at the most specific requested level.
 if(c.tile!==null){
  if(!detectorMatches(c,c.tile))return '0';
  parts.push(`hit.det == ${Number(c.tile)}`);
 }else if(c.system==='showermax')parts.push(smSelection());
 else if(c.system!=='all')parts.push(mainSelection());
 else if(c.region!=='all')parts.push(`((${mainSelection()}) || (${smSelection()}))`);
 if(c.particle!=='all')parts.push(`hit.pid == ${Number(c.particle)}`);
 if(c.mode!=='raw'&&c.system!=='showermax'&&c.system!=='all'){
  // An explicit e−, e+ or proton selection already guarantees nonzero charge.
  if(!['11','-11','2212'].includes(String(c.particle)))parts.push('charged(hit.pid)');
  parts.push('hit.k >= 1');
 }
 if(c.mode==='secondary'||c.ancestry==='secondary')parts.push('hit.mtrid > 0');
 else if(c.ancestry==='primary')parts.push('hit.mtrid == 0');
 if(c.direction!=='all')parts.push(`hit.pz ${c.direction==='forward'?'>':'<'} 0`);
 return parts.join(' &&\n')||'1';
}
function validatePlotRange(low,high,log){
 if(!Number.isFinite(low)||!Number.isFinite(high)||low>=high||(log&&low<=0))throw Error('Choose increasing axis limits; log axes require positive limits.');
 if(Math.abs(low)>1e15||Math.abs(high)>1e15||!Number.isFinite(high-low)||(log&&!Number.isFinite(high/low)))throw Error('Axis range is too large. Use closer limits.');
 return [low,high];
}
function signalAxisValues(values,previous,next){
 const copy={...values},factor=previous==='pe'&&next==='current'?1.602176634e-10:previous==='current'&&next==='pe'?1/1.602176634e-10:null;
 for(const key of ['plot-y-min','plot-y-max'])if(copy[key]!==undefined&&copy[key]!=='')copy[key]=factor===null?'':String(Number(copy[key])*factor);
 return copy;
}
const mainCut=buildCut({mode:'hits',system:'main',region:'all',ancestry:'all',direction:'all',particle:'all',tile:null});
if(typeof module!=='undefined'&&module.exports){module.exports={buildCut,mainCut,validatePlotRange,signalAxisValues};}
else {
const $=id=>document.getElementById(id);
let config={mode:'hits',system:'main',region:'all',ancestry:'all',direction:'all',particle:'all',tile:null};
let signalSelection='hits',rawHitData=null;
let geometryFaces=[],hitData=null,revision=0,initialized=false,stopped=false,busy=false,choosing=false,requesting=false,hitKey='',applied=null,editing=false;
const surface=new HitMap3D.SurfaceView($('surface-map'),{onSelect:selectTile});
// Keep the scale outside the canvas and consolidate camera/navigation controls.
const toolbar=$('surface-map').querySelector('.surfaceToolbar');
const navigation=$('surface-map').querySelector('.surfaceNavigation');
for(const group of [...navigation.children])if(group.tagName==='DIV')toolbar.append(group);
navigation.remove();$('surface-map').append($('surface-map').querySelector('.surfaceLegend'));
// Fit the actual projected extent to this workspace's larger canvas, with a margin.
const baseFit=surface.fit.bind(surface);
surface.fit=()=>{
 baseFit();if(surface.unwrapped||!surface.pickFaces.length)return;
 const r=surface.canvas.getBoundingClientRect(),camera={yaw:surface.yaw,pitch:surface.pitch,zoom:1,pan:[0,0]};
 if(!r.width||!r.height)return;
 const points=surface.pickFaces.flatMap(f=>f.points.map(p=>HitMap3D.projectPoint(p,camera,r.width,r.height)));
 const xs=points.map(p=>p[0]),ys=points.map(p=>p[1]);
 const x0=Math.min(...xs),x1=Math.max(...xs),y0=Math.min(...ys),y1=Math.max(...ys);
 surface.zoom=Math.min(.90*r.width/(x1-x0),.90*r.height/(y1-y0));
 surface.pan=[-(x0+x1)/2*surface.zoom,-(y0+y1)/2*surface.zoom];surface.draw();
};
toolbar.querySelectorAll('[data-camera]').forEach(b=>b.addEventListener('click',()=>surface.fit()));
function filePaths(){return $('hit-file').value.split('\n').map(s=>s.trim()).filter(Boolean);}

function selectTile(det){
 if(busy||requesting){surface.selectedDetector=config.tile;return;}
 config.tile=det;
 if(det!==null){config.region=detectorRegion(det);if(det>=73001&&det<=73028){config.system='showermax';if(config.mode!=='raw')config.direction='forward';}else if(config.system==='showermax'||config.system.startsWith('ring'))config.system='main';}
 surface.selectedDetector=det;
 selectAndRun();
}
function syncControls(){
 for(const id of ['system','region','ancestry','particle'])$(id).value=config[id];
 document.querySelectorAll('[data-mode]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.mode===config.mode&&!editing)));
 document.querySelectorAll('[data-direction]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.direction===config.direction)));
 $('ancestry-row').hidden=config.mode==='secondary';
 const charged=config.mode!=='raw'&&config.system!=='showermax';
 for(const option of $('particle').options)option.disabled=charged&&['22','2112'].includes(option.value);
 $('system').querySelector('[value="all"]').disabled=config.mode!=='raw';
 $('tile-clear').hidden=config.tile===null;
 $('selection-label').textContent=config.tile===null?'All tiles':`Detector ${config.tile}`;
 $('cut-state').textContent=editing?'Custom':'Preset';
}
function clearTile(){config.tile=null;surface.selectedDetector=null;surface.focusedDetector=null;surface.fit();}
function selectAndRun(){
 editing=false;syncControls();$('hit-cut').value=buildCut(config);
 sessionStorage.setItem('detector-selection',JSON.stringify({config,cut:$('hit-cut').value}));
 if($('hit-file').value.trim())loadHits();else{$('hit-status').textContent='Choose ROOT files.';showApplied();}
}
function setGeometry(){
 if(!geometryFaces.length){surface.clear('Load geometry first.');return;}
 const cells=hitData?.cells||[],shown=cells.reduce((n,c)=>n+c[3],0);
 surface.setData({faces:geometryFaces,cells,cell_mm:10,total:shown,plane:'all_detectors',unit:hitData?.unit||'hits',signalDisplay:'pe',linear:config.mode==='secondary'});
 if(!hitData)surface.status('');
 else if(!hitData.positions_available)surface.status('Hit positions unavailable');
 else if(!hitData.surfaces_available)surface.status('Load geometry, then apply cuts');
 else if(hitData.hits>0&&!shown)surface.status('Selected hits do not match these surfaces');
 $('hit-legend').textContent=!hitData?'Click a tile to select it.':`${HitMap3D.formatNumber(shown)} ${hitData.unit||'hits'} mapped`+(hitData.unsupported_pe?` · ${hitData.unsupported_pe.toLocaleString()} outside PE model support`:'')+
  (hitData.unmatched_positions?` · ${hitData.unmatched_positions.toLocaleString()} outside surfaces`:'')+
  (hitData.other_detectors?` · ${hitData.other_detectors.toLocaleString()} on other detectors`:'')+
  (hitData.invalid_positions?` · ${hitData.invalid_positions.toLocaleString()} invalid positions`:'');
}
function showApplied(){
 const pending=applied&&($('hit-cut').value.trim()!==applied.cut||$('hit-file').value.trim()!==applied.path);
 $('applied-cut').textContent=applied?(pending?`Currently applied: ${applied.cut}`:'Applied to 3D + distributions · '+(hitData?.unit||'hits')):'No selection applied';
 $('shared-cuts').textContent='Cuts · same selection as 3D';
 $('charge-help').hidden=!$('hit-cut').value.includes('charged(');
 $('shared-cut-state').textContent=pending?'Edits not applied':!applied?'No selection applied':'';
}
const svgNS='http://www.w3.org/2000/svg';
function svgNode(name,attrs={},text){const n=document.createElementNS(svgNS,name);for(const [k,v] of Object.entries(attrs))n.setAttribute(k,v);if(text!==undefined)n.textContent=text;return n;}
// Display controls never alter selection, weighting, or ROOT histogram bins.
let plotDimension=1,plotBinning={},binDraft=false;
const plotViewInputs=['plot-x-min','plot-x-max','plot-y-min','plot-y-max','plot-title','plot-x-label','plot-y-label'];
const plotViewChecks=['plot-log-x','plot-log-y','plot-log-color','plot-grid-x','plot-grid-y'];
const plotViews=new Map();let activePlotKey='1:radius';
function plotKey(){return plotDimension+':'+$(plotDimension===2?'plot-pair':'plot-kind').value;}
function rememberPlotView(){plotViews.set(activePlotKey,{values:Object.fromEntries(plotViewInputs.map(k=>[k,$(k).value])),checks:Object.fromEntries(plotViewChecks.map(k=>[k,$(k).checked]))});}
function changePlot(){
 rememberPlotView();const next=plotKey();if(next===activePlotKey){renderPlot();return;}
 activePlotKey=next;binDraft=false;const saved=plotViews.get(next);
 if(saved){for(const [k,v] of Object.entries(saved.values))$(k).value=v;for(const [k,v] of Object.entries(saved.checks))$(k).checked=v;syncBinControls();renderPlot();}
 else resetPlot();
}
const plotNames={radius:'Radius',energy:'Kinetic energy',total_energy:'Total energy',edep:'Deposited energy',x:'Hit x',y:'Hit y',z:'Hit z',vertex:'Creation z',vertex_r:'Creation radius',phi:'Position azimuth',theta:'Momentum polar angle',time:'Hit time',momentum:'Momentum',px:'Momentum x',py:'Momentum y',pz:'Momentum z',particles:'Particle type',detectors:'Detector ID',asymmetry:'Generator asymmetry',q2:'Generator Q²',w:'Generator W',beam_p:'Beam momentum',rate:'Stored event rate value',spin_x:'Polarization x',spin_y:'Polarization y',spin_z:'Polarization z'};
const pairNames={xy:'Hit Y vs X',energy_radius:'Energy vs radius',creation_rz:'Creation radius vs z',radius_z:'Radius vs z',radius_phi:'Radius vs azimuth',energy_theta:'Energy vs polar angle',energy_time:'Energy vs time',energy_creation:'Energy vs creation z',momentum_radius:'Momentum vs radius',px_py:'Momentum Y vs X',edep_energy:'Deposited energy vs kinetic energy',q2_w:'Q² vs W',asymmetry_q2:'Asymmetry vs Q²'};
for(const [id, names] of [['plot-kind',plotNames],['plot-pair',pairNames]])for(const [key,label] of Object.entries(names)){const o=document.createElement('option');o.value=key;o.textContent=label;$(id).append(o);}
const unsupported=document.createElement('optgroup');unsupported.label='Need response / calibration / normalization';
for(const label of ['ADC / voltage / gain','Acceptance / efficiency','Noise / resolution / MC uncertainty','Beam current / charge','Dose / temperature']){const o=document.createElement('option');o.textContent=label;o.disabled=true;unsupported.append(o);}$('plot-kind').append(unsupported);
function resetPlot(){
 binDraft=false;
 for(const id of ['plot-x-min','plot-x-max','plot-y-min','plot-y-max'])$(id).value='';
 const h=hitData?.histograms, pair=hitData?.histograms2d?.[$('plot-pair').value];
 $('plot-log-x').checked=!!h?.[plotDimension===2?pair?.x:$('plot-kind').value]?.log_bins;
 $('plot-log-y').checked=plotDimension===2&&!!h?.[pair?.y]?.log_bins;
 $('plot-grid-x').checked=true;$('plot-grid-y').checked=true;
 for(const id of ['plot-title','plot-x-label','plot-y-label'])$(id).value='';
 syncBinControls();$('plot-log-color').checked=true;renderPlot();
}
function renderPlot(){
 const svg=$('plot'),empty=$('plot-empty'),kind=$('plot-kind').value,two=plotDimension===2;
 svg.replaceChildren();
 const pair=hitData?.histograms2d?.[$('plot-pair').value],h=hitData?.histograms?.[two?pair?.x:kind],hy=two?hitData?.histograms?.[pair?.y]:null;
 const categorical=!two&&['particles','detectors'].includes(kind);
 const entries=categorical?Object.entries(hitData?.[kind]||{}).sort((a,b)=>Number(a[0])-Number(b[0])):null;
 const available=two?pair?.available:kind==='particles'?hitData?.particles_available:kind==='detectors'?!!hitData:h?.available;
 $('plot-log-x').disabled=categorical;for(const id of ['plot-x-min','plot-x-max'])$(id).disabled=categorical;
 $('axis-status').textContent=categorical?'X uses detector or particle categories.':'';
 const fail=message=>{svg.setAttribute('hidden','');empty.hidden=false;empty.textContent=message;};
 showApplied();
 $('plot-note').textContent='Shared ROOT cuts · Unweighted hit counts';
 if(!hitData){fail('Load hits to see the distribution.');return;}
 if(!available){fail('Unavailable in one or more files: '+(two?[h?.expression,hy?.expression].filter(Boolean).join(' / '):h?.expression||'hit.pid'));return;}
 svg.removeAttribute('hidden');empty.hidden=true;
 const bounds=svg.parentElement.getBoundingClientRect(),W=Math.max(240,Math.floor(bounds.width-6)),H=Math.max(100,Math.floor(bounds.height-6)),left=74,right=two?86:16,top=$('plot-title').value.trim()?34:12,bottom=39,pw=W-left-right,ph=H-top-bottom;
 svg.setAttribute('viewBox',`0 0 ${W} ${H}`);
 const counts=entries?entries.map(e=>e[1]):two?pair.cells.map(c=>c[2]):h.counts;
 const max=Math.max(...counts,0)||1,logX=!categorical&&$('plot-log-x').checked,logY=$('plot-log-y').checked;
 const getRange=(axis,edges,log)=>{
  let low=edges[0],high=edges.at(-1);
  if(log&&low<=0)low=edges.find(v=>v>0);
  $('plot-'+axis+'-min').placeholder=Number.isFinite(low)?Number(low.toPrecision(5)).toString():'Auto';
  $('plot-'+axis+'-max').placeholder=Number.isFinite(high)?Number(high.toPrecision(5)).toString():'Auto';
  if($("plot-"+axis+"-min").value!=='')low=Number($("plot-"+axis+"-min").value);
  if($("plot-"+axis+"-max").value!=='')high=Number($("plot-"+axis+"-max").value);
  return validatePlotRange(low,high,log);
 };
 let xr,yr;try{xr=categorical?[0,Math.max(1,entries.length)]:getRange('x',h.edges,logX);yr=getRange('y',two?hy.edges:[logY?Math.min(...counts.filter(v=>v>0),max)/10:0,max*1.1],logY);}catch(e){$('axis-status').textContent=e.message;fail(e.message);return;}
 svg.dataset.xMin=xr[0];svg.dataset.xMax=xr[1];svg.dataset.yMin=yr[0];svg.dataset.yMax=yr[1];
 const transform=(v,range,log)=>log?Math.log(v/range[0])/Math.log(range[1]/range[0]):(v-range[0])/(range[1]-range[0]);
 const xpos=v=>left+pw*transform(v,xr,logX),ypos=v=>top+ph-ph*transform(v,yr,logY);
 const fmt=HitMap3D.formatNumber;
 const text=(x,y,t,attrs={})=>svg.append(svgNode('text',{x,y,fill:'#526879','font-size':10,'font-family':'system-ui, sans-serif',...attrs},t));
 const ticks=(r,log,n)=>{
  if(log){const a=Math.ceil(Math.log10(r[0])),b=Math.floor(Math.log10(r[1]));if(b>a){const step=Math.max(1,Math.ceil((b-a)/n));return Array.from({length:Math.floor((b-a)/step)+1},(_,i)=>10**(a+i*step));}}
  return Array.from({length:n+1},(_,i)=>log?r[0]*Math.pow(r[1]/r[0],i/n):r[0]+(r[1]-r[0])*i/n);
 };
 const defs=svgNode('defs'),clip=svgNode('clipPath',{id:'distribution-clip'});clip.append(svgNode('rect',{x:left,y:top,width:pw,height:ph}));defs.append(clip);svg.append(defs);
 const bars=svgNode('g',{'clip-path':'url(#distribution-clip)'});svg.append(bars);
 let visible=0,hidden=0;
 const rect=(x0,x1,y0,y1,n,color,title)=>{
  // A log axis cannot display zero/negative bins; count them in the note.
  if(x1<=xr[0]||x0>=xr[1]||(two&&(y1<=yr[0]||y0>=yr[1]))){hidden+=n;return;}
  const a=xpos(Math.max(x0,xr[0])),b=xpos(Math.min(x1,xr[1]));
  const c=ypos(Math.min(y1,yr[1])),d=ypos(Math.max(y0,yr[0]));
  if(![a,b,c,d].every(Number.isFinite)||d<c){hidden+=n;return;}
  const r=svgNode('rect',{x:a,y:c,width:Math.max(0,b-a),height:Math.max(0,d-c),fill:color});r.append(svgNode('title',{},title));bars.append(r);visible+=n;
 };
 if(two){
  for(const [i,j,n] of pair.cells){const t=HitMap3D.colorFraction(n,max,$('plot-log-color').checked);
   const color=`rgb(${HitMap3D.viridis(t).map(v=>Math.round(v*255)).join(',')})`;
   rect(h.edges[i],h.edges[i+1],hy.edges[j],hy.edges[j+1],n,color,`${h.label}: ${fmt(h.edges[i])} – ${fmt(h.edges[i+1])}\n${hy.label}: ${fmt(hy.edges[j])} – ${fmt(hy.edges[j+1])}\n${n.toLocaleString()} ${hitData.unit||'hits'}`);
  }
  for(let i=0;i<40;i++)svg.append(svgNode('rect',{x:W-right+12,y:top+ph*(1-(i+1)/40),width:9,height:ph/40+.3,fill:`rgb(${HitMap3D.viridis(i/39).map(v=>Math.round(v*255)).join(',')})`}));
  text(W-right+24,top+8,fmt(max),{'font-size':9});text(W-right+24,top+ph,'0',{'font-size':9});text(W-right+10,top-3,hitData.unit||'Hits',{'font-size':9});
 }else counts.forEach((n,i)=>{if(!n)return;rect(categorical?i:h.edges[i],categorical?i+1:h.edges[i+1],yr[0],n,n,'#17809a',`${categorical?entries[i][0]:fmt(h.edges[i])+' – '+fmt(h.edges[i+1])}: ${n.toLocaleString()} ${hitData.unit||'hits'}`);});
 for(const v of ticks(yr,logY,Math.max(2,Math.min(4,Math.floor(ph/40))))){const y=ypos(v);if($('plot-grid-y').checked)svg.append(svgNode('line',{x1:left,x2:W-right,y1:y,y2:y,stroke:'#a9bdca','stroke-opacity':.35}));text(left-6,y+3,fmt(v),{'text-anchor':'end'});}
 if(categorical){
  const names={'11':'e−','-11':'e+','22':'γ','2212':'p','2112':'n'},stride=Math.max(1,Math.ceil(entries.length/(pw/46)));
  entries.forEach(([id],i)=>{if(i%stride===0)text(xpos(i+.5),top+ph+14,kind==='particles'?names[id]||id:id,{'text-anchor':'middle','font-size':9});});
 }else for(const v of ticks(xr,logX,Math.max(2,Math.min(5,Math.floor(pw/70))))){const x=xpos(v);if($('plot-grid-x').checked)svg.append(svgNode('line',{x1:x,x2:x,y1:top,y2:top+ph,stroke:'#a9bdca','stroke-opacity':.35}));text(x,top+ph+14,fmt(v),{'text-anchor':'middle'});}
 svg.append(svgNode('path',{d:`M ${left} ${top} V ${top+ph} H ${W-right}`,fill:'none',stroke:'#7a8c99'}));
 if($('plot-title').value.trim())text(left+pw/2,18,$('plot-title').value.trim(),{'text-anchor':'middle','font-size':14,'font-weight':600});
 text(left+pw/2,H-6,$('plot-x-label').value.trim()||(categorical?plotNames[kind]:h.label),{'text-anchor':'middle'});
 text(12,top+ph/2,$('plot-y-label').value.trim()||(two?hy.label:(hitData.unit||'Hits')+' per bin'),{transform:`rotate(-90 12 ${top+ph/2})`,'text-anchor':'middle'});
 if(!counts.some(n=>n>0))text(left+pw/2,top+ph/2,'No hits in stored range',{'text-anchor':'middle'});
 const invalid=categorical?(kind==='particles'?hitData.invalid_particles:0):two?pair.invalid:h.invalid;
 const excluded=categorical?0:two?pair.outside:h.underflow+h.overflow;
 const excludedUnit=two?'hits':hitData.unit||'hits';
 const scope=h?.event_value||hy?.event_value?' · Event value per selected hit':'';
 $('plot-note').textContent=(hitData.unit==='nA'?'Cathode current · 65 µA':hitData.weighting==='pe'?'PE rate · 65 µA':'Unweighted hits')+scope+(two?' · '+h.expression+' × '+hy.expression:categorical?'': ' · '+h.expression)+(excluded?` · ${fmt(excluded)} ${excludedUnit} outside stored range`:'')+(invalid?` · ${invalid.toLocaleString()} invalid hits`:'')+(hidden?` · ${fmt(hidden)} ${hitData.unit||'hits'} outside view`:'')+((logX||logY)?' · Log axes':'')+(hitData.weighting==='pe'?' · PE model: main forward e± >1.1 MeV; ShowerMax forward':'');
}
document.querySelectorAll('[data-dimension]').forEach(b=>b.onclick=()=>{
 plotDimension=Number(b.dataset.dimension);document.querySelectorAll('[data-dimension]').forEach(x=>x.setAttribute('aria-pressed',String(x===b)));
 $('plot-kind').hidden=plotDimension===2;$('plot-pair').hidden=plotDimension===1;$('plot-color-option').hidden=plotDimension===1;changePlot();
});
let settingsLayout=null;
function showPlotSettings(open){
 if(open&&!binDraft)syncBinControls();
 // A tiny window cannot split a short plot again without crushing both panes.
 // Temporarily use the existing focus layout; closing settings restores it.
 if(open&&document.querySelector('.plotpanel').getBoundingClientRect().height<300&&focusedPanel!=='distribution'){
  settingsLayout={focused:focusedPanel,visible:{...panelVisible}};
  focusedPanel='distribution';for(const key in panelVisible)panelVisible[key]=key==='distribution';updatePanels();
 }else if(!open&&settingsLayout){
  if(focusedPanel==='distribution'){focusedPanel=settingsLayout.focused;Object.assign(panelVisible,settingsLayout.visible);updatePanels();}
  settingsLayout=null;
 }
 $('plot-settings').hidden=!open;document.querySelector('.plotpanel').classList.toggle('settings-open',open);
 $('plot-settings-toggle').setAttribute('aria-expanded',String(open));
}
$('plot-settings-toggle').onclick=()=>showPlotSettings($('plot-settings').hidden);
$('plot-settings').onkeydown=e=>{if(e.key==='Escape'){showPlotSettings(false);$('plot-settings-toggle').focus();}};
for(const id of [...plotViewInputs,...plotViewChecks]){
 const update=()=>{rememberPlotView();renderPlot();};$(id).oninput=update;$(id).onchange=update;
}
$('axis-auto').onclick=()=>{for(const id of plotViewInputs.slice(0,4))$(id).value='';rememberPlotView();renderPlot();};
$('plot-pair').onchange=changePlot;$('plot-reset').onclick=()=>{resetPlot();rememberPlotView();};
$('bin-axes').oninput=()=>{binDraft=true;};$('bin-axes').onchange=()=>{binDraft=true;};
function selectedPlotAxes(){
 const pair=hitData?.histograms2d?.[$('plot-pair').value];
 return plotDimension===2?[['X',pair?.x],['Y',pair?.y]]:[['X',$('plot-kind').value]];
}
function syncBinControls(){
 if(binDraft)return;
 $('bin-axes').replaceChildren();$('bin-status').textContent='';
 let count=0;
 for(const [axis,key] of selectedPlotAxes()){
  const h=hitData?.histograms?.[key];if(!h?.available)continue;count++;
  const field=document.createElement('fieldset');field.className='bin-axis';field.dataset.key=key;
  const legend=document.createElement('legend');legend.textContent=axis+' · '+h.label;field.append(legend);
  const grid=document.createElement('div');grid.className='bin-fields';field.append(grid);
  const values={min:h.edges[0],max:h.edges.at(-1),bins:h.counts.length,width:(h.edges.at(-1)-h.edges[0])/h.counts.length,spacing:h.log_bins?'log':'linear'};
  for(const [name,label] of [['min','Minimum'],['max','Maximum'],['bins','Number of bins'],['width','Bin width'],['spacing','Spacing']]){
   const wrap=document.createElement('label');wrap.textContent=label;
   const input=document.createElement(name==='spacing'?'select':'input');input.dataset.bin=name;input.setAttribute('aria-label',axis+' '+(['min','max'].includes(name)?'bin ':'')+label.toLowerCase());
   if(name==='spacing')for(const [v,t] of [['linear','Equal width'],['log','Logarithmic']]){const o=document.createElement('option');o.value=v;o.textContent=t;input.append(o);}
   else {input.type='number';input.step=name==='bins'?'1':'any';if(name==='bins'){input.min='1';input.max='256';}}
   input.value=values[name];wrap.append(input);grid.append(wrap);
  }
  const get=name=>grid.querySelector(`[data-bin="${name}"]`);
  const updateWidth=()=>{get('width').disabled=get('spacing').value==='log';get('width').value=get('width').disabled?'':Number(((Number(get('max').value)-Number(get('min').value))/Number(get('bins').value)).toPrecision(10));};
  for(const name of ['min','max','bins','spacing'])get(name).onchange=updateWidth;
  get('width').onchange=()=>{const width=Number(get('width').value),low=Number(get('min').value),high=Number(get('max').value);if(width>0&&high>low){const n=Math.ceil((high-low)/width-1e-10);get('bins').value=n;get('max').value=Number((low+n*width).toPrecision(12));$('bin-status').textContent='Maximum aligned to a whole number of bins.';}};
  updateWidth();$('bin-axes').append(field);
 }
 $('bin-help').textContent=count?'Changing bins rereads the applied selection in ROOT.':'No numeric histogram available.';
 $('apply-bins').disabled=$('reset-bins').disabled=!count||busy||requesting;
}
async function applyBins(reset=false){
 if(!applied||busy||requesting)return;
 const next={...plotBinning};
 try{
  for(const field of $('bin-axes').children){
   const key=field.dataset.key;if(reset){delete next[key];continue;}
   const read=name=>field.querySelector(`[data-bin="${name}"]`);
   if(['min','max','bins'].some(n=>read(n).value===''))throw Error('Set the limits and number of bins.');
   if(read('spacing').value==='linear'&&(!Number.isFinite(Number(read('width').value))||Number(read('width').value)<=0))throw Error('Bin width must be positive.');
   const v={min:Number(read('min').value),max:Number(read('max').value),bins:Number(read('bins').value),spacing:read('spacing').value};
   if(!Number.isInteger(v.bins)||v.bins<1||v.bins>256)throw Error('Use 1–256 bins per axis.');
   if(!Number.isFinite(v.min)||!Number.isFinite(v.max)||Math.abs(v.min)>1e15||Math.abs(v.max)>1e15||v.min>=v.max||(v.spacing==='log'&&v.min<=0))throw Error('Use increasing limits; log bins need positive limits.');
   next[key]=v;
  }
  $('bin-status').textContent='Reading ROOT hits…';
  await loadHits({rebin:true,binning:next});
 }catch(error){$('bin-status').textContent=error.message;}
}
$('apply-bins').onclick=()=>applyBins();$('reset-bins').onclick=()=>applyBins(true);
$('plot-export').onclick=()=>{if($('plot').hasAttribute('hidden'))return;const copy=$('plot').cloneNode(true);copy.prepend(svgNode('desc',{},$('plot-note').textContent+'; ROOT cut: '+(applied?.cut||'none')));copy.setAttribute('xmlns',svgNS);copy.setAttribute('width',$('plot').viewBox.baseVal.width);copy.setAttribute('height',$('plot').viewBox.baseVal.height);const url=URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(copy)],{type:'image/svg+xml'}));const a=document.createElement('a');a.href=url;a.download='distribution.svg';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};

async function api(path,data){const response=await fetch(path,data===undefined?{}:{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(data)});const result=await response.json();if(!response.ok)throw Error(result.error||'Request failed');return result;}
function lockControls(){
 const locked=busy||choosing||requesting;
 for(const id of ['load','choose','load-hits','choose-hits','apply-cut','clear-cut','tile-clear','system','region','ancestry','particle','apply-bins','reset-bins'])$(id).disabled=locked;
 document.querySelectorAll('[data-mode],[data-direction],[data-signal]').forEach(b=>b.disabled=locked);
 $('apply-bins').disabled=$('reset-bins').disabled=locked||!applied||!$('bin-axes').children.length;
}
async function loadHits(options={}){
 if(busy||requesting)return false;
 const previousData=hitData,previousApplied=applied;
 const paths=options.rebin?applied.path.split('\n'):filePaths(),cut=options.rebin?applied.cut:$('hit-cut').value,binning=options.binning||plotBinning;
 requesting=true;lockControls();hitData=null;applied=null;showApplied();setGeometry();renderPlot();$('hit-status').textContent='Reading selected hits…';
 try{await api('/load-hits',{paths,cut,binning,weighting:signalSelection==='hits'?'hits':'pe'});rawHitData=null;plotBinning=binning;binDraft=false;busy=true;return true;}
 catch(error){hitData=previousData;signalSelection=hitData?.weighting==='pe'?(hitData.unit==='nA'?'current':'pe'):'hits';syncSignalControls();applied=previousApplied;showApplied();setGeometry();renderPlot();$('hit-status').textContent=error.message;$('bin-status').textContent=error.message;return false;}
 finally{requesting=false;lockControls();}
}
async function poll(){
 if(stopped)return;
 try{
  const state=await api('/status');busy=state.busy;lockControls();document.dispatchEvent(new CustomEvent('workspace-status',{detail:state}));
  if(!initialized){
   signalSelection=state.hits.weighting==='pe'?'pe':'hits';syncSignalControls();
   $('root').value=state.default_root;$('hit-file').value=(state.hits.paths||[state.hits.path||'']).join('\n');
   plotBinning=state.hits.binning||{};
   if(state.hits.path){$('hit-cut').value=state.hits.cut||'1';editing=true;try{const saved=JSON.parse(sessionStorage.getItem('detector-selection'));if(saved?.cut===state.hits.cut){config=saved.config;editing=false;}}catch{} }
   else $('hit-cut').value=buildCut(config);
   syncControls();initialized=true;
  }
  $('status').textContent=state.message;
  if(state.ready&&state.revision!==revision){const geometry=await api('/geometry');geometryFaces=geometry.faces;document.dispatchEvent(new CustomEvent('workspace-geometry',{detail:geometry}));revision=state.revision;setGeometry();surface.fit();}
  const key=JSON.stringify(state.hits);
  if(key!==hitKey){
   hitKey=key;$('hit-status').textContent=state.hits.message;
   applied=state.hits.result?{cut:state.hits.cut||'1',path:(state.hits.paths||[state.hits.path]).join('\n')}:null;showApplied();
   if(!state.hits.result){rawHitData=null;hitData=null;setGeometry();renderPlot();}
   else{
    const full=await api('/hits');if(key===hitKey){rawHitData=full;displaySignalData();plotBinning=state.hits.binning||{};setGeometry();syncBinControls();renderPlot();}
    $('hit-details').hidden=false;$('hit-fields').replaceChildren();$('hit-counts').replaceChildren();
    for(const field of state.hits.result.fields){const b=document.createElement('button');b.textContent=field;b.onclick=()=>{const input=$('hit-cut');if(input.value.trim()==='1')input.value='';input.focus();input.setRangeText(field,input.selectionStart,input.selectionEnd,'end');input.dispatchEvent(new Event('input'));};$('hit-fields').append(b);}
    for(const [det,count] of Object.entries(state.hits.result.detectors)){const b=document.createElement('button');b.textContent=`${det}: ${count.toLocaleString()}`;b.onclick=()=>selectTile(Number(det));$('hit-counts').append(b);}
   }
  }
 }catch(error){$('status').textContent=error.message;}
 if(!stopped)setTimeout(poll,500);
}
for(const id of ['system','region','ancestry','particle'])$(id).onchange=()=>{
 clearTile();config[id]=$(id).value;
 if(id==='system'){config.direction=config.system==='showermax'&&config.mode!=='raw'?'forward':'all';if(config.mode!=='raw'&&config.system!=='showermax'&&['22','2112'].includes(config.particle))config.particle='all';}
 selectAndRun();
};
document.querySelectorAll('[data-direction]').forEach(b=>b.onclick=()=>{config.direction=b.dataset.direction;selectAndRun();});
document.querySelectorAll('[data-mode]').forEach(b=>b.onclick=()=>{
 config.mode=b.dataset.mode;clearTile();config.ancestry='all';config.particle='all';
 if(config.mode==='raw'){config.system='all';config.direction='all';config.region='all';}
 else if(config.system==='all')config.system='main';
 $('plot-kind').value=config.mode==='secondary'?'energy':'radius';changePlot();selectAndRun();
});
$('tile-clear').onclick=()=>{clearTile();selectAndRun();};
$('plot-kind').onchange=changePlot;
$('hit-cut').oninput=()=>{editing=true;if(config.tile!==null)clearTile();syncControls();showApplied();$('hit-status').textContent='Cut edited. Apply to update both views.';};
$('hit-cut').onkeydown=e=>{if(e.key==='Enter'&&(e.ctrlKey||e.metaKey)){e.preventDefault();loadHits();}};
$('apply-cut').onclick=$('load-hits').onclick=loadHits;
$('clear-cut').onclick=()=>{config={mode:'raw',system:'all',region:'all',ancestry:'all',direction:'all',particle:'all',tile:null};clearTile();selectAndRun();};
$('hit-file').oninput=showApplied;
$('choose-hits').onclick=async()=>{choosing=true;lockControls();try{const result=await api('/choose-hits',{});if(result.paths?.length){$('hit-file').value=[...new Set([...filePaths(),...result.paths])].join('\n');$('hit-file').dispatchEvent(new Event('input'));$('hit-status').textContent=filePaths().length+' files selected · Click Load';}}catch(error){$('hit-status').textContent=error.message;}finally{choosing=false;lockControls();}};
$('choose').onclick=async()=>{try{const result=await api('/choose',{});if(result.path)$('root').value=result.path;}catch(error){$('status').textContent=error.message;}};
$('load').onclick=async()=>{try{await api('/load',{root:$('root').value});geometryFaces=[];hitData=null;setGeometry();renderPlot();}catch(error){$('status').textContent=error.message;}};
$('quit').onclick=async()=>{try{await api('/quit',{});stopped=true;$('hit-status').textContent='Stopped. Relaunch with python3 geometry/run.py.';document.querySelectorAll('button').forEach(b=>b.disabled=true);}catch(error){$('hit-status').textContent=error.message;}};
// Layout controls only: no scan, selection or weighting changes.
const workspace=document.querySelector('.workspace'),divider=$('view-divider'),cutsPanel=$('cuts-panel');
const panelVisible={scene:true,selection:true,distribution:true};
let focusedPanel=null;
function updatePanels(){
 for(const [name,visible] of Object.entries(panelVisible))workspace.classList.toggle('hide-'+name,!visible);
 workspace.classList.toggle('focus-scene',focusedPanel==='scene');workspace.classList.toggle('focus-distribution',focusedPanel==='distribution');
 document.querySelectorAll('[data-panel]').forEach(b=>b.setAttribute('aria-pressed',String(panelVisible[b.dataset.panel])));
 document.querySelectorAll('[data-panel-focus]').forEach(b=>{const active=focusedPanel===b.dataset.panelFocus;b.textContent=active?'↙':'⛶';b.setAttribute('aria-label',active?'Restore panels':b.dataset.panelFocus==='scene'?'Expand 3D viewer':'Expand Distribution');b.title=b.getAttribute('aria-label');});
 requestAnimationFrame(()=>{renderPlot();surface.draw();});
}
function restorePanels(){for(const key in panelVisible)panelVisible[key]=true;focusedPanel=null;updatePanels();}
function togglePanel(name){
 focusedPanel=null;panelVisible[name]=!panelVisible[name];
 if(!Object.values(panelVisible).some(Boolean))panelVisible[name==='distribution'?'scene':'distribution']=true;
 updatePanels();
}
document.querySelectorAll('[data-panel]').forEach(b=>b.onclick=()=>togglePanel(b.dataset.panel));
document.querySelectorAll('[data-minimize]').forEach(b=>b.onclick=()=>togglePanel(b.dataset.minimize));
document.querySelectorAll('[data-panel-focus]').forEach(b=>b.onclick=()=>{if(focusedPanel===b.dataset.panelFocus){restorePanels();return;}focusedPanel=b.dataset.panelFocus;for(const key in panelVisible)panelVisible[key]=key===focusedPanel;updatePanels();});
$('restore-panels').onclick=restorePanels;
let sceneShare=matchMedia('(max-width:850px)').matches?50:58,resizeFrame=0;
function setSplit(value){
 sceneShare=Math.max(30,Math.min(72,value));
 workspace.style.setProperty('--scene-share',`${sceneShare}fr`);
 workspace.style.setProperty('--plot-share',`${100-sceneShare}fr`);
 divider.setAttribute('aria-valuenow',String(Math.round(sceneShare)));
 divider.setAttribute('aria-valuetext',`${Math.round(sceneShare)} percent detector view`);
}
setSplit(sceneShare);
divider.onpointerdown=e=>{if(e.button!==0)return;e.preventDefault();divider.setPointerCapture(e.pointerId);document.body.classList.add('resizing');};
divider.onpointermove=e=>{
 if(!divider.hasPointerCapture(e.pointerId))return;
 const scene=workspace.querySelector('.scene').getBoundingClientRect(),analysis=workspace.querySelector('.analysis').getBoundingClientRect();
 if(matchMedia('(max-width:850px)').matches){const plot=document.querySelector('.plotpanel').getBoundingClientRect();setSplit(100*(e.clientY-scene.top-divider.clientHeight/2)/(scene.height+plot.height));}
 else setSplit(100*(e.clientX-scene.left-divider.clientWidth/2)/(scene.width+analysis.width));
};
function endResize(e){if(divider.hasPointerCapture(e.pointerId))divider.releasePointerCapture(e.pointerId);document.body.classList.remove('resizing');}
divider.onpointerup=endResize;divider.onpointercancel=endResize;divider.onlostpointercapture=()=>document.body.classList.remove('resizing');
divider.ondblclick=()=>setSplit(matchMedia('(max-width:850px)').matches?50:58);
divider.onkeydown=e=>{if(!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','Home'].includes(e.key))return;e.preventDefault();setSplit(e.key==='Home'?(matchMedia('(max-width:850px)').matches?50:58):sceneShare+(['ArrowRight','ArrowDown'].includes(e.key)?2:-2));};
// Vertical split trades selection-panel space for plot height; no data reload.
const heightDivider=$('plot-height-divider'),analysisPanel=document.querySelector('.analysis');
let selectionShare=35;
try{const saved=Number(sessionStorage.getItem('distribution-height-share'));if(saved>=15&&saved<=70)selectionShare=saved;}catch{}
function verticalSpace(){return Math.max(0,analysisPanel.clientHeight-analysisPanel.querySelector('.modes').getBoundingClientRect().height-24);}
function setPlotHeight(value){
 selectionShare=Math.max(15,Math.min(70,value));
 if(!matchMedia('(max-width:850px)').matches){
  const available=verticalSpace(),selection=Math.max(0,Math.min(available*selectionShare/100,available-140));
  analysisPanel.style.setProperty('--selection-height',`${selection}px`);
 }
 heightDivider.setAttribute('aria-valuenow',String(Math.round(selectionShare)));
 heightDivider.setAttribute('aria-valuetext',`${Math.round(100-selectionShare)} percent distribution height`);
 try{sessionStorage.setItem('distribution-height-share',String(selectionShare));}catch{}
}
heightDivider.onpointerdown=e=>{if(e.button!==0)return;e.preventDefault();heightDivider.setPointerCapture(e.pointerId);document.body.classList.add('resizing-height');};
heightDivider.onpointermove=e=>{
 if(!heightDivider.hasPointerCapture(e.pointerId))return;
 const selection=analysisPanel.querySelector('.subcuts').getBoundingClientRect();
 setPlotHeight(100*(e.clientY-selection.top-10)/Math.max(1,verticalSpace()));
};
function endHeightResize(e){if(heightDivider.hasPointerCapture(e.pointerId))heightDivider.releasePointerCapture(e.pointerId);document.body.classList.remove('resizing-height');}
heightDivider.onpointerup=endHeightResize;heightDivider.onpointercancel=endHeightResize;
heightDivider.onlostpointercapture=()=>document.body.classList.remove('resizing-height');
heightDivider.ondblclick=()=>setPlotHeight(35);
heightDivider.onkeydown=e=>{if(!['ArrowUp','ArrowDown','Home'].includes(e.key))return;e.preventDefault();setPlotHeight(e.key==='Home'?35:selectionShare+(e.key==='ArrowUp'?-3:3));};
new ResizeObserver(()=>setPlotHeight(selectionShare)).observe(analysisPanel);
new ResizeObserver(()=>setPlotHeight(selectionShare)).observe(analysisPanel.querySelector('.modes'));
setPlotHeight(selectionShare);
function setCutsOpen(open){
 cutsPanel.classList.toggle('is-collapsed',!open);workspace.classList.toggle('cuts-open',open);
 $('shared-cuts').setAttribute('aria-expanded',String(open));
 $('toggle-cuts').setAttribute('aria-expanded',String(open));$('toggle-cuts').setAttribute('aria-label',open?'Hide ROOT cuts':'Show ROOT cuts');
 $('toggle-cuts').title=open?'Hide ROOT cuts':'Show ROOT cuts';$('toggle-cuts').innerHTML=open?'Hide <span aria-hidden="true">−</span>':'Cuts <span aria-hidden="true">＋</span>';
}
$('toggle-cuts').onclick=$('shared-cuts').onclick=()=>setCutsOpen(cutsPanel.classList.contains('is-collapsed'));
new ResizeObserver(()=>{const narrow=matchMedia('(max-width:850px)').matches;divider.setAttribute('aria-orientation',narrow?'horizontal':'vertical');cancelAnimationFrame(resizeFrame);resizeFrame=requestAnimationFrame(renderPlot);}).observe($('plot').parentElement);
poll();

function syncSignalControls(){document.querySelectorAll('[data-signal]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.signal===signalSelection)));}
function displaySignalData(){hitData=HitMap3D.signalData(rawHitData,signalSelection);}
function switchSignalAxes(previous,next){
 rememberPlotView();
 for(const [key,view] of plotViews)if(key.startsWith('1:'))view.values=signalAxisValues(view.values,previous,next);
 const active=plotViews.get(activePlotKey);if(active)for(const key of ['plot-y-min','plot-y-max'])$(key).value=active.values[key];
}

document.querySelectorAll('[data-signal]').forEach(b=>b.onclick=()=>{
 const previous=signalSelection;if(previous===b.dataset.signal)return;
 if(!filePaths().length&&!rawHitData){$('hit-status').textContent='Choose ROOT files.';return;}
 if(!filePaths().length&&(previous==='hits'||b.dataset.signal==='hits')){$('hit-status').textContent='Choose ROOT files.';return;}
 signalSelection=b.dataset.signal;syncSignalControls();
 if(previous!=='hits'&&signalSelection!=='hits'&&rawHitData?.weighting==='pe'){switchSignalAxes(previous,signalSelection);displaySignalData();setGeometry();renderPlot();}
 else if(filePaths().length)loadHits().then(ok=>{if(ok)switchSignalAxes(previous,signalSelection);});else $('hit-status').textContent='Choose ROOT files.';
});
document.addEventListener('workspace-page',event=>{if(event.detail==='detector'){surface.fit();renderPlot();}});
// Same selected response files feed the hit weights and ROOT deconvolution.

function responseFields(){
 $('response-local').hidden=$('response-mode').value!=='local';
 document.querySelectorAll('[data-response]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.response===$('response-mode').value)));
 if(!$('response-remoll').value)$('response-remoll').value=$('root').value;
}
function showResponse(data){
 const r=data.resources;
 $('response-status').textContent=r.mode==='bundled'
  ? 'Bundled · MOLLER DocDB · July 14, 2026 · MainDet_ResponseFunction.zip. ShowerMax: remoll develop.'
  : 'Local · 12 main tables, 8 ShowerMax tables and main fit code verified.';
 $('response-status').title='SHA256: '+r.fingerprint+' · '+r.evaluator;
}
$('response-mode').onchange=()=>{responseFields();$('response-status').textContent='Selection not applied.';};
for(const id of ['response-remoll','response-main','response-sm'])$(id).oninput=()=>{$('response-status').textContent='Selection not applied.';};
document.querySelectorAll('[data-response-settings]').forEach(b=>b.onclick=()=>{responseFields();$('response-options').showModal();});
$('response-close').onclick=()=>$('response-options').close();
document.querySelectorAll('[data-response]').forEach(b=>b.onclick=()=>{$('response-mode').value=b.dataset.response;responseFields();$('response-status').textContent='Selection not applied.';});
$('response-save').onclick=async()=>{
 const button=$('response-save');button.disabled=true;
 try{const data=await api('/responses',{mode:$('response-mode').value,remoll:$('response-remoll').value,main_tables:$('response-main').value,showermax_tables:$('response-sm').value});showResponse(data);document.dispatchEvent(new Event('response-changed'));if(signalSelection!=='hits'&&filePaths().length)loadHits();}
 catch(error){$('response-status').textContent=error.message;}
 finally{button.disabled=false;}
};
api('/responses').then(data=>{
 $('response-mode').value=data.config.mode;
 $('response-remoll').value=data.config.remoll||'';
 $('response-main').value=data.config.main_tables||'';
 $('response-sm').value=data.config.showermax_tables||'';
 responseFields();showResponse(data);
}).catch(error=>{$('response-status').textContent=error.message;});

// Capture the displayed result, not un-applied edits in the cut/input fields.
document.addEventListener('save-current-view', event=>{
 if(!$('fit-page').hidden)return;
 if(!hitData||!applied||busy||requesting){event.detail.error='Load a view first';return;}
 const pair=hitData.histograms2d?.[$('plot-pair').value];
 const keys=plotDimension===2?[pair?.x,pair?.y]:[$('plot-kind').value];
 const histograms=Object.fromEntries(keys.filter(k=>hitData.histograms?.[k]).map(k=>[k,hitData.histograms[k]]));
 event.detail.snapshot={page:'detector',selection:{paths:applied.path.split('\n'),cut:applied.cut,weighting:signalSelection},
  plot:{dimension:plotDimension,quantity:$(plotDimension===2?'plot-pair':'plot-kind').value,
   values:Object.fromEntries(plotViewInputs.map(k=>[k,$(k).value])),checks:Object.fromEntries(plotViewChecks.map(k=>[k,$(k).checked]))},
  camera:{yaw:surface.yaw,pitch:surface.pitch,zoom:surface.zoom,pan:surface.pan,unwrapped:surface.unwrapped,selectedDetector:surface.selectedDetector},
  geometry:{faces:geometryFaces},
  result:{...hitData,histograms,histograms2d:plotDimension===2&&pair?{[$('plot-pair').value]:pair}:{}},
  layout:{panels:panelVisible,focus:focusedPanel}};
});

}
