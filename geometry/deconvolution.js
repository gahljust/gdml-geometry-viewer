/* Adapted from the published campaign correlation renderer (SHA256
 * 2c9a1e9a62e82a7b8737c1cb09a997cd179b6d61ab346198bcda7b076b046fe7).
 * Scoped to this panel; inputs come from the current ROOT result and loaded GDML.
 * Bird palette, sqrt display scaling, picking, focus and unwrapping are retained.
 */
function createDetectorCorrelations(root,D){
 const events=new AbortController();
const $=id=>root.querySelector('#'+id),C=D.C,rows=D.rows,cv=$('scene'),ctx=cv.getContext('2d');
let focusGroup='all',correlationActive=true;
const treatments={neutral:true,flat:false,outlines:false,glow:true,unwrapped:false};
let selected=Math.max(0,rows.findIndex(r=>r.detector===150252)),mode='correlation',yaw=.93,pitch=.32,zoom=1,pan=[0,0],W=1,H=1,hitFaces=[],hover=-1,drag=null,maxValue=1,values=[],partners=[];
const label=r=>r.system==='showermax'?`ShowerMax ${r.detector}`:`Ring ${r.system.slice(4)} · ${r.detector}`;
const group=r=>r.system==='showermax'?'ShowerMax':`Ring ${r.system.slice(4)}`;
const world=p=>[p[0]/1000,p[1]/1000,(p[2]-22850)/1000];
const centers=Array.from({length:rows.length},()=>[0,0,0]),counts=Array(rows.length).fill(0),faces=[];
for(const m of D.meshes){const center=world(m.center);centers[m.index]=centers[m.index].map((v,i)=>v+center[i]);counts[m.index]++;for(const f of m.faces){const p=f.map(world),a=p[1].map((v,i)=>v-p[0][i]),b=p[2].map((v,i)=>v-p[0][i]),n=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],len=Math.hypot(...n)||1;faces.push({index:m.index,p,n:n.map(x=>x/len)})}}
centers.forEach((p,i)=>centers[i]=p.map(v=>v/counts[i]));
// Radial and longitudinal envelopes of the existing geometry, in metres.
const bounds=[...new Set(rows.map(r=>r.system))].map(system=>{const points=faces.filter(f=>rows[f.index].system===system).flatMap(f=>f.p);return {system,r0:Math.min(...points.map(p=>Math.hypot(p[0],p[1]))),r1:Math.max(...points.map(p=>Math.hypot(p[0],p[1]))),z0:Math.min(...points.map(p=>p[2])),z1:Math.max(...points.map(p=>p[2]))}});
function transformed(p,system){return p;}
function envelope(b,z,r){const points=[];for(let i=0;i<=120;i++){const a=i*2*Math.PI/120;points.push(project(transformed([r*Math.cos(a),r*Math.sin(a),z],b.system)))}return points}
function drawBounds(){const system=focusGroup==='all'?(treatments.outlines&&correlationActive?rows[selected].system:null):focusGroup;if(!system)return;const b=bounds.find(b=>b.system===system);ctx.fillStyle='#bcecff';ctx.globalAlpha=.055;for(const z of [b.z0-.012,b.z1+.012]){ctx.beginPath();for(const r of [b.r1+.015,Math.max(0,b.r0-.015)]){envelope(b,z,r).forEach((p,i)=>i?ctx.lineTo(p[0],p[1]):ctx.moveTo(p[0],p[1]));ctx.closePath()}ctx.fill('evenodd')}ctx.globalAlpha=1;if(treatments.outlines){ctx.strokeStyle='#bcecff66';ctx.lineWidth=1;for(const r of [b.r0,b.r1]){polygon(envelope(b,(b.z0+b.z1)/2,r));ctx.stroke()}}}
function strength(i){return Math.sqrt(Math.max(0,Math.min(1,values[i]/maxValue)))}
function displayColor(t,light=1){if(!treatments.neutral)return color(t,light);const c=bird[Math.round(t*(bird.length-1))],mix=Math.min(1,t*3);return `rgb(${c.map((v,j)=>Math.round(([79,91,105][j]*(1-mix)+v*mix)*light)).join(',')})`}
function tileColor(i,light){return displayColor(strength(i),light)}
function updateLegend(){root.querySelector('.ramp').style.background='linear-gradient(90deg,'+bird.map((c,i)=>`${displayColor(i/(bird.length-1))} ${100*(i/(bird.length-1))**2}%`).join(',')+')'}
function paintTile(f,light=1){const sel=correlationActive&&f.index===selected,high=f.index===hover,t=strength(f.index);ctx.globalAlpha=focusGroup==='all'||rows[f.index].system===focusGroup?1:.25+.75*t;polygon(f.screen);ctx.fillStyle=tileColor(f.index,light);if(treatments.glow&&t>.15){ctx.shadowColor=color(t);ctx.shadowBlur=8*t}ctx.fill();ctx.shadowBlur=0;ctx.strokeStyle=sel?'#bcf2ff':high?'#ffffff':'#90b4c13d';ctx.lineWidth=sel?(treatments.outlines?3:1.6):high?1.8:.55;ctx.stroke();if(focusGroup==='all'||rows[f.index].system===focusGroup)hitFaces.push(f)}
// Unroll by moving and rotating whole meshes, retaining their local geometry.
// Ring spacing is for presentation; lengths and within-ring z offsets remain metres.
let unwrapAmount=0,unwrapFrame=0,unwrapMotion=null,unwrapLayout=null;
const TAU=2*Math.PI,unwrapPitch=.60;
const wrapAngle=a=>((a+Math.PI)%TAU+TAU)%TAU-Math.PI;
const stripGroups=bounds.map((b,j)=>({...b,j,radius:(b.r0+b.r1)/2,zmid:(b.z0+b.z1)/2}));
let stripHeight=stripGroups.reduce((s,b)=>s+b.r1-b.r0+.35,0)-.35,stripY=stripHeight/2;
for(const b of stripGroups){b.y=stripY-(b.r1-b.r0)/2;stripY-=b.r1-b.r0+.35;}
// Cut halfway between neighboring channel centers near phi=0.
// Each ring uses its own circumference; no seam-duplicate detectors.
const ring5Angles=centers.filter((_,i)=>rows[i].system==='ring5').map(p=>(Math.atan2(p[1],p[0])+TAU)%TAU).sort((a,b)=>a-b);
const unwrapSeam=(ring5Angles[0]+ring5Angles.at(-1)-TAU)/2;
const tileUnwrap=centers.map((c,i)=>{
 const b=stripGroups.find(b=>b.system===rows[i].system),phi=Math.atan2(c[1],c[0]);
 const theta=((phi-unwrapSeam)%TAU+TAU)%TAU;
 return {center:c,angle:wrapAngle(-phi-Math.PI/2),group:b,
   target:[b.radius*(theta-Math.PI),b.y-(Math.hypot(c[0],c[1])-b.radius),c[2]-b.zmid]};
});
function unfoldPoint(p,index,t=unwrapAmount){
 const u=tileUnwrap[index],start=transformed(u.center,rows[index].system);
 const scale=1,a=u.angle*t,c=Math.cos(a),s=Math.sin(a);
 const d=p.map((v,j)=>(v-u.center[j])*scale);
 return [c*d[0]-s*d[1],s*d[0]+c*d[1],d[2]].map((v,j)=>v+(1-t)*start[j]+t*u.target[j]);
}
function unfoldedNormal(n,index){const a=tileUnwrap[index].angle*unwrapAmount,c=Math.cos(a),s=Math.sin(a);return[c*n[0]-s*n[1],s*n[0]+c*n[1],n[2]];}
function stripCamera([x,y,z]){const c=Math.cos(unwrapPitch),s=Math.sin(unwrapPitch);return[x,c*y-s*z,s*y+c*z];}
function fitUnwrap(){
 const points=faces.flatMap(f=>f.p.map(p=>stripCamera(unfoldPoint(p,f.index,1))));
 const xs=points.map(p=>p[0]),ys=points.map(p=>p[1]);
 const x0=Math.min(...xs),x1=Math.max(...xs),y0=Math.min(...ys),y1=Math.max(...ys);
 unwrapLayout={scale:Math.min(Math.max(30,W-Math.min(100,W*.20))/(x1-x0),Math.max(30,H-Math.min(70,H*.25))/(y1-y0)),cx:(x0+x1)/2,cy:(y0+y1)/2};
}
function animateUnwrap(on){
 if(unwrapFrame)cancelAnimationFrame(unwrapFrame);
 unwrapFrame=0;unwrapMotion=null;
 const start=unwrapAmount,target=on?1:0,duration=1300*Math.abs(target-start);
 if(!duration)return;
 unwrapMotion={start,target,time:performance.now(),duration};
 const tick=now=>{const m=unwrapMotion,q=Math.min(1,(now-m.time)/m.duration),ease=q*q*(3-2*q);
  unwrapAmount=m.start+(m.target-m.start)*ease;draw();
  if(q<1)unwrapFrame=requestAnimationFrame(tick);else{unwrapFrame=0;unwrapMotion=null;}
 };unwrapFrame=requestAnimationFrame(tick);
}
function stripAnglePoint(b,phi,y=b.y){
 const theta=((phi-unwrapSeam)%TAU+TAU)%TAU;
 return [b.radius*(theta-Math.PI),y,0];
}
function drawStripAngles(){
 const alpha=Math.max(0,(unwrapAmount-.72)/.28);if(!alpha)return;
 ctx.save();ctx.globalAlpha=alpha;ctx.fillStyle='#9ab4ca';ctx.font='11px system-ui';
 ctx.textAlign='center';ctx.strokeStyle='#809ab15c';ctx.lineWidth=.8;ctx.setLineDash([2,5]);
 const first=stripGroups[0],tick0=project(stripAnglePoint(first,0)),tick1=project(stripAnglePoint(first,TAU/7));
 const labelStride=Math.max(1,Math.ceil(44/Math.max(1,Math.abs(tick1[0]-tick0[0]))));
 for(let septant=0;septant<7;septant++){
  const phi=septant*TAU/7,first=stripGroups[0],last=stripGroups.at(-1);
  const start=project(stripAnglePoint(first,phi,first.y+(first.r1-first.r0)/2+.16));
  ctx.beginPath();ctx.moveTo(start[0],start[1]);
  // Follow the same azimuth across unequal ring circumferences.
  for(const b of stripGroups){const p=project(stripAnglePoint(b,phi));ctx.lineTo(p[0],p[1]);}
  const end=project(stripAnglePoint(last,phi,last.y-(last.r1-last.r0)/2-.08));
  ctx.lineTo(end[0],end[1]);ctx.stroke();
  if(septant%labelStride===0)ctx.fillText(`${(septant*360/7).toFixed(1)}°`,start[0],start[1]-8);
 }
 ctx.restore();
}
function drawStripLabels(){
 const alpha=Math.max(0,(unwrapAmount-.72)/.28);if(!alpha)return;
 ctx.save();ctx.globalAlpha=alpha;ctx.fillStyle='#9ab4ca';ctx.font='12px system-ui';ctx.textAlign='right';
 for(const b of stripGroups){
  const anchor=project([-Math.PI*b.radius-.12,b.y,0]);
  ctx.fillText(b.system==='showermax'?'ShowerMax':`Ring ${b.j+1}`,anchor[0],anchor[1]+4);
 }
 ctx.restore();
}
// Coalesce input bursts into one frame; no background rendering while idle.
let drawFrame=0,projectionKey='',projectedFaces=[],pickGrid=new Map();
const pickCell=48;
function requestDraw(){if(!drawFrame)drawFrame=requestAnimationFrame(()=>{drawFrame=0;draw()});}
function projectedGeometry(){
 const key=[W,H,yaw,pitch,zoom,...pan,unwrapAmount].join(',');
 if(key===projectionKey)return projectedFaces;
 projectionKey=key;
 const t=unwrapAmount,l=unwrapLayout,cy=Math.cos(yaw*(1-t)),sy=Math.sin(yaw*(1-t));
 const cp=Math.cos(pitch*(1-t)+unwrapPitch*t),sp=Math.sin(pitch*(1-t)+unwrapPitch*t);
 const scale=(Math.min(W/4.5,H/3.9)*(1-t)+l.scale*t)*zoom;
 const ox=W*.51*(1-t)+(W+Math.min(65,W*.15))/2*t+pan[0],oy=H*.51*(1-t)+(H-Math.min(25,H*.12))/2*t+pan[1];
 const localScale=1;
 const transforms=tileUnwrap.map((u,i)=>{
  const a=u.angle*t,c=Math.cos(a),s=Math.sin(a),start=transformed(u.center,rows[i].system);
  return {u,c,s,center:start.map((v,j)=>(1-t)*v+t*u.target[j])};
 });
 projectedFaces=faces.map(f=>{
  const {u,c,s,center}=transforms[f.index];
  const screen=f.p.map(p=>{
   const dx=(p[0]-u.center[0])*localScale,dy=(p[1]-u.center[1])*localScale;
   const x=c*dx-s*dy+center[0],y=s*dx+c*dy+center[1],z=(p[2]-u.center[2])*localScale+center[2];
   const rx=cy*x+sy*z,rz=-sy*x+cy*z,ry=cp*y-sp*rz,depth=sp*y+cp*rz;
   const perspective=(1-t)*6/(6-depth)+t;
   return [ox+(rx-l.cx*t)*scale*perspective,oy-(ry-l.cy*t)*scale*perspective,depth];
  });
  const nx=c*f.n[0]-s*f.n[1],ny=s*f.n[0]+c*f.n[1],nz=f.n[2];
  const rx=cy*nx+sy*nz,rz=-sy*nx+cy*nz;
  const shade=Math.abs(rx*.3+(cp*ny-sp*rz)*.5+(sp*ny+cp*rz)*.81);
  const xs=screen.map(p=>p[0]),ys=screen.map(p=>p[1]);
  return {index:f.index,screen,shade,z:screen.reduce((s,p)=>s+p[2],0)/screen.length,
   x0:Math.min(...xs),x1:Math.max(...xs),y0:Math.min(...ys),y1:Math.max(...ys)};
 }).sort((a,b)=>a.z-b.z);
 // Screen-space bins retain painter order for exactly the same frontmost picking.
 pickGrid=new Map();
 for(const f of projectedFaces){
  const x0=Math.max(0,Math.floor(f.x0/pickCell)),x1=Math.min(Math.floor(W/pickCell),Math.floor(f.x1/pickCell));
  const y0=Math.max(0,Math.floor(f.y0/pickCell)),y1=Math.min(Math.floor(H/pickCell),Math.floor(f.y1/pickCell));
  for(let y=y0;y<=y1;y++)for(let x=x0;x<=x1;x++){
   const k=x+','+y;if(!pickGrid.has(k))pickGrid.set(k,[]);pickGrid.get(k).push(f);
  }
 }
 return projectedFaces;
}
function rotate([x,y,z]){const y0=yaw*(1-unwrapAmount),p0=pitch*(1-unwrapAmount)+unwrapPitch*unwrapAmount;
 const a=Math.cos(y0)*x+Math.sin(y0)*z,b=-Math.sin(y0)*x+Math.cos(y0)*z;return[a,Math.cos(p0)*y-Math.sin(p0)*b,Math.sin(p0)*y+Math.cos(p0)*b];}
function project(p){
 const [x,y,z]=rotate(p),t=unwrapAmount,l=unwrapLayout;
 const base=Math.min(W/4.5,H/3.9),scale=(base*(1-t)+l.scale*t)*zoom;
 const perspective=(1-t)*6/(6-z)+t;
 return[W*.51*(1-t)+(W+Math.min(65,W*.15))/2*t+(x-l.cx*t)*scale*perspective+pan[0],
        H*.51*(1-t)+(H-Math.min(25,H*.12))/2*t-(y-l.cy*t)*scale*perspective+pan[1],z];
}
// ROOT kBird, sampled from the installed ROOT palette. Square-root scaling reveals weak correlations.
const bird=[[53,42,135],[52,44,138],[51,46,141],[49,47,143],[48,49,146],[47,50,149],[46,52,152],[45,54,154],[43,55,157],[42,57,160],[41,58,163],[40,60,166],[38,62,168],[37,63,171],[36,65,174],[35,66,177],[33,68,180],[32,69,182],[31,71,185],[30,73,188],[29,74,191],[27,76,194],[26,77,196],[25,79,199],[24,81,202],[22,82,205],[21,84,207],[20,85,210],[19,87,213],[18,89,216],[16,90,219],[15,92,221],[15,93,221],[15,94,221],[16,95,221],[16,96,220],[16,98,220],[16,99,220],[16,100,220],[16,101,220],[16,102,219],[17,103,219],[17,104,219],[17,106,219],[17,107,218],[17,108,218],[17,109,218],[17,110,218],[18,111,217],[18,112,217],[18,114,217],[18,115,217],[18,116,216],[18,117,216],[19,118,216],[19,119,216],[19,121,215],[19,122,215],[19,123,215],[19,124,215],[19,125,215],[20,126,214],[20,127,214],[20,129,214],[19,130,213],[19,131,213],[19,132,213],[18,133,212],[18,134,212],[17,135,212],[17,136,211],[16,137,211],[16,138,210],[16,140,210],[15,141,210],[15,142,209],[14,143,209],[14,144,209],[13,145,208],[13,146,208],[12,147,207],[12,148,207],[12,149,207],[11,151,206],[11,152,206],[10,153,206],[10,154,205],[9,155,205],[9,156,204],[9,157,204],[8,158,204],[8,159,203],[7,160,203],[7,161,203],[6,163,202],[6,164,202],[7,164,201],[8,165,199],[10,165,198],[11,166,197],[12,167,196],[13,167,195],[15,168,194],[16,169,192],[17,169,191],[18,170,190],[20,170,189],[21,171,188],[22,172,186],[23,172,185],[25,173,184],[26,173,183],[27,174,182],[28,175,180],[30,175,179],[31,176,178],[32,176,177],[33,177,176],[35,178,175],[36,178,173],[37,179,172],[38,179,171],[40,180,170],[41,181,169],[42,181,167],[43,182,166],[45,182,165],[46,183,164],[49,183,162],[52,184,161],[54,184,160],[57,184,158],[60,184,157],[63,185,155],[65,185,154],[68,185,153],[71,185,151],[74,186,150],[77,186,148],[79,186,147],[82,186,146],[85,187,144],[88,187,143],[91,187,141],[93,187,140],[96,188,139],[99,188,137],[102,188,136],[105,188,134],[107,189,133],[110,189,132],[113,189,130],[116,189,129],[118,190,127],[121,190,126],[124,190,125],[127,190,123],[130,191,122],[132,191,120],[135,191,119],[137,191,118],[140,191,117],[142,191,116],[144,191,115],[147,190,114],[149,190,113],[151,190,112],[154,190,111],[156,190,111],[158,190,110],[160,190,109],[163,189,108],[165,189,107],[167,189,106],[170,189,105],[172,189,104],[174,189,103],[177,189,102],[179,189,101],[181,188,100],[183,188,99],[186,188,98],[188,188,98],[190,188,97],[193,188,96],[195,188,95],[197,188,94],[200,187,93],[202,187,92],[204,187,91],[206,187,90],[209,187,89],[210,187,88],[212,188,87],[213,188,86],[214,189,84],[216,189,83],[217,189,82],[219,190,81],[220,190,79],[221,191,78],[223,191,77],[224,192,76],[226,192,75],[227,192,73],[228,193,72],[230,193,71],[231,194,70],[233,194,68],[234,195,67],[236,195,66],[237,195,65],[238,196,64],[240,196,62],[241,197,61],[243,197,60],[244,198,59],[245,198,58],[247,198,56],[248,199,55],[250,199,54],[251,200,53],[252,200,51],[254,200,50],[254,202,49],[254,204,48],[253,205,47],[253,207,46],[253,208,44],[253,210,43],[253,211,42],[253,213,41],[253,215,40],[252,216,39],[252,218,38],[252,219,37],[252,221,35],[252,222,34],[252,224,33],[251,226,32],[251,227,31],[251,229,30],[251,230,29],[251,232,27],[251,233,26],[251,235,25],[250,237,24],[250,238,23],[250,240,22],[250,241,21],[250,243,19],[250,244,18],[249,246,17],[249,248,16],[249,249,15]];
function color(t,light=1){const c=bird[Math.round(Math.max(0,Math.min(1,t))*(bird.length-1))];return `rgb(${c.map(v=>Math.round(v*light)).join(',')})`}
root.querySelector('.ramp').style.background='linear-gradient(90deg,'+bird.map((c,i)=>`rgb(${c.join(',')}) ${100*(i/(bird.length-1))**2}%`).join(',')+')';
const value=(a,b)=>mode==='correlation'?D.correlation[a][b]:C[a][b];
const fmt=v=>mode==='correlation'?v.toFixed(3):v.toPrecision(3)+' ppb²';
function update(){values=rows.map((_,i)=>correlationActive?value(selected,i):0);maxValue=mode==='correlation'?1:(Math.max(...values.filter((v,i)=>i!==selected))||1);partners=rows.map((r,i)=>({i,v:values[i]})).filter(p=>p.i!==selected).sort((a,b)=>b.v-a.v).slice(0,8);$('detector').textContent=label(rows[selected]);root.querySelectorAll('.detectorChoice').forEach(b=>b.setAttribute('aria-pressed',String(correlationActive&&+b.dataset.i===selected)));$('selectedName').textContent=label(rows[selected]);$('selectedRegion').textContent=rows[selected].region+' region';$('legendTitle').textContent=mode==='correlation'?'Correlation with selected tile':'Covariance with selected tile';$('legendMax').textContent=fmt(maxValue);$('partnersTitle').textContent=mode==='correlation'?'Strongest correlations':'Largest covariances';$('partners').innerHTML=partners.map(p=>`<button class="partner" data-i="${p.i}" title="${label(rows[p.i])} · ${rows[p.i].region} · ${fmt(p.v)}"><div>${label(rows[p.i])}<span>${rows[p.i].region}</span></div><strong>${fmt(p.v)}</strong></button>`).join('');root.querySelectorAll('.partner').forEach(b=>{b.onclick=()=>choose(+b.dataset.i);b.onpointerenter=()=>{hover=+b.dataset.i;draw()};b.onpointerleave=()=>{hover=-1;draw()}});$('detector').textContent=correlationActive?label(rows[selected]):'Select a detector';$('selectedName').textContent=correlationActive?label(rows[selected]):'Select a tile';$('selectedRegion').textContent=correlationActive?rows[selected].region+' region':'';for(const el of [root.querySelector('.legend'),$('partnersTitle'),$('partners'),root.querySelector('.key')])el.hidden=!correlationActive;draw()}
function choose(i){correlationActive=!(correlationActive&&i===selected);selected=i;hover=-1;update()}
function polygon(p){ctx.beginPath();p.forEach((v,i)=>i?ctx.lineTo(v[0],v[1]):ctx.moveTo(v[0],v[1]));ctx.closePath()}
function draw(){if(drawFrame){cancelAnimationFrame(drawFrame);drawFrame=0;}ctx.clearRect(0,0,W,H);hitFaces=[];if(!unwrapLayout)fitUnwrap();
// A restrained beam axis provides depth without apparatus clutter.
ctx.globalAlpha=1-unwrapAmount;let a=project([0,0,-1.3]),b=project([0,0,1.65]);ctx.strokeStyle='#809ab133';ctx.setLineDash([5,9]);ctx.beginPath();ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.stroke();ctx.setLineDash([]);ctx.globalAlpha=1;
drawStripAngles();
const sorted=projectedGeometry();
for(const f of sorted){const focused=focusGroup===rows[f.index].system,light=treatments.flat?1:(focused?.88:.78)+(focused?.12:.22)*f.shade;paintTile(f,light)}
ctx.globalAlpha=1;if(unwrapAmount===0)drawBounds();drawStripLabels();


ctx.globalAlpha=1-unwrapAmount;ctx.font='12px system-ui';ctx.fillStyle='#9ab4ca';for(const [s,p]of [['Main quartz',[0,-1.3,-.35]],['ShowerMax',[0,-1.3,1.05]]]){const v=project(p);ctx.fillText(s,v[0]-35,v[1]+18)}ctx.globalAlpha=1;
}
function inside(x,y,p){let c=false;for(let i=0,j=p.length-1;i<p.length;j=i++){if(((p[i][1]>y)!==(p[j][1]>y))&&(x<(p[j][0]-p[i][0])*(y-p[i][1])/(p[j][1]-p[i][1])+p[i][0]))c=!c}return c}
function pick(x,y){const candidates=pickGrid.get(Math.floor(x/pickCell)+','+Math.floor(y/pickCell))||[];
 for(let k=candidates.length-1;k>=0;k--){const f=candidates[k];if(focusGroup!=='all'&&rows[f.index].system!==focusGroup)continue;
  if(x>=f.x0&&x<=f.x1&&y>=f.y0&&y<=f.y1&&inside(x,y,f.screen))return f.index;
 }return -1;
}
function coords(e){const r=cv.getBoundingClientRect();return[e.clientX-r.left,e.clientY-r.top]}
cv.onpointerdown=e=>{const p=coords(e);drag={p,start:p,shift:e.shiftKey,dist:0};cv.setPointerCapture(e.pointerId);$('tip').style.display='none'};
cv.onpointermove=e=>{const p=coords(e);if(drag){const dx=p[0]-drag.p[0],dy=p[1]-drag.p[1];drag.dist+=Math.hypot(dx,dy);if(drag.shift||treatments.unwrapped){pan[0]+=dx;pan[1]+=dy}else{yaw+=dx*.006;pitch=Math.max(-1.5,Math.min(1.5,pitch+dy*.006))}drag.p=p;requestDraw();return}const i=pick(...p);cv.style.cursor=i<0?'grab':'pointer';const tip=$('tip');if(i<0){tip.style.display='none';return}const text=label(rows[i])+' · '+rows[i].region+(correlationActive?' · '+(i===selected?'selected; self '+fmt(values[i]):fmt(values[i])):'');if(tip.textContent!==text)tip.textContent=text;tip.style.display='block';tip.style.left='0';tip.style.top='0';tip.style.transform=`translate(${Math.max(0,Math.min(W-295,p[0]+16))}px,${Math.max(0,Math.min(H-70,p[1]+16))}px)`};
cv.onpointerup=e=>{if(drag&&drag.dist<5){const i=pick(...coords(e));if(i>=0)choose(i)}drag=null};cv.onpointercancel=()=>drag=null;cv.onpointerleave=()=>{$('tip').style.display='none'};
cv.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.max(.4,Math.min(4,zoom*Math.exp(-e.deltaY*.001)));requestDraw()},{passive:false});
cv.onkeydown=e=>{if(e.key.startsWith('Arrow')){e.preventDefault();if(treatments.unwrapped){if(e.key==='ArrowLeft')pan[0]-=20;if(e.key==='ArrowRight')pan[0]+=20;if(e.key==='ArrowUp')pan[1]-=20;if(e.key==='ArrowDown')pan[1]+=20;draw();return;}if(e.key==='ArrowLeft')yaw-=.08;if(e.key==='ArrowRight')yaw+=.08;if(e.key==='ArrowUp')pitch-=.08;if(e.key==='ArrowDown')pitch+=.08;draw()}};
for(const system of [...new Set(rows.map(r=>r.system))]){
 const section=document.createElement('section'),heading=document.createElement('h2');heading.textContent=group(rows.find(r=>r.system===system));section.append(heading);
 const members=rows.filter(r=>r.system===system),regions=[...new Set(members.map(r=>r.region))].sort((a,b)=>['open','closed','transition'].indexOf(a)-['open','closed','transition'].indexOf(b));
 for(const region of regions){const h=document.createElement('h3');h.textContent=region[0].toUpperCase()+region.slice(1)+' Region';section.append(h);const grid=document.createElement('div');grid.className='detectorChoices';for(const r of members.filter(r=>r.region===region)){const button=document.createElement('button');button.className='detectorChoice';button.dataset.i=r.index;button.textContent=r.detector;button.setAttribute('aria-label',label(r)+' · '+region);button.onclick=()=>{choose(r.index);$('detectorMenu').open=false;$('detector').focus()};grid.append(button)}section.append(grid)}$('detectorOptions').append(section)
}
document.addEventListener('pointerdown',e=>{if(!$('detectorMenu').contains(e.target))$('detectorMenu').open=false},{signal:events.signal});
$('detectorMenu').addEventListener('keydown',e=>{if(e.key==='Escape'){$('detectorMenu').open=false;$('detector').focus()}});
root.querySelectorAll('[data-mode]').forEach(b=>b.onclick=()=>{mode=b.dataset.mode;root.querySelectorAll('[data-mode]').forEach(x=>x.setAttribute('aria-selected',String(x===b)));update()});
root.querySelectorAll('[data-treatment]').forEach(input=>input.onchange=()=>{treatments[input.dataset.treatment]=input.checked;$('layoutNote').hidden=!treatments.unwrapped;root.querySelector('.hint').innerHTML=treatments.unwrapped?'Click a tile to select · drag to pan<br>Scroll to zoom':'Click a tile to select · drag to rotate<br>Shift-drag to pan · scroll to zoom';if(input.dataset.treatment==='unwrapped')animateUnwrap(input.checked);updateLegend();draw()});
root.querySelectorAll('[data-focus]').forEach(b=>b.onclick=()=>{focusGroup=b.dataset.focus;root.querySelectorAll('[data-focus]').forEach(x=>x.setAttribute('aria-pressed',String(x===b)));draw()});
function view(y,p){yaw=y;pitch=p;zoom=1;pan=[0,0];draw()}
$('reset').onclick=()=>view(.93,.32);$('front').onclick=()=>view(0,0);$('side').onclick=()=>view(Math.PI,0);
const observer=new ResizeObserver(()=>{const r=cv.getBoundingClientRect();W=r.width;H=r.height;fitUnwrap();const dpr=Math.min(2,devicePixelRatio||1);cv.width=W*dpr;cv.height=H*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);draw()});observer.observe(cv);updateLegend();update();


return {snapshot(){return {selectedDetector:rows[selected]?.detector,correlationActive,focusGroup,yaw,pitch,zoom,pan,unwrapped:treatments.unwrapped};},destroy(){events.abort();observer.disconnect();cancelAnimationFrame(drawFrame);cancelAnimationFrame(unwrapFrame);},refresh(){requestDraw();}};
}

/* ROOT result presentation. No fit or covariance calculation in JavaScript. */
(()=>{
const $=id=>document.getElementById(id);
let page='detector',revision=-1,result=null,geometry=null,correlationView=null,renderedResult=null,renderedGeometry=null,matrix='uncertainties',busy=false,inputVersion=0,startedVersion=0;
document.addEventListener('save-current-view', event=>{
 if(page!=='deconvolution')return;
 if(!result?.matrices||busy){event.detail.error='Run deconvolution first';return;}
 event.detail.snapshot={page:'deconvolution',matrix,result,geometry,camera:correlationView?.snapshot()||null};
});
function selectPage(next){
 page=next;document.querySelector('.workspace').hidden=next!=='detector';document.querySelector('.statusrow').hidden=next!=='detector';$('fit-page').hidden=next!=='deconvolution';
 document.querySelectorAll('[data-page]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.page===next)));
 document.dispatchEvent(new CustomEvent('workspace-page',{detail:next}));
 if(next==='deconvolution')requestAnimationFrame(drawCorrelations);
}
document.querySelectorAll('[data-page]').forEach(b=>b.onclick=()=>selectPage(b.dataset.page));
document.addEventListener('workspace-geometry',e=>{geometry=e.detail;drawCorrelations();});
function drawCorrelations(){
 const host=$('fit-map');
 if(!geometry||!result?.matrices){correlationView?.destroy();correlationView=null;renderedResult=null;renderedGeometry=null;$('fit-correlation-details').replaceChildren();host.textContent=geometry?'Run deconvolution to view correlations.':'Load detector geometry to view correlations.';return;}
 if(renderedResult===result&&renderedGeometry===geometry){correlationView?.refresh();return;}
 correlationView?.destroy();correlationView=null;
 const rows=result.detectors.map((detector,index)=>({detector,index,region:detectorRegion(detector),system:detector>=73001&&detector<=73028?'showermax':Object.keys(ringEndings).find(r=>ringEndings[r].includes(detector%100))}));
 const lookup=new Map(rows.map(r=>[r.detector,r.index]));
 const meshes=geometry.meshes.filter(m=>m.hit_surface&&lookup.has(m.geometry_detector)).map(m=>({index:lookup.get(m.geometry_detector),center:[0,1,2].map(j=>m.vertices.reduce((s,p)=>s+p[j],0)/m.vertices.length),faces:m.faces.map(f=>f.map(i=>m.vertices[i]))}));
 if(new Set(meshes.map(m=>m.index)).size!==rows.length){host.textContent='Loaded geometry does not contain every detector in this result.';return;}
 host.replaceChildren($('correlation-template').content.cloneNode(true));
 $('fit-correlation-details').replaceChildren(host.querySelector('aside'));
 correlationView=createDetectorCorrelations($('fit-page'),{rows,C:result.matrices.C,correlation:result.matrices.correlation,meshes});
 renderedResult=result;renderedGeometry=geometry;
}
function el(tag,text,cls){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n;}
function table(headers,rows,caption){
 const t=el('table');if(caption)t.append(el('caption',caption));const head=el('thead'),tr=el('tr');headers.forEach(h=>tr.append(el('th',h)));head.append(tr);t.append(head);
 const body=el('tbody');for(const row of rows){const tr=el('tr');row.forEach(v=>tr.append(el('td',v)));body.append(tr);}t.append(body);return t;
}
const fmt=n=>Number.isFinite(n)?Number(n.toPrecision(4)).toString():'Unavailable';
function render(){
 const host=$('fit-matrix');host.replaceChildren();if(!result?.matrices)return;
 document.querySelectorAll('[data-matrix]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.matrix===matrix)));
 const m=result.matrices,labels=result.labels;
 if(matrix==='uncertainties'){
  host.append(table(['Interaction','Experimental σ [ppb]','MC σ [ppb]','Combined σ [ppb]'],labels.map((name,i)=>[name,...['experimental_covariance','mc_covariance','combined_covariance'].map(k=>fmt(Math.sqrt(Math.max(0,m[k][i][i]))))]),'Projected statistical uncertainty; MC uncertainty propagated through H.'));
 }else if(matrix==='C'||matrix==='correlation'){
  const a=m[matrix],canvas=el('canvas');canvas.width=canvas.height=a.length;canvas.setAttribute('aria-label',matrix==='C'?'Detector asymmetry covariance':'Detector correlation matrix');
  const ctx=canvas.getContext('2d'),data=ctx.createImageData(a.length,a.length);let max=matrix==='correlation'?1:Math.max(...a.map(row=>Math.max(...row.map(Math.abs))));
  for(let i=0;i<a.length;i++)for(let j=0;j<a.length;j++){
   const rgb=HitMap3D.viridis(Math.max(0,Math.min(1,a[i][j]/(max||1))));let p=(i*a.length+j)*4;rgb.forEach((v,k)=>data.data[p+k]=Math.round(v*255));data.data[p+3]=255;
  }
  ctx.putImageData(data,0,0);host.append(el('p',matrix==='C'?'C · ppb² · both axes follow detector ID order':'Correlation · 0 to 1 · both axes follow detector ID order'),canvas);
  canvas.onmousemove=e=>{const r=canvas.getBoundingClientRect(),j=Math.min(a.length-1,Math.floor((e.clientX-r.left)/r.width*a.length)),i=Math.min(a.length-1,Math.floor((e.clientY-r.top)/r.height*a.length));canvas.title=`${result.detectors[i]} × ${result.detectors[j]}: ${fmt(a[i][j])}`;};
 }else host.append(table([matrix==='M'?'Interaction':'Detector',...labels],m[matrix].map((row,i)=>[matrix==='M'?labels[i]:result.detectors[i],...row.map(fmt)]),matrix==='M'?'M = Hᵀ C⁻¹ H · ppb⁻²':matrix==='H'?'H = dilution × regional asymmetry / reference asymmetry':'F · PE signal fractions'));
}
document.querySelectorAll('[data-matrix]').forEach(b=>b.onclick=()=>{matrix=b.dataset.matrix;render();});
function request(){return {paths:$('hit-file').value.split('\n').map(s=>s.trim()).filter(Boolean),extra:$('fit-extra').value.split(/[\s,]+/).filter(Boolean),days:Number($('fit-days').value),polarization:Number($('fit-polarization').value),current:Number($('fit-current').value)};}
function invalidate(){inputVersion++;result=null;$('fit-results').hidden=true;$('fit-readiness').replaceChildren();$('fit-status').textContent='Inputs changed. Check files before running.';}
for(const id of ['hit-file','fit-extra','fit-days','fit-polarization','fit-current'])$(id).addEventListener('input',invalidate);
document.addEventListener('response-changed',invalidate);
async function start(action){
 if(busy)return;startedVersion=inputVersion;busy=true;$('fit-check').disabled=$('fit-run').disabled=true;$('fit-results').hidden=true;$('fit-readiness').replaceChildren();
 $('fit-status').textContent=action==='check'?'Checking ROOT files…':'Running ROOT deconvolution…';
 try{const r=await fetch('/deconvolution/'+action,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(request())});const data=await r.json();if(!r.ok)throw Error(data.error);}
 catch(e){$('fit-status').textContent=e.message;busy=false;$('fit-check').disabled=$('fit-run').disabled=false;}
}
$('fit-check').onclick=()=>start('check');$('fit-run').onclick=()=>start('run');
document.addEventListener('workspace-status',async event=>{
 const state=event.detail.fit;if(!state)return;busy=event.detail.busy;$('fit-check').disabled=$('fit-run').disabled=busy;
 if(state.busy)$('fit-status').textContent=state.message;
 if(state.revision===revision)return;revision=state.revision;
 if(revision===0)return;
 $('fit-status').textContent=state.message;
 try{
  const response=await fetch('/deconvolution');const data=await response.json();if(inputVersion!==startedVersion){$('fit-status').textContent='Inputs changed. Run again for this selection.';return;}if(data){const source=await(await fetch('/responses')).json();if(source.resources?.fingerprint!==data.response_source?.fingerprint){$('fit-status').textContent='Response changed. Run again with the selected source.';return;}}if(data&&!inputVersion&&!$('hit-file').value.trim()){$('hit-file').value=data.paths.join('\n');$('fit-days').value=data.forecast.days;$('fit-polarization').value=data.forecast.polarization;$('fit-current').value=data.forecast.current;$('fit-extra').value=data.extra.join(', ');}if(data&&JSON.stringify(request())!==JSON.stringify({paths:data.paths,extra:data.extra,days:data.forecast.days,polarization:data.forecast.polarization,current:data.forecast.current})){$('fit-status').textContent='Inputs changed. Run again for this selection.';return;}result=data;$('fit-readiness').replaceChildren();$('fit-results').hidden=!data?.matrices;
  if(!data)return;
  if(data.files){for(const file of data.files){const row=el('div',undefined,'fit-readiness-row');row.append(el('span',file.path.split(/[\\/]/).at(-1)),el('span',file.process||'Unknown'),el('span',file.status==='Ready'?`${file.events.toLocaleString()} events · Ready`:file.status));$('fit-readiness').append(row);}if(data.missing.length)$('fit-readiness').append(el('p','Missing: '+data.missing.join(', '),'fit-missing'));}
  if(data.matrices){$('fit-summary').textContent=`${data.paths.length} files · ${data.detectors.length} detectors · ${data.forecast.days} days · ${100*data.forecast.polarization}% polarization · ${data.forecast.current} µA`;render();drawCorrelations();}
 }catch(e){$('fit-status').textContent=e.message;}
});
})();
