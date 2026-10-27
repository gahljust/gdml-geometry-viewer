/* Detector-resolved surface heat maps. WebGL draws binned signal, never tracks. */
(function(root){
'use strict';
// Current magnitude at the photocathode, before PMT multiplication. Exact SI charge.
function displaySignal(unit,display='current'){if(unit==='PE/s'&&display==='pe')return {factor:1,unit,label:'PE rate'};return unit==='PE/s'?{factor:1.602176634e-10,unit:'nA',label:'Cathode current'}:{factor:1,unit,label:'Rate'}}
// Compact scientific notation keeps axis/legend labels inside their reserved space.
function formatNumber(value){
 if(!Number.isFinite(value))return 'Unavailable';
 if(value===0)return '0';
 if(Math.abs(value)>=1e4||Math.abs(value)<1e-3){const [m,e]=value.toExponential(2).split('e');return `${Number(m)}e${Number(e)}`;}
 return Number(value.toPrecision(3)).toString();
}
// A dimensionless scale: changing PE/s to nA must not change colors.
function colorFraction(value,maximum,log=false){
 const fraction=maximum>0?Math.max(0,Math.min(1,value/maximum)):0;
 return log?Math.log1p(999*fraction)/Math.log(1000):fraction;
}
function signalData(raw,display){
 const data=raw?structuredClone(raw):null;
 if(!data||data.weighting!=='pe')return data;
 const factor=display==='current'?1.602176634e-10:1;
 data.unit=display==='current'?'nA':'PE/s';
 for(const cell of data.cells)cell[3]*=factor;
 for(const h of Object.values(data.histograms)){h.counts=h.counts.map(x=>x*factor);h.underflow*=factor;h.overflow*=factor;}
 for(const h of Object.values(data.histograms2d))for(const cell of h.cells)cell[2]*=factor;
 data.signal_total*=factor;
 for(const field of ['detector_signal','particle_signal'])data[field]=Object.fromEntries(Object.entries(data[field]).map(([k,v])=>[k,v*factor]));
 data.detectors={...data.detector_signal};data.particles={...data.particle_signal};
 return data;
}
// Face the upstream entry surfaces in both the view and its fit calculation.
const unwrapCamera={yaw:Math.PI,pitch:.60};
const regions=['closed','transition','open'],directions=['forward','backward','zero','uncertain'];
const sum=(a,b)=>a.map((v,i)=>v+b[i]),scale=(a,s)=>a.map(v=>v*s);
function viridis(t){const s=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];t=Math.max(0,Math.min(.999999,t))*4;const i=Math.floor(t),f=t-i;return s[i].map((v,k)=>(v+(s[i+1][k]-v)*f)/255)}
function aggregate(products,{secondary=false,source=null,region='all',momentum='all',mode='0',plane='main_detector'}={}){
 const first=products[0];if(!first)return null;
 const cells=new Map(),selected=Array.isArray(source)?new Set(source):source==null?null:new Set([source]);
 let total=0;
 for(const data of products){
  if(data.faces.length!==first.faces.length||data.cell_mm!==first.cell_mm)throw Error('Surface geometry mismatch');
  const column=mode==='2'?7:6;
  for(const c of data.cells){
   if(secondary&&c[3]<0)continue;
   if(selected&&(c[3]<0||!selected.has(data.sources[c[3]])))continue;
   if(region!=='all'&&regions[c[4]]!==region)continue;
   if(momentum!=='all'&&directions[c[5]]!==momentum)continue;
   const f=data.faces[c[0]];
   if(plane.endsWith('_bf')&&f.detector>=150000||plane.endsWith('_ff')&&f.detector<150000)continue;
   const value=c[column];if(value<=0)continue;
   const key=c.slice(0,3).join(',');const old=cells.get(key);
   if(old)old[3]+=value;else cells.set(key,[c[0],c[1],c[2],value]);total+=value;
  }
 }
 return {faces:first.faces,rows:first.rows,cell_mm:first.cell_mm,cells:[...cells.values()],total,plane,
  unit:plane==='circle'||mode==='2'?'PE/s':'crossings/s',linear:secondary};
}
function clipCell(face,x,y,size){
 let p=[[x*size,y*size],[(x+1)*size,y*size],[(x+1)*size,(y+1)*size],[x*size,(y+1)*size]];
 const boundary=face.polygon;
 let area=0;for(let i=0;i<boundary.length;i++){const a=boundary[i],b=boundary[(i+1)%boundary.length];area+=a[0]*b[1]-b[0]*a[1]}
 const sign=area>=0?1:-1;
 for(let i=0;i<boundary.length&&p.length;i++){
  const a=boundary[i],b=boundary[(i+1)%boundary.length],side=q=>sign*((b[0]-a[0])*(q[1]-a[1])-(b[1]-a[1])*(q[0]-a[0]));
  const out=[];
  for(let j=0;j<p.length;j++){const u=p[j],v=p[(j+1)%p.length],du=side(u),dv=side(v),inside=du>=-1e-8;
   if(inside)out.push(u);if(inside!==(dv>=-1e-8)){const t=du/(du-dv);out.push([u[0]+t*(v[0]-u[0]),u[1]+t*(v[1]-u[1])])}}
  p=out;
 }
 return p;
}
function point(face,uv,outset=0){return face.origin.map((v,i)=>v+face.u[i]*uv[0]+face.v[i]*uv[1]+face.normal[i]*outset)}
// Screen-space anchoring keeps the point under the pointer stationary while zooming.
function zoomCamera(zoom,pan,next,anchor=[0,0],center=false){
 const value=Math.max(.1,Math.min(32,next)),ratio=value/zoom;
 return {zoom:value,pan:anchor.map((a,i)=>(center?0:a)-(a-pan[i])*ratio)};
}
// Shared camera math for click picking and tile focus; mirrors the vertex shader.
function rotatePoint(p,yaw,pitch){const cy=Math.cos(yaw),sy=Math.sin(yaw),cp=Math.cos(pitch),sp=Math.sin(pitch),z=-sy*p[0]+cy*p[2];return [cy*p[0]+sy*p[2],cp*p[1]-sp*z,sp*p[1]+cp*z]}
function projectPoint(p,camera,width,height){const [x,y,z]=rotatePoint(p,camera.yaw,camera.pitch),s=Math.min(width,height)/3.05*camera.zoom*6/(6-z);return [x*s+camera.pan[0],-y*s+camera.pan[1],z]}
function pickDetector(faces,x,y,camera,width,height){
 const scale=Math.min(width,height)/3.05*camera.zoom;
 const ray=[(x-camera.pan[0])/scale,-(y-camera.pan[1])/scale,-6],origin=[0,0,6];
 const sub=(a,b)=>a.map((v,i)=>v-b[i]),dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
 let nearest=Infinity,detector=null;
 for(const f of faces){const p=f.points.map(p=>rotatePoint(p,camera.yaw,camera.pitch));
  for(let j=1;j+1<p.length;j++){
   const a=p[0],e1=sub(p[j],a),e2=sub(p[j+1],a),h=cross(ray,e2),det=dot(e1,h);if(Math.abs(det)<1e-12)continue;
   const s=sub(origin,a),u=dot(s,h)/det;if(u<0||u>1)continue;
   const q=cross(s,e1),v=dot(ray,q)/det;if(v<0||u+v>1)continue;
   const t=dot(e2,q)/det;if(t>0&&t<nearest){nearest=t;detector=f.detector;}
  }
 }return detector;
}
// Same rigid-tile cylindrical unwrapping used by the correlation view.
function unwrapGeometry(faces,centerZ){
 const tau=2*Math.PI,tiles=new Map();
 for(const f of faces){if(!tiles.has(f.detector))tiles.set(f.detector,{system:f.system,points:[]});tiles.get(f.detector).points.push(...f.polygon.map(uv=>point(f,uv).map((v,i)=>(v-(i===2?centerZ:0))/1000)));}
 const systems=[...new Set(faces.map(f=>f.system))].sort((a,b)=>a==='showermax'?1:b==='showermax'?-1:a.localeCompare(b));
 for(const t of tiles.values()){const lo=[0,1,2].map(i=>Math.min(...t.points.map(p=>p[i]))),hi=[0,1,2].map(i=>Math.max(...t.points.map(p=>p[i])));t.center=lo.map((v,i)=>(v+hi[i])/2);}
 const groups=systems.map(system=>{const ts=[...tiles.values()].filter(t=>t.system===system),points=ts.flatMap(t=>t.points),r=points.map(p=>Math.hypot(p[0],p[1])),z=points.map(p=>p[2]);return {system,r0:Math.min(...r),r1:Math.max(...r),zmid:(Math.min(...z)+Math.max(...z))/2};});
 let height=groups.reduce((a,b)=>a+b.r1-b.r0+.35,0)-.35,top=height/2;
 for(const g of groups){g.radius=(g.r0+g.r1)/2;g.y=top-(g.r1-g.r0)/2;top-=g.r1-g.r0+.35;}
 const seamTiles=[...tiles.values()].filter(t=>t.system===(systems.includes('ring5')?'ring5':systems[0]));
 const angles=seamTiles.map(t=>(Math.atan2(t.center[1],t.center[0])+tau)%tau).sort((a,b)=>a-b),seam=(angles[0]+angles.at(-1)-tau)/2;
 for(const t of tiles.values()){const g=groups.find(g=>g.system===t.system),phi=Math.atan2(t.center[1],t.center[0]),theta=((phi-seam)%tau+tau)%tau;t.angle=((-phi-Math.PI/2+Math.PI)%tau+tau)%tau-Math.PI;t.target=[g.radius*(theta-Math.PI),g.y-(Math.hypot(t.center[0],t.center[1])-g.radius),t.center[2]-g.zmid];}
 return {tiles,groups,seam};
}
function unwrapPoint(p,t,amount){const a=t.angle*amount,c=Math.cos(a),s=Math.sin(a),d=p.map((v,i)=>v-t.center[i]);return [c*d[0]-s*d[1],s*d[0]+c*d[1],d[2]].map((v,i)=>v+t.center[i]*(1-amount)+t.target[i]*amount);}
class SurfaceView{
 constructor(element,options={}){
  this.element=element;this.onSelect=options.onSelect;this.selectedDetector=null;this.focusedDetector=null;this.pickFaces=[];this.yaw=Math.PI+.45;this.pitch=.22;this.zoom=1;this.pan=[0,0];this.unwrapped=false;this.unwrapAmount=0;this.meshes=[];this.data=null;this.pending=false;this.dragMode='rotate';this.camera={yaw:this.yaw,pitch:this.pitch,zoom:this.zoom,pan:[...this.pan]};
  element.innerHTML=`<div class="surfaceToolbar"><div role="group" aria-label="Camera view"><button data-camera="oblique" aria-pressed="true">Oblique</button><button data-camera="upstream" aria-pressed="false">Upstream</button><button data-camera="downstream" aria-pressed="false">Downstream</button></div><label><input type="checkbox" class="surfaceUnwrap"> Unwrapped</label></div><div class="surfaceNavigation"><div role="group" aria-label="Drag action"><button data-drag="rotate" aria-pressed="true">Rotate</button><button data-drag="pan" aria-pressed="false">Pan</button></div><span>Double click to focus</span><div role="group" aria-label="Zoom"><button data-zoom="out" aria-label="Zoom out">−</button><button data-zoom="in" aria-label="Zoom in">+</button><button data-zoom="fit">Fit</button></div></div><div class="surfaceStage"><canvas tabindex="0" aria-label="3D detector hit map. Drag to rotate or choose Pan. Right drag or shift drag pans. Scroll zooms at the pointer. Double click focuses. Arrow keys rotate; shift arrows pan; Home fits."></canvas><canvas class="surfaceAngleOverlay" aria-hidden="true"></canvas><small class="surfaceLayoutNote" hidden>Tile spacing is illustrative only</small><div class="surfaceMessage" role="status"></div><div class="surfaceLegend"><div class="surfaceRamp"></div><div class="surfaceEnds"><span>0</span><span class="surfaceMaximum"></span></div><span class="surfaceUnit"></span></div><span class="surfaceHelp">Scroll to zoom · Shift drag to pan</span></div>`;
  this.canvas=element.querySelector('canvas');this.message=element.querySelector('.surfaceMessage');
  this.gl=this.canvas.getContext('webgl2',{antialias:true,alpha:false,preserveDrawingBuffer:true});
  if(!this.gl){this.message.textContent='This browser does not support WebGL 2.';return;}
  const gl=this.gl;
  const shader=(type,source)=>{const s=gl.createShader(type);gl.shaderSource(s,source);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))throw Error(gl.getShaderInfoLog(s));return s};
  const vertex=shader(gl.VERTEX_SHADER,`#version 300 es
   in vec3 position; in vec3 color; in vec3 tileCenter; in vec3 tileTarget; in float tileAngle; out vec3 rgb;
   uniform float unwrap;
   uniform vec2 angle; uniform vec2 screen; uniform vec2 pan; uniform float zoom;
   void main(){float cy=cos(angle.x),sy=sin(angle.x),cp=cos(angle.y),sp=sin(angle.y);
    float a=tileAngle*unwrap;vec3 d0=position-tileCenter;
    vec3 q=vec3(cos(a)*d0.x-sin(a)*d0.y,sin(a)*d0.x+cos(a)*d0.y,d0.z)+mix(tileCenter,tileTarget,unwrap);
    vec3 p=vec3(cy*q.x+sy*q.z,q.y,-sy*q.x+cy*q.z);
    p=vec3(p.x,cp*p.y-sp*p.z,sp*p.y+cp*p.z);float d=6.-p.z;float s=min(screen.x,screen.y)/3.05*zoom;
    gl_Position=vec4(12.*s*p.x/screen.x+2.*pan.x/screen.x*d,12.*s*p.y/screen.y-2.*pan.y/screen.y*d,1.004008*d-.2004008,d);rgb=color;}`);
  const fragment=shader(gl.FRAGMENT_SHADER,`#version 300 es
   precision highp float;in vec3 rgb;out vec4 outputColor;void main(){outputColor=vec4(rgb,1.);}`);
  this.program=gl.createProgram();gl.attachShader(this.program,vertex);gl.attachShader(this.program,fragment);gl.linkProgram(this.program);
  if(!gl.getProgramParameter(this.program,gl.LINK_STATUS))throw Error(gl.getProgramInfoLog(this.program));
  gl.deleteShader(vertex);gl.deleteShader(fragment);
  this.uniforms=Object.fromEntries(['angle','screen','pan','zoom','unwrap'].map(n=>[n,gl.getUniformLocation(this.program,n)]));
  element.querySelectorAll('[data-camera]').forEach(b=>b.onclick=()=>{this.unwrapped=false;element.querySelector('.surfaceUnwrap').checked=false;element.querySelector('.surfaceLayoutNote').hidden=true;this.focusedDetector=null;if(this.data)this.setData(this.data);this.yaw=b.dataset.camera==='oblique'?Math.PI+.45:b.dataset.camera==='upstream'?0:Math.PI;this.pitch=b.dataset.camera==='oblique'?.22:0;this.zoom=1;this.pan=[0,0];element.querySelectorAll('[data-camera]').forEach(x=>x.setAttribute('aria-pressed',String(x===b)));this.draw()});
  element.querySelector('.surfaceUnwrap').onchange=e=>{this.unwrapped=e.target.checked;this.dragMode=this.unwrapped?'pan':'rotate';element.querySelectorAll('[data-drag]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.drag===this.dragMode)));element.querySelector('.surfaceLayoutNote').hidden=!this.unwrapped;this.yaw=this.unwrapped?unwrapCamera.yaw:Math.PI+.45;this.pitch=this.unwrapped?unwrapCamera.pitch:.22;this.fit();element.querySelectorAll('[data-camera]').forEach(b=>b.setAttribute('aria-pressed',String(!this.unwrapped&&b.dataset.camera==='oblique')));};
  element.querySelectorAll('[data-drag]').forEach(b=>b.onclick=()=>{this.dragMode=b.dataset.drag;element.querySelectorAll('[data-drag]').forEach(x=>x.setAttribute('aria-pressed',String(x===b)));});
  element.querySelectorAll('[data-zoom]').forEach(b=>b.onclick=()=>{if(b.dataset.zoom==='fit')this.fit();else this.zoomAt(this.zoom*(b.dataset.zoom==='in'?1.5:1/1.5));});
  const position=e=>{const r=this.canvas.getBoundingClientRect();return [e.clientX-r.left-r.width/2,e.clientY-r.top-r.height/2]};
  const pointers=new Map();let drag=null;
  const pair=()=>{const p=[...pointers.values()];return {center:[(p[0][0]+p[1][0])/2,(p[0][1]+p[1][1])/2],distance:Math.hypot(p[0][0]-p[1][0],p[0][1]-p[1][1])}};
  this.canvas.oncontextmenu=e=>e.preventDefault();
  this.canvas.onpointerdown=e=>{if(e.button>2)return;this.clickStart={x:e.clientX,y:e.clientY,button:e.button,moved:false};if(pointers.size)this.clickStart.moved=true;e.preventDefault();this.canvas.focus({preventScroll:true});pointers.set(e.pointerId,position(e));drag={x:e.clientX,y:e.clientY,pan:e.shiftKey||e.button===1||e.button===2||this.dragMode==='pan'||this.unwrapped};this.canvas.classList.add('dragging');this.canvas.setPointerCapture(e.pointerId)};
  this.canvas.onpointermove=e=>{if(this.clickStart&&Math.hypot(e.clientX-this.clickStart.x,e.clientY-this.clickStart.y)>4)this.clickStart.moved=true;if(!pointers.has(e.pointerId))return;
   if(pointers.size===2){const before=pair();pointers.set(e.pointerId,position(e));const after=pair();if(before.distance>0)this.zoomAt(this.zoom*after.distance/before.distance,before.center);this.pan=this.pan.map((v,i)=>v+after.center[i]-before.center[i]);this.draw();return;}
   pointers.set(e.pointerId,position(e));if(!drag)return;const dx=e.clientX-drag.x,dy=e.clientY-drag.y;
   if(drag.pan||e.shiftKey){this.pan[0]+=dx;this.pan[1]+=dy}else{this.yaw+=dx*.006;this.pitch=Math.max(-1.5,Math.min(1.5,this.pitch+dy*.006));}
   drag.x=e.clientX;drag.y=e.clientY;this.draw();};
  const releasePointer=e=>{pointers.delete(e.pointerId);drag=null;this.clickStart=null;if(!pointers.size)this.canvas.classList.remove('dragging');};
  this.canvas.onpointerup=e=>{if(this.onSelect&&this.clickStart&&!this.clickStart.moved&&this.clickStart.button===0){const p=position(e),r=this.canvas.getBoundingClientRect(),det=pickDetector(this.currentPickFaces(),...p,this.camera,r.width,r.height);if(det!==null){if(det===this.selectedDetector)this.clearSelection();else this.selectDetector(det);}}releasePointer(e);};
  this.canvas.onpointercancel=this.canvas.onlostpointercapture=releasePointer;
  this.canvas.addEventListener('wheel',e=>{e.preventDefault();const delta=e.deltaY*(e.deltaMode===1?16:e.deltaMode===2?this.canvas.clientHeight:1);this.zoomAt(this.zoom*Math.exp(-Math.max(-500,Math.min(500,delta))*.0018),position(e));},{passive:false});
  this.canvas.ondblclick=e=>{e.preventDefault();if(!this.onSelect)this.zoomAt(this.zoom*2,position(e),true)};
  this.canvas.onkeydown=e=>{if(e.key==='Escape'&&this.onSelect){e.preventDefault();this.clearSelection();return}if(e.key==='Home'||e.key==='0'){e.preventDefault();this.fit();return}if(['+','=','-','_'].includes(e.key)){e.preventDefault();this.zoomAt(this.zoom*(['-','_'].includes(e.key)?1/1.4:1.4));return}if(!e.key.startsWith('Arrow'))return;e.preventDefault();
   if(e.shiftKey||this.unwrapped){this.pan[0]+=e.key==='ArrowLeft'?-30:e.key==='ArrowRight'?30:0;this.pan[1]+=e.key==='ArrowUp'?-30:e.key==='ArrowDown'?30:0;}
   else{if(e.key==='ArrowLeft')this.yaw-=.08;if(e.key==='ArrowRight')this.yaw+=.08;if(e.key==='ArrowUp')this.pitch-=.08;if(e.key==='ArrowDown')this.pitch+=.08;}this.draw();};
  this.observer=new ResizeObserver(()=>this.draw());this.observer.observe(this.canvas);
 }
 clearSelection(){this.selectedDetector=null;this.fit();this.onSelect?.(null);}
 selectDetector(detector,focus=true){
  this.selectedDetector=detector;this.focusedDetector=focus?detector:null;if(this.data)this.setData(this.data);
  if(focus){const points=this.currentPickFaces().filter(f=>f.detector===detector).flatMap(f=>f.points),r=this.canvas.getBoundingClientRect();
   if(points.length&&r.width&&r.height){const c={yaw:this.yaw,pitch:this.pitch,zoom:1,pan:[0,0]},p=points.map(v=>projectPoint(v,c,r.width,r.height)),xs=p.map(v=>v[0]),ys=p.map(v=>v[1]);
    const x0=Math.min(...xs),x1=Math.max(...xs),y0=Math.min(...ys),y1=Math.max(...ys);
    this.zoom=Math.max(1,Math.min(18,.56*Math.min(r.width/Math.max(1,x1-x0),r.height/Math.max(1,y1-y0))));
    this.pan=[-(x0+x1)/2*this.zoom,-(y0+y1)/2*this.zoom];this.draw();
   }
  }this.onSelect?.(detector);
 }
 zoomAt(next,anchor=[0,0],center=false){const view=zoomCamera(this.zoom,this.pan,next,anchor,center);this.zoom=view.zoom;this.pan=view.pan;this.draw();}
 currentPickFaces(){return this.pickFaces.map(f=>({detector:f.detector,points:f.points.map(p=>unwrapPoint(p,this.unwrapLayout.tiles.get(f.detector),this.unwrapAmount))}));}
 fit(){this.focusedDetector=null;this.zoom=1;this.pan=[0,0];if(this.data)this.setData(this.data);else this.draw();if(this.unwrapped&&this.pickFaces.length){
  const r=this.canvas.getBoundingClientRect(),camera={...unwrapCamera,zoom:1,pan:[0,0]},points=this.pickFaces.flatMap(f=>f.points.map(p=>projectPoint(unwrapPoint(p,this.unwrapLayout.tiles.get(f.detector),1),camera,r.width,r.height)));
  const xs=points.map(p=>p[0]),ys=points.map(p=>p[1]),x0=Math.min(...xs),x1=Math.max(...xs),y0=Math.min(...ys),y1=Math.max(...ys);
  this.zoom=Math.max(.1,Math.min(.80*r.width/(x1-x0),.70*r.height/(y1-y0)));this.pan=[-(x0+x1)/2*this.zoom,-(y0+y1)/2*this.zoom];this.draw();
 }}
 status(text){this.message.textContent=text;this.message.hidden=!text;}
 clear(text='Loading 3D hit map…'){this.data=null;this.release();this.status(text);this.element.querySelector('.surfaceLegend').hidden=true;this.draw();}
 release(){if(!this.gl)return;for(const m of this.meshes){this.gl.deleteBuffer(m.buffer);this.gl.deleteVertexArray(m.vao)}this.meshes=[];}
 dispose(){this.disposed=true;this.release();this.observer?.disconnect();this.gl?.getExtension('WEBGL_lose_context')?.loseContext();}
 buffer(values,mode){const gl=this.gl,vao=gl.createVertexArray(),buffer=gl.createBuffer();gl.bindVertexArray(vao);gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(values),gl.STATIC_DRAW);for(const [name,offset] of [['position',0],['color',12],['tileCenter',24],['tileTarget',36],['tileAngle',48]]){const a=gl.getAttribLocation(this.program,name);gl.enableVertexAttribArray(a);gl.vertexAttribPointer(a,name==='tileAngle'?1:3,gl.FLOAT,false,52,offset)}this.meshes.push({vao,buffer,mode,count:values.length/13});}
 setData(data){
  this.data=data;if(!this.gl)return;this.release();
  const plane=data.plane,allowed=f=>plane==='all_detectors'?true:plane==='main_detector'?f.system!=='showermax':plane==='circle'?f.system==='showermax':f.system===plane.slice(0,5)&&(!plane.endsWith('_bf')||f.detector<150000)&&(!plane.endsWith('_ff')||f.detector>=150000);
  const selectedFaces=data.faces.filter(allowed);if(!selectedFaces.length){this.pickFaces=[];this.status('No detector geometry');return;}const zs=selectedFaces.map(f=>f.origin[2]),centerZ=(Math.min(...zs)+Math.max(...zs))/2;
  this.unwrapLayout=unwrapGeometry(selectedFaces,centerZ);
  const world=p=>[p[0]/1000,p[1]/1000,(p[2]-centerZ)/1000];
  const vertex=(p,color,f)=>{const t=this.unwrapLayout.tiles.get(f.detector);return [...world(p),...color,...t.center,...t.target,t.angle];};
  const background=[],edges=[],signal=[],selectedEdges=[];this.pickFaces=[];
  const tileMetric=data.tileValues,metricValues=tileMetric?Object.values(tileMetric).filter(Number.isFinite):[];
  let max=tileMetric?Math.max(0,...metricValues):0;
  const triangle=(dst,a,b,c,color,f)=>{for(const p of [a,b,c])dst.push(...vertex(p,color,f))};
  const polygon=(dst,p,color,sys)=>{for(let j=1;j+1<p.length;j++)triangle(dst,p[0],p[j],p[j+1],color,sys)};
  for(const f of selectedFaces){
   const p=f.polygon.map(uv=>point(f,uv)),worldPoints=p.map(v=>world(v,f.system));if(this.focusedDetector===null||f.detector===this.focusedDetector)this.pickFaces.push({detector:f.detector,points:worldPoints});
   const shade=.20+.08*Math.abs(f.normal[0]*.3+f.normal[1]*.4+f.normal[2]*.86),metric=tileMetric?.[f.detector];
   const t=max>0&&Number.isFinite(metric)?(colorFraction(metric,max,!data.linear)):0;
   const color=tileMetric&&Number.isFinite(metric)?viridis(t).map(v=>v*(.82+.18*Math.abs(f.normal[2]))):[shade*.8,shade*.95,shade*1.12];
   const focused=this.focusedDetector===null||f.detector===this.focusedDetector;
   if(focused)polygon(background,p,color,f);
   const edgeColor=focused?[.34,.43,.51]:[.10,.16,.22];
   for(let j=0;j<p.length;j++)edges.push(...vertex(p[j],edgeColor,f),...vertex(p[(j+1)%p.length],edgeColor,f));
   if(f.detector===this.selectedDetector){const pp=f.polygon.map(uv=>point(f,uv,.06));for(let j=0;j<pp.length;j++)selectedEdges.push(...vertex(pp[j],[.77,.96,1],f),...vertex(pp[(j+1)%pp.length],[.77,.96,1],f));}
  }
  if(!tileMetric){for(const c of data.cells)max=Math.max(max,c[3]);
   for(const [fi,x,y,w] of data.cells){const f=data.faces[fi];if(!allowed(f))continue;const p=clipCell(f,x,y,data.cell_mm).map(uv=>point(f,uv,.035));if(p.length<3)continue;const t=colorFraction(w,max,!data.linear);polygon(signal,p,viridis(t),f);}}
  this.buffer(background,this.gl.TRIANGLES);this.buffer(edges,this.gl.LINES);this.buffer(signal,this.gl.TRIANGLES);if(selectedEdges.length)this.buffer(selectedEdges,this.gl.LINES);
  this.status(tileMetric?'':data.total>0?'':'No hits for this selection');
  const legend=this.element.querySelector('.surfaceLegend');legend.hidden=!(max>0);
  const display=displaySignal(data.unit,data.signalDisplay);
  this.element.querySelector('.surfaceMaximum').textContent=formatNumber(max*display.factor);
  this.element.querySelector('.surfaceUnit').textContent=data.legendLabel||`${display.unit} / surface bin${data.linear?'':' · log color'}`;
  this.element.querySelector('.surfaceRamp').style.background='linear-gradient(90deg,'+Array.from({length:12},(_,i)=>`rgb(${viridis(i/11).map(x=>Math.round(x*255)).join(',')})`).join(',')+')';
  this.draw();
 }
 drawAngles(r,dpr){
 const canvas=this.element.querySelector('.surfaceAngleOverlay');if(canvas.width!==this.canvas.width||canvas.height!==this.canvas.height){canvas.width=this.canvas.width;canvas.height=this.canvas.height;}
 const ctx=canvas.getContext('2d');ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,r.width,r.height);if(this.unwrapAmount<.95||!this.unwrapLayout)return;
 const {groups,seam}=this.unwrapLayout,tau=2*Math.PI;const project=p=>{const q=projectPoint(p,this.camera,r.width,r.height);return [q[0]+r.width/2,q[1]+r.height/2];};
 const at=(g,phi,y=g.y)=>project([g.radius*((((phi-seam)%tau+tau)%tau)-Math.PI),y,0]);ctx.strokeStyle='#809ab15c';ctx.fillStyle='#9ab4ca';ctx.font='11px system-ui';ctx.setLineDash([2,5]);ctx.textAlign='center';
 const first=groups[0],last=groups.at(-1),stride=Math.max(1,Math.ceil(43/Math.abs(at(first,tau/7)[0]-at(first,0)[0])));
 for(let i=0;i<7;i++){const phi=i*tau/7,start=at(first,phi,first.y+(first.r1-first.r0)/2+.16);ctx.beginPath();ctx.moveTo(...start);for(const g of groups)ctx.lineTo(...at(g,phi));ctx.lineTo(...at(last,phi,last.y-(last.r1-last.r0)/2-.08));ctx.stroke();if(i%stride===0)ctx.fillText((i*360/7).toFixed(1)+'°',start[0],start[1]-7);}
 ctx.textAlign='right';for(const g of groups){const q=project([Math.PI*g.radius+.09,g.y,0]);ctx.fillText(g.system==='showermax'?'ShowerMax':'Ring '+g.system.slice(4),q[0],q[1]+4);}
 }

 draw(){if(!this.gl||this.pending||this.disposed)return;this.pending=true;requestAnimationFrame(now=>{this.pending=false;if(this.disposed)return;
 const dt=this.lastDraw?Math.min(50,now-this.lastDraw):16;this.lastDraw=now;
 const alpha=1-Math.exp(-dt/65),c=this.camera;let moving=false;const delta=(this.unwrapped?1:0)-this.unwrapAmount;if(Math.abs(delta)>1e-4){this.unwrapAmount+=delta*(1-Math.exp(-dt/190));moving=true}else this.unwrapAmount=this.unwrapped?1:0;
 for(const key of ['yaw','pitch','zoom']){const d=this[key]-c[key];if(Math.abs(d)>1e-5){c[key]+=d*alpha;moving=true}else c[key]=this[key];}
 for(let i=0;i<2;i++){const d=this.pan[i]-c.pan[i];if(Math.abs(d)>.02){c.pan[i]+=d*alpha;moving=true}else c.pan[i]=this.pan[i];}
const gl=this.gl,r=this.canvas.getBoundingClientRect();if(!r.width||!r.height)return;const dpr=Math.min(2,root.devicePixelRatio||1),w=Math.round(r.width*dpr),h=Math.round(r.height*dpr);if(this.canvas.width!==w||this.canvas.height!==h){this.canvas.width=w;this.canvas.height=h}gl.viewport(0,0,w,h);gl.clearColor(.047,.078,.125,1);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);gl.enable(gl.DEPTH_TEST);gl.depthFunc(gl.LEQUAL);gl.useProgram(this.program);gl.uniform2f(this.uniforms.angle,c.yaw,c.pitch);gl.uniform2f(this.uniforms.screen,r.width,r.height);gl.uniform2f(this.uniforms.pan,...c.pan);gl.uniform1f(this.uniforms.zoom,c.zoom);gl.uniform1f(this.uniforms.unwrap,this.unwrapAmount);for(const m of this.meshes){gl.bindVertexArray(m.vao);gl.drawArrays(m.mode,0,m.count)}gl.bindVertexArray(null);this.drawAngles(r,dpr);if(moving)this.draw();});}
}
const api={formatNumber,colorFraction,signalData,SurfaceView,aggregate,clipCell,viridis,displaySignal,zoomCamera,rotatePoint,projectPoint,pickDetector,unwrapGeometry,unwrapPoint,unwrapCamera};
if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.HitMap3D=api;
})(typeof window!=='undefined'?window:globalThis);
