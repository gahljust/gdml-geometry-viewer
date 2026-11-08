/* Website gestures and descriptive fits adapted to the ROOT histogram SVG.
 * Display transforms never alter bins, weights, or analysis results. */
window.PlotInteractions={create({svg,redraw,reset,beforeFit,afterFit}){
 const $=id=>document.getElementById(id),math=PlotAxesMath,core=CurveFit;
 const node=(tag,attrs,text)=>{const n=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const[k,v]of Object.entries(attrs))n.setAttribute(k,v);if(text!==undefined)n.textContent=text;return n;};
 const panel=document.createElement('aside');panel.id='curve-fit';panel.hidden=true;
 panel.innerHTML=`<header><b>Curve fit</b><button aria-label="Close curve fit">×</button></header><select aria-label="Fit model"><option value="gaussian">Gaussian</option><option value="power">Power</option><option value="constant">Constant</option></select><label>From <output id="curve-lo-value"></output></label><input type="range" id="curve-lo" aria-label="Fit from bin"><label>To <output id="curve-hi-value"></output></label><input type="range" id="curve-hi" aria-label="Fit to bin"><p id="curve-fit-status" role="status"></p><small>Descriptive fit · selected bins</small><p><button id="curve-clear">Clear fit</button></p>`;
 svg.closest('.plotbody').append(panel);
 let plot=null,source=null,key=null,result=null,drag=null,timer=null,gesture=null,remembered={};
 const close=()=>{panel.hidden=true;svg.closest('.plotpanel').classList.remove('fit-open');afterFit();};
 panel.querySelector('header button').onclick=()=>{close();redraw();};
 $('curve-clear').onclick=()=>{result=null;redraw();};
 const state=axis=>({domain:[...(axis==='x'?plot.xr:plot.yr)],mode:(axis==='x'?plot.logX:plot.logY)?'log':'linear'});
 const set=(axis,range)=>{range.forEach((v,i)=>$('plot-'+axis+(i?'-max':'-min')).value=Number(v.toPrecision(12)));redraw();};
 const position=e=>{const b=svg.getBoundingClientRect();return [(e.clientX-b.left)*plot.W/b.width,(e.clientY-b.top)*plot.H/b.height];};
 const anchor=(e,axis)=>{const[x,y]=position(e);return axis==='x'?(x-plot.left)/plot.pw:1-(y-plot.top)/plot.ph;};
 const toggle=axis=>{
  const s=state(axis),next=s.mode==='log'?'linear':'log';
  const range=remembered[axis+next]||(next==='log'?math.logRange(s.domain,axis==='x'?plot.xValues:plot.yValues):s.domain);
  if(!range)return;remembered[axis+s.mode]=s.domain;$('plot-log-'+axis).checked=next==='log';set(axis,range);
 };
 const resetAxes=()=>{remembered={};reset();};
 const fit=()=>{
  if(!plot?.h)return;const h=plot.h;
  let lo=Number($('curve-lo').value),hi=Number($('curve-hi').value);
  if(lo>hi){if(document.activeElement===$('curve-lo'))hi=lo;else lo=hi;}
  $('curve-lo').value=lo;$('curve-hi').value=hi;
  $('curve-lo-value').textContent=math.label(h.edges[lo]);$('curve-hi-value').textContent=math.label(h.edges[hi+1]);
  const kind=panel.querySelector('select').value,points=h.counts.map((v,i)=>[(h.edges[i]+h.edges[i+1])/2,v,h.edges[i],h.edges[i+1]]);
  try{result=core.fit({kind,histogram:kind==='gaussian',points:points.slice(lo,hi+1),range:[h.edges[lo],h.edges[hi+1]]});$('curve-fit-status').textContent=result.boundary?'Fit reached its search limit.':'';}
  catch(e){result=null;$('curve-fit-status').textContent=e.message.replace('positive N','positive x');}
  redraw();
 };
 panel.querySelector('select').onchange=fit;$('curve-lo').oninput=fit;$('curve-hi').oninput=fit;
 const open=()=>{
  if(!plot?.h)return;beforeFit();panel.hidden=false;svg.closest('.plotpanel').classList.add('fit-open');
  const h=plot.h,n=h.counts.length;
  for(const id of ['curve-lo','curve-hi']){$(id).min=0;$(id).max=n-1;$(id).step=1;}
  $('curve-lo').value=Math.max(0,h.edges.findIndex((v,i)=>i<n&&h.edges[i+1]>plot.xr[0]));
  $('curve-hi').value=Math.max(0,h.edges.findLastIndex((v,i)=>i<n&&v<plot.xr[1]));fit();
 };
 svg.addEventListener('dblclick',e=>{if(!plot)return;clearTimeout(timer);e.preventDefault();if(e.target.dataset.color){$('plot-log-color').checked=!$('plot-log-color').checked;redraw();return;}const axis=e.target.dataset.axis;if(axis)toggle(axis);else resetAxes();});
 svg.addEventListener('click',e=>{if(!plot||e.target.dataset.axis||e.detail>1)return;const[x,y]=position(e);if(x<plot.left||x>plot.left+plot.pw||y<plot.top||y>plot.top+plot.ph)return;clearTimeout(timer);timer=setTimeout(open,260);});
 svg.addEventListener('pointerdown',e=>{const axis=e.target.dataset.axis;if(!plot||!axis||e.button!==0)return;e.preventDefault();const b=svg.getBoundingClientRect();drag={axis,start:state(axis),origin:axis==='x'?e.clientX:e.clientY,pointerId:e.pointerId,captured:false,extent:axis==='x'?plot.pw*b.width/plot.W:plot.ph*b.height/plot.H};});
 svg.addEventListener('pointermove',e=>{if(!drag)return;const d=drag.axis==='x'?e.clientX-drag.origin:drag.origin-e.clientY;if(!drag.captured){if(Math.abs(d)<3)return;svg.setPointerCapture(drag.pointerId);drag.captured=true;}set(drag.axis,math.scaledRange(drag.start,d,drag.extent));});
 for(const name of ['pointerup','pointercancel','lostpointercapture'])svg.addEventListener(name,()=>drag=null);
 svg.addEventListener('wheel',e=>{const axis=e.target.dataset.axis;if(!plot||!axis||gesture)return;e.preventDefault();e.stopPropagation();set(axis,math.zoomRange(state(axis),anchor(e,axis),math.wheelFactor(e,axis==='x'?plot.pw:plot.ph)));},{passive:false});
 svg.addEventListener('gesturestart',e=>{const axis=e.target.dataset.axis;if(!plot||!axis)return;e.preventDefault();gesture={axis,start:state(axis),anchor:anchor(e,axis)};},{passive:false});
 svg.addEventListener('gesturechange',e=>{if(!gesture||!(e.scale>0))return;e.preventDefault();set(gesture.axis,math.zoomRange(gesture.start,gesture.anchor,1/e.scale));},{passive:false});
 svg.addEventListener('gestureend',()=>gesture=null);
 svg.addEventListener('keydown',e=>{const axis=e.target.dataset.axis;if(!axis||!plot)return;if(e.key==='Enter'){e.preventDefault();toggle(axis);}else if(e.key==='Escape'){e.preventDefault();resetAxes();}});
 svg.closest('.plotpanel').querySelector('h2').ondblclick=resetAxes;
 function update(p){
  plot=p;if(source!==p.data||key!==p.key){source=p.data;key=p.key;result=null;remembered={};close();}
  if(result&&p.h){
   const g=node('g',{'clip-path':'url(#distribution-clip)'}),points=[];
   p.h.counts.forEach((v,i)=>{const a=p.h.edges[i],b=p.h.edges[i+1],x=(a+b)/2;if(x<result.range[0]||x>result.range[1])return;
    // Integrate the Gaussian across actual bin boundaries, including negative x.
    const y=result.kind==='gaussian'?result.params[0]*core.binShape('gaussian',a,b,result.params.slice(1)):core.predict(result,x);
    const px=p.xpos(x),py=p.ypos(y);if(Number.isFinite(px)&&Number.isFinite(py))points.push([px,py]);
   });
   if(points.length)g.append(node('path',{d:points.map(([x,y],i)=>(i?'L':'M')+x+' '+y).join(' '),fill:'none',stroke:'#b52d32','stroke-width':2}));svg.append(g);
   const f=v=>Number(v.toPrecision(4)),q=result.params;
   const lines=result.kind==='gaussian'?['Gaussian',`μ = ${f(q[1])}   σ = ${f(q[2])}`,`Area = ${f(q[0])}`]:result.kind==='power'?['Power',`y = ${f(q[0])} (x/${f(q[2])})^${f(-q[1])}`]:['Constant',`y = ${f(q[0])}`];
   const width=Math.min(p.pw-10,Math.max(...lines.map(l=>l.length))*5.6+12),height=lines.length*13+8;
   const boxes=p.h.counts.flatMap((v,i)=>{const a=p.xpos(Math.max(p.xr[0],p.h.edges[i])),b=p.xpos(Math.min(p.xr[1],p.h.edges[i+1])),y=p.ypos(v);return b>a&&Number.isFinite(y)?[[a-p.left,Math.max(0,y-p.top),b-a,p.ph-Math.max(0,y-p.top)]]:[];});
   const box=core.popoverBox({width:p.pw,height:p.ph,size:[width,height],boxes,anchor:[p.pw,0]});
   if(box){const[x,y]=box;lines.forEach((t,i)=>svg.append(node('text',{x:p.left+x+5,y:p.top+y+14+i*13,fill:'#b52d32','font-size':10,'font-family':'system-ui',stroke:'white','stroke-width':3,'paint-order':'stroke'},t)));}
  }
  if(!p.h&&!p.categorical){const color=node('rect',{x:p.left+p.pw+8,y:p.top,width:40,height:p.ph,fill:'transparent','data-color':'true',role:'button',tabindex:0,'aria-label':'Color scale: double click for log or linear'});color.append(node('title',{},'Double click for log / linear color'));svg.append(color);}
  for(const axis of ['x','y']){
   if(axis==='x'&&p.categorical)continue;
   const attrs=axis==='x'?{x:p.left,y:p.top+p.ph-3,width:p.pw,height:38}:{x:25,y:p.top,width:p.left-22,height:p.ph-3};
   const handle=node('rect',{...attrs,fill:'transparent','data-axis':axis,tabindex:0,role:'button','aria-label':`${axis.toUpperCase()} axis: drag or scroll to scale; double click or Enter for log; Escape to reset`,style:`cursor:${axis==='x'?'ew':'ns'}-resize`});
   handle.append(node('title',{},'Drag / scroll to scale · double click for log'));svg.append(handle);
  }
 }
 return {update,close,clear(){plot=null;result=null;close();}};
}};
