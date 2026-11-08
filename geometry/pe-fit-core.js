/* Descriptive curve fits only. No Poisson count assumption, error bars, p-values,
 * or independent-checkpoint likelihood. Fit the saved ordinates in their native
 * units, regardless of the display's linear/log settings. */
(function(root){
'use strict';
const finite=v=>typeof v==='number'&&Number.isFinite(v);
const clamp=(v,a,b)=>Math.max(a,Math.min(b,v));
const sum=a=>a.reduce((s,v)=>s+v,0);
function erf(x){const sign=x<0?-1:1;x=Math.abs(x);const t=1/(1+.3275911*x);return sign*(1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-.284496736)*t+.254829592)*t*Math.exp(-x*x));}
const cdf=z=>.5*(1+erf(z/Math.SQRT2));
function binShape(kind,lo,hi,p){
 if(kind==='gaussian')return Math.max(0,cdf((hi-p[0])/p[1])-cdf((lo-p[0])/p[1]));
 if(kind==='lognormal')return hi<=0?0:Math.max(0,cdf((Math.log(hi)-p[0])/p[1])-(lo<=0?0:cdf((Math.log(lo)-p[0])/p[1])));
 throw Error('Unknown peak model');
}
function amplitude(shape,y){const d=sum(shape.map(v=>v*v));return d>1e-30?Math.max(0,sum(shape.map((v,i)=>v*y[i]))/d):0;}
const loss=(shape,y,a)=>sum(shape.map((v,i)=>(a*v-y[i])**2))/y.length;
function golden(f,lo,hi){const g=(Math.sqrt(5)-1)/2;let a=hi-g*(hi-lo),b=lo+g*(hi-lo),fa=f(a),fb=f(b);
 for(let i=0;i<75;i++){if(fa<fb){hi=b;b=a;fb=fa;a=hi-g*(hi-lo);fa=f(a);}else{lo=a;a=b;fa=fb;b=lo+g*(hi-lo);fb=f(b);}}return fa<fb?a:b;}
function simplex(f,start,bounds){
 const bound=p=>p.map((v,i)=>clamp(v,...bounds[i]));
 let s=[start,start.map((v,i)=>v+(i===0?(bounds[0][1]-bounds[0][0])*.06:0)),start.map((v,i)=>v+(i===1?.18:0))].map(p=>({p:bound(p),v:0}));s.forEach(q=>q.v=f(q.p));
 const evalAt=p=>{p=bound(p);return {p,v:f(p)}};
 for(let k=0;k<200;k++){
  s.sort((a,b)=>a.v-b.v);const c=s[0].p.map((v,i)=>(v+s[1].p[i])/2),r=evalAt(c.map((v,i)=>2*v-s[2].p[i]));
  if(r.v<s[0].v){const e=evalAt(c.map((v,i)=>v+2*(r.p[i]-v)));s[2]=e.v<r.v?e:r;}
  else if(r.v<s[1].v)s[2]=r;
  else{const o=r.v<s[2].v,co=evalAt(c.map((v,i)=>v+.5*((o?r.p[i]:s[2].p[i])-v)));
   if(co.v<(o?r.v:s[2].v))s[2]=co;else s=s.map((q,i)=>i?evalAt(q.p.map((v,j)=>(v+s[0].p[j])/2)):q);
  }
  if(k>30&&Math.max(...s.map(q=>q.v))-Math.min(...s.map(q=>q.v))<1e-14)break;
 }
 return s.sort((a,b)=>a.v-b.v)[0];
}
function trend(points,kind){
 const norm=Math.max(...points.map(p=>Math.abs(p[1])));if(!(norm>0))throw Error('No nonzero values');
 const x=points.map(p=>p[0]),y=points.map(p=>p[1]/norm);
 if(kind==='constant'){const c=sum(y)/y.length;return {kind,params:[c*norm],loss:sum(y.map(v=>(v-c)**2))/y.length};}
 if(kind!=='power')throw Error('Unknown trend model');
 if(x.some(v=>v<=0)||y.some(v=>v<0))throw Error('Power fit needs positive N and nonnegative values');
 const pivot=Math.exp(sum(x.map(Math.log))/x.length),shape=p=>x.map(v=>(v/pivot)**(-p));
 const score=p=>{const g=shape(p);return loss(g,y,amplitude(g,y));};
 let best={p:0,s:Infinity};
 for(let p=-3;p<=3;p+=.1){const s=score(p);if(s<best.s)best={p,s};}
 const p=golden(score,Math.max(-3,best.p-.15),Math.min(3,best.p+.15));
 const a=amplitude(shape(p),y)*norm;
 return {kind,params:[a,p,pivot],loss:score(p),boundary:Math.abs(p)>2.99};
}
function peak(points,kind){
 const max=Math.max(...points.map(p=>p[1]));if(!(max>0))throw Error('No accepted entries');
 const y=points.map(p=>p[1]/max),x=points.map(p=>p[0]),total=sum(y),mean=sum(x.map((v,i)=>v*y[i]))/total;
 const sd=Math.sqrt(sum(x.map((v,i)=>(v-mean)**2*y[i]))/total),lo=points[0][2],hi=points.at(-1)[3],span=hi-lo;
 const bw=Math.min(...points.map(p=>p[3]-p[2]));
 if(!(span>bw))throw Error('Choose a wider interval');
 let bounds,initial;
 if(kind==='gaussian'){bounds=[[lo-span,hi+span],[Math.log(bw*.25),Math.log(span*2)]];initial=[mean,Math.log(Math.max(sd,bw))];}
 else{
  const lx=x.map(v=>Math.log(Math.max(v,bw*.1))),m=sum(lx.map((v,i)=>v*y[i]))/total,s=Math.sqrt(sum(lx.map((v,i)=>(v-m)**2*y[i]))/total);
  bounds=[[Math.log(bw*.01),Math.log(Math.max(hi,1)*20)],[Math.log(.02),Math.log(3)]];initial=[m,Math.log(Math.max(.1,s))];
 }
 const physical=p=>[p[0],Math.exp(p[1])],shape=p=>points.map(q=>binShape(kind,q[2],q[3],physical(p)));
 const objective=p=>{const g=shape(p);return loss(g,y,amplitude(g,y));};
 const starts=[initial,[initial[0],initial[1]-.5],[initial[0],initial[1]+.5]];
 let best=starts.map(s=>simplex(objective,s,bounds)).sort((a,b)=>a.v-b.v)[0];
 const p=physical(best.p),area=amplitude(shape(best.p),y)*max;
 const meanFit=kind==='gaussian'?p[0]:Math.exp(p[0]+p[1]**2/2);
 const rms=kind==='gaussian'?p[1]:Math.sqrt((Math.exp(p[1]**2)-1)*Math.exp(2*p[0]+p[1]**2));
 return {kind,params:[area,...p],mean:meanFit,rms,loss:best.v,
  boundary:best.p.some((v,i)=>Math.min(v-bounds[i][0],bounds[i][1]-v)<1e-4)};
}
function histogramStats(input){
 const {range}=input;
 if(!Array.isArray(range)||range.length!==2||!range.every(finite)||range[1]<=range[0])throw Error('Choose a PE interval');
 const all=input.points.map(p=>[...p]).sort((a,b)=>a[2]-b[2]);
 if(!all.length||all.some(p=>p.length!==4||!p.every(finite)||p[1]<0||p[3]<=p[2]))throw Error('Invalid histogram bins');
 if(all.some((p,i)=>i&&p[2]<all[i-1][3]))throw Error('Overlapping histogram bins');
 // Values are already integrated entry rates in each bin, not densities.
 // Select whole bins: no interpolation or assumption about their interior.
 const points=all.filter(p=>p[3]>range[0]&&p[2]<range[1]);
 if(!points.length)throw Error('Choose a PE interval');
 if(points.some(p=>p[2]<range[0]||p[3]>range[1]))throw Error('Select whole PE bins');
 const rate=sum(points.map(p=>p[1])),maximum=Math.max(...points.map(p=>p[1]));
 let mean=null,rms=null;
 if(maximum>0){
  const weights=points.map(p=>p[1]/maximum),total=sum(weights);
  mean=sum(points.map((p,i)=>(p[2]+p[3])/2*weights[i]))/total;
  rms=Math.sqrt(sum(points.map((p,i)=>(((p[2]+p[3])/2)-mean)**2*weights[i]))/total);
 }
 let binned=true;
 const full=points.length===all.length&&input.fullMoments;
 if(full){
  if(![full.rate,full.mean,full.rms].every(finite)||full.rate<=0||full.rms<0||Math.abs(full.rate-rate)>Math.max(1e-12,rate*1e-9))throw Error('Saved PE moments do not match bins');
  mean=full.mean;rms=full.rms;binned=false;
 }
 return {kind:'summary',rate,mean,rms,rmsOverMean:mean>0?rms/mean:null,binned,n:points.length,
  range:[points[0][2],points.at(-1)[3]]};
}
function fit(input){
 if(input.kind==='summary')return histogramStats(input);
 const {kind,range}=input;let points=input.points.filter(p=>p.every(finite)&&p[0]>=range[0]&&p[0]<=range[1]).map(p=>[...p]);
 if(points.length<(input.histogram?5:3))throw Error('Choose a wider interval');
 points.sort((a,b)=>a[0]-b[0]);
 if(input.histogram&&points.some(p=>p.length<4||p[1]<0||p[3]<=p[2]))throw Error('Invalid histogram bins');
 if(kind==='auto'){
  const choices=input.histogram?['gaussian','lognormal']:['constant','power'];
  const results=choices.flatMap(k=>{try{return [fit({...input,kind:k})]}catch{return []}});
  if(!results.length)throw Error('No fit for this interval');
  results.sort((a,b)=>a.loss-b.loss);let chosen=results[0];
  // Prefer the simpler constant unless a trend reduces squared residuals by
  // at least 10%. This is a UI heuristic, not AIC or a hypothesis test.
  const constant=results.find(r=>r.kind==='constant');if(constant&&chosen.loss>=.9*constant.loss)chosen=constant;
  return {...chosen,suggested:true};
 }
 const result=input.histogram?peak(points,kind):trend(points,kind);
 return {...result,n:points.length,range:[points[0][2]??points[0][0],points.at(-1)[3]??points.at(-1)[0]],descriptive:true};
}
function predict(result,x,width=1){
 const p=result.params;
 if(result.kind==='constant')return p[0];
 if(result.kind==='power')return x>0?p[0]*(x/p[2])**(-p[1]):NaN;
 return p[0]*binShape(result.kind,Math.max(0,x-width/2),x+width/2,p.slice(1));
}
// Screen-space placement only; this never changes plotted observations.
// A small occupancy grid protects entire line segments and label rectangles.
function popoverBox({width,height,size,segments=[],boxes=[],anchor=[0,0],previous=null}){
 const [w,h]=size,g=6,pad=8;
 if(width<w+2*pad||height<h+2*pad)return null;
 const nx=Math.ceil(width/g),ny=Math.ceil(height/g),grid=new Uint8Array(nx*ny);
 const mark=(x,y,rw,rh)=>{
  const l=Math.max(0,Math.floor((x-5)/g)),r=Math.min(nx,Math.ceil((x+rw+5)/g));
  const t=Math.max(0,Math.floor((y-5)/g)),b=Math.min(ny,Math.ceil((y+rh+5)/g));
  for(let j=t;j<b;j++)for(let i=l;i<r;i++)grid[j*nx+i]=1;
 };
 for(const box of boxes)if(box.every(finite))mark(...box);
 for(const [x,y,xx,yy]of segments){
  if(![x,y,xx,yy].every(finite))continue;
  const dx=xx-x,dy=yy-y;let lo=0,hi=1,valid=true;
  for(const [p,q]of [[-dx,x],[dx,width-x],[-dy,y],[dy,height-y]]){
   if(p===0){if(q<0){valid=false;break;}}else{const t=q/p;if(p<0)lo=Math.max(lo,t);else hi=Math.min(hi,t);}
  }
  if(!valid||hi<lo)continue;
  const steps=Math.max(1,Math.ceil(Math.max(Math.abs(dx),Math.abs(dy))*(hi-lo)/g));
  for(let i=0;i<=steps;i++){const t=lo+(hi-lo)*i/steps;mark(x+t*dx,y+t*dy,1,1);}
 }
 const stride=nx+1,integral=new Uint32Array((nx+1)*(ny+1));
 for(let j=0;j<ny;j++){let row=0;for(let i=0;i<nx;i++){row+=grid[j*nx+i];integral[(j+1)*stride+i+1]=integral[j*stride+i+1]+row;}}
 const covered=(x,y)=>{
  if(x<pad||y<pad||x+w>width-pad||y+h>height-pad)return Infinity;
  const l=Math.floor(x/g),r=Math.ceil((x+w)/g),t=Math.floor(y/g),b=Math.ceil((y+h)/g);
  return integral[b*stride+r]-integral[t*stride+r]-integral[b*stride+l]+integral[t*stride+l];
 };
 if(previous&&covered(previous[0],previous[1])===0)return [previous[0],previous[1],w,h];
 const xs=[pad,width-pad-w],ys=[pad,height-pad-h];
 for(let x=pad+g;x<width-pad-w;x+=g)xs.push(x);
 for(let y=pad+g;y<height-pad-h;y+=g)ys.push(y);
 let best=null,score=Infinity,overlap=Infinity;
 for(const y of ys)for(const x of xs){
  const n=covered(x,y);
  const d=(x+w/2-anchor[0])**2+(y+h/2-anchor[1])**2;
  if(n<overlap||(n===overlap&&d<score)){overlap=n;score=d;best=[x,y,w,h];}
 }
 if(previous&&covered(previous[0],previous[1])<=overlap)return [previous[0],previous[1],w,h];
 return best;
}
const api={fit,predict,binShape,histogramStats,popoverBox};
if(typeof module!=='undefined'&&module.exports)module.exports=api;
else if(typeof document==='undefined')root.onmessage=e=>{try{root.postMessage({id:e.data.id,result:fit(e.data.input)})}catch(error){root.postMessage({id:e.data.id,error:error.message})}};
else root.CurveFit=api;
})(typeof self!=='undefined'?self:globalThis);
