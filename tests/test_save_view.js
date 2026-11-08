const assert=require('node:assert/strict');
const fs=require('node:fs');
// Plot gestures preserve the value under the cursor in either display mode.
const axes=require('../geometry/plot-axes.js'),curves=require('../geometry/pe-fit-core.js');
for(const state of [{domain:[-10,30],mode:'linear'},{domain:[.1,1000],mode:'log'}]){
 const anchor=.37,range=axes.zoomRange(state,anchor,.4),f=state.mode==='log'?Math.log:x=>x;
 const pivot=f(state.domain[0])+anchor*(f(state.domain[1])-f(state.domain[0]));
 assert.ok(Math.abs(f(range[0])+anchor*(f(range[1])-f(range[0]))-pivot)<1e-12);
 assert.ok(axes.scaledRange(state,30,300)[1]<state.domain[1]);
}
assert.equal(axes.logRange([-10,0],[-8,0]),null);
// Fit integrated bins, including negative coordinates and unequal widths.
const edges=[-10,-8,-6,-5,-4,-3,-2,-1,0,1,2,4,7];
const bins=edges.slice(0,-1).map((lo,i)=>[(lo+edges[i+1])/2,1200*curves.binShape('gaussian',lo,edges[i+1],[-2,1.5]),lo,edges[i+1]]);
const original=JSON.stringify(bins);
for(const scale of [1,1.602176634e-10]){
 const result=curves.fit({kind:'gaussian',histogram:true,points:bins.map(p=>[p[0],p[1]*scale,p[2],p[3]]),range:[-10,7]});
 assert.ok(Math.abs(result.params[1]+2)<.001);assert.ok(Math.abs(result.params[2]-1.5)<.001);
 assert.ok(Math.abs(result.params[0]/scale/1200-1)<.001);
}
assert.equal(JSON.stringify(bins),original);
const {signalAxisValues}=require('../geometry/view.js');
const bounds={'plot-y-min':'0','plot-y-max':'1e12','plot-x-max':'10'};
const amps=signalAxisValues(bounds,'pe','current');
assert.equal(Number(amps['plot-y-max']),160.2176634);assert.equal(amps['plot-x-max'],'10');
assert.ok(Math.abs(Number(signalAxisValues(amps,'current','pe')['plot-y-max'])/1e12-1)<1e-12);
assert.equal(signalAxisValues(bounds,'pe','hits')['plot-y-max'],'');
const vm=require('node:vm');
// Enlarging a plot must not relabel the correlation viewer's ring buttons.
const source=fs.readFileSync(require.resolve('../geometry/view.js'),'utf8');
const control=(text,dataset)=>({textContent:text,dataset,setAttribute(k,v){this[k]=v;},getAttribute(k){return this[k];}});
const panel=control('',{panelFocus:'distribution'}),ring=control('Ring 5',{focus:'ring5'});
const layout={panelVisible:{scene:true,selection:true,distribution:true},focusedPanel:'distribution',
 workspace:{classList:{toggle(){}}},requestAnimationFrame(){},
 document:{querySelectorAll:s=>s==='[data-panel-focus]'?[panel]:s==='[data-focus]'?[panel,ring]:[]}};
vm.runInNewContext(source.slice(source.indexOf('function updatePanels(){'),source.indexOf('function restorePanels(){'))+'\nupdatePanels();',layout);
assert.equal(panel.textContent,'↙');assert.equal(ring.textContent,'Ring 5');assert.equal(ring['aria-label'],undefined);
const {formatNumber,colorFraction,signalData}=require('../geometry/hit-map-3d.js');
assert.equal(formatNumber(1.16e13),'1.16e13');
assert.equal(formatNumber(1.73e10),'1.73e10');
assert.equal(formatNumber(.05),'0.05');
assert.equal(formatNumber(1e-6),'1e-6');
assert.equal(formatNumber(0),'0');
for(const log of [false,true]) for(const f of [0,1e-6,.001,.1,1]) assert.ok(Math.abs(colorFraction(f*1e12,1e12,log)-colorFraction(f*160.2176634,160.2176634,log))<1e-14);
const raw={weighting:'pe',unit:'PE/s',cells:[[0,0,0,1e9]],histograms:{x:{counts:[1e9],underflow:2e9,overflow:3e9}},histograms2d:{xy:{cells:[[0,0,1e9]]}},signal_total:1e13,detector_signal:{1:1e9},particle_signal:{11:1e9}};
const before=JSON.stringify(raw),current=signalData(raw,'current');
const close=(a,b)=>assert.ok(Math.abs(a-b)<Math.max(1,Math.abs(b))*1e-12,`${a} != ${b}`);
close(current.signal_total,1602.176634);
for(const n of [current.cells[0][3],current.histograms.x.counts[0],current.histograms2d.xy.cells[0][2],current.detectors[1],current.particles[11],current.detector_signal[1],current.particle_signal[11]])close(n,.1602176634);
close(current.histograms.x.underflow,.3204353268);
assert.equal(current.unit,'nA');assert.equal(JSON.stringify(raw),before);
assert.equal(signalData(raw,'pe').signal_total,1e13);
assert.deepEqual(signalData({weighting:'hits',cells:[[0,0,0,9]]},'current'),{weighting:'hits',cells:[[0,0,0,9]]});
let snapshot=null,blob=null,downloads=0,filename='';
const button={textContent:'Save view'};
const context={Blob,Date,JSON,URL:{createObjectURL:b=>(blob=b,'blob:test'),revokeObjectURL:()=>{}},setTimeout:()=>{},CustomEvent:class{constructor(type,args){this.detail=args.detail;}},document:{getElementById:()=>button,dispatchEvent:e=>{e.detail.snapshot=snapshot;},createElement:()=>({set download(v){filename=v;},click(){downloads++;}})}};
vm.runInNewContext(fs.readFileSync(require.resolve('../geometry/save-view.js'),'utf8'),context);
assert.equal(downloads,0);button.onclick();assert.equal(downloads,0);
(async()=>{for(const page of ['detector','deconvolution']){snapshot={page,result:{value:123},selection:{cut:'hit.pid == 11'}};button.onclick();const saved=JSON.parse(await blob.text());assert.equal(saved.schema,'moller-view-v1');assert.deepEqual(saved.result,snapshot.result);assert.deepEqual(saved.selection,snapshot.selection);assert.equal(filename,`moller-${page}-view.json`);}assert.equal(downloads,2);console.log('Save view and signal unit tests passed');})().catch(e=>{console.error(e);process.exitCode=1;});
