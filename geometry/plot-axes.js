/* Offline plot gestures. Coordinates are exported from the existing figures. */
(function (root) {
  'use strict';
  const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
  const valid = v => typeof v === 'number' && Number.isFinite(v);
  const equal = (a, b) => a.length === b.length && a.every((v, i) => v === b[i]);

  function unit(v, state) {
    if (!valid(v) || (state.mode === 'log' && v <= 0)) return NaN;
    const f = state.mode === 'log' ? Math.log : v => v;
    return (f(v) - f(state.domain[0])) / (f(state.domain[1]) - f(state.domain[0]));
  }
  function scaledRange(state, dx, width) {
    const factor = Math.exp(clamp(-2.5 * dx / Math.max(40, width), -7, 7));
    const [lo, hi] = state.domain;
    if (state.mode === 'log') {
      const decades = clamp(Math.log10(hi / lo) * factor, .035, 16);
      return [lo, lo * Math.pow(10, decades)];
    }
    return [lo, lo + clamp((hi - lo) * factor, 1e-8, 1e100)];
  }
  // Keep the data coordinate beneath the cursor fixed, in linear or log space.
  function zoomRange(state, anchor, factor) {
    anchor = clamp(anchor, 0, 1);
    factor = clamp(factor, .05, 20);
    const log = state.mode === 'log';
    const [lo, hi] = state.domain.map(v => log ? Math.log(v) : v);
    const pivot = lo + anchor * (hi - lo);
    const span = clamp((hi - lo) * factor, log ? .035*Math.LN10 : 1e-8,
                       log ? 16*Math.LN10 : 1e100);
    const next = [pivot-anchor*span, pivot+(1-anchor)*span];
    if (log && (next[0] < -690 || next[1] > 690)) return [...state.domain];
    // Avoid turning an original zero into a tiny positive lower bound after
    // zooming out again; that would create a spurious many-decade log range.
    const epsilon=8*Number.EPSILON*Math.max(Math.abs(pivot),span);
    return next.map(v => log ? Math.exp(v) : Math.abs(v)<epsilon ? 0 : v);
  }
  function wheelFactor(event, extent) {
    const pixels = event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? extent : 1;
    const delta = Math.abs(event.deltaY) >= Math.abs(event.deltaX) ? event.deltaY : event.deltaX;
    return Math.exp(clamp(delta * pixels * (event.ctrlKey ? .01 : .0025), -1, 1));
  }
  function tagPoint(points, fraction, box, padding) {
    const [x,y,w,h] = box, finite = points.filter(([a,b]) => Number.isFinite(a) && Number.isFinite(b));
    if (!finite.length) return null;
    const left = Math.max(x,Math.min(...finite.map(p=>p[0])));
    const right = Math.min(x+w,Math.max(...finite.map(p=>p[0])));
    const target = left + fraction*(right-left);
    const inside = finite.filter(([a,b]) => a >= x+padding && a <= x+w-padding && b >= y+padding/3 && b <= y+h-padding/3);
    return inside.reduce((best,p) => !best || Math.abs(p[0]-target) < Math.abs(best[0]-target) ? p : best,null);
  }
  function logRange(domain, values) {
    const positive = values.filter(v => valid(v) && v > 0);
    if (!positive.length) return null;
    let lo = domain[0] > 0 ? domain[0] : positive.reduce((a,b)=>Math.min(a,b),Infinity);
    let hi = domain[1];
    // If zoom excluded every positive sample, start from the actual samples.
    if (hi <= lo) hi = positive.reduce((a,b)=>Math.max(a,b),-Infinity);
    if (hi <= lo) hi = lo * 10;
    return [lo, hi];
  }
  function ticks(state, count) {
    const [lo, hi] = state.domain;
    if (state.mode === 'log') {
      const span = Math.log10(hi / lo), result = [];
      const start = Math.floor(Math.log10(lo)), end = Math.ceil(Math.log10(hi));
      const step = Math.max(1, Math.ceil((end - start) / (count + 2)));
      const multiples = span < 1.5 ? [1, 2, 5] : [1];
      for (let p = start; p <= end; p += step)
        for (const m of multiples) {
          const v = m * Math.pow(10, p);
          if (v >= lo * (1 - 1e-12) && v <= hi * (1 + 1e-12)) result.push(v);
        }
      if (result.length < 2) return [lo, hi];
      return result;
    }
    const rough = (hi - lo) / count, power = Math.pow(10, Math.floor(Math.log10(rough)));
    const step = [1, 2, 2.5, 5, 10].find(v => v * power >= rough) * power;
    const result = [];
    for (let i = Math.ceil(lo / step - 1e-10); i <= Math.floor(hi / step + 1e-10); i++)
      result.push(Number((i * step).toPrecision(12)));
    return result;
  }
  function label(v) {
    if (v === 0) return '0';
    if (Math.abs(v) >= 1e5 || Math.abs(v) < .001) return v.toExponential(1).replace('.0e', 'e').replace('e+', 'e');
    return Number(v.toPrecision(4)).toString();
  }
  function makeStates(model) {
    const states = new Map();
    for (const axis of model.axes) for (const direction of ['x', 'y']) {
      const spec = axis[direction];
      if (!states.has(spec.group)) states.set(spec.group, {
        domain: [...spec.domain], mode: spec.mode, original: spec,
        remembered: {}, touched: false,
      });
    }
    return states;
  }
  function reset(state) {
    state.domain = [...state.original.domain]; state.mode = state.original.mode;
    state.remembered = {}; state.touched = false;
  }
  function toggle(state, values) {
    const next = state.mode === 'log' ? 'linear' : 'log';
    const domain = state.remembered[next] || (next === 'log' ? logRange(state.domain, values) : state.original.domain);
    if (!domain) return false;
    state.remembered[state.mode] = [...state.domain];
    state.mode = next; state.domain = [...domain];
    state.touched = next !== state.original.mode || !equal(state.domain, state.original.domain);
    return true;
  }
  const math = {unit, scaledRange, zoomRange, wheelFactor, tagPoint, logRange, ticks, label, makeStates, reset, toggle};
  if (typeof module !== 'undefined' && module.exports) module.exports = math;
  if (!root.document) return;
  root.PlotAxesMath = math;

  const catalog = root.CONVERGENCE_PLOTS;
  if (!catalog) return;
  const figures = new Map(), cache = new Map();
  const rgba = c => `rgba(${c.slice(0,3).map(v => Math.round(255 * v)).join(',')},${c[3]})`;

  class Figure {
    constructor(element) {
      this.element = element; this.img = element.querySelector('img');
      this.canvas = document.createElement('canvas'); this.canvas.className = 'plot-overlay';
      this.canvas.setAttribute('aria-hidden', 'true'); element.append(this.canvas);
      this.ctx = this.canvas.getContext('2d'); this.handles = [];
      this.resize = new ResizeObserver(() => this.requestDraw()); this.resize.observe(this.img);
      this.img.addEventListener('load', () => { this.mount(); this.requestDraw(); });
      element.addEventListener('dblclick', e => {
        const axis = this.axisAt(e);
        if (!axis) {
          const bounds = this.img.getBoundingClientRect();
          const y = (e.clientY-bounds.top) * this.model.height / bounds.height;
          if (y >= 0 && y < Math.min(...this.model.axes.map(a=>a.rect[1]))-10) {
            e.preventDefault();for(const state of this.states.values())reset(state);this.requestDraw();
          }
          return;
        }
        e.preventDefault(); e.stopPropagation();
        reset(this.states.get(axis.x.group)); reset(this.states.get(axis.y.group));
        this.requestDraw();
      });
      this.mount();
    }
    axisAt(event) {
      if (!this.model) return null;
      const bounds = this.img.getBoundingClientRect();
      if (!bounds.width || !bounds.height) return null;
      const x = (event.clientX - bounds.left) * this.model.width / bounds.width;
      const y = (event.clientY - bounds.top) * this.model.height / bounds.height;
      return this.model.axes.find(a => {
        const [left,top,width,height] = a.rect;
        return x > left && x < left+width && y > top && y < top+height;
      });
    }
    mount() {
      const name = this.img.getAttribute('src').split(/[?#]/)[0].split('/').pop().replace('.png', '');
      if (name === this.name) return;
      const model = catalog.figures[name]; if (!model) return;
      this.name = name; this.model = model;
      if (!cache.has(name)) cache.set(name, makeStates(model));
      this.states = cache.get(name);
      this.handles.forEach(h => h.remove()); this.handles = [];
      model.axes.forEach((axis, index) => ['x', 'y'].forEach(direction => {
        const h = document.createElement('button'); h.type = 'button'; h.className = 'axis-handle';
        h.dataset.axis = direction; h.dataset.panel = index;
        h.setAttribute('aria-label', `${axis.label || 'Panel '+(index+1)} ${direction} axis: wheel or pinch to zoom at cursor, drag ${direction === 'x' ? 'left/right' : 'up/down'} to rescale, double click or Enter for log scale, Escape to reset`);
        h.addEventListener('click', e => { e.preventDefault(); e.stopPropagation(); });
        h.addEventListener('dblclick', e => {
          e.preventDefault(); e.stopPropagation();
          const state = this.states.get(axis[direction].group);
          if (e.shiftKey) reset(state); else toggle(state, this.values(axis[direction].group, direction));
          this.requestDraw();
        });
        h.addEventListener('pointerdown', e => this.startDrag(e, h, axis, direction));
        h.addEventListener('wheel', e => {
          e.preventDefault(); e.stopPropagation();
          if (this.gesture) return;
          const state = this.states.get(axis[direction].group);
          state.domain = zoomRange(state,this.anchorAt(e,axis,direction),wheelFactor(e,direction==='x'?this.img.clientWidth:this.img.clientHeight));
          state.touched = true; this.requestDraw();
        }, {passive:false});
        // WebKit exposes native Mac trackpad pinch as gesture events instead
        // of the ctrl+wheel events used by Chromium and Firefox.
        h.addEventListener('gesturestart', e => {
          e.preventDefault();e.stopPropagation();
          const state = this.states.get(axis[direction].group);
          this.gesture={state,start:{mode:state.mode,domain:[...state.domain]},anchor:this.anchorAt(e,axis,direction)};
        }, {passive:false});
        h.addEventListener('gesturechange', e => {
          if (!this.gesture || !valid(e.scale) || e.scale<=0) return;
          e.preventDefault();e.stopPropagation();
          const g=this.gesture;g.state.domain=zoomRange(g.start,g.anchor,1/e.scale);g.state.touched=true;this.requestDraw();
        }, {passive:false});
        h.addEventListener('gestureend', e => {e.preventDefault();e.stopPropagation();this.gesture=null;}, {passive:false});
        h.addEventListener('keydown', e => {
          h.classList.remove('pointer-focus');
          const state = this.states.get(axis[direction].group);
          if (e.key === 'Escape' || e.key === 'Home') reset(state);
          else if (e.key === 'Enter' || e.key === ' ') toggle(state, this.values(axis[direction].group, direction));
          else if (e.key === (direction === 'x' ? 'ArrowRight' : 'ArrowUp') || e.key === (direction === 'x' ? 'ArrowLeft' : 'ArrowDown')) {
            state.domain = scaledRange(state, e.key === 'ArrowRight' || e.key === 'ArrowUp' ? 8 : -8, 200); state.touched = true;
          } else return;
          e.preventDefault(); e.stopPropagation(); this.requestDraw();
        });
        this.element.append(h); this.handles.push(h);
      }));
      this.requestDraw();
    }
    anchorAt(event,axis,direction) {
      const bounds=this.img.getBoundingClientRect(),[x,y,w,h]=axis.rect;
      return direction==='x' ? ((event.clientX-bounds.left)*this.model.width/bounds.width-x)/w
        : 1-((event.clientY-bounds.top)*this.model.height/bounds.height-y)/h;
    }
    values(group, direction) {
      const arrays = new Set();
      for (const axis of this.model.axes) if (axis[direction].group === group)
        for (const line of axis.lines) if (line.space[direction === 'x' ? 0 : 1] === 'data') arrays.add(line[direction]);
      // Reduce each array first: avoid spread argument limits on large plots.
      return [...arrays].flatMap(id => catalog.arrays[id].filter(v => valid(v) && v > 0));
    }
    startDrag(e, handle, axis, direction) {
      if (e.button !== 0) return;
      e.preventDefault(); e.stopPropagation(); handle.classList.add('pointer-focus'); handle.focus({preventScroll:true});
      handle.setPointerCapture(e.pointerId);
      const state = this.states.get(axis[direction].group), start = {mode:state.mode, domain:[...state.domain]};
      const vertical = direction === 'y', origin = vertical ? e.clientY : e.clientX;
      const extent = vertical ? axis.rect[3] * this.img.clientHeight / this.model.height
                              : axis.rect[2] * this.img.clientWidth / this.model.width;
      const move = event => {
        const distance = vertical ? origin - event.clientY : event.clientX - origin;
        if (Math.abs(distance) < 3 && !this.dragging) return;
        this.dragging = true; this.element.classList.add('axis-dragging');
        state.domain = scaledRange(start, distance, extent); state.touched = true; this.requestDraw();
      };
      const end = () => {
        this.dragging = false; this.element.classList.remove('axis-dragging');
        handle.removeEventListener('pointermove', move); handle.removeEventListener('pointerup', end);
        handle.removeEventListener('pointercancel', end); handle.removeEventListener('lostpointercapture', end);
      };
      handle.addEventListener('pointermove', move); handle.addEventListener('pointerup', end);
      handle.addEventListener('pointercancel', end); handle.addEventListener('lostpointercapture', end);
    }
    requestDraw() {
      if (!this.frame) this.frame = requestAnimationFrame(() => {this.frame = null; this.draw();});
    }
    positionHandles() {
      const w = this.img.clientWidth, h = this.img.clientHeight;
      const sx = w / this.model.width, sy = h / this.model.height;
      this.handles.forEach(handle => {
        const axis = this.model.axes[Number(handle.dataset.panel)];
        const state = this.states.get(axis[handle.dataset.axis].group);
        handle.dataset.mode = state.mode;
        handle.dataset.range = JSON.stringify(state.domain);
        handle.setAttribute('aria-pressed', String(state.mode === 'log'));
        const [x, y, pw, ph] = axis.rect;
        const box = handle.dataset.axis === 'x'
          ? [x*sx+9, (y+ph)*sy-5, pw*sx-9, Math.max(18, 30*sy)]
          : [x*sx-Math.max(19, 42*sx), y*sy, Math.max(25, 49*sx), ph*sy-5];
        ['left','top','width','height'].forEach((p,i) => handle.style[p] = box[i]+'px');
      });
    }
    draw() {
      if (!this.model || !this.img.complete) return;
      const w = this.img.clientWidth, h = this.img.clientHeight;
      if (!w || !h) {this.canvas.width = this.canvas.height = 1; return;}
      this.positionHandles();
      const dpr = Math.min(root.devicePixelRatio || 1, 2);
      const width = Math.round(w*dpr), height = Math.round(h*dpr);
      if (this.canvas.width !== width || this.canvas.height !== height) {this.canvas.width=width;this.canvas.height=height;}
      const ctx = this.ctx;
      ctx.setTransform(width/this.model.width, 0, 0, height/this.model.height, 0, 0);
      ctx.clearRect(0,0,this.model.width,this.model.height);
      let changed = false;
      for (const axis of this.model.axes) {
        const xs=this.states.get(axis.x.group), ys=this.states.get(axis.y.group);
        if (!xs.touched && !ys.touched) continue;
        changed = true; this.drawAxis(axis,xs,ys);
      }
      if (changed) {ctx.fillStyle='white';for(const b of this.model.range_notes || []) ctx.fillRect(...b);}
      root.FitLens?.drawFigure(this);
    }
    drawAxis(axis, xs, ys) {
      const ctx=this.ctx,[x,y,w,h]=axis.rect;
      ctx.save(); ctx.fillStyle='white';
      ctx.fillRect(x-3,y-3,w+6,h+6);
      for (const direction of ['x','y']) if(axis[direction].mask)ctx.fillRect(...axis[direction].mask);
      ctx.fillStyle=rgba(axis.face);ctx.fillRect(x,y,w,h);
      const xt=ticks(xs,5),yt=ticks(ys,4);
      const point=(v,space,direction)=> {
        const fraction = space === 'axes' ? v : unit(v,direction===0?xs:ys);
        return direction===0 ? x+fraction*w : y+h-fraction*h;
      };
      ctx.save();ctx.beginPath();ctx.rect(x,y,w,h);ctx.clip();
      ctx.lineWidth=1.5;ctx.strokeStyle='rgba(0,0,0,.16)';ctx.beginPath();
      for(const v of xt){const p=point(v,'data',0);ctx.moveTo(p,y);ctx.lineTo(p,y+h);}
      for(const v of yt){const p=point(v,'data',1);ctx.moveTo(x,p);ctx.lineTo(x+w,p);}
      ctx.stroke();ctx.restore();
      const tags=[];
      for(const line of axis.lines){
        ctx.save();if(line.clip){ctx.beginPath();ctx.rect(x,y,w,h);ctx.clip();}
        ctx.globalAlpha=line.alpha;ctx.strokeStyle=rgba(line.color);ctx.lineWidth=line.width;
        ctx.setLineDash(line.dash);ctx.lineJoin='round';ctx.lineCap='butt';
        const xv=catalog.arrays[line.x],yv=catalog.arrays[line.y],n=Math.max(xv.length,yv.length);
        const points=[];
        for(let i=0;i<n;i++){
          const a=xv[xv.length===1?0:i],b=yv[yv.length===1?0:i];
          points.push(valid(a)&&valid(b)?[point(a,line.space[0],0),point(b,line.space[1],1)]:[NaN,NaN]);
        }
        if(line.stroke){
          ctx.beginPath();let pen=false;
          for(const [a,b]of points){if(!Number.isFinite(a)||!Number.isFinite(b)){pen=false;continue;}
            if(pen)ctx.lineTo(a,b);else ctx.moveTo(a,b);pen=true;}
          if(line.white_halo){ctx.save();ctx.strokeStyle='white';ctx.lineWidth+=line.white_halo;ctx.stroke();ctx.restore();}
          ctx.stroke();
        }
        if(['o','s','^'].includes(line.marker)){
          ctx.setLineDash([]);ctx.strokeStyle=rgba(line.marker_edge);ctx.fillStyle=rgba(line.marker_face);ctx.lineWidth=line.marker_width;
          const radius=line.marker_size/2;ctx.beginPath();
          for(const[a,b]of points){
            if(!Number.isFinite(a)||!Number.isFinite(b)||a<x-radius||a>x+w+radius||b<y-radius||b>y+h+radius)continue;
            if(line.marker==='o'){ctx.moveTo(a+radius,b);ctx.arc(a,b,radius,0,2*Math.PI);}
            else if(line.marker==='^'){ctx.moveTo(a,b-radius);ctx.lineTo(a+radius,b+radius);ctx.lineTo(a-radius,b+radius);ctx.closePath();}
            else ctx.rect(a-radius,b-radius,2*radius,2*radius);
          }
          if(line.marker_face[3])ctx.fill();if(line.marker_edge[3])ctx.stroke();
        }
        if(line.tag)tags.push({line,points});
        ctx.restore();
      }
      // Labels belong to their curves, so re-anchor them within the visible
      // samples after zooming. A small white halo cuts them into the line.
      ctx.save();ctx.beginPath();ctx.rect(x,y,w,h);ctx.clip();
      const zoom = state => {
        const f=state.mode==='log'?Math.log:v=>v,base=state.original.domain;
        const span=f(state.domain[1])-f(state.domain[0]);
        return state.mode===state.original.mode ? (f(base[1])-f(base[0]))/span : 1;
      };
      const growth=clamp(Math.sqrt(Math.max(zoom(xs),zoom(ys))),1,1.8);
      for(const {line,points} of tags){
        const size=line.tag.size*growth;
        ctx.font=`600 ${size}px 'DejaVu Sans',Arial,sans-serif`;
        const p=tagPoint(points,line.tag.fraction,axis.rect,ctx.measureText(line.tag.text).width/2+5);
        if(!p)continue;
        ctx.textAlign='center';ctx.textBaseline='middle';ctx.lineJoin='round';ctx.setLineDash([]);
        ctx.strokeStyle='white';ctx.lineWidth=4;ctx.strokeText(line.tag.text,...p);
        ctx.fillStyle=rgba(line.color);ctx.fillText(line.tag.text,...p);
      }
      ctx.restore();
      for(const[direction,values,state]of [['x',xt,xs],['y',yt,ys]]){
        const spec=axis[direction];ctx.strokeStyle=rgba(spec.color);ctx.fillStyle='#202020';ctx.lineWidth=spec.width;ctx.setLineDash([]);
        if(spec.spine){ctx.beginPath();ctx.moveTo(x,y+(direction==='x'?h:0));ctx.lineTo(direction==='x'?x+w:x,y+h);ctx.stroke();}
        for(const v of values){
          const p=point(v,'data',direction==='x'?0:1);
          if(spec.ticks){ctx.beginPath();if(direction==='x'){ctx.moveTo(p,y+h);ctx.lineTo(p,y+h+8);}else{ctx.moveTo(x,p);ctx.lineTo(x-8,p);}ctx.stroke();}
          if(spec.labels){
            const text=label(v);ctx.font=`${spec.size}px 'DejaVu Sans',Arial,sans-serif`;
            if(direction==='x'){ctx.textAlign='center';ctx.textBaseline='top';ctx.fillText(text,p,y+h+13,Math.max(30,w/values.length-6));}
            else{ctx.textAlign='right';ctx.textBaseline='middle';const room=spec.mask?x-13-spec.mask[0]:60;ctx.fillText(text,x-13,p,Math.max(22,room));}
          }
        }
      }
      for(const t of axis.texts){
        ctx.font=`${t.size}px 'DejaVu Sans',Arial,sans-serif`;ctx.fillStyle=rgba(t.color);
        ctx.textAlign={left:'left',right:'right',center:'center'}[t.ha]||'left';
        ctx.textBaseline={top:'top',bottom:'bottom',center:'middle',baseline:'alphabetic'}[t.va]||'alphabetic';
        ctx.fillText(t.text,point(t.position[0],t.space[0],0),point(t.position[1],t.space[1],1));
      }
      ctx.restore();
    }
  }
  function refresh(){document.querySelectorAll('.figure').forEach(el=>{
    if(!figures.has(el))figures.set(el,new Figure(el));
    else{figures.get(el).mount();figures.get(el).requestDraw();}
  });}
  document.querySelector('[data-reset-plots]')?.addEventListener('dblclick',e=>{
    e.preventDefault();
    for(const states of cache.values())for(const state of states.values())reset(state);
    for(const figure of figures.values())figure.requestDraw();
  });
  root.ScalePlots={refresh,figures,catalog};refresh();
})(typeof window === 'undefined' ? globalThis : window);
