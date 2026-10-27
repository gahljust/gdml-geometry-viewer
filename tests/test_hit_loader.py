"""Known-answer ROOT fixtures and the geometry server's hit-loading path.

python3 tests/test_hit_loader.py
Uses develop's hit schema; all scratch files are deleted on exit.
"""
import collections
import importlib.util
import json
import math
import os
import re
import subprocess
import tempfile
import threading
import time
from pathlib import Path
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parents[1]
REMOLL = Path(os.environ['REMOLL_DIR']).expanduser().resolve() if os.environ.get('REMOLL_DIR') else next((p for p in HERE.parents if (p / 'include/remolltypes.hh').is_file()), None)
if REMOLL is None or not (REMOLL / 'include/remolltypes.hh').is_file():
    raise SystemExit('Set REMOLL_DIR to your remoll source directory.')
spec = importlib.util.spec_from_file_location('geometry_server', HERE / 'geometry/run.py')
app = importlib.util.module_from_spec(spec)
spec.loader.exec_module(app)


def run_root(macro, cwd):
    result = subprocess.run(['root', '-n', '-l', '-b', '-q', str(macro)], cwd=cwd,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr


# Compare every preset combination against the previous selection semantics.
# Includes range edges, all suffixes, neutral particles, threshold equality,
# backwards/zero momentum and incompatible tile/region combinations.
subprocess.run(['node','-'],input=r"""
const assert=require('node:assert/strict');
const {buildCut,validatePlotRange}=require('./geometry/view.js');
assert.deepEqual(validatePlotRange(-100,100,false),[-100,100]);
assert.deepEqual(validatePlotRange(.001,20000,true),[.001,20000]);
for(const args of [[-1e308,1e308,false],[1e308,1.1e308,false],[5e-324,1e308,true],[0,10,true],[2,1,false],[1,Infinity,false]])assert.throws(()=>validatePlotRange(...args));
const endings=[10,20,30,40,51,52,53,60], charged=pid=>[11,-11,2212].includes(pid);
const ids=[-1,0,73000,73029,110000,111999,112000,150000,151999,152000,70001,150254];
for(const base of [110000,150000])for(let s=1;s<=14;s++)for(const end of endings)ids.push(base+100*s+end);
for(let id=73001;id<=73028;id++)ids.push(id);
const hits=ids.map((det,i)=>({det,pid:[11,-11,22,2212,2112][i%5],k:[.5,1,2][i%3],mtrid:i%3,pz:[-1,0,1][i%3]}));
function previous(c,h){
 const main=((h.det>=110000&&h.det<112000)||(h.det>=150000&&h.det<152000))&&endings.includes(h.det%100);
 const sm=h.det>=73001&&h.det<=73028;
 if(c.system==='showermax'&&!sm)return false;
 if(c.system!=='all'&&c.system!=='showermax'){
  if(!main)return false;
  if(c.mode!=='raw'&&(!charged(h.pid)||h.k<1))return false;
  if(c.system.startsWith('ring')){const r=Number(c.system.slice(4));if(!(r===5?[51,52,53].includes(h.det%100):h.det%100===r*10))return false;}
 }
 if((c.mode==='secondary'||c.ancestry==='secondary')?h.mtrid<=0:c.ancestry==='primary'&&h.mtrid!==0)return false;
 if(c.direction!=='all'&&!(c.direction==='forward'?h.pz>0:h.pz<0))return false;
 if(c.particle!=='all'&&h.pid!==Number(c.particle))return false;
 if(c.region!=='all'){
  const mr=c.region==='transition'?h.det>=110000&&h.det<112000:h.det>=150000&&h.det<152000&&(c.region==='open'?h.det%200<100:h.det%200>=100);
  const sr=c.region==='transition'?h.det%2===0:h.det%4===(c.region==='open'?3:1);
  if(!(c.system==='showermax'?sr:c.system==='all'?(main&&mr)||(sm&&sr):mr))return false;
 }
 return c.tile===null||h.det===c.tile;
}
let configurations=0,checks=0;
for(const system of ['all','main','ring1','ring2','ring3','ring4','ring5','ring6','showermax'])
for(const mode of ['raw','hits','secondary'])for(const region of ['all','open','closed','transition'])
for(const ancestry of ['all','primary','secondary'])for(const direction of ['all','forward','backward'])
for(const particle of ['all','11','-11','22','2212','2112'])for(const tile of [null,150252,110451,73003]){
 const c={system,mode,region,ancestry,direction,particle,tile},cut=buildCut(c);
 const evaluate=new Function('hit','charged','return Boolean('+cut+');');
 for(const h of hits){assert.equal(evaluate(h,charged),previous(c,h),JSON.stringify({c,h,cut}));checks++;}
 configurations++;
 if(tile!==null&&cut!=='0')assert.ok(!cut.includes('hit.det%')&&!cut.includes('hit.det >='));
 if(system==='ring5'&&tile===null)assert.ok(!cut.includes('hit.det%100 == 10'));
}
const c={mode:'raw',system:'ring5',region:'open',ancestry:'all',direction:'all',particle:'all',tile:null};
assert.equal(buildCut(c),'hit.det >= 150000 && hit.det < 152000 && hit.det%200 < 100 && (hit.det%100 == 51 || hit.det%100 == 52 || hit.det%100 == 53)');
console.log(`PASS ${configurations} preset combinations, ${checks} hit decisions: simplified cuts match previous selections`);
""",cwd=HERE,text=True,check=True)

app.OUTPUT.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(prefix='hit-loader-test-', dir=app.OUTPUT) as temporary:
    scratch = Path(temporary)
    app.OUTPUT = scratch
    # No custom remoll library and no saved macro, rate or normalization metadata.
    header = REMOLL / 'include/remolltypes.hh'
    try:
        stock = subprocess.check_output(['git', '-C', str(REMOLL), 'show', 'develop:include/remolltypes.hh'], text=True, stderr=subprocess.DEVNULL)
    except (OSError, subprocess.CalledProcessError):
        stock = header.read_text()
    (scratch / 'remolltypes.hh').write_text(stock)
    fixture = scratch / 'fixtures.C'
    fixture.write_text(r'''
#include "remolltypes.hh"
#include <TFile.h>
#include <TTree.h>
#include <TInterpreter.h>
#include <limits>
void fixtures() {
  {TFile f("no_positions.root","RECREATE");TTree t("T","IDs only");
   int det[2]={150252,73003};t.Branch("hit.det",det,"hit.det[2]/I");t.Fill();t.Write();}
  for(int kind=0;kind<3;++kind){
    const char* names[]={"short_cut","fractional_pid","control_name"};
    TFile f((std::string(names[kind])+".root").c_str(),"RECREATE");TTree t("T","adversarial arrays");
    int n=2, shortN=1,det[2]={73001,73002},pid[2]={11,22},trid[2]={1,0};double badPid[2]={11.9,22.7},extra=7;
    t.Branch("n",&n,"n/I");t.Branch("shortN",&shortN,"shortN/I");t.Branch("hit.det",det,"hit.det[n]/I");
    if(kind==1)t.Branch("hit.pid",badPid,"hit.pid[n]/D");else t.Branch("hit.pid",pid,"hit.pid[n]/I");
    if(kind==0)t.Branch("hit.trid",trid,"hit.trid[shortN]/I");
    if(kind==2)t.Branch("bad\nname",&extra,"bad\nname/D");
    t.Fill();t.Fill();t.Write();
  }
  gInterpreter->GenerateDictionary("vector<remollGenericDetectorHit_t>","remolltypes.hh;vector");
  for (std::string mode : {"valid", "empty", "zero_hits", "missing_tree", "missing_det", "wrong_det", "short_field", "scalar", "other_array", "views", "bad_position", "hist_edges", "plot_known"}) {
    TFile f((mode+".root").c_str(),"RECREATE");
    TTree t(mode=="missing_tree"?"wrong":"T","Test");
    std::vector<remollGenericDetectorHit_t> hits;
    std::vector<double> wrong;
    int dummy=0;
    double rate=0, q2=10000, w2=1000000, asym=-30;
    if(mode=="plot_known"){t.Branch("ev.Q2",&q2);t.Branch("ev.W2",&w2);t.Branch("ev.A",&asym);}
    std::vector<double> other;
    if(mode=="scalar") t.Branch("rate", &rate);
    if(mode=="other_array") t.Branch("other", &other);
    if(mode=="missing_det") t.Branch("dummy",&dummy);
    else if(mode=="wrong_det") t.Branch("hit.det",&wrong);
    else t.Branch("hit",&hits,32000,99);
    for (int e=0; e<(mode=="empty"?0:3); ++e) {
      hits.clear(); wrong.clear(); other.assign(6,1); rate=10*(e+1);
      if (e!=1 && mode!="zero_hits")
        for(int id : {150252,150252,73003,70001,999,-1}) {
          remollGenericDetectorHit_t h{};
          h.det=id; h.pid=e==0?22:11; h.k=.2; h.pz=-1; h.x=1.5;
          if(mode=="views") {h.k=2;h.mtrid=hits.size()%2;h.pz=id==73003?1:-1;}
          if(mode=="hist_edges") {double values[]={.0001,.001,20000,std::numeric_limits<double>::quiet_NaN(),1,-1};h.k=values[hits.size()];}
          if(mode=="bad_position") h.x=std::numeric_limits<double>::quiet_NaN();
          hits.push_back(h); wrong.push_back(id);
        }
      if(mode=="plot_known") {
        for(size_t i=0;i<hits.size();++i){auto& h=hits[i];h.x=100+100*i;h.y=-100-200*double(i);h.vx=3;h.vy=4;h.vz=1000;h.t=10*i;h.edep=.5*i;h.sx=.25;h.k=10+i;}
        q2=e==0?10000:20000;w2=e==0?1000000:-1;
      }
      t.Fill();
    }
    if(mode=="short_field") t.GetBranch("hit.x")->SetEntries(1);
    t.Write();
  }
}
''')
    run_root(fixture, scratch)
    (scratch / 'corrupt.root').write_text('not a ROOT file')
    # Quoting and spaces in paths must work without building a shell command.
    (scratch / 'quoted " input.root').write_bytes((scratch / 'valid.root').read_bytes())
    # Two known square faces, normal +z; hits at (1.5,0,0) occupy cell (0,0).
    (scratch/'surface-faces.tsv').write_text(''.join(f'{i} {det} 0 0 0 1 0 0 0 1 0 0 0 1 4 0 0 10 0 10 10 0 10\n' for i,det in enumerate([150252,73003])))
    app.state['ready']=True
    server = app.ThreadingHTTPServer(('127.0.0.1', 0), app.Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    url = f'http://127.0.0.1:{server.server_port}'

    def api(path, data=None):
        request = Request(url + path, data=None if data is None else json.dumps(data).encode(),
                          headers={'Content-Type': 'application/json'})
        with urlopen(request, timeout=10) as response:
            return json.load(response)

    def load(name, cut='1', binning=None):
        api('/load-hits', {**({'paths': [str(scratch / n) for n in name]} if isinstance(name,list) else {'path': str(scratch / name)}), 'cut': cut, 'binning': binning or {}})
        deadline = time.monotonic() + 30
        while time.monotonic() < deadline:
            state = api('/status')
            if not state['busy']:
                return state['hits']
            time.sleep(.03)
        raise AssertionError('Hit load timed out')

    try:
        from urllib.error import HTTPError
        malformed=[('/load-hits',x) for x in ([],7,'text')]
        malformed += [('/load-hits',{'paths':[str(scratch/'valid.root')],'cut':x}) for x in (None,7,[],{})]
        malformed += [('/load',{'root':x}) for x in (None,7,[],{})]
        for endpoint,payload in malformed:
            try: api(endpoint,payload)
            except HTTPError as error:
                assert error.code==400 and json.load(error)['error']
            else: raise AssertionError(('Malformed request accepted',payload))
            assert not api('/status')['busy']
        request=Request(url+'/load-hits',data=b'null',headers={'Content-Type':'application/json'})
        try: urlopen(request,timeout=10)
        except HTTPError as error: assert error.code==400 and json.load(error)['error']
        else: raise AssertionError('Null request accepted')
        for value in [b'{"x":NaN}',b'{"x":Infinity}',b'{"x":-Infinity}']:
            try: urlopen(Request(url+'/load-hits',data=value,headers={'Content-Type':'application/json'}),timeout=10)
            except HTTPError as error: assert error.code==400 and json.load(error)['error']
            else: raise AssertionError('Nonfinite JSON accepted')
        for path in ['bad\nfile.root','bad\rfile.root','bad\0file.root']:
            try: api('/load-hits',{'paths':[path]})
            except HTTPError as error: assert error.code==400
            else: raise AssertionError('Ambiguous path accepted')
        print('PASS malformed request types return JSON errors and server stays responsive',flush=True)
        for endpoint in ['/status','/hits','/geometry','/responses']:
            for headers in [{'Host':'untrusted.invalid:'+str(server.server_port)},
                            {'Origin':'https://untrusted.invalid'}, {'Sec-Fetch-Site':'cross-site'}]:
                try: urlopen(Request(url+endpoint,headers=headers),timeout=10)
                except HTTPError as error: assert error.code==403
                else: raise AssertionError('Foreign host/origin accepted')
        try: urlopen(Request(url+'/load-hits',data=b'['*2000+b']'*2000,headers={'Content-Type':'application/json'}),timeout=10)
        except HTTPError as error: assert error.code==400 and json.load(error)['error']
        else: raise AssertionError('Deep JSON accepted')
        assert api('/status')['busy'] is False
        for name,cut,error in [('short_cut.root','hit.trid==1','Cut field length mismatch'),('fractional_pid.root','1','integer particle IDs')]:
            rejected=load(name,cut);assert rejected['result'] is None and error in rejected['message'],rejected
        assert load('control_name.root')['result']['hits']==4
        assert 'bad\nname' in api('/hits')['fields']
        print('PASS adversarial array lengths, fractional PID rejection, JSON escaping and local request boundary',flush=True)
        expected = dict(collections.Counter(str(i) for i in [150252,150252,73003,70001,999,-1]*2))
        for name in ['valid.root', 'quoted " input.root']:
            result = load(name)
            assert result['result'] is not None, result
            assert all(result['result'][k] == v for k,v in {'events': 3, 'hits': 12, 'detectors': expected, 'stored_hits':12, 'selected_events':2}.items()), result
            assert 'hit.det' in result['result']['fields'], result
            print('PASS all detectors, repeats, gammas, low energy, backward hits, no rate:', name)
        for name, events in [('empty.root', 0), ('zero_hits.root', 3)]:
            result = load(name)
            assert result['result'] is not None, result
            assert all(result['result'][k] == v for k,v in {'events':events,'hits':0,'detectors':{},'stored_hits':0,'selected_events':0}.items()), result
            print('PASS', name)
        for name, error in [('missing_tree.root','Missing event tree'), ('missing_det.root','Missing hit.det'),
                            ('wrong_det.root','integer detector IDs'), ('short_field.root','Incomplete hit branch'),
                            ('corrupt.root','Cannot open ROOT file')]:
            result = load(name)
            assert result['result'] is None and error in result['message'], result
            print('PASS rejected:', name)
        for cut, hits, detectors in [
            ('hit.det == 73003', 2, {'73003':2}),
            ('hit.det == 150252 && hit.pid == 11', 2, {'150252':2}),
            ('hit.det == 73003 || hit.det == 150252', 6, {'73003':2, '150252':4}),
            ('hit.k >= 1', 0, {}),
            ('charged(hit.pid)', 6, {k:v//2 for k,v in expected.items()}),
            ('hit.mtrid == 0', 12, expected),
            ('hit.mtrid > 0', 0, {}),
            ('hit.pz < 0 && abs(hit.pid) == 11', 6, {k:v//2 for k,v in expected.items()}),
            ('sqrt(hit.x*hit.x) > 1 && Entry$ == 2', 6, {k:v//2 for k,v in expected.items()}),
            ('Entry$ == 0', 6, {k:v//2 for k,v in expected.items()}),
            ('Sum$(hit.pid == 22) > 0', 6, {k:v//2 for k,v in expected.items()}),
            ('hit.det[0] == 150252', 12, expected),
            ('2', 12, expected), # Nonzero means selected; never weights counts.
            ('0', 0, {}),
        ]:
            result = load('valid.root', cut)
            assert result['result'] is not None, result
            assert result['result']['hits'] == hits and result['result']['detectors'] == detectors, result
            assert result['result']['stored_hits'] == 12, result
            assert result['cut'] == cut, result
            print('PASS ROOT cut:', cut)
        assert load('scalar.root','rate > 20 && hit.det == 73003')['result']['hits'] == 1
        assert 'different particle array' in load('other_array.root','other > 0')['message']
        print('PASS scalar event field; unrelated array rejected')
        for cut, error in [('hit.no_such_field > 0','Invalid ROOT cut'),
                           ('hit.det ==','Invalid ROOT cut'),
                           ('TMath::Sqrt(-abs(hit.x))','Nonfinite cut value'),
                           ('hit.det[99] > 0','one value per hit')]:
            result=load('valid.root', cut)
            assert result['result'] is None and error in result['message'], result
            print('PASS rejected cut:', cut)
        result=load('no_positions.root')['result']
        assert result['hits']==2 and not result['positions_available']
        assert api('/hits')['cells']==[]
        assert not api('/hits')['histograms']['energy']['available']
        assert not api('/hits')['histograms']['radius']['available']
        print('PASS missing positions: counts available, plots unavailable')
        main=subprocess.check_output(['node','-e',"console.log(require('./geometry/view.js').mainCut)"],cwd=HERE,text=True).strip()
        for name,cut,hits in [('main',main,2),('primary',main+' && hit.mtrid == 0',1),
                              ('secondary',main+' && hit.mtrid > 0',1),('backward',main+' && hit.pz < 0',2),
                              ('showermax','hit.det >= 73001 && hit.det <= 73028 && hit.pz > 0',2)]:
            result=load('views.root',cut)
            assert result['result'] is not None and result['result']['hits']==hits,result
            assert sum(b[3] for b in api('/hits')['cells'])==hits
            print('PASS preset:',name,hits)
        assert load('bad_position.root')['result']['hits']==12
        assert api('/hits')['cells']==[] and api('/hits')['invalid_positions']==12
        print('PASS invalid coordinates counted, explicitly omitted from display')
        assert load('valid.root')['result']['hits'] == 12
        spatial=api('/hits')
        assert spatial['positions_available'] and spatial['invalid_positions']==0
        assert spatial['cells']==[[0,0,0,4],[1,0,0,2]], spatial
        assert spatial['other_detectors']==6 and spatial['unmatched_positions']==0
        assert 'cells' not in api('/status')['hits']['result']
        print('PASS exact surface cells, detector separation, other-detector counts and compact polling')
        h=spatial['histograms']
        assert sum(h['energy']['counts'])==sum(h['radius']['counts'])==sum(h['vertex']['counts'])==12
        assert h['radius']['counts'][0]==12 and h['vertex']['counts'][2]==12
        assert spatial['particles']=={'11':6,'22':6}
        assert load('valid.root','hit.det == 150252 && hit.pid == 11')['result']['hits']==2
        selected=api('/hits');assert sum(selected['histograms']['energy']['counts'])==sum(c[3] for c in selected['cells'])==2
        assert selected['particles']=={'11':2}
        load('hist_edges.root');energy=api('/hits')['histograms']['energy']
        assert sum(energy['counts'])==4 and energy['counts'][0]==2
        assert (energy['underflow'],energy['overflow'],energy['invalid'])==(4,2,2), energy
        print('PASS linked map/histogram cuts, particle counts, exact histogram edges and invalid values')
        load('plot_known.root')
        plotted=api('/hits');hs=plotted['histograms'];js=plotted['histograms2d']
        assert len(js['xy']['cells'])==6 and all(c[2]==2 for c in js['xy']['cells']),js['xy']
        assert sum(c[2] for c in js['xy']['cells'])==12
        # Six anti-correlated x/y pairs, not 36 combinations of marginal bins.
        for ix,iy,n in js['xy']['cells']:
            assert hs['x']['edges'][ix] in [100,200,300,400,500,600]
            assert hs['y']['edges'][iy]==100-2*hs['x']['edges'][ix]
        assert sum(hs['q2']['counts'])==12 and hs['q2']['event_value']
        assert sum(hs['w']['counts'])==6 and hs['w']['invalid']==6
        assert sum(c[2] for c in js['q2_w']['cells'])==6 and js['q2_w']['invalid']==6
        assert sum(hs['vertex_r']['counts'])==12 and hs['vertex_r']['counts'][0]==12
        assert 'histograms2d' not in api('/status')['hits']['result']
        load('plot_known.root','hit.det==73003')
        assert sum(c[2] for c in api('/hits')['histograms2d']['xy']['cells'])==2
        load('hist_edges.root');joint=api('/hits')['histograms2d']['energy_radius']
        assert sum(c[2] for c in joint['cells'])==4 and joint['outside']==6 and joint['invalid']==2
        load('bad_position.root');assert api('/hits')['histograms2d']['xy']['invalid']==12
        print('PASS 2D known joint pairs, scalar units/per-hit counting, invalid/overflow handling, linked cuts')
        custom={'x':{'bins':3,'min':0,'max':600,'spacing':'linear'},
                'y':{'bins':3,'min':-1200,'max':0,'spacing':'linear'},
                'energy':{'bins':3,'min':1,'max':1000,'spacing':'log'}}
        load('plot_known.root',binning=custom);data=api('/hits')
        assert data['histograms']['x']['edges']==[0,200,400,600]
        assert data['histograms']['x']['counts']==[2,4,4]
        assert data['histograms']['x']['overflow']==2
        assert all(math.isclose(a,b,rel_tol=1e-14) for a,b in zip(data['histograms']['energy']['edges'],[1,10,100,1000]))
        assert data['histograms']['energy']['counts']==[0,12,0]
        assert data['histograms2d']['xy']['cells']==[[0,2,2],[1,1,2],[1,2,2],[2,0,2],[2,1,2]]
        assert data['histograms2d']['xy']['outside']==2
        assert data['hits']==12 and api('/status')['hits']['binning']==custom
        load(['valid.root','plot_known.root'],binning=custom);data=api('/hits')
        assert sum(data['histograms']['x']['counts'])==22
        assert data['histograms']['x']['overflow']==2
        assert sum(c[2] for c in data['histograms2d']['xy']['cells'])+data['histograms2d']['xy']['outside']==24
        from urllib.error import HTTPError
        invalid=[{'bins':0,'min':0,'max':1,'spacing':'linear'},
                 {'bins':257,'min':0,'max':1,'spacing':'linear'},
                 {'bins':2.5,'min':0,'max':1,'spacing':'linear'},
                 {'bins':2,'min':1,'max':1,'spacing':'linear'},
                 {'bins':2,'min':0,'max':1,'spacing':'log'},
                 {'bins':2,'min':0,'max':float('nan'),'spacing':'linear'},
                 {'bins':2,'min':0,'max':1,'spacing':'unknown'},
                 {'bins':80,'min':-(10**400),'max':2000,'spacing':'linear'},
                 {'bins':80,'min':0,'max':10**400,'spacing':'linear'}]
        for bad in invalid:
            try:load('valid.root',binning={'x':bad})
            except HTTPError as error:assert error.code==400
            else:raise AssertionError('Invalid binning accepted')
        print('PASS custom ROOT bin edges, logarithmic bins, exact 2D pairs, overflow, multi-file sums and invalid settings')
        # Execute the browser's composed sub-cuts in ROOT, including a clicked tile.
        script="const {buildCut}=require('./geometry/view.js');console.log(JSON.stringify(['all','open','closed','transition'].map(region=>buildCut({mode:'raw',system:'main',ancestry:'all',direction:'all',particle:'all',tile:150252,region}))));"
        cuts=json.loads(subprocess.check_output(['node','-e',script],cwd=HERE,text=True))
        for cut,expected_hits in zip(cuts,[4,4,0,0]):
            selected=load('valid.root',cut)['result'];assert selected is not None and selected['hits']==expected_hits,(cut,selected)
            assert sum(api('/hits')['histograms']['radius']['counts'])==expected_hits
        # The fixture deliberately puts open tile 150252 at phi=0: ID must win.
        script="const {buildCut}=require('./geometry/view.js');console.log(JSON.stringify(['all','open','closed','transition'].map(region=>buildCut({mode:'raw',system:'showermax',ancestry:'all',direction:'all',particle:'all',tile:null,region}))));"
        smcuts=json.loads(subprocess.check_output(['node','-e',script],cwd=HERE,text=True))
        for cut,n in zip(smcuts,[2,2,0,0]):
            assert 'atan2' not in cut and load('valid.root',cut)['result']['hits']==n
        assert load('no_positions.root',smcuts[1])['result']['hits']==1
        print('PASS physical ID region cuts in ROOT: both detectors at contrary phi; no positions required')
        # Multiple files are independent ROOT scans, then additive counts.
        combined=load(['valid.root','quoted " input.root'])['result']
        assert combined['hits']==24 and combined['events']==6 and combined['file_count']==2
        assert combined['detectors']=={k:2*v for k,v in expected.items()}
        combined=api('/hits')
        assert sum(combined['histograms']['energy']['counts'])==24
        assert sum(c[2] for c in combined['histograms2d']['xy']['cells'])==24
        assert sum(c[3] for c in combined['cells'])==12
        assert load(['valid.root','quoted " input.root'],'hit.pid==11')['result']['hits']==12
        mixed=load(['valid.root','no_positions.root'])['result']
        assert mixed['hits']==14 and not mixed['positions_available']
        assert api('/hits')['cells']==[] and not api('/hits')['histograms']['energy']['available']
        assert not api('/hits')['histograms2d']['xy']['available']
        assert load(['valid.root','corrupt.root'])['result'] is None
        assert load(['valid.root','no_positions.root'],'hit.pid==11')['result'] is None
        from urllib.error import HTTPError
        for names in [['valid.root','valid.root'],['missing.root']]:
            try: load(names)
            except HTTPError as error: assert error.code==400
            else: raise AssertionError('Duplicate/missing file accepted')
        assert load('valid.root')['result']['hits']==12
        print('PASS multiple-file sums, shared cuts, missing fields, no partial results and duplicate rejection')
        assert not list(scratch.glob('hit-counts*.json'))
        print('PASS reload after errors; temporary receipt removed')
        # Geometry is independent: no ROOT load should claim to replace it.
        app.state['ready']=False
        assert app.state['ready'] is False and app.state['revision'] == 0
    finally:
        server.shutdown()
        server.server_close()
        thread.join()
        app.stop()
