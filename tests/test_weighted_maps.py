"""Known PE rates through real ROOT files, compiled metadata and multi-file pooling."""
import csv, importlib.util, json, math, shutil, subprocess, sys, tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'app'))
import responses
spec=importlib.util.spec_from_file_location('geometry_run',HERE/'geometry/run.py')
run=importlib.util.module_from_spec(spec);spec.loader.exec_module(run)

def main():
 with tempfile.TemporaryDirectory(dir=run.OUTPUT.parent,prefix='weighted-test-') as tmp:
  base=Path(tmp);run.OUTPUT=base/'output';run.OUTPUT.mkdir()
  code=base/'remoll/analysis';shutil.copytree(HERE/'vendor/main_response',code)
  tables=base/'tables';shutil.copytree(HERE/'data/main_pe',tables)
  sm=code/'showermax/response_lookup/data';shutil.copytree(HERE/'data/response',sm)
  for path in tables.glob('*.csv'):
   rows=list(csv.reader(path.open()))
   for row in rows[1:]:
    if 'NewLookup' in path.name: row[8:13]=['10','0','0','0','0']
    else: row[1]='1'
   with path.open('w',newline='') as f:csv.writer(f).writerows(rows)
  for path in sm.glob('*.csv'):
   rows=list(csv.reader(path.open()))
   for row in rows[1:]:row[1:]=['100','0','100','0','0']
   with path.open('w',newline='') as f:csv.writer(f).writerows(rows)
  selected=responses.resolve(dict(mode='local',remoll=str(base/'remoll'),main_tables=str(tables)))
  maker=base/'Make.C'
  maker.write_text('#include '+json.dumps(str(HERE/'Deconvolution.C'))+'\n'+r'''
void Make(){
 for(int seed=1;seed<=10;++seed){
  std::string path=std::string(gSystem->Getenv("FIXTURE_DIR"))+"/"+std::to_string(seed)+".root";
  TFile f(path.c_str(),"RECREATE");TTree t("T","known weighted hits");
  int n=0,det[4],pid[4];double rate,x[4],y[4],z[4],xl[4],yl[4],zl[4],px[4],py[4],pz[4],k[4],e[4];
  t.Branch("n",&n,"n/I");t.Branch("rate",&rate);
  t.Branch("hit.det",det,"hit.det[n]/I");t.Branch("hit.pid",pid,"hit.pid[n]/I");
  #define FIELD(a) t.Branch("hit." #a,a,"hit." #a "[n]/D")
  FIELD(x);FIELD(y);FIELD(z);FIELD(xl);FIELD(yl);FIELD(zl);FIELD(px);FIELD(py);FIELD(pz);FIELD(k);FIELD(e);
  for(int j=0;j<4;j++){pid[j]=11;x[j]=1100;y[j]=0;z[j]=28000;xl[j]=yl[j]=zl[j]=px[j]=py[j]=0;pz[j]=100;k[j]=100;e[j]=1000;}
  n=3;rate=2;det[0]=det[1]=150252;det[2]=73001;if(seed==7)det[0]=det[1]=150010;t.Fill();
  n=1;rate=5;t.Fill();n=0;rate=0;t.Fill();t.Write();
  remollRunData d;d.fNthrown=3;d.fSeed=seed;d.fGitInfo="known-rate-test";
  std::string macro="/remoll/evgen/set moller\n/remoll/beamcurr 65 uA\n/run/beamOn 3\n";
  if(seed>=8)macro="/remoll/evgen/set elastic\n/remoll/beamcurr 65 uA\n"+std::string(seed==9?"/remoll/rasx 2 mm\n":"")+"/run/beamOn 3\n";
  macro.insert(macro.find("/run/beamOn"),std::string("/remoll/SD/enable_all\n")+(seed==10?"/remoll/SD/detect all 150252\n/remoll/SD/detect all 73001\n":"/remoll/SD/detect_range boundaryhits 110000 152000\n/remoll/SD/detect_range surfacehits 73001 73028\n"));
  d.fMacro.fBufferSize=macro.size()+1;d.fMacro.fBuffer=new char[d.fMacro.fBufferSize];std::strcpy(d.fMacro.fBuffer,macro.c_str());
  if(seed==3)d.fGitInfo="";
  if(seed!=4){
   remollTextFile g;std::string name=gSystem->Getenv("FIXTURE_GEOMETRY_NAME"),contents=root_deconv::text(gSystem->Getenv("FIXTURE_GEOMETRY"));
   if(seed==5)contents="<gdml/>";
   g.fFilenameSize=name.size()+1;g.fFilename=new char[g.fFilenameSize];std::strcpy(g.fFilename,name.c_str());
   g.fBufferSize=contents.size()+1;g.fBuffer=new char[g.fBufferSize];std::strcpy(g.fBuffer,contents.c_str());d.fGDMLFiles.push_back(g);
  }
  d.Write("run_data");
 }
}
''')
  import os
  remoll=Path(os.environ['REMOLL_DIR']).expanduser().resolve() if os.environ.get('REMOLL_DIR') else next((p for p in HERE.parents if (p/'geometry/mollerMother.gdml').is_file()),None)
  if remoll is None or not (remoll/'geometry/mollerMother.gdml').is_file():raise RuntimeError('Set REMOLL_DIR to your remoll source directory.')
  geometry_name=(HERE/'data/geometry.md5').read_text().splitlines()[0].split()[1]
  p=subprocess.run(['root' ,'-n','-l','-b','-q',str(maker)+'+'],env={**os.environ,'FIXTURE_DIR':str(base),'FIXTURE_GEOMETRY_NAME':geometry_name,'FIXTURE_GEOMETRY':str(remoll/'geometry'/geometry_name)},capture_output=True,text=True)
  assert p.returncode==0,p.stdout+p.stderr
  for files in [[base/'1.root'],[base/'1.root',base/'2.root']]:
   run.load_hit_file(files,'charged(hit.pid) && hit.k >= 1',response=selected)
   data=run.state['hits']['result'];assert data is not None,run.state['hits']['message']
   assert math.isclose(data['signal_total'],156,rel_tol=1e-12),data['signal_total']
   assert math.isclose(data['detector_signal']['150252'],90,rel_tol=1e-12)
   assert math.isclose(data['detector_signal']['73001'],66,rel_tol=1e-12)
   print('PASS',len(files),'file(s): main 90 + ShowerMax 66 = 156 PE/s')
  run.load_hit_file([base/'1.root',base/'1.root'],'1',response=selected)
  assert run.state['hits']['result'] is None
  assert 'Repeated seed' in run.state['hits']['message'],run.state['hits']['message']
  print('PASS duplicate seed rejected; no partial output')
  run.load_hit_file([base/'1.root',base/'8.root'],'1',response=selected)
  assert math.isclose(run.state['hits']['result']['signal_total'],312,rel_tol=1e-12),run.state['hits']['message']
  run.load_hit_file([base/'1.root',base/'9.root'],'1',response=selected)
  assert run.state['hits']['result'] is None
  assert 'Incompatible transport' in run.state['hits']['message'],run.state['hits']['message']
  print('PASS different interactions sum; mismatched transport rejected across interactions')
  run.load_hit_file([base/'10.root'],'1',response=selected)
  assert run.state['hits']['result'] is None
  assert 'entry recording' in run.state['hits']['message'],run.state['hits']['message']
  run.load_hit_file([base/'10.root'],'1')
  assert run.state['hits']['result']['hits']==4
  print('PASS all-step PE weighting rejected; raw hit counts remain available')
  (base/'6.neyman.tsv').write_text('independent malicious fixture: sampled history normalization')
  for seed,message in [(3,'source identity'),(4,'Saved geometry'),(5,'does not match'),(6,'stratified normalization')]:
   run.load_hit_file([base/f'{seed}.root'],'1',response=selected)
   assert run.state['hits']['result'] is None, seed
   assert message in run.state['hits']['message'],run.state['hits']['message']
   run.load_hit_file([base/f'{seed}.root'],'1')
   assert run.state['hits']['result']['hits']==4
  run.load_hit_file([base/'7.root'],'1',response=selected)
  assert run.state['hits']['result']['unsupported_pe']==3
  assert math.isclose(run.state['hits']['result']['signal_total'],66,rel_tol=1e-12)
  run.load_hit_file([base/'1.root'],'1',response=selected)
  assert math.isclose(run.state['hits']['result']['signal_total'],156,rel_tol=1e-12)
  print('PASS missing source/geometry, mismatched PE geometry, sampling sidecar, unknown detector; hit counts and recovery remain available')
if __name__=='__main__':main()
