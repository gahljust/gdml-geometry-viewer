// Stock remoll split ROOT branches only; no remoll headers or custom track fields.
#include "detector_selection.hh"
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TDatabasePDG.h>
#include <vector>
#include <map>
#include <array>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdexcept>
namespace portable_deconv {
struct Moments {long long entries=0,unsupported_sm=0; std::vector<double> y,cross,noise,coverage;};
using Coeff=std::array<double,5>;
using Fits=std::map<int,std::map<double,Coeff>>;
template<size_t N> std::vector<std::array<double,N>> numeric_table(const std::string& path){
 std::ifstream in(path);if(!in)throw std::runtime_error("Missing PE table: "+path);
 std::string line;std::getline(in,line);
 if(std::count(line.begin(),line.end(),',')!=N-1)throw std::runtime_error("Invalid PE header: "+path);
 std::vector<std::array<double,N>> result;int row=1;
 while(std::getline(in,line)){
  ++row;if(line.find_first_not_of(" \t\r")==std::string::npos)continue;
  auto invalid=[&](){return std::runtime_error("Invalid PE table row "+std::to_string(row)+": "+path);};
  if(std::count(line.begin(),line.end(),',')!=N-1)throw invalid();
  std::istringstream columns(line);std::array<double,N> values;std::string cell;
  for(double& value:values){
   if(!std::getline(columns,cell,','))throw invalid();
   auto first=cell.find_first_not_of(" \t\r"),last=cell.find_last_not_of(" \t\r");
   if(first==std::string::npos)throw invalid();cell=cell.substr(first,last-first+1);
   if(cell.size()>1&&cell.front()=='"'&&cell.back()=='"')cell=cell.substr(1,cell.size()-2);
   size_t used=0;try{value=std::stod(cell,&used);}catch(...){throw invalid();}
   if(cell.find_first_not_of(" \t\r",used)!=std::string::npos||!std::isfinite(value))throw invalid();
  }
  result.push_back(values);
 }
 if(result.empty())throw std::runtime_error("Empty PE table: "+path);
 return result;
}
Fits load(const std::string& folder){
 Fits fits;const int pid[]={11,-11,22,13,-13,-211,211,2112};
 const char* names[]={"e-","e+","gamma","mu-","mu+","pi-","pi+","neutron"};
 for(int i=0;i<8;i++){
  const std::string path=folder+"/fit_param_xy_"+names[i]+"_ifarm.csv";
  for(const auto& values:numeric_table<6>(path)){
   const double e=values[0];Coeff c;std::copy(values.begin()+1,values.end(),c.begin());
   if(e<=0||(fits[pid[i]].count(e)&&fits[pid[i]].at(e)!=c))throw std::runtime_error("Invalid energy or conflicting PE coefficients: "+path);
   fits[pid[i]][e]=c; // Identical duplicate rows occur in the upstream tables.
  }
 }
 return fits;
}
bool charged(int pid){long long n=std::abs((long long)pid);if(n>=1000000000LL)return (n/10000)%1000>0;if(!pid)return false;auto* p=TDatabasePDG::Instance()->GetParticle(pid);if(!p)throw std::runtime_error("Unknown PDG charge: "+std::to_string(pid));return p->Charge()!=0;}
double pe(const Fits& fits,int det,int pid,double e,double x,double y){auto f=fits.find(pid);if(f==fits.end()||e<f->second.begin()->first||e>f->second.rbegin()->first)return -1;const auto& tab=f->second;auto hi=tab.lower_bound(e),lo=hi;if(hi!=tab.begin()&&hi->first!=e)--lo;auto local=detector_selection::showermax_xy(det,x,y);double lx=local.first,ly=local.second;auto eval=[&](const Coeff& c){return .5*(c[0]+lx*c[1])+ .5*(c[2]+ly*c[3]+ly*ly*c[4]);};double v=eval(lo->second);if(hi!=lo)v+=(eval(hi->second)-v)*(e-lo->first)/(hi->first-lo->first);return std::max(0.,v*((det-73000)%4==3?.22:.33));}
Moments scan(const std::string& path,const std::vector<int>& ids,int channel,const std::string& response,double kmin){
 TFile file(path.c_str(),"READ");auto* tree=dynamic_cast<TTree*>(file.Get("T"));if(!tree)throw std::runtime_error("Missing T tree");
 const char* required[]={"rate","ev.A","ev.W2","hit.det","hit.pid","hit.k","hit.e","hit.x","hit.y","hit.pz"};for(auto s:required)if(!tree->GetBranch(s)&&!tree->GetLeaf(s))throw std::runtime_error(std::string("Missing stock branch ")+s);
 TTreeReader reader(tree);TTreeReaderValue<double> rate(reader,"rate"),asym(reader,"ev.A"),w2(reader,"ev.W2");TTreeReaderArray<int> det(reader,"hit.det"),pid(reader,"hit.pid");TTreeReaderArray<double> k(reader,"hit.k"),e(reader,"hit.e"),x(reader,"hit.x"),y(reader,"hit.y"),pz(reader,"hit.pz");
 const int n=ids.size(),d=n+1,nc=channel==2?3:1,dim=2*d;Moments out;out.y.resize(nc*dim);out.cross.resize(nc*dim*dim);out.noise.resize(n*n);out.coverage.resize(n);std::map<int,int> index;bool sm=false;for(int i=0;i<n;i++){index[ids[i]]=i;sm|=ids[i]>=73001&&ids[i]<=73028;}Fits fits;if(sm)fits=load(response);
 while(reader.Next()){out.entries++;if(!std::isfinite(*rate)||*rate<0||!std::isfinite(*asym))throw std::runtime_error("Nonfinite/negative rate or asymmetry");if(*rate==0)continue;const auto size=det.GetSize();if(pid.GetSize()!=size||k.GetSize()!=size||e.GetSize()!=size||x.GetSize()!=size||y.GetSize()!=size||pz.GetSize()!=size)throw std::runtime_error("Inconsistent hit arrays");
 std::vector<double> q(d,0);for(size_t j=0;j<size;j++){auto it=index.find(det[j]);if(it==index.end())continue;int i=it->second;out.coverage[i]++;double v=0;if(!std::isfinite(k[j])||!std::isfinite(e[j])||!std::isfinite(x[j])||!std::isfinite(y[j])||!std::isfinite(pz[j]))throw std::runtime_error("Nonfinite hit");if(det[j]>=73001&&det[j]<=73028){if(pz[j]>0){v=pe(fits,det[j],pid[j],e[j],x[j],y[j]);if(v<0){out.unsupported_sm++;v=0;}}}else if(charged(pid[j])&&k[j]>=kmin){v=1;const int suffix=det[j]%100;double phi=fmod(atan2(y[j],x[j])+2*M_PI,2*M_PI/7);if(suffix>=51&&suffix<=53&&phi>=3*M_PI/28&&phi<5*M_PI/28)q[n]+=1;}q[i]+=v;}
 int c=0;if(channel==2){double w=*w2>0?sqrt(*w2)*.001:-1;c=w>=1&&w<1.4?0:w>=1.4&&w<2.5?1:w>=2.5&&w<6?2:-1;}if(c<0){if(std::any_of(q.begin(),q.end(),[](double v){return v>0;}))throw std::runtime_error("Signal outside configured inelastic W bins [1,6) GeV");continue;}
 std::vector<int> nz;for(int i=0;i<d;i++)if(q[i])nz.push_back(i);for(int i:nz){double yi=*rate*q[i],zi=yi*(*asym);out.y[c*dim+i]+=yi;out.y[c*dim+d+i]+=zi;for(int j:nz){double yy=yi*(*rate)*q[j];out.cross[c*dim*dim+i*dim+j]+=yy;out.cross[c*dim*dim+i*dim+d+j]+=yy*(*asym);out.cross[c*dim*dim+(d+i)*dim+j]+=yy*(*asym);out.cross[c*dim*dim+(d+i)*dim+d+j]+=yy*(*asym)*(*asym);if(i<n&&j<n)out.noise[i*n+j]+=*rate*q[i]*q[j];}}}
 if(out.entries!=tree->GetEntries())throw std::runtime_error("Incomplete ROOT read");return out;
}
}
