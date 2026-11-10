// Correlated deconvolution from raw remoll ROOT files.
// Edit fnms below: select Moller, ep-elastic and ep-inelastic files.
// Individual filenames and wildcards are supported.
// Inputs need the required detector entries and saved run_data metadata.
//
// Default: main-detector crossings + ShowerMax expected PE.
//   root -l -b -q 'deconvolution_correlated.C+()'
// For PE in both detectors, use 'deconvolution_correlated.C+(true)'.
//
// Arguments, in order:
// 1. usePE           Main detector: false = crossings, true = expected PE.
// 2. tableDirectory  PE CSV folder; "" searches automatically.
// 3. output          Output ROOT filename; existing files are not overwritten.
// 4. beamDays        Beam exposure in days (344).
// 5. polarization    Fraction, not percent (0.8).
// 6. current_uA      Beam current in microamps (65).
// 7. minimumK        Crossing kinetic-energy cut in MeV (1).
// 8. useShowerMaxPE  ShowerMax: true = expected PE, false = charged crossings.
//
// Default mode needs ShowerMax tables. Main PE also needs its tables,
// minimumK=1 and ShowerMax PE enabled. PE models apply additional selections.
// Set both PE switches to false for crossing rates without tables.
// Turning PE off does not remove detectors. See README.md for examples.
//
// Sum each event's signal per detector; accumulate detector-pair products
// with one event rate weight to calculate covariance. Finite-MC error is
// calculated separately. Additional optical/electronics/beam noise is excluded.
//
// Adapted directly from gahljust/gdml-geometry-viewer, commit 5077001171756458963314261943a0dddd9bfa40:
// root/Deconvolution.C and root/root_solver.hh. The numerical solver and PE
// evaluators below are copied without changing their equations. Adaptations:
// one-file packaging, editable input patterns, PE discovery, and the rate-only option.
// Real-cohort support follows experiments/tile_convergence/src/scan.cc and
// convergence_stats.py: physical Ring 5 open reference and per-plan-stratum MC.
// develop:analysis/deconvolution.C at ae0c5fabc802cb2c2139c458d9d849677e6d98d3
// supplies the editable-file-list style, not the correlated estimator.

#include <TChain.h>
#include <TChainElement.h>
#include <TDatabasePDG.h>
#include <TDecompChol.h>
#include <TDecompSVD.h>
#include <TFile.h>
#include <TH2D.h>
#include <TMD5.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TObject.h>
#include <TROOT.h>
#include <TRotation.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <unordered_map>
#include <vector>

// ------------------------- INPUT FILES / OPTIONS -------------------------
std::vector<std::string> fnms = {
    "moller*.root",
    "ep_elastic*.root",
    "ep_inelastic*.root"
};
// Set true only for an explicitly selected mixed-source cohort, such as the
// accepted 9M set. Generator domains are still checked within each interaction.
bool allowMixedSources = false;
// When enabled, every source/geometry/physical-command setup must match an
// independently reviewed cohort manifest. The normal compatibility checks
// remain the default. The transport MD5 is printed in a failed preflight.
std::set<std::string> trustedMixedTransportMD5;
// Fixed-quota inputs only: independently verified generation-plan MD5 values.
// Do not approve a hash merely because it was found beside a ROOT file.
std::set<std::string> trustedSamplingPlanMD5;

// PE search: explicit folder, otherwise MOLLER_PE_TABLES, otherwise common
// data/tables folders near this macro, the working directory and input files.
// A supplied folder/environment variable is authoritative: missing tables fail.
// A data folder may contain main_pe/ and response/ subdirectories.
// Default needs fit_param_xy_{e-,e+,gamma,mu-,mu+,pi-,pi+,neutron}_ifarm.csv.
// usePE=true additionally needs r1..r6_NewLookupTable.csv and r1..r6_EnergyDepTable.csv.
void deconvolution_correlated(bool usePE=false, const char* tableDirectory="",
                   const char* output="deconvolution_correlated.root", double beamDays=344,
                   double polarization=.8, double current_uA=65, double minimumK=1,
                   bool useShowerMaxPE=true);
// ------------------------------------------------------------------------

// Source: root/metadata_schema.hh | SHA256 8df2e8e577a05c162a3f41239c100b4346fe6df11290c020258e46556aee5fc1
// Read-only adaptation of remollRunData v2 and remollTextFile v1. No remoll
// filesystem methods are needed. Fixed-width run counters avoid Windows long.
struct filedata_t {
    char filename[255];
    char hashsum[255];
    TTimeStamp timestamp;
};
class remollTextFile:public TObject {
public:
    unsigned int fFilenameSize=0;
    char* fFilename=nullptr;
    unsigned long long int fBufferSize=0;
    char* fBuffer=nullptr;
    remollTextFile()=default;
    remollTextFile(const remollTextFile& other):TObject(other) {assign(other);}
    remollTextFile& operator=(const remollTextFile& other) {
        if(this!=&other){TObject::operator=(other);delete[] fFilename;delete[] fBuffer;assign(other);}return *this;
    }
    ~remollTextFile() override {delete[] fFilename;delete[] fBuffer;}
private:
    void assign(const remollTextFile& other) {
        fFilenameSize=other.fFilenameSize;fBufferSize=other.fBufferSize;
        fFilename=other.fFilename?new char[fFilenameSize]:nullptr;
        fBuffer=other.fBuffer?new char[fBufferSize]:nullptr;
        if(fFilename)std::memcpy(fFilename,other.fFilename,fFilenameSize);
        if(fBuffer)std::memcpy(fBuffer,other.fBuffer,fBufferSize);
    }
    ClassDefOverride(remollTextFile,1);
};
class remollRunData:public TObject {
public:
    std::vector<filedata_t> fMagData;
    remollTextFile fMacro;
    std::vector<remollTextFile> fGDMLFiles;
    Long64_t fNthrown=0;
    Long64_t fSeed=0;
    std::string fGitInfo;
    TTimeStamp fRunTime;
    std::string fRunPath;
    std::string fHostName;
    ClassDefOverride(remollRunData,3);
};


// Source: root/LinkDef.h | SHA256 7d58a21902948b52618a3ea74c0ac3995b6f850e11046d7d5b5b442765214a16
#ifdef __CLING__
#pragma link C++ class filedata_t+;
#pragma link C++ class std::vector<filedata_t>+;
#pragma link C++ class remollTextFile+;
#pragma link C++ class std::vector<remollTextFile>+;
#pragma link C++ class remollRunData+;
#endif


// Source: vendor/portable/detector_selection.hh | SHA256 0c4bb572add57655ddd7b6590b12acc8da4c202d91ef6d48652e1795182d534a
#ifndef REMOLL_ANALYSIS_DETECTOR_SELECTION_HH
#define REMOLL_ANALYSIS_DETECTOR_SELECTION_HH

// Current ThinQuartz and 28 ShowerMax entrance-plane numbering.
// Region order matches campaign products: closed, transition, open.
namespace standalone_detector_selection {
inline bool main_quartz(int det) {
  const int suffix=det%100, sector=(det%10000)/100;
  const bool tile=suffix==10||suffix==20||suffix==30||suffix==40||
                  suffix==51||suffix==52||suffix==53||suffix==60;
  return tile && sector>=1 && sector<=14 &&
         ((det>=110000&&det<112000)||(det>=150000&&det<152000));
}
inline int region(int det) {
  if (det>=73001 && det<=73028) return det%4==3 ? 2 : det%4==1 ? 0 : 1;
  if (!main_quartz(det)) throw std::invalid_argument("Unknown detector region: "+std::to_string(det));
  if (det<112000) return 1;
  return ((det%10000)/100)%2==0 ? 2 : 0;
}
inline bool ring5_open(int det) {
  return det>=150000 && det<152000 && det%100>=51 && det%100<=53 && region(det)==2;
}
inline std::pair<double,double> showermax_xy(int det,double x,double y) {
  if(det<73001||det>73028)throw std::invalid_argument("Not a ShowerMax entrance plane");
  // Preserve the lookup's nominal Qsim axes: x radial, y tangential, mm;
  // radius 1100 mm and 28 equally spaced modules. Only the frame identity changes.
  const double a=(det-73001)*(2*std::acos(-1.0)/28);
  const double dx=x-1100*std::cos(a),dy=y-1100*std::sin(a);
  return {dx*std::cos(a)+dy*std::sin(a),-dx*std::sin(a)+dy*std::cos(a)};
}
}
#endif


// Source: vendor/portable/scan.C | SHA256 99e49d7a03474ee3ad27a9df60a6f3dcf4ff0af10ade71a926b991a1eb1aa837
// Stock remoll split ROOT branches only; no remoll headers or custom track fields.
namespace standalone_portable {
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
double pe(const Fits& fits,int det,int pid,double e,double x,double y){auto f=fits.find(pid);if(f==fits.end()||e<f->second.begin()->first||e>f->second.rbegin()->first)return -1;const auto& tab=f->second;auto hi=tab.lower_bound(e),lo=hi;if(hi!=tab.begin()&&hi->first!=e)--lo;auto local=standalone_detector_selection::showermax_xy(det,x,y);double lx=local.first,ly=local.second;auto eval=[&](const Coeff& c){return .5*(c[0]+lx*c[1])+ .5*(c[2]+ly*c[3]+ly*ly*c[4]);};double v=eval(lo->second);if(hi!=lo)v+=(eval(hi->second)-v)*(e-lo->first)/(hi->first-lo->first);return std::max(0.,v*((det-73000)%4==3?.22:.33));}
}


// Source: root/root_solver.hh | SHA256 346d5f34c2ff8de5717641e19e28e718ea4a84fe5603bfc36e1bd0059cc90f74
// ROOT-only equivalent of the launcher's validated GLS solver.
// Inputs: Y=[rate, rate*asymmetry] per process, joint finite-MC covariance V,
// and the physical second moment B (one rate weight, not rate squared).
namespace standalone_deconv {
// Center scalar event scores without subtracting two almost equal raw moments.
// Only nonzero scores need updates; skipped histories are merged as a zero group.
struct CenteredScore {
    long long count=0;
    long double origin=0,mean=0,m2=0;
    long double average()const{return origin+mean;}
    void add(double value){
        if(!std::isfinite(value))throw std::runtime_error("Nonfinite MC score");
        if(value==0)return;
        if(!count)origin=value;
        ++count;long double delta=(static_cast<long double>(value)-origin)-mean;
        m2+=delta*delta*(count-1)/count;mean+=delta/count;
    }
    double variance_of_sum(long long histories)const{
        if(histories<2||count>histories)throw std::runtime_error("Invalid MC history count");
        const long double centered=m2+average()*average()*count*(histories-count)/histories;
        return static_cast<double>(centered*histories/(histories-1));
    }
};
// Center the histories where both scores are nonzero, then merge the one-score
// and zero-score groups. This keeps sparse scans sparse, including implicit zeros.
struct CenteredPair {
    long long count=0;
    long double origin_x=0,origin_y=0,mean_x=0,mean_y=0,cross=0;
    void add(double x,double y){
        if(x==0||y==0)return;
        if(!count){origin_x=x;origin_y=y;}
        ++count;long double dx=(static_cast<long double>(x)-origin_x)-mean_x;
        long double dy=(static_cast<long double>(y)-origin_y)-mean_y;
        cross+=dx*dy*(count-1)/count;mean_x+=dx/count;mean_y+=dy/count;
    }
    double covariance_of_sum(const CenteredScore& x,const CenteredScore& y,long long histories)const{
        if(histories<2||count>x.count||count>y.count||x.count+y.count-count>histories)
            throw std::runtime_error("Invalid joint MC history count");
        long long used=count;
        const long double joint_x=origin_x+mean_x,joint_y=origin_y+mean_y;
        long double mx=joint_x,my=joint_y,centered=cross;
        auto merge=[&](long long n,long double ax,long double ay){
            if(!n)return;
            long double dx=ax-mx,dy=ay-my,fraction=static_cast<long double>(n)/(used+n);
            centered+=dx*dy*used*fraction;mx+=dx*fraction;my+=dy*fraction;used+=n;
        };
        long long only_x=x.count-count,only_y=y.count-count;
        if(only_x)merge(only_x,x.average()+count*(x.average()-joint_x)/only_x,0);
        if(only_y)merge(only_y,0,y.average()+count*(y.average()-joint_y)/only_y);
        merge(histories-x.count-y.count+count,0,0);
        return static_cast<double>(centered*histories/(histories-1));
    }
};
struct Fit {
    TMatrixD rates,F,F_mc_sigma,H,B,C,correlation,M;
    TMatrixD experimental_covariance,mc_covariance,combined_covariance;
    TVectorD reference_asymmetry_ppb,extracted_asymmetry_ppb;
    double scaled_condition=0;
};
inline TMatrixD transpose(const TMatrixD& m){return TMatrixD(TMatrixD::kTransposed,m);}
inline Fit solve(const TMatrixD& y,const TMatrixD& v,const TMatrixD& b,
                 double days=344,double polarization=.8){
    const int n=b.GetNrows(),d=n+1,dim=2*d,components=y.GetNrows();
    if(components<1||n<components||y.GetNcols()!=dim||v.GetNrows()!=components*dim||v.GetNcols()!=components*dim||b.GetNcols()!=n)
        throw std::runtime_error("Incompatible moment dimensions");
    if(!(std::isfinite(days)&&days>0&&std::isfinite(polarization)&&polarization>0&&polarization<=1))throw std::runtime_error("Invalid exposure or polarization");
    for(int c=0;c<components;c++)for(int i=0;i<dim;i++){
        if(!std::isfinite(y(c,i)))throw std::runtime_error("Nonfinite rate/asymmetry moment");
        if(i<d&&y(c,i)<0)throw std::runtime_error("Negative detector rate");
        if(i<d&&y(c,i)==0&&y(c,d+i)!=0)throw std::runtime_error("Asymmetry sum without detector rate");
    }
    auto validate_symmetric=[](const TMatrixD& m,const std::string& name,bool nonnegative){
        for(int i=0;i<m.GetNrows();i++)for(int j=0;j<m.GetNcols();j++){
            if(!std::isfinite(m(i,j)))throw std::runtime_error("Nonfinite "+name);
            if((nonnegative||i==j)&&m(i,j)<0)throw std::runtime_error("Negative "+name);
            if(std::abs(m(i,j)-m(j,i))>1e-10*std::max({1.,std::abs(m(i,j)),std::abs(m(j,i))}))
                throw std::runtime_error("Asymmetric "+name);
        }
    };
    validate_symmetric(b,"physical second moment",true);
    validate_symmetric(v,"MC covariance",false);
    // Allow singular/zero MC covariance, but reject negative-variance directions.
    // Normalize heterogeneous rate/asymmetry units; jitter is for validation only.
    std::vector<int> active;
    for(int i=0;i<v.GetNrows();i++){
        if(v(i,i)>0)active.push_back(i);
        else for(int j=0;j<v.GetNcols();j++)if(v(i,j)!=0)
            throw std::runtime_error("MC covariance with zero variance and nonzero covariance");
    }
    if(!active.empty()){
        int size=active.size();TMatrixD normalized(size,size);
        double tolerance=128*std::numeric_limits<double>::epsilon()*size;
        for(int i=0;i<size;i++)for(int j=0;j<size;j++){
            int a=active[i],b=active[j];double scale=std::sqrt(v(a,a))*std::sqrt(v(b,b));
            double x=v(a,b)/scale,z=v(b,a)/scale;
            if(!std::isfinite(x)||std::abs(x-z)>1e-10*std::max({1.,std::abs(x),std::abs(z)}))
                throw std::runtime_error("Asymmetric/nonfinite normalized MC covariance");
            normalized(i,j)=.5*(x+z)+(i==j?tolerance:0);
        }
        TDecompChol check(normalized);
        if(!check.Decompose())throw std::runtime_error("MC covariance is not positive semidefinite");
    }
    Fit f;f.rates.ResizeTo(n,components);f.F.ResizeTo(n,components);f.H.ResizeTo(n,components);f.B.ResizeTo(b);f.B=b;
    f.C.ResizeTo(n,n);f.correlation.ResizeTo(n,n);f.reference_asymmetry_ppb.ResizeTo(components);
    f.extracted_asymmetry_ppb.ResizeTo(components);
    TVectorD r(n),sd(n),observed(n);
    for(int c=0;c<components;c++){
        if(!(y(c,n)>0))throw std::runtime_error("Missing Ring 5 open reference for component "+std::to_string(c));
        f.reference_asymmetry_ppb[c]=y(c,d+n)/y(c,n);
        if(!std::isfinite(f.reference_asymmetry_ppb[c])||f.reference_asymmetry_ppb[c]==0)throw std::runtime_error("Invalid reference asymmetry");
    }
    for(int i=0;i<n;i++){
        for(int c=0;c<components;c++){r[i]+=y(c,i);observed[i]+=y(c,d+i);f.rates(i,c)=y(c,i);}
        if(!(r[i]>0))throw std::runtime_error("Empty detector row "+std::to_string(i)+"; collect more events or enable missing detectors");
        observed[i]/=r[i];
        for(int c=0;c<components;c++){f.F(i,c)=y(c,i)/r[i];f.H(i,c)=y(c,d+i)/(r[i]*f.reference_asymmetry_ppb[c]);}
    }
    // Forecast covariance in ppb^2. It is distinct from finite-MC covariance V.
    for(int i=0;i<n;i++)for(int j=0;j<n;j++){
        f.C(i,j)=b(i,j)/(polarization*polarization*days*86400*r[i]*r[j])*1e18;
        if(!std::isfinite(f.C(i,j)))throw std::runtime_error("Nonfinite projected covariance");
    }
    for(int i=0;i<n;i++){
        sd[i]=std::sqrt(f.C(i,i));if(!(sd[i]>0))throw std::runtime_error("Zero detector variance");
    }
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)f.correlation(i,j)=f.C(i,j)/(sd[i]*sd[j]);
    TDecompChol chol(f.correlation);
    if(!chol.Decompose())throw std::runtime_error("Singular detector covariance; insufficient statistics or redundant channels");
    const TMatrixD& u=chol.GetU(); // rho = U^T U; L = U^T.
    TMatrixD linv(n,n);
    for(int col=0;col<n;col++)for(int i=0;i<n;i++){
        double value=(i==col?1.:0.);for(int k=0;k<i;k++)value-=u(k,i)*linv(k,col);
        linv(i,col)=value/u(i,i);
    }
    TMatrixD weighted_h=f.H;
    for(int i=0;i<n;i++)for(int c=0;c<components;c++)weighted_h(i,c)/=sd[i];
    TMatrixD a=linv*weighted_h,scaled=a;
    TVectorD scales(components);
    for(int c=0;c<components;c++){
        for(int i=0;i<n;i++)scales[c]+=a(i,c)*a(i,c);
        scales[c]=std::sqrt(scales[c]);if(!(scales[c]>0))throw std::runtime_error("Missing fitted component");
        for(int i=0;i<n;i++)scaled(i,c)/=scales[c];
    }
    TDecompSVD svd(scaled);
    if(!svd.Decompose())throw std::runtime_error("SVD failed");
    const TVectorD& singular=svd.GetSig();
    if(singular[components-1]<=singular[0]*std::max(n,components)*std::numeric_limits<double>::epsilon())throw std::runtime_error("Rank deficient component fit");
    f.scaled_condition=singular[0]/singular[components-1];
    const TMatrixD& su=svd.GetU();const TMatrixD& sv=svd.GetV();
    TMatrixD pinv(components,n);
    for(int c=0;c<components;c++)for(int i=0;i<n;i++)for(int k=0;k<components;k++)pinv(c,i)+=sv(c,k)*su(i,k)/(singular[k]*scales[c]);
    TMatrixD g=pinv*linv;
    for(int c=0;c<components;c++)for(int i=0;i<n;i++)g(c,i)/=sd[i];
    for(int c=0;c<components;c++){
        double fitted=0;for(int i=0;i<n;i++)fitted+=g(c,i)*observed[i];
        f.extracted_asymmetry_ppb[c]=fitted;
        if(std::abs(fitted-f.reference_asymmetry_ppb[c])>1e-8*(1+std::abs(f.reference_asymmetry_ppb[c])))throw std::runtime_error("Reference closure failed");
    }
    f.M.ResizeTo(components,components);f.M=transpose(a)*a;
    f.experimental_covariance.ResizeTo(components,components);f.experimental_covariance=pinv*transpose(pinv);
    // Propagate correlated dilution, asymmetry-shape and reference errors in H.
    TMatrixD j(n,components*dim);
    for(int c=0;c<components;c++)for(int i=0;i<n;i++){
        j(i,c*dim+i)=-observed[i]/r[i];j(i,c*dim+d+i)=1/r[i];
        j(i,c*dim+n)=f.H(i,c)*f.reference_asymmetry_ppb[c]/y(c,n);
        j(i,c*dim+d+n)=-f.H(i,c)/y(c,n);
    }
    TMatrixD derivative=g*j;
    f.mc_covariance.ResizeTo(components,components);f.mc_covariance=derivative*v*transpose(derivative);
    double max_mc=1;
    for(int i=0;i<components;i++)for(int k=i;k<components;k++){
        double mean=.5*(f.mc_covariance(i,k)+f.mc_covariance(k,i));f.mc_covariance(i,k)=f.mc_covariance(k,i)=mean;max_mc=std::max(max_mc,std::abs(mean));
    }
    for(int c=0;c<components;c++)if(f.mc_covariance(c,c)<-1e-9*max_mc)throw std::runtime_error("Negative propagated MC variance");
    f.combined_covariance.ResizeTo(components,components);f.combined_covariance=f.experimental_covariance+f.mc_covariance;
    f.F_mc_sigma.ResizeTo(n,components);
    for(int i=0;i<n;i++)for(int c=0;c<components;c++){
        double variance=0;
        for(int k=0;k<components;k++)for(int l=0;l<components;l++)variance+=((c==k?1.:0.)-f.F(i,c))*v(k*dim+i,l*dim+i)*((c==l?1.:0.)-f.F(i,c))/(r[i]*r[i]);
        f.F_mc_sigma(i,c)=std::sqrt(std::max(0.,variance));
    }
    for(const auto* matrix:{&f.F,&f.H,&f.C,&f.correlation,&f.M,&f.experimental_covariance,&f.mc_covariance,&f.combined_covariance,&f.F_mc_sigma})
        for(int i=0;i<matrix->GetNrows();i++)for(int j=0;j<matrix->GetNcols();j++)
            if(!std::isfinite((*matrix)(i,j)))throw std::runtime_error("Nonfinite fit output");
    return f;
}
}


namespace standalone_main_pe {
struct remollGenericDetectorHit_t {int det=0;double x=0,y=0,z=0,xl=0,yl=0,zl=0,px=0,py=0,pz=0;};
using std::vector;
// Source: vendor/main_response/lookuptable.C | SHA256 66ff13be9912a7e092a1a6cb33851cda13b492ec3e00a154f889bda446c8b0b0
/* 	
    This file provides utility functions to load lookup tables for electron position and energy response, 
    and to calculate the mean photoelectron (PE) yield based on hit position and energy for the thin quartz detectors.

    
    - There are five functions you need to use to calculate the PE yield:
        1. LoadLookupTable_electron_position: Loads the position-based lookup table from a CSV file.
        2. LoadLookupTable_electron_energy: Loads the energy-based lookup table from a CSV file.
        3. GetMeanPE(this is just for testing): Returns the mean PE yield for a given hit position (h, v) using the position lookup table. 
        4. cal_PE: Calculates the PE yield for a given hit position (h, v), polar angle (theta), and azimuthal angle (phi) using the position lookup table.
        5. cal_PE_energy: Returns the energy scaling factor for a given hit energy using the energy lookup table.
    How to use:
    -----------
    1. Include this file in your ROOT macro or C++ analysis script.
    2. load the lookup tables for the desired ring number:
           std::vector<tableEntry> lookup = LoadLookupTable_electron_position(Form("table/r%d_LookupTable.csv", ringnumber));
           std::vector<tableEntry_energy> lookup_energy = LoadLookupTable_electron_energy(Form("table/r%d_EnergyDepTable.csv", ringnumber));
    3. For each hit (of type RemollHit), call:
           double pe = cal_PE(hit.xl, hit.yl, theta, phi, lookup);
           double energy_scale = cal_PE_energy(hit.e, lookup_energy);
    4. Use the returned values for further analysis or histogramming.       

    All functions are documented below.
*/

// Structure to hold position-based lookup table entries
struct tableEntry {
    double hmin, hmax;
    double vmin, vmax;
    double meanPE;
    double langauPE;
    double RMS, Resolution;
    double p0, p1, p2, p3, p4; // fit parameters
};
// Structure to hold energy-based lookup table entries
struct tableEntry_energy {
    double energy, scale;
};

// Function to load position-based lookup table from CSV file
std::vector<tableEntry> LoadLookupTable_electron_position(const std::string& filename) {
    std::vector<tableEntry> lookupTable;
    lookupTable.clear();
    std::ifstream file(filename);
    std::string line;
    getline(file, line); // skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        tableEntry e;
        std::string tmp;
        std::getline(ss, tmp, ','); e.hmin = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.hmax = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.vmin = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.vmax = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.meanPE = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.langauPE = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.RMS = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.Resolution = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.p0 = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.p1 = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.p2 = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.p3 = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.p4 = atof(tmp.c_str());
        lookupTable.push_back(e);
    }
    std::cout << "Loaded " << lookupTable.size() << " New table entries.\n";
    return lookupTable;
}
// Function to load energy-based lookup table from CSV file
std::vector<tableEntry_energy> LoadLookupTable_electron_energy(const std::string& filename) {
    std::vector<tableEntry_energy> lookupTable;
    lookupTable.clear();
    std::ifstream file(filename);
    std::string line;
    getline(file, line); // skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        tableEntry_energy e;
        std::string tmp;
        std::getline(ss, tmp, ','); e.energy = atof(tmp.c_str());
        std::getline(ss, tmp, ','); e.scale = atof(tmp.c_str());
        lookupTable.push_back(e);
    }
    std::cout << "Loaded " << lookupTable.size() << " Energy table entries.\n";
    return lookupTable;
}

// Function to get mean PE yield for a given hit position (h, v) using the position lookup table
double GetMeanPE(double h, double v,vector<tableEntry>& lookupTable) {
    for (auto& e : lookupTable) {
        if (h >= e.hmin && h < e.hmax && v >= e.vmin && v < e.vmax)
            return e.meanPE;
    }
    return -1; // not found
}
// Function to calculate PE yield for a given hit position (h, v), polar angle (theta), and azimuthal angle (phi) using the position lookup table
double cal_PE(double h, double v,double theta,double phi,vector<tableEntry>& lookupTable) {
    for (auto& e : lookupTable) {
        if (h > e.hmin-1e-9 && h < e.hmax+1e-9 && v >= e.vmin-1e-9 && v < e.vmax+1e-9) {
            double theta1 = theta;
            double phi1 = TMath::DegToRad() * phi;
            double px = e.p0;
            double py = e.p0+e.p1*theta1+e.p2*theta1*theta1+e.p3*theta1*theta1*theta1+e.p4*theta1*theta1*theta1*theta1;
            if(py == 0){
                return 0;
            }
            double PE = px*py/sqrt(pow(px*sin(phi1),2)+pow(py*cos(phi1),2));
            return PE;
        }
    }
    return -1; // not found if position is not in the table, return -1 not to be used for further analysis
}
//  Function to return the energy scaling factor for a given hit energy using the energy lookup table
double cal_PE_energy(double energy,vector<tableEntry_energy>& lookupTable) {
    for (auto& e : lookupTable) {
        if (energy >= e.energy-1 && energy < e.energy+1) {
            return e.scale;
        }
    }
    return 1; // not found if energy is not in the table, return 1 (no scaling)
}

// Source: vendor/main_response/convert_labframe_to_local_thinquartz.C | SHA256 20853195e7d9e27aa3e757d8b69b750a4e0e3a4e4eb7f742592f4528b80b6814
/* 	
    This file provides utility functions to convert detector hit data from the lab frame
    to the local coordinate system of the main detectors (thin quartz).
    
    - The main function is `convert_labframe_to_local_thinquartz`, which takes a single hit
      and returns a `mainquartz_local_info` struct with local coordinates, rotated momentum,
      local angles, wedge information, and detector ID.
    - Several helper functions are provided for coordinate and cut calculations.
    - The detector geometry parameters are stored in `maindetector_detid_map`.
	---Important Note---
	This function can only be used for the thin quartz detectors which my naming system and the hits on the front plane or wedge plane of the detector!
	if you want to use it for other names or other detectors, please modify the map and the cut function accordingly.
	if you want to include the back face of the detector, please comment out the is_edge part in the main function.

    How to use:
    -----------
    1. Include this file in your ROOT macro or C++ analysis script.
    2. For each hit (of type RemollHit), call:
           mainquartz_local_info info = convert_labframe_to_local_thinquartz(hit);
    3. Use the returned struct for further analysis or histogramming.

    All functions are documented below.
*/

bool equal(double a, double b, double epsilon = 1e-4) {
	return std::fabs(a - b) < epsilon;
}

// Type definitions for detector hit and hit list
typedef remollGenericDetectorHit_t RemollHit;
typedef std::vector<RemollHit> hit_list;

// Structure to hold processed hit information
struct mainquartz_local_info {
    double lpx;     // Rotated px
    double lpy;     // Rotated py
    double lpz;     // Rotated pz
    double lthx;    // Local theta x
    double lthy;    // Local theta y
    double lazimuthal;  // Local azimuthal angle momentum direction in xy plane.
	double lpolar;  // Local polar angle momentum direction from z axis which is theta consider the beam line
	bool is_wedge; // Is hit on wedge 
	bool is_edge;  // Is hit on front face
    double detid;   // Detector ID
};

// Detector parameter map: key is detector ID, value is tuple of parameters
// you can change the parameters here to match your detector setup
// Currently we have angleX, angleZ, height, thickness
std::map <int,std::tuple<double,double,double,double,double>> maindetector_detid_map{
	{150151, {0.0000, -2.9976, 90.0000,70 ,8.5 }},
	{150153, {0.0000, -2.9976, 90.0000,70 ,8.5 }},
	{150152, {0.0000, -2.9976, 90.0000,70 ,8.5 }},
	{150140, {0.0000, -2.9971, 90.0000,60 ,10 }},
	{150130, {0.0000, -3.0022, 90.0000,30 ,10 }},
	{150120, {0.0000, -3.0022, 90.0000,30 ,10 }},
	{150110, {0.0000, -2.9923, 90.0000,5 ,10 }},
	{150160, {0.0000, -2.9997, 90.0000,50 ,10 }},
	{110152, {0.6687, -2.9257, 77.1193,70 ,8.5 }},
	{110153, {0.6687, -2.9256, 77.1193,70 ,8.5 }},
	{110151, {0.6688, -2.9259, 77.1193,70 ,8.5 }},
	{110140, {0.6686, -2.9259, 77.1224,60 ,10 }},
	{110130, {0.6672, -2.9197, 77.1229,30 ,10 }},
	{110120, {0.6672, -2.9198, 77.1230,30 ,10 }},
	{110110, {0.6692, -2.9282, 77.1209,5 ,10 }},
	{110160, {0.6683, -2.9245, 77.1226,50 ,10 }},
	{150251, {1.3015, -2.7005, 64.2563,70 ,8.5 }},
	{150253, {1.3015, -2.7004, 64.2563,70 ,8.5 }},
	{150252, {1.3015, -2.7005, 64.2563,70 ,8.5 }},
	{150240, {1.3015, -2.7001, 64.2536,60 ,10 }},
	{150230, {1.3035, -2.7046, 64.2551,30 ,10 }},
	{150220, {1.3035, -2.7047, 64.2560,30 ,10 }},
	{150210, {1.2992, -2.6956, 64.2550,5 ,10 }},
	{150260, {1.3025, -2.7026, 64.2557,50 ,10 }},
	{110252, {1.8723, -2.3461, 51.3913,70 ,8.5 }},
	{110253, {1.8726, -2.3456, 51.3813,70 ,8.5 }},
	{110251, {1.8724, -2.3459, 51.3869,70 ,8.5 }},
	{110240, {1.8724, -2.3459, 51.3871,60 ,10 }},
	{110230, {1.8685, -2.3410, 51.3873,30 ,10 }},
	{110220, {1.8688, -2.3414, 51.3875,30 ,10 }},
	{110210, {1.8739, -2.3478, 51.3876,5 ,10 }},
	{110260, {1.8715, -2.3448, 51.3872,50 ,10 }},
	{150351, {2.3444, -1.8683, 38.5321,70 ,8.5 }},
	{150353, {2.3443, -1.8686, 38.5366,70 ,8.5 }},
	{150352, {2.3443, -1.8686, 38.5366,70 ,8.5 }},
	{150340, {2.3441, -1.8683, 38.5346,60 ,10 }},
	{150330, {2.3481, -1.8714, 38.5339,30 ,10 }},
	{150320, {2.3481, -1.8714, 38.5340,30 ,10 }},
	{150310, {2.3402, -1.8651, 38.5333,5 ,10 }},
	{150360, {2.3462, -1.8698, 38.5331,50 ,10 }},
	{110352, {2.7046, -1.3013, 25.6759,70 ,8.5 }},
	{110353, {2.7045, -1.3016, 25.6823,70 ,8.5 }},
	{110351, {2.7045, -1.3016, 25.6823,70 ,8.5 }},
	{110340, {2.7044, -1.3015, 25.6814,60 ,10 }},
	{110330, {2.6989, -1.2988, 25.6808,30 ,10 }},
	{110320, {2.6989, -1.2989, 25.6813,30 ,10 }},
	{110310, {2.7072, -1.3028, 25.6806,5 ,10 }},
	{110360, {2.7033, -1.3010, 25.6809,50 ,10 }},
	{150451, {2.9226, -0.6667, 12.8396,70 ,8.5 }},
	{150453, {2.9223, -0.6670, 12.8466,70 ,8.5 }},
	{150452, {2.9226, -0.6666, 12.8380,70 ,8.5 }},
	{150440, {2.9219, -0.6666, 12.8409,60 ,10 }},
	{150430, {2.9269, -0.6679, 12.8431,30 ,10 }},
	{150420, {2.9273, -0.6678, 12.8398,30 ,10 }},
	{150410, {2.9173, -0.6656, 12.8410,5 ,10 }},
	{150460, {2.9248, -0.6672, 12.8390,50 ,10 }},
	{110452, {3.0011, 0.0000, 0.0000,70 ,8.5 }},
	{110453, {3.0013, 0.0000, 0.0000,70 ,8.5 }},
	{110451, {3.0011, 0.0000, 0.0000,70 ,8.5 }},
	{110440, {3.0012, 0.0001, -0.0027,60 ,10 }},
	{110430, {2.9951, 0.0000, 0.0000,30 ,10 }},
	{110420, {2.9951, 0.0000, 0.0000,30 ,10 }},
	{110410, {3.0038, 0.0002, -0.0034,5 ,10 }},
	{110460, {2.9997, 0.0000, 0.0000,50 ,10 }},
	{150551, {2.9226, 0.6667, -12.8396,70 ,8.5 }},
	{150553, {2.9226, 0.6666, -12.8380,70 ,8.5 }},
	{150552, {2.9225, 0.6666, -12.8381,70 ,8.5 }},
	{150540, {2.9222, 0.6667, -12.8409,60 ,10 }},
	{150530, {2.9271, 0.6678, -12.8401,30 ,10 }},
	{150520, {2.9270, 0.6676, -12.8367,30 ,10 }},
	{150510, {2.9173, 0.6655, -12.8384,5 ,10 }},
	{150560, {2.9248, 0.6672, -12.8390,50 ,10 }},
	{110552, {2.7043, 1.3017, -25.6854,70 ,8.5 }},
	{110553, {2.7043, 1.3017, -25.6854,70 ,8.5 }},
	{110551, {2.7044, 1.3020, -25.6888,70 ,8.5 }},
	{110540, {2.7044, 1.3017, -25.6850,60 ,10 }},
	{110530, {2.6989, 1.2991, -25.6848,30 ,10 }},
	{110520, {2.6988, 1.2991, -25.6855,30 ,10 }},
	{110510, {2.7066, 1.3027, -25.6837,5 ,10 }},
	{110560, {2.7033, 1.3012, -25.6838,50 ,10 }},
	{150651, {2.3446, 1.8682, -38.5265,70 ,8.5 }},
	{150653, {2.3443, 1.8686, -38.5366,70 ,8.5 }},
	{150652, {2.3445, 1.8684, -38.5321,70 ,8.5 }},
	{150640, {2.3441, 1.8682, -38.5325,60 ,10 }},
	{150630, {2.3481, 1.8713, -38.5320,30 ,10 }},
	{150620, {2.3482, 1.8713, -38.5315,30 ,10 }},
	{150610, {2.3402, 1.8651, -38.5333,5 ,10 }},
	{150660, {2.3462, 1.8698, -38.5313,50 ,10 }},
	{110652, {1.8723, 2.3458, -51.3869,70 ,8.5 }},
	{110653, {1.8722, 2.3461, -51.3913,70 ,8.5 }},
	{110651, {1.8721, 2.3460, -51.3913,70 ,8.5 }},
	{110640, {1.8721, 2.3460, -51.3926,60 ,10 }},
	{110630, {1.8684, 2.3413, -51.3915,30 ,10 }},
	{110620, {1.8684, 2.3413, -51.3920,30 ,10 }},
	{110610, {1.8737, 2.3480, -51.3923,5 ,10 }},
	{110660, {1.8713, 2.3450, -51.3921,50 ,10 }},
	{150751, {1.3015, 2.7004, -64.2563,70 ,8.5 }},
	{150753, {1.3018, 2.7004, -64.2498,70 ,8.5 }},
	{150752, {1.3015, 2.7006, -64.2563,70 ,8.5 }},
	{150740, {1.3015, 2.7000, -64.2524,60 ,10 }},
	{150730, {1.3038, 2.7048, -64.2511,30 ,10 }},
	{150720, {1.3038, 2.7050, -64.2531,30 ,10 }},
	{150710, {1.2992, 2.6953, -64.2535,5 ,10 }},
	{150760, {1.3026, 2.7024, -64.2528,50 ,10 }},
	{110752, {0.6683, 2.9260, -77.1279,70 ,8.5 }},
	{110753, {0.6683, 2.9260, -77.1279,70 ,8.5 }},
	{110751, {0.6684, 2.9259, -77.1263,70 ,8.5 }},
	{110740, {0.6683, 2.9260, -77.1276,60 ,10 }},
	{110730, {0.6670, 2.9202, -77.1287,30 ,10 }},
	{110720, {0.6669, 2.9203, -77.1292,30 ,10 }},
	{110710, {0.6692, 2.9288, -77.1242,5 ,10 }},
	{110760, {0.6681, 2.9246, -77.1269,50 ,10 }},
	{150851, {0.0000, 2.9976, -90.0000,70 ,8.5 }},
	{150853, {0.0000, 2.9974, -90.0000,70 ,8.5 }},
	{150852, {0.0000, 2.9976, -90.0000,70 ,8.5 }},
	{150840, {0.0001, 2.9971, -89.9973,60 ,10 }},
	{150830, {0.0002, 3.0022, -89.9970,30 ,10 }},
	{150820, {0.0000, 3.0022, -90.0000,30 ,10 }},
	{150810, {0.0000, 2.9917, -90.0000,5 ,10 }},
	{150860, {0.0001, 2.9997, -89.9978,50 ,10 }},
	{110852, {-0.6684, 2.9259, -102.8737,70 ,8.5 }},
	{110853, {-0.6686, 2.9257, -102.8791,70 ,8.5 }},
	{110851, {-0.6684, 2.9259, -102.8737,70 ,8.5 }},
	{110840, {-0.6685, 2.9259, -102.8750,60 ,10 }},
	{110830, {-0.6670, 2.9198, -102.8742,30 ,10 }},
	{110820, {-0.6672, 2.9201, -102.8770,30 ,10 }},
	{110810, {-0.6690, 2.9283, -102.8758,5 ,10 }},
	{110860, {-0.6682, 2.9246, -102.8752,50 ,10 }},
	{150951, {-1.3015, 2.7005, -115.7437,70 ,8.5 }},
	{150953, {-1.3016, 2.7004, -115.7468,70 ,8.5 }},
	{150952, {-1.3012, 2.7007, -115.7373,70 ,8.5 }},
	{150940, {-1.3013, 2.7002, -115.7428,60 ,10 }},
	{150930, {-1.3035, 2.7048, -115.7421,30 ,10 }},
	{150920, {-1.3035, 2.7047, -115.7426,30 ,10 }},
	{150910, {-1.2991, 2.6956, -115.7435,5 ,10 }},
	{150960, {-1.3025, 2.7026, -115.7433,50 ,10 }},
	{110952, {-1.8724, 2.3459, -128.6131,70 ,8.5 }},
	{110953, {-1.8722, 2.3461, -128.6087,70 ,8.5 }},
	{110951, {-1.8722, 2.3461, -128.6087,70 ,8.5 }},
	{110940, {-1.8724, 2.3461, -128.6112,60 ,10 }},
	{110930, {-1.8685, 2.3412, -128.6104,30 ,10 }},
	{110920, {-1.8685, 2.3413, -128.6100,30 ,10 }},
	{110910, {-1.8736, 2.3476, -128.6098,5 ,10 }},
	{110960, {-1.8714, 2.3449, -128.6093,50 ,10 }},
	{151051, {-2.3444, 1.8683, -141.4679,70 ,8.5 }},
	{151053, {-2.3443, 1.8686, -141.4634,70 ,8.5 }},
	{151052, {-2.3442, 1.8685, -141.4623,70 ,8.5 }},
	{151040, {-2.3439, 1.8683, -141.4637,60 ,10 }},
	{151030, {-2.3480, 1.8715, -141.4638,30 ,10 }},
	{151020, {-2.3481, 1.8714, -141.4660,30 ,10 }},
	{151010, {-2.3402, 1.8652, -141.4646,5 ,10 }},
	{151060, {-2.3461, 1.8700, -141.4638,50 ,10 }},
	{111052, {-2.7044, 1.3016, -154.3177,70 ,8.5 }},
	{111053, {-2.7045, 1.3016, -154.3177,70 ,8.5 }},
	{111051, {-2.7044, 1.3016, -154.3177,70 ,8.5 }},
	{111040, {-2.7046, 1.3016, -154.3186,60 ,10 }},
	{111030, {-2.6990, 1.2990, -154.3179,30 ,10 }},
	{111020, {-2.6989, 1.2991, -154.3145,30 ,10 }},
	{111010, {-2.7064, 1.3025, -154.3193,5 ,10 }},
	{111060, {-2.7032, 1.3009, -154.3191,50 ,10 }},
	{151151, {-2.9226, 0.6667, -167.1604,70 ,8.5 }},
	{151153, {-2.9225, 0.6667, -167.1604,70 ,8.5 }},
	{151152, {-2.9226, 0.6667, -167.1604,70 ,8.5 }},
	{151140, {-2.9219, 0.6668, -167.1565,60 ,10 }},
	{151130, {-2.9273, 0.6680, -167.1563,30 ,10 }},
	{151120, {-2.9269, 0.6679, -167.1570,30 ,10 }},
	{151110, {-2.9172, 0.6656, -167.1583,5 ,10 }},
	{151160, {-2.9247, 0.6674, -167.1567,50 ,10 }},
	{111152, {-3.0013, 0.0000, -180.0000,70 ,8.5 }},
	{111153, {-3.0013, 0.0000, -180.0000,70 ,8.5 }},
	{111151, {-3.0013, 0.0000, -180.0000,70 ,8.5 }},
	{111140, {-3.0012, 0.0000, -180.0000,60 ,10 }},
	{111130, {-2.9951, 0.0000, -180.0000,30 ,10 }},
	{111120, {-2.9951, 0.0000, -180.0000,30 ,10 }},
	{111110, {-3.0038, 0.0000, -180.0000,5 ,10 }},
	{111160, {-3.0000, 0.0000, -180.0000,50 ,10 }},
	{151251, {-2.9225, -0.6662, 167.1689,70 ,8.5 }},
	{151253, {-2.9226, -0.6666, 167.1620,70 ,8.5 }},
	{151252, {-2.9225, -0.6663, 167.1673,70 ,8.5 }},
	{151240, {-2.9222, -0.6665, 167.1623,60 ,10 }},
	{151230, {-2.9271, -0.6676, 167.1628,30 ,10 }},
	{151220, {-2.9270, -0.6677, 167.1609,30 ,10 }},
	{151210, {-2.9175, -0.6654, 167.1623,5 ,10 }},
	{151260, {-2.9249, -0.6671, 167.1636,50 ,10 }},
	{111252, {-2.7043, -1.3017, 154.3146,70 ,8.5 }},
	{111253, {-2.7043, -1.3017, 154.3146,70 ,8.5 }},
	{111251, {-2.7045, -1.3016, 154.3177,70 ,8.5 }},
	{111240, {-2.7046, -1.3017, 154.3186,60 ,10 }},
	{111230, {-2.6990, -1.2990, 154.3179,30 ,10 }},
	{111220, {-2.6989, -1.2991, 154.3145,30 ,10 }},
	{111210, {-2.7069, -1.3028, 154.3179,5 ,10 }},
	{111260, {-2.7033, -1.3010, 154.3182,50 ,10 }},
	{151351, {-2.3447, -1.8682, 141.4735,70 ,8.5 }},
	{151353, {-2.3447, -1.8682, 141.4735,70 ,8.5 }},
	{151352, {-2.3446, -1.8681, 141.4735,70 ,8.5 }},
	{151340, {-2.3440, -1.8679, 141.4691,60 ,10 }},
	{151330, {-2.3483, -1.8711, 141.4723,30 ,10 }},
	{151320, {-2.3480, -1.8710, 141.4710,30 ,10 }},
	{151310, {-2.3404, -1.8649, 141.4715,5 ,10 }},
	{151360, {-2.3464, -1.8698, 141.4700,50 ,10 }},
	{111352, {-1.8722, -2.3460, 128.6087,70 ,8.5 }},
	{111353, {-1.8724, -2.3459, 128.6131,70 ,8.5 }},
	{111351, {-1.8724, -2.3459, 128.6131,70 ,8.5 }},
	{111340, {-1.8723, -2.3461, 128.6091,60 ,10 }},
	{111330, {-1.8685, -2.3412, 128.6109,30 ,10 }},
	{111320, {-1.8685, -2.3412, 128.6105,30 ,10 }},
	{111310, {-1.8740, -2.3478, 128.6146,5 ,10 }},
	{111360, {-1.8715, -2.3451, 128.6097,50 ,10 }},
	{151451, {-1.3017, -2.7002, 115.7502,70 ,8.5 }},
	{151453, {-1.3017, -2.7002, 115.7502,70 ,8.5 }},
	{151452, {-1.3017, -2.7002, 115.7502,70 ,8.5 }},
	{151440, {-1.3016, -2.7000, 115.7488,60 ,10 }},
	{151430, {-1.3038, -2.7046, 115.7489,30 ,10 }},
	{151420, {-1.3038, -2.7044, 115.7512,30 ,10 }},
	{151410, {-1.2996, -2.6958, 115.7496,5 ,10 }},
	{151460, {-1.3026, -2.7024, 115.7472,50 ,10 }},
	{111452, {-0.6688, -2.9259, 102.8807,70 ,8.5 }},
	{111453, {-0.6687, -2.9259, 102.8791,70 ,8.5 }},
	{111451, {-0.6683, -2.9258, 102.8721,70 ,8.5 }},
	{111440, {-0.6685, -2.9259, 102.8750,60 ,10 }},
	{111430, {-0.6671, -2.9202, 102.8742,30 ,10 }},
	{111420, {-0.6671, -2.9202, 102.8739,30 ,10 }},
	{111410, {-0.6689, -2.9283, 102.8725,5 ,10 }},
	{111460, {-0.6682, -2.9246, 102.8752,50 ,10 }}
};





/**
 * @brief Read all lines from a text file.
 * @param filename Path to the file.
 * @return Vector of strings, each string is a line from the file.
 *
 * Usage: auto files = readlines("filelist.txt");
 */
std::vector<std::string> readlines(std::string filename){
	std::vector<std::string> all_lines;
	std::string line; std::ifstream fileobj(filename);
	while(std::getline(fileobj,line)){
		all_lines.push_back(std::move(line));
	}
	return all_lines;
}


/**
 * @brief Apply wedge cut to a hit. defined if hits hitting on the wedge of the main detector or not. 
 * @param height Detector height.
 * @param hit Detector hit object.
 * @return True if hit is on the wedge, false otherwise.
 */
bool def_wedge(RemollHit& hit,double height){
	 
	if((hit.det/10)%10 == 1){
		if(hit.yl < 5){ // this is specific for ring 1 detector, becasue the the mean part is from -15 to 5 mm based on the new design.
			return false;
		}else{
			return true;
		}
	}else{
		if (hit.yl < height){
			return false;
		}else{
			return true;
		}
	}
}


/**
 * @brief Apply front cut to a hit. we are only interested in hits hitting on the front of the main detector.
 * @param thickness Detector thickness.
 * @param height Detector height.
 * @param hit Detector hit object.
 * @return True if hit passes front cut, false otherwise.
 */
bool cut_front(RemollHit& hit,double height,double thickness){ //dx means the width of the detector
	if(def_wedge(hit,height)){
		double b = height + thickness - 0.5;
		if(hit.yl > hit.zl + b){
			return true;
		}else{
			return false;
		}
	}else{
		if(hit.zl < 0.5 - thickness ){
			return true;
		}else{
			return false;
		}
	}
		
}

/**
 * @brief Rotate a 3D vector by given Z and X angles. basicly rotate the coordinate system for our detector.
 * @param point Input vector (TVector3).
 * @param rotateZ Rotation angle around Z axis (degrees).
 * @param rotateX Rotation angle around X axis (degrees).
 * @return Rotated vector (TVector3).
 */
TVector3 cal_newP(TVector3 point, double rotateX, double rotateY, double rotateZ){
	TRotation rotation;
	
   

	// Rotate around the X-axis
    rotation.RotateX(rotateX * TMath::DegToRad());

	// Rotate around the Y-axis
	rotation.RotateY(rotateY * TMath::DegToRad());

	// Rotat around the Z-axis
	rotation.RotateZ(rotateZ * TMath::DegToRad());


	// Apply the rotation to the point
	return rotation * point;
}

/**
 * @brief Convert a single hit from lab frame to local detector coordinates.
 *        Returns a struct with all relevant local and global information.
 * @param hit The detector hit (RemollHit).
 * @return mainquartz_local_info struct with local coordinate angle, rotated momentum, wedge information and detector ID.
 *
 * Usage example:
 *     mainquartz_local_info info = convert_labframe_to_local_thinquartz(hit);
 */
mainquartz_local_info convert_labframe_to_local_thinquartz(RemollHit hit){
	mainquartz_local_info info;
	auto [angleX, angleY, angleZ, height, thickness] = maindetector_detid_map[hit.det];
	double rotate_z = angleZ;
	double rotate_x = angleX;
	double rotate_y = angleY;

	bool is_edge = cut_front(hit,height,thickness);
	//if you want to include the back face, just comment out the is_edge part
	if(!is_edge){
		//std::cerr << "Error: hit not in front face or wedge face" << std::endl;
        //return mainquartz_local_info{};
	}
	//bool is_wedge = def_wedge(hit,height);
	//if (info.lx != hit.xl) std::cout << "lx incorrect" << std::endl;

	//wstd::cout << hit.yl << " " << info.ly << std::endl;
	TVector3 point(hit.px, hit.py, hit.pz);
	TVector3 rotatedMomentum = cal_newP(point,rotate_x,rotate_y,rotate_z);
	//all the angle here is the angle between the momentum and the z axis

	info.lpx = rotatedMomentum.X(); 
	info.lpy = rotatedMomentum.Y();
	info.lpz = rotatedMomentum.Z();
	info.lthx = atan2(rotatedMomentum.Z(),rotatedMomentum.X())*180.0/M_PI;
	info.lthy = atan2(rotatedMomentum.Z(),rotatedMomentum.Y())*180.0/M_PI;
	info.lpolar = atan2(sqrt(rotatedMomentum.X()*rotatedMomentum.X()+rotatedMomentum.Y()*rotatedMomentum.Y()),rotatedMomentum.Z())*180.0/M_PI;
	info.lazimuthal = atan2(rotatedMomentum.Y(),rotatedMomentum.X())*180.0/M_PI;
	info.is_wedge = def_wedge(hit,height);
	info.is_edge = is_edge;
	info.detid = hit.det;
	return info;
}	


}
// Source: root/pe_response.hh | SHA256 9df5aa91b671acfc1fcaaa9db5b1a878674fed0db9f7bcb1909f7c0d0a0fdb0b
namespace standalone_deconv {
namespace fs=std::filesystem;
struct MainModel {
    std::array<std::vector<standalone_main_pe::tableEntry>,6> position;
    std::array<std::vector<standalone_main_pe::tableEntry_energy>,6> energy;
    explicit MainModel(const fs::path& folder){for(int r=0;r<6;r++){
        auto stem=folder/("r"+std::to_string(r+1));
        // Strict parsing; keep the upstream fit functions unchanged.
        const auto positions=standalone_portable::numeric_table<13>(stem.string()+"_NewLookupTable.csv");
        for(const auto& p:positions)if(p[0]>=p[1]||p[2]>=p[3])throw std::runtime_error("Invalid main PE position bounds: "+stem.string());
        const auto energies=standalone_portable::numeric_table<2>(stem.string()+"_EnergyDepTable.csv");
        if(energies.size()!=500)throw std::runtime_error("Main PE energy table needs 500 rows: "+stem.string());
        for(size_t i=0;i<energies.size();i++)if(std::abs(energies[i][0]-(2.511+2*i))>.001||energies[i][1]<0)
            throw std::runtime_error("Invalid main PE energy grid/scaling: "+stem.string());
        for(const auto& p:positions)position[r].push_back({p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12]});
        for(const auto& e:energies)energy[r].push_back({e[0],e[1]});
    }}
    // The vendor energy table takes total energy (hit.e), not kinetic energy.
    double pe(standalone_main_pe::remollGenericDetectorHit_t h,double totalEnergy){
        if(!standalone_detector_selection::main_quartz(h.det))throw std::runtime_error("Unknown main quartz detector: "+std::to_string(h.det));
        int suffix=h.det%100,ring=suffix>=51&&suffix<=53?4:suffix/10-1;
        if(ring<0||ring>=6)throw std::runtime_error("Unknown main quartz detector");
        auto angle=standalone_main_pe::convert_labframe_to_local_thinquartz(h);
        const double offset[]={15,30,30,60,70,50};
        double value=standalone_main_pe::cal_PE(h.xl,h.yl+offset[ring],std::min(angle.lpolar,4.),angle.lazimuthal,position[ring]);
        return value<0?-1:value*standalone_main_pe::cal_PE_energy(totalEnergy,energy[ring]);
    }
};
}


// Source: root/Deconvolution.C | SHA256 afddcf1dfdd6c1f36a9ecf3d6a7b489ee84beed400c3ed058c1507b985d9cbdc
namespace standalone_deconv {
namespace fs=std::filesystem;
const char* component_names[]={"Moller","ep elastic","ep inelastic Delta","ep inelastic resonance","ep inelastic continuum"};
inline std::string text(const fs::path& p){std::ifstream f(p);if(!f)throw std::runtime_error("Cannot open "+p.string());return std::string(std::istreambuf_iterator<char>(f),{});}
inline std::string trim(std::string s){auto a=s.find_first_not_of(" \t\r\n");return a==s.npos?"":s.substr(a,s.find_last_not_of(" \t\r\n")-a+1);}
inline std::string digest(const fs::path& p){std::unique_ptr<TMD5> h(TMD5::FileChecksum(p.string().c_str()));if(!h)throw std::runtime_error("Cannot fingerprint "+p.string());return h->AsString();}
inline std::string digest_text(const std::string& s){TMD5 h;h.Update(reinterpret_cast<const UChar_t*>(s.data()),s.size());h.Final();return h.AsString();}
// Shared provenance checks for both PE maps and the deconvolution.
inline fs::path check_sampling(TFile& file,TTree* tree,const fs::path& path,std::string macro){
    // Comments are not sampling instructions. Conservatively reject recorded quota plans.
    std::istringstream lines(macro);std::string line,commands;
    while(std::getline(lines,line))commands+=line.substr(0,line.find('#'))+"\n";
    macro=commands;
    std::transform(macro.begin(),macro.end(),macro.begin(),[](unsigned char c){return std::tolower(c);});
    if(tree->GetBranch("hit_weight")||tree->GetBranch("sampling_weight"))
        throw std::runtime_error("Per-hit transport weights need an independently validated physical second moment; unsupported by this macro");
    auto plan=path;plan.replace_extension(".neyman.tsv");
    const bool fixed=macro.find("neyman")!=std::string::npos||macro.find("fixed_quota")!=std::string::npos||file.Get("sampling_plan")||fs::exists(plan);
    if(!fixed)return {};
    if(!fs::is_regular_file(plan)){
        auto canonical=fs::canonical(path);canonical.replace_extension(".neyman.tsv");
        if(fs::is_regular_file(canonical))plan=canonical;
        else throw std::runtime_error("Fixed-quota input needs its matching .neyman.tsv sampling plan: "+path.string());
    }
    const auto fingerprint=digest(plan);
    if(!trustedSamplingPlanMD5.count(fingerprint))
        throw std::runtime_error("Sampling plan identity is unconfirmed. Verify it against the generation receipt, then add its MD5 to trustedSamplingPlanMD5: "+plan.string()+" MD5 "+fingerprint);
    return plan;
}
inline std::map<std::string,std::string> saved_geometry(const remollRunData& data){
    if(data.fGitInfo.empty())throw std::runtime_error("Simulation source identity is missing.");
    std::map<std::string,std::string> gdml;
    for(const auto& g:data.fGDMLFiles){
        if(!g.fFilename||!g.fBuffer)continue;
        std::string name=g.fFilename;std::replace(name.begin(),name.end(),'\\','/');
        auto pos=name.find("/geometry/");
        if(pos!=name.npos)name=name.substr(pos+10);else if(name.rfind("geometry/",0)==0)name=name.substr(9);
        const auto hash=digest_text(g.fBuffer);
        if(gdml.count(name)&&gdml.at(name)!=hash)throw std::runtime_error("Conflicting saved geometry: "+name);
        gdml[name]=hash;
    }
    if(gdml.empty())throw std::runtime_error("Saved geometry is missing.");
    return gdml;
}
inline void check_main_geometry(const std::map<std::string,std::string>& gdml,const std::map<std::string,std::string>& expected){
    bool thin=false;
    for(const auto& item:gdml)if(item.first.rfind("detector/ThinQuartz/",0)==0){
        thin=true;
        if(!expected.count(item.first)||expected.at(item.first)!=item.second)
            throw std::runtime_error("Main quartz geometry does not match the PE model.");
    }
    if(!thin)throw std::runtime_error("Main PE needs saved thin-quartz geometry.");
}
// Known rate-corrected proposal controls; other bias commands remain significant.
inline bool proposal_setting(const std::string& command){
    return command=="/remoll/bias/beamp/mode"||command=="/remoll/bias/beamp/min"||
           command=="/remoll/bias/beamp/max"||command=="/remoll/bias/beamp/physicalFraction";
}
// Compare transport across interactions, not only generator settings within one.
// Preserve ordered commands: multiple field scale/current commands are meaningful.
inline std::vector<std::string> transport_details(const remollRunData& data,const std::string& macro){
    std::vector<std::string> identity={"source="+digest_text(data.fGitInfo)};
    for(const auto& g:saved_geometry(data))identity.push_back("geometry "+g.first+"="+g.second);
    std::vector<std::string> fields;
    for(const auto& field:data.fMagData){
        auto bounded=[](const char* value,size_t size){const auto* end=static_cast<const char*>(std::memchr(value,0,size));
            if(!end||end==value)throw std::runtime_error("Invalid saved magnetic-field identity");return std::string(value,end);};
        auto name=bounded(field.filename,sizeof(field.filename));
        std::replace(name.begin(),name.end(),'\\','/');name=name.substr(name.find_last_of('/')+1);
        fields.push_back(name+"="+bounded(field.hashsum,sizeof(field.hashsum)));
    }
    std::sort(fields.begin(),fields.end());for(const auto& f:fields)identity.push_back("field "+f);
    std::istringstream lines(macro);std::string line;
    while(std::getline(lines,line)){
        line=line.substr(0,line.find('#'));std::istringstream words(line);std::string command,word,value;words>>command;
        if(command=="/run/beamOn")break;
        if(command=="/remoll/addfield"){
            if(fields.empty())throw std::runtime_error("Saved magnetic-field fingerprints are missing");continue;
        }
        if(command.rfind("/remoll/evgen/",0)==0||command.rfind("/remoll/SD/",0)==0||
           command=="/remoll/seed"||command=="/remoll/beamcurr"||command=="/remoll/filename"||
           command=="/remoll/interval"||
           // These proposal settings change sampling, with correction already in rate.
           proposal_setting(command)||
           command=="/remoll/setgeofile"||command=="/remoll/setgeofileParallel"||
           command=="/run/initialize"||command=="/run/numberOfThreads"||command=="/run/verbose"||command=="/run/printProgress")continue;
        // Includes may move geometry after GDML loading. Compare their ordered calls;
        // never read today's local macro as evidence of what a past run executed.
        if(command.rfind("/remoll/",0)!=0&&command.rfind("/process/",0)!=0&&command.rfind("/physics/",0)!=0&&command.rfind("/run/",0)!=0&&command.rfind("/control/",0)!=0)continue;
        if(command.find("verbose")!=std::string::npos||command=="/remoll/print"||command=="/remoll/target/print")continue;
        while(words>>word)value+=(value.empty()?"":" ")+word;
        identity.push_back(command+" "+value);
    }
    return identity;
}
inline std::string transport_identity(const remollRunData& data,const std::string& macro){
    std::string text;for(const auto& item:transport_details(data,macro))text+=item+"\n";
    return digest_text(text);
}
struct Input {
    fs::path path,samplingPlan;long long events=0,generated=0,seed=0;int process=-1;double current=0;
    std::string commit,geometry,macro,transport;
    std::vector<std::string> transportDetails;
    std::map<std::string,std::string> settings;
};
inline void check_transport(const Input& a,const Input& b){
    if(a.transport==b.transport)return;
    size_t i=0;while(i<a.transportDetails.size()&&i<b.transportDetails.size()&&a.transportDetails[i]==b.transportDetails[i])++i;
    auto item=[&](const Input& r){return i<r.transportDetails.size()?r.transportDetails[i]:"<no corresponding record>";};
    const bool control=item(a).rfind("/control/",0)==0||item(b).rfind("/control/",0)==0;
    throw std::runtime_error("Incompatible transport: first differing record in "+a.path.filename().string()+" ["+item(a)+"] versus "+b.path.filename().string()+" ["+item(b)+"]."+
        (control?" Included macro contents are not stored; differing calls cannot be assumed equivalent.":" Analyze separately."));
}
inline void check_generation_settings(const Input& a,const Input& b){
    if(a.process!=b.process)return;
    std::map<std::string,std::string> left,right;
    auto physical=[](const std::string& key){std::string lower=key;std::transform(lower.begin(),lower.end(),lower.begin(),[](unsigned char c){return std::tolower(c);});return lower.rfind("/remoll/evgen/",0)==0&&lower.find("neyman")==std::string::npos;};
    for(const auto& setting:a.settings)if(physical(setting.first))left.insert(setting);
    for(const auto& setting:b.settings)if(physical(setting.first))right.insert(setting);
    if(left!=right){
        std::set<std::string> keys;for(const auto& s:left)keys.insert(s.first);for(const auto& s:right)keys.insert(s.first);
        std::string changed;
        for(const auto& key:keys)if(!left.count(key)||!right.count(key)||left.at(key)!=right.at(key))changed+=(changed.empty()?"":", ")+key;
        throw std::runtime_error("Incompatible generator settings ("+changed+") within one interaction; analyze separately. Mixed generation domains cannot be averaged.");
    }
}
// Replay remoll recording commands in order; shared by fit and PE map readers.
inline void detector_command(const std::string& command,const std::string& value,
                             const std::vector<int>& ids,std::set<int>& enabled,std::set<int>& entries){
    if(command=="/remoll/SD/disable_all")enabled.clear();
    if(command=="/remoll/SD/enable_all")enabled=std::set<int>(ids.begin(),ids.end());
    const bool enable=command=="/remoll/SD/enable"||command=="/remoll/SD/enable_range";
    const bool disable=command=="/remoll/SD/disable"||command=="/remoll/SD/disable_range";
    const bool detect=command=="/remoll/SD/detect"||command=="/remoll/SD/detect_range";
    const bool undetect=command=="/remoll/SD/undetect"||command=="/remoll/SD/undetect_range";
    if(enable||disable||detect||undetect){
        std::istringstream args(value);std::string kind,extra;double low=0,high=0;
        if(detect||undetect)args>>kind;
        args>>low;high=low;
        if(command.find("_range")!=std::string::npos)args>>high;
        if(!args||!std::isfinite(low)||!std::isfinite(high)||args>>extra)
            throw std::runtime_error("Invalid detector setting: "+command);
        std::transform(kind.begin(),kind.end(),kind.begin(),[](unsigned char c){return std::tolower(c);});
        // remoll range commands include both endpoints. Recording flags accumulate.
        for(int det:ids)if(det>=low&&det<=high){
            if(enable)enabled.insert(det);if(disable)enabled.erase(det);
            if(kind=="boundaryhits"||kind=="surfacehits"){
                if(detect)entries.insert(det);if(undetect)entries.erase(det);
            }
        }
    }
}
inline bool entry_recorded(const std::string& macro,int detector){
    std::set<int> enabled,entries;std::istringstream lines(macro);std::string line;
    while(std::getline(lines,line)){
        line=line.substr(0,line.find('#'));std::istringstream words(line);std::string command,value;words>>command;
        if(command=="/run/beamOn")break;
        if(command=="/control/execute"||command=="/control/loop"||command=="/control/foreach"){enabled.clear();entries.clear();continue;}
        std::getline(words,value);detector_command(command,trim(value),{detector},enabled,entries);
    }
    return enabled.count(detector)&&entries.count(detector);
}
inline Input inspect(const fs::path& path,const std::vector<int>& ids,bool mainPE,const std::map<std::string,std::string>& thin_md5,const std::vector<std::string>& extra={}){
    TFile file(path.string().c_str(),"READ");auto* tree=file.Get<TTree>("T");auto* saved=file.Get<remollRunData>("run_data");
    if(saved)file.GetList()->Remove(saved);
    std::unique_ptr<remollRunData> data(saved);
    if(file.IsZombie()||file.TestBit(TFile::kRecovered))throw std::runtime_error("Unreadable or recovered ROOT file: "+path.string());
    if(!tree||!data)throw std::runtime_error("Missing T or run_data in "+path.string());
    Input r;r.path=path;r.generated=data->fNthrown;r.events=tree->GetEntries();r.seed=data->fSeed;r.commit=data->fGitInfo;r.macro=data->fMacro.fBuffer?data->fMacro.fBuffer:"";
    r.samplingPlan=check_sampling(file,tree,path,r.macro);
    if(r.events<2||r.generated<r.events||(r.samplingPlan.empty()&&r.generated!=r.events))throw std::runtime_error("Generated/stored event count mismatch: run_data="+std::to_string(r.generated)+", T="+std::to_string(r.events)+" in "+path.filename().string()+". Ordinary input requires exact agreement; fixed-quota input additionally requires its full transport plan.");
    if(r.commit.empty())throw std::runtime_error("Simulation source identity is missing");
    std::set<int> enabled,entries;int runs=0;std::istringstream lines(r.macro);std::string line;
    while(std::getline(lines,line)){
        line=line.substr(0,line.find('#'));std::istringstream words(line);std::string command,value;words>>command;
        if(command=="/run/beamOn"){runs++;continue;}if(runs)continue;
        if(command=="/control/execute"||command=="/control/loop"||command=="/control/foreach"){r.settings.clear();enabled.clear();entries.clear();continue;}
        std::getline(words,value);value=trim(value);
        if(command.rfind("/remoll/",0)==0&&!value.empty())r.settings[command]=value;
        detector_command(command,value,ids,enabled,entries);
    }
    if(runs!=1)throw std::runtime_error("Cannot establish one effective beamOn");
    std::string proc=r.settings["/remoll/evgen/set"];
    r.process=proc=="moller"?0:proc=="elastic"?1:proc=="inelastic"?2:-1;
    if(r.process<0){auto found=std::find(extra.begin(),extra.end(),proc);if(found!=extra.end())r.process=3+std::distance(extra.begin(),found);}
    if(r.process<0)throw std::runtime_error("Generator is missing or not in the configured background list: "+proc);
    std::istringstream beam(r.settings["/remoll/beamcurr"]);std::string unit;beam>>r.current>>unit;
    if(unit=="mA")r.current*=1000;else if(unit=="ampere")r.current*=1e6;else if(unit=="nA")r.current*=.001;else if(unit!="uA"&&unit!="microampere")r.current=0;
    if(!(r.current>0)||!std::isfinite(r.current))throw std::runtime_error("Saved beam current is missing/invalid");
    for(int det:ids)if(!enabled.count(det)||!entries.count(det))throw std::runtime_error("Enablement and entry recording are not established for detector "+std::to_string(det));
    const auto gdml=saved_geometry(*data);
    if(mainPE)check_main_geometry(gdml,thin_md5);
    for(const auto& g:gdml)r.geometry+=g.first+" "+g.second+"\n";
    r.transport=transport_identity(*data,r.macro);
    r.transportDetails=transport_details(*data,r.macro);
    return r;
}
// Physical beam energy remains an invariant even for a reviewed source mixture.
inline double beam_energy_MeV(const Input& input){
    const auto found=input.settings.find("/remoll/beamene");
    if(found==input.settings.end())throw std::runtime_error("Saved beam energy is missing");
    std::istringstream words(found->second);double value=0;std::string unit,extra;
    words>>value>>unit;
    if(!words||words>>extra||!std::isfinite(value)||value<=0)throw std::runtime_error("Invalid saved beam energy");
    if(unit=="GeV")value*=1000;else if(unit=="keV")value*=.001;else if(unit=="eV")value*=.000001;
    else if(unit!="MeV")throw std::runtime_error("Unknown saved beam energy unit: "+unit);
    return value;
}
// Same source-stratum definition as acceptance_strata/analyze_neyman_batch.C.
// Stored rates already contain the fixed-quota design correction. Only the MC
// centering changes: center within each original stratum, then add variances.
struct SamplingPlan {
    int pBins=0,thetaBins=0,invalid=-1;
    std::vector<double> pEdges,thetaEdges;
    std::vector<int> cells;
    std::map<int,std::pair<double,long long>> groups;
    explicit SamplingPlan(const fs::path& path){
        std::ifstream file(path);if(!file)throw std::runtime_error("Cannot read sampling plan: "+path.string());
        struct Row{int p,t,g;double p0,p1,t0,t1,probability;long long quota;};
        std::vector<Row> rows;std::string line;std::getline(file,line);
        while(std::getline(file,line)){
            if(trim(line).empty())continue;
            std::istringstream fields(line);std::vector<std::string> cell;std::string value;
            while(std::getline(fields,value,'\t'))cell.push_back(value);
            if(cell.size()!=9)throw std::runtime_error("Invalid sampling plan columns: "+path.string());
            Row r{std::stoi(cell[0]),std::stoi(cell[1]),std::stoi(cell[6]),std::stod(cell[2]),std::stod(cell[3]),std::stod(cell[4]),std::stod(cell[5]),std::stod(cell[7]),std::stoll(cell[8])};
            if(r.g<0||!std::isfinite(r.probability)||r.probability<0||r.quota<2)
                throw std::runtime_error("Invalid sampling plan row: "+path.string());
            const auto known=groups.find(r.g);
            if(known!=groups.end()&&(std::abs(known->second.first-r.probability)>1e-12||known->second.second!=r.quota))
                throw std::runtime_error("Inconsistent sampling group: "+path.string());
            groups[r.g]={r.probability,r.quota};rows.push_back(r);
            if(r.p>=0&&r.t>=0){pBins=std::max(pBins,r.p+1);thetaBins=std::max(thetaBins,r.t+1);}
        }
        if(rows.empty()||pBins<=0||thetaBins<=0)throw std::runtime_error("Empty sampling plan");
        pEdges.resize(pBins+1);thetaEdges.resize(thetaBins+1);cells.assign(pBins*thetaBins,-1);
        for(const auto& r:rows){
            if(r.p<0||r.t<0){invalid=r.g;continue;}
            pEdges[r.p]=r.p0;pEdges[r.p+1]=r.p1;thetaEdges[r.t]=r.t0;thetaEdges[r.t+1]=r.t1;cells[r.p*thetaBins+r.t]=r.g;
        }
        double sum=0;for(const auto& g:groups)sum+=g.second.first;
        if(std::abs(sum-1)>1e-6||!std::is_sorted(pEdges.begin(),pEdges.end())||!std::is_sorted(thetaEdges.begin(),thetaEdges.end())||
           std::find(cells.begin(),cells.end(),-1)!=cells.end())throw std::runtime_error("Incomplete sampling plan");
    }
    int stratum(double p,double theta)const{
        if(!std::isfinite(p)||!std::isfinite(theta)||p<pEdges.front()||p>=pEdges.back()||theta<thetaEdges.front()||theta>=thetaEdges.back())return invalid;
        const int ip=std::upper_bound(pEdges.begin(),pEdges.end(),p)-pEdges.begin()-1;
        const int it=std::upper_bound(thetaEdges.begin(),thetaEdges.end(),theta)-thetaEdges.begin()-1;
        return cells[ip*thetaBins+it];
    }
    long long transports()const{long long n=0;for(const auto& g:groups)n+=g.second.second;return n;}
};
struct StratumMoments {
    long long seen=0;
    std::map<int,CenteredScore> scores;
    std::unordered_map<std::uint64_t,CenteredPair> pairs;
    static std::uint64_t key(int a,int b){if(a>b)std::swap(a,b);return(std::uint64_t(a)<<32)|std::uint32_t(b);}
    void add(const std::vector<std::pair<int,double>>& values){
        for(const auto& v:values)scores[v.first].add(v.second);
        for(size_t i=0;i<values.size();i++)for(size_t j=i+1;j<values.size();j++){
            if(values[i].first<values[j].first)pairs[key(values[i].first,values[j].first)].add(values[i].second,values[j].second);
            else pairs[key(values[i].first,values[j].first)].add(values[j].second,values[i].second);
        }
    }
    double covariance(int a,int b)const{
        if(a>b)std::swap(a,b);
        const auto x=scores.find(a),y=scores.find(b);if(x==scores.end()||y==scores.end())return 0;
        const double vx=x->second.variance_of_sum(seen),vy=y->second.variance_of_sum(seen);
        if(vx==0||vy==0)return 0;if(a==b)return vx;
        const auto found=pairs.find(key(a,b));const CenteredPair empty;
        return (found==pairs.end()?empty:found->second).covariance_of_sum(x->second,y->second,seen);
    }
};
struct ScannedMoments:standalone_portable::Moments {
    std::map<int,StratumMoments> strata;
    double reference_noise=0;
    long long main_pe_candidates=0,unmatched_main=0,sm_pe_candidates=0;
    double mc_covariance(int c,int i,int other,int j,int dim,long long)const{
        double sum=0;for(const auto& group:strata)sum+=group.second.covariance(c*dim+i,other*dim+j);return sum;
    }
    void add_covariance(TMatrixD& target,int first,int dim,double scale)const{
        for(const auto& group:strata){
            const auto& state=group.second;
            std::vector<std::pair<int,const CenteredScore*>> active;
            std::vector<double> variance;
            for(const auto& s:state.scores){
                const double v=s.second.variance_of_sum(state.seen);
                if(v>0){active.push_back({s.first,&s.second});variance.push_back(v);}
            }
            for(size_t i=0;i<active.size();i++)for(size_t j=i;j<active.size();j++){
                const int a=active[i].first,b=active[j].first;
                double value=variance[i];
                if(i!=j){const auto p=state.pairs.find(StratumMoments::key(a,b));const CenteredPair zero;
                    value=(p==state.pairs.end()?zero:p->second).covariance_of_sum(*active[i].second,*active[j].second,state.seen);}
                const int row=first*dim+a,col=first*dim+b;target(row,col)+=scale*value;
                if(row!=col)target(col,row)+=scale*value;
            }
        }
    }
};
inline ScannedMoments scan(const Input& input,const std::vector<int>& ids,const standalone_portable::Fits& sm,
                                    MainModel* main,double minimumK,bool useShowerMaxPE=true){
    TFile file(input.path.string().c_str(),"READ");auto* tree=file.Get<TTree>("T");
    auto has=[&](const char* key){return tree->GetLeaf(key)||tree->GetBranch(key);};
    for(auto key:{"rate","ev.A","ev.W2","hit.det","hit.pid","hit.k","hit.e","hit.x","hit.y","hit.pz"})if(!has(key))throw std::runtime_error(std::string("Missing ROOT field: ")+key);
    if(main)for(auto key:{"hit.z","hit.xl","hit.yl","hit.zl","hit.px","hit.py"})if(!has(key))throw std::runtime_error(std::string("Main PE needs ROOT field: ")+key);
    TTreeReader reader(tree);TTreeReaderValue<double> rate(reader,"rate"),asym(reader,"ev.A"),w2(reader,"ev.W2");
    TTreeReaderArray<int> det(reader,"hit.det"),pid(reader,"hit.pid");
    TTreeReaderArray<double> k(reader,"hit.k"),e(reader,"hit.e"),x(reader,"hit.x"),y(reader,"hit.y"),pz(reader,"hit.pz");
    auto optional=[&](const char* key){return main?std::make_unique<TTreeReaderArray<double>>(reader,key):nullptr;};
    auto z=optional("hit.z"),xl=optional("hit.xl"),yl=optional("hit.yl"),zl=optional("hit.zl"),px=optional("hit.px"),py=optional("hit.py");
    std::unique_ptr<SamplingPlan> plan;
    std::unique_ptr<TTreeReaderArray<int>> partPid;
    std::unique_ptr<TTreeReaderArray<double>> partP,partX,partY,partZ;
    if(!input.samplingPlan.empty()){
        plan=std::make_unique<SamplingPlan>(input.samplingPlan);
        if(plan->transports()!=input.events)throw std::runtime_error("Sampling-plan transport count differs from the raw T tree");
        for(const char* key:{"part.pid","part.p","part.px","part.py","part.pz"})if(!has(key))throw std::runtime_error(std::string("Sampling plan needs ")+key);
        partPid=std::make_unique<TTreeReaderArray<int>>(reader,"part.pid");partP=std::make_unique<TTreeReaderArray<double>>(reader,"part.p");
        partX=std::make_unique<TTreeReaderArray<double>>(reader,"part.px");partY=std::make_unique<TTreeReaderArray<double>>(reader,"part.py");partZ=std::make_unique<TTreeReaderArray<double>>(reader,"part.pz");
    }
    const int n=ids.size(),d=n+1,dim=2*d,nc=input.process==2?3:1;
    ScannedMoments out;out.y.resize(nc*dim);out.noise.resize(n*n);out.coverage.resize(n);
    std::map<int,int> index;for(int i=0;i<n;i++)index[ids[i]]=i;
    while(reader.Next()){
        ++out.entries;
        if(!std::isfinite(*rate)||*rate<0||!std::isfinite(*asym))throw std::runtime_error("Invalid event rate/asymmetry");
        int group=0;
        if(plan){
            const size_t count=partPid->GetSize();
            for(size_t length:{partP->GetSize(),partX->GetSize(),partY->GetSize(),partZ->GetSize()})if(length!=count)throw std::runtime_error("Inconsistent primary-particle arrays");
            int best=-1,any=-1;
            for(size_t i=0;i<count;i++)if((*partPid)[i]==11){
                if(any<0||(*partP)[i]>(*partP)[any])any=i;
                if((*partZ)[i]>0&&(best<0||(*partP)[i]>(*partP)[best]))best=i;
            }
            if(best<0)best=any;
            group=plan->stratum(best<0?-1:(*partP)[best],best<0?-1:std::atan2(std::hypot((*partX)[best],(*partY)[best]),(*partZ)[best])*1000);
            if(!plan->groups.count(group))throw std::runtime_error("Raw history is outside the sampling-plan groups");
        }
        auto& stratum=out.strata[group];++stratum.seen;
        if(*rate==0)continue;size_t size=det.GetSize();
        for(size_t length:{pid.GetSize(),k.GetSize(),e.GetSize(),x.GetSize(),y.GetSize(),pz.GetSize()})if(length!=size)throw std::runtime_error("Inconsistent hit arrays");
        if(main)for(size_t length:{z->GetSize(),xl->GetSize(),yl->GetSize(),zl->GetSize(),px->GetSize(),py->GetSize()})if(length!=size)throw std::runtime_error("Inconsistent PE hit arrays");
        std::vector<double> q(d);
        for(size_t j=0;j<size;j++){
            auto found=index.find(det[j]);if(found==index.end())continue;int i=found->second;out.coverage[i]++;
            for(double v:{k[j],e[j],x[j],y[j],pz[j]})if(!std::isfinite(v))throw std::runtime_error("Nonfinite hit state");
            if(main)for(double v:{(*z)[j],(*xl)[j],(*yl)[j],(*zl)[j],(*px)[j],(*py)[j]})
                if(!std::isfinite(v))throw std::runtime_error("Nonfinite PE hit state");
            if(det[j]>=73001&&det[j]<=73028){
                // ShowerMax cut, followed by supported expected-PE lookup.
                if(pz[j]<=0)continue;
                if(!useShowerMaxPE){
                    if(standalone_portable::charged(pid[j])&&k[j]>=minimumK)q[i]+=1;
                    continue;
                }
                ++out.sm_pe_candidates;
                double pe=standalone_portable::pe(sm,det[j],pid[j],e[j],x[j],y[j]);
                if(pe<0){out.unsupported_sm++;continue;}q[i]+=pe;
            }else{
                // Main rate selection: charged particles, kinetic energy >= minimumK.
                if(!standalone_portable::charged(pid[j])||k[j]<minimumK)continue;
                // Reference stays rate-weighted in both PE and rate modes.
                if(standalone_detector_selection::ring5_open(det[j]))q[n]+=1;
                if(!main){q[i]+=1;continue;}
                // Upstream main-PE model support restrictions.
                if(std::abs(pid[j])!=11||pz[j]<=0||k[j]<=1.1)continue;
                standalone_main_pe::remollGenericDetectorHit_t h;
                h.det=det[j];h.x=x[j];h.y=y[j];h.z=(*z)[j];h.xl=(*xl)[j];h.yl=(*yl)[j];h.zl=(*zl)[j];h.px=(*px)[j];h.py=(*py)[j];h.pz=pz[j];
                ++out.main_pe_candidates;
                const double response=main->pe(h,e[j]);
                if(response<0){++out.unmatched_main;continue;}
                q[i]+=response;
            }
        }
        int c=0;
        if(input.process==2){double w=*w2>0?std::sqrt(*w2)*.001:-1;c=w>=1&&w<1.4?0:w>=1.4&&w<2.5?1:w>=2.5&&w<6?2:-1;}
        if(c<0){if(std::any_of(q.begin(),q.end(),[](double v){return v>0;}))throw std::runtime_error("Signal outside inelastic W bins [1,6) GeV");continue;}
        std::vector<int> nz;for(int i=0;i<d;i++)if(q[i])nz.push_back(i);
        std::vector<std::pair<int,double>> values;values.reserve(2*nz.size());
        for(int i:nz){
            double yi=*rate*q[i];out.y[c*dim+i]+=yi;out.y[c*dim+d+i]+=yi*(*asym);
            if(yi)values.push_back({c*dim+i,yi});
            if(yi*(*asym))values.push_back({c*dim+d+i,yi*(*asym)});
            for(int j:nz){
                if(i<n&&j<n)out.noise[i*n+j]+=*rate*q[i]*q[j];
            }
        }
        stratum.add(values);out.reference_noise+=*rate*q[n]*q[n];
    }
    if(out.entries!=input.events)throw std::runtime_error("Incomplete ROOT read");
    if(plan)for(const auto& g:plan->groups)if(out.strata[g.first].seen!=g.second.second)
        throw std::runtime_error("Raw history count does not match sampling stratum "+std::to_string(g.first));
    return out;
}
inline void write_fit(TFile& file,Fit& fit,const std::vector<int>& ids,const TMatrixD& sums,const std::string& provenance,const std::vector<std::string>& labels={}){
    file.cd();
    for(auto item:std::initializer_list<std::pair<const char*,TMatrixD*>>{
        {"rates",&fit.rates},{"F",&fit.F},{"F_mc_sigma",&fit.F_mc_sigma},{"H",&fit.H},{"B",&fit.B},{"C",&fit.C},{"correlation",&fit.correlation},{"M",&fit.M},
        {"experimental_covariance",&fit.experimental_covariance},{"mc_covariance",&fit.mc_covariance},{"combined_covariance",&fit.combined_covariance}})item.second->Write(item.first);
    fit.reference_asymmetry_ppb.Write("reference_asymmetry_ppb");sums.Write("rate_asymmetry_sums");
    fit.extracted_asymmetry_ppb.Write("extracted_asymmetry_ppb");
    TVectorD detectors(ids.size());for(size_t i=0;i<ids.size();i++)detectors[i]=ids[i];detectors.Write("detector_order");
    const int components=fit.M.GetNrows();
    TH2D uncertainty("uncertainties","Extracted asymmetry uncertainty [ppb];Component;Uncertainty source",components,0,components,3,0,3);
    for(int c=0;c<components;c++){
        uncertainty.GetXaxis()->SetBinLabel(c+1,labels.empty()?component_names[c]:labels.at(c).c_str());
        uncertainty.SetBinContent(c+1,1,std::sqrt(std::max(0.,fit.experimental_covariance(c,c))));
        uncertainty.SetBinContent(c+1,2,std::sqrt(std::max(0.,fit.mc_covariance(c,c))));
        uncertainty.SetBinContent(c+1,3,std::sqrt(std::max(0.,fit.combined_covariance(c,c))));
    }
    uncertainty.GetYaxis()->SetBinLabel(1,"Projected statistical");uncertainty.GetYaxis()->SetBinLabel(2,"Simulation MC");uncertainty.GetYaxis()->SetBinLabel(3,"Combined");uncertainty.Write();
    TObjString(provenance.c_str()).Write("provenance");
}
}


namespace standalone_deconv {
inline std::vector<int> channel_ids(){return {
    110110,110210,110310,110410,110510,110610,110710,110810,110910,111010,111110,111210,111310,111410,
    150110,150210,150310,150410,150510,150610,150710,150810,150910,151010,151110,151210,151310,151410,
    110120,110220,110320,110420,110520,110620,110720,110820,110920,111020,111120,111220,111320,111420,
    150120,150220,150320,150420,150520,150620,150720,150820,150920,151020,151120,151220,151320,151420,
    110130,110230,110330,110430,110530,110630,110730,110830,110930,111030,111130,111230,111330,111430,
    150130,150230,150330,150430,150530,150630,150730,150830,150930,151030,151130,151230,151330,151430,
    110140,110240,110340,110440,110540,110640,110740,110840,110940,111040,111140,111240,111340,111440,
    150140,150240,150340,150440,150540,150640,150740,150840,150940,151040,151140,151240,151340,151440,
    110151,110152,110153,110251,110252,110253,110351,110352,110353,110451,110452,110453,110551,110552,
    110553,110651,110652,110653,110751,110752,110753,110851,110852,110853,110951,110952,110953,111051,
    111052,111053,111151,111152,111153,111251,111252,111253,111351,111352,111353,111451,111452,111453,
    150151,150152,150153,150251,150252,150253,150351,150352,150353,150451,150452,150453,150551,150552,
    150553,150651,150652,150653,150751,150752,150753,150851,150852,150853,150951,150952,150953,151051,
    151052,151053,151151,151152,151153,151251,151252,151253,151351,151352,151353,151451,151452,151453,
    110160,110260,110360,110460,110560,110660,110760,110860,110960,111060,111160,111260,111360,111460,
    150160,150260,150360,150460,150560,150660,150760,150860,150960,151060,151160,151260,151360,151460,
    73001,73002,73003,73004,73005,73006,73007,73008,73009,73010,73011,73012,73013,73014,
    73015,73016,73017,73018,73019,73020,73021,73022,73023,73024,73025,73026,73027,73028
};}
// Geometry identities for the copied thin-quartz coordinate model; no PE tables.
inline std::map<std::string,std::string> thin_geometry(){return {
    {"detector/ThinQuartz/DetectorArray/BFSeg10_R1-BF-Quartz.gdml","5b599c2f891c373d7ba7934f9030e192"},
    {"detector/ThinQuartz/DetectorArray/BFSeg10_R2-BF-Quartz.gdml","15e466615f3ccf0d65c7ec7c19f75b18"},
    {"detector/ThinQuartz/DetectorArray/BFSeg10_R3-BF-Quartz.gdml","9cebde3578f53a508375caacf738f82e"},
    {"detector/ThinQuartz/DetectorArray/BFSeg10_R4-BF-Quartz.gdml","09d2afcd123a5674a103cb3edd39d27f"},
    {"detector/ThinQuartz/DetectorArray/BFSeg10_R5-BF-Quartz_L.gdml","89a32aba20bb5fa9b8e427b514b45363"},
    {"detector/ThinQuartz/DetectorArray/BFSeg10_R5-BF-Quartz_R.gdml","7b9b95b809f93c8acd46e272a721ed10"},
    {"detector/ThinQuartz/DetectorArray/BFSeg10_R5-FF-Quartz_C.gdml","941b39dbcdd050a6ae386d03161a14c6"},
    {"detector/ThinQuartz/DetectorArray/BFSeg10_R6-BF-Quartz.gdml","136b67466dc41bef4307fbe01097334b"},
    {"detector/ThinQuartz/DetectorArray/BFSeg11_R1-BF-Quartz.gdml","2cf20af468d5681cd22df60e52c19dd2"},
    {"detector/ThinQuartz/DetectorArray/BFSeg11_R2-BF-Quartz.gdml","a6bee9e56a9760d04af762280a445207"},
    {"detector/ThinQuartz/DetectorArray/BFSeg11_R3-BF-Quartz.gdml","bebde9c6ef362f0d3f6f8affe415cdca"},
    {"detector/ThinQuartz/DetectorArray/BFSeg11_R4-BF-Quartz.gdml","7f6e4f02caf4a514c39e78f218b94d8c"},
    {"detector/ThinQuartz/DetectorArray/BFSeg11_R5-BF-Quartz_L.gdml","cae5806cc56524c5dcc48de0035e41a3"},
    {"detector/ThinQuartz/DetectorArray/BFSeg11_R5-BF-Quartz_R.gdml","a6a8b8d00670808b0e962b9160498822"},
    {"detector/ThinQuartz/DetectorArray/BFSeg11_R5-FF-Quartz_C.gdml","f69c88191836a3fe1f85c0defc6fe05b"},
    {"detector/ThinQuartz/DetectorArray/BFSeg11_R6-BF-Quartz.gdml","776fb216188866628e81d3ba425bd042"},
    {"detector/ThinQuartz/DetectorArray/BFSeg12_R1-BF-Quartz.gdml","1d2b32bc5025c64cd91d42f52eee4e31"},
    {"detector/ThinQuartz/DetectorArray/BFSeg12_R2-BF-Quartz.gdml","d34adc8a9c556e1b4ce91a79078d5ee9"},
    {"detector/ThinQuartz/DetectorArray/BFSeg12_R3-BF-Quartz.gdml","71b57c6038326307318488d662404c4b"},
    {"detector/ThinQuartz/DetectorArray/BFSeg12_R4-BF-Quartz.gdml","e27cadcbe16d0f5628e0bf34b6c6a34a"},
    {"detector/ThinQuartz/DetectorArray/BFSeg12_R5-BF-Quartz_L.gdml","b084a8cfb17f5f09901048a8180a94a3"},
    {"detector/ThinQuartz/DetectorArray/BFSeg12_R5-BF-Quartz_R.gdml","b0b762b4ddbdbf36c7258cac19551f2e"},
    {"detector/ThinQuartz/DetectorArray/BFSeg12_R5-FF-Quartz_C.gdml","9ec1b256fcfccc66ce6a117433da3409"},
    {"detector/ThinQuartz/DetectorArray/BFSeg12_R6-BF-Quartz.gdml","c29435240f445fa8e29fe62b6d7e2d92"},
    {"detector/ThinQuartz/DetectorArray/BFSeg13_R1-BF-Quartz.gdml","935f14991ce44181c13962c6e7d0acda"},
    {"detector/ThinQuartz/DetectorArray/BFSeg13_R2-BF-Quartz.gdml","02dc264c926481d8ee15c0c0ce6fa265"},
    {"detector/ThinQuartz/DetectorArray/BFSeg13_R3-BF-Quartz.gdml","89f4949e856b2b0467016a9c6f650532"},
    {"detector/ThinQuartz/DetectorArray/BFSeg13_R4-BF-Quartz.gdml","4d45afde18b67af4e021c71e3b12e72c"},
    {"detector/ThinQuartz/DetectorArray/BFSeg13_R5-BF-Quartz_L.gdml","6c1155cfd54adfeca254255cb38108ff"},
    {"detector/ThinQuartz/DetectorArray/BFSeg13_R5-BF-Quartz_R.gdml","8df6d8b3ed874ca1bf404858652d5f7c"},
    {"detector/ThinQuartz/DetectorArray/BFSeg13_R5-FF-Quartz_C.gdml","d6b850f7ea6a04a64a332bd48886006d"},
    {"detector/ThinQuartz/DetectorArray/BFSeg13_R6-BF-Quartz.gdml","b07840c5d366a81a1972909b84e2213d"},
    {"detector/ThinQuartz/DetectorArray/BFSeg14_R1-BF-Quartz.gdml","3d18bf622c1c92e4d1d938a38ef93fa2"},
    {"detector/ThinQuartz/DetectorArray/BFSeg14_R2-BF-Quartz.gdml","b3a66c4ffe7dc2c81abfaa4eb071d420"},
    {"detector/ThinQuartz/DetectorArray/BFSeg14_R3-BF-Quartz.gdml","9cb1f63717dfd50afeb019e8c3331581"},
    {"detector/ThinQuartz/DetectorArray/BFSeg14_R4-BF-Quartz.gdml","811bf5ed5432fb72778e51dc06c6b1a3"},
    {"detector/ThinQuartz/DetectorArray/BFSeg14_R5-BF-Quartz_L.gdml","8ef2759d4ecb664103cae8cf4dae5cd6"},
    {"detector/ThinQuartz/DetectorArray/BFSeg14_R5-BF-Quartz_R.gdml","a2e8269b99c837c51b911f41ed2d64ef"},
    {"detector/ThinQuartz/DetectorArray/BFSeg14_R5-FF-Quartz_C.gdml","adfab7ef6ed16c76ec6aed76e7da4164"},
    {"detector/ThinQuartz/DetectorArray/BFSeg14_R6-BF-Quartz.gdml","e494b0e0e5149525361cce62e418bc14"},
    {"detector/ThinQuartz/DetectorArray/BFSeg1_R1-BF-Quartz.gdml","8acdfe4be4ae6068cfd635b96193e7c1"},
    {"detector/ThinQuartz/DetectorArray/BFSeg1_R2-BF-Quartz.gdml","35d4bed86bc0f63859256afc8fcf4c56"},
    {"detector/ThinQuartz/DetectorArray/BFSeg1_R3-BF-Quartz.gdml","124dc1552c2280198421dedc31979898"},
    {"detector/ThinQuartz/DetectorArray/BFSeg1_R4-BF-Quartz.gdml","cbcf9923ddd10079ac307eccce765ef0"},
    {"detector/ThinQuartz/DetectorArray/BFSeg1_R5-BF-Quartz_L.gdml","c3d43a1c857ee63de12ff81d23d83f03"},
    {"detector/ThinQuartz/DetectorArray/BFSeg1_R5-BF-Quartz_R.gdml","33869f22a2deebded055dda3d34c2789"},
    {"detector/ThinQuartz/DetectorArray/BFSeg1_R5-FF-Quartz_C.gdml","f4756daeee3a1dc27b30244e99e9eb06"},
    {"detector/ThinQuartz/DetectorArray/BFSeg1_R6-BF-Quartz.gdml","9cd8343b5b34336b1e4ee7b2ac0a1f8e"},
    {"detector/ThinQuartz/DetectorArray/BFSeg2_R1-BF-Quartz.gdml","176937e1e2bd0a23b72a39407ea0fec6"},
    {"detector/ThinQuartz/DetectorArray/BFSeg2_R2-BF-Quartz.gdml","cb32fe89e8e2d65c13be93240423ff07"},
    {"detector/ThinQuartz/DetectorArray/BFSeg2_R3-BF-Quartz.gdml","8cad757bffb14217d579d93a8185d0ff"},
    {"detector/ThinQuartz/DetectorArray/BFSeg2_R4-BF-Quartz.gdml","b34b0bead28ba8c4ad20fe463b9f4ba1"},
    {"detector/ThinQuartz/DetectorArray/BFSeg2_R5-BF-Quartz_L.gdml","43b06f465f1350abfe179861c963f509"},
    {"detector/ThinQuartz/DetectorArray/BFSeg2_R5-BF-Quartz_R.gdml","829d8bc0cc05a3d0d2466027b89f46f4"},
    {"detector/ThinQuartz/DetectorArray/BFSeg2_R5-FF-Quartz_C.gdml","9e4275b7ae18387d2b80451fa3e446f5"},
    {"detector/ThinQuartz/DetectorArray/BFSeg2_R6-BF-Quartz.gdml","f7bcb326dab5c62fa1f8983e048f7298"},
    {"detector/ThinQuartz/DetectorArray/BFSeg3_R1-BF-Quartz.gdml","a19ff8070c1daa251e43c7d76023f458"},
    {"detector/ThinQuartz/DetectorArray/BFSeg3_R2-BF-Quartz.gdml","eda774e0a5a9bfb9c2340bccf2454fe0"},
    {"detector/ThinQuartz/DetectorArray/BFSeg3_R3-BF-Quartz.gdml","06e4616f7b62b021e36f055628850474"},
    {"detector/ThinQuartz/DetectorArray/BFSeg3_R4-BF-Quartz.gdml","627227746b203495f39efb70fd3284bb"},
    {"detector/ThinQuartz/DetectorArray/BFSeg3_R5-BF-Quartz_L.gdml","506dfe8bc80fab5af70353f3529a9d53"},
    {"detector/ThinQuartz/DetectorArray/BFSeg3_R5-BF-Quartz_R.gdml","f55f085430e1a5fbe90325159dd69d7c"},
    {"detector/ThinQuartz/DetectorArray/BFSeg3_R5-FF-Quartz_C.gdml","e05523aaa7e1b3a26646ed63fb97c79b"},
    {"detector/ThinQuartz/DetectorArray/BFSeg3_R6-BF-Quartz.gdml","3615ca4dce0824476359a771e01a7b9c"},
    {"detector/ThinQuartz/DetectorArray/BFSeg4_R1-BF-Quartz.gdml","b2405dfd39fc6d96f3414153fe90764f"},
    {"detector/ThinQuartz/DetectorArray/BFSeg4_R2-BF-Quartz.gdml","92df1d13f2ecae116f2b448f8dc14103"},
    {"detector/ThinQuartz/DetectorArray/BFSeg4_R3-BF-Quartz.gdml","0680a6ed1140fa009adb916ebce5e1f4"},
    {"detector/ThinQuartz/DetectorArray/BFSeg4_R4-BF-Quartz.gdml","cde51e38e3474fc46ed7415d7d587b5c"},
    {"detector/ThinQuartz/DetectorArray/BFSeg4_R5-BF-Quartz_L.gdml","efebaf9e20dd8a9b79d3a713795a522b"},
    {"detector/ThinQuartz/DetectorArray/BFSeg4_R5-BF-Quartz_R.gdml","1c13d0e626f9dbd84f262704bb3fca86"},
    {"detector/ThinQuartz/DetectorArray/BFSeg4_R5-FF-Quartz_C.gdml","2dda5e9b7bcf90ec0a0f94c75885ecf0"},
    {"detector/ThinQuartz/DetectorArray/BFSeg4_R6-BF-Quartz.gdml","cf5eb3aef053f4d27053fc07ebd9362f"},
    {"detector/ThinQuartz/DetectorArray/BFSeg5_R1-BF-Quartz.gdml","821e2ad101733c146c9b4bd7a20dd6ea"},
    {"detector/ThinQuartz/DetectorArray/BFSeg5_R2-BF-Quartz.gdml","5f0c5a5bb29ddf3a318181522f7fd281"},
    {"detector/ThinQuartz/DetectorArray/BFSeg5_R3-BF-Quartz.gdml","087266e48d362e46a9085c2220656485"},
    {"detector/ThinQuartz/DetectorArray/BFSeg5_R4-BF-Quartz.gdml","8f4fe5435b35b2f805b409da81792984"},
    {"detector/ThinQuartz/DetectorArray/BFSeg5_R5-BF-Quartz_L.gdml","f095e1429da6067f9e37d49172aa4adf"},
    {"detector/ThinQuartz/DetectorArray/BFSeg5_R5-BF-Quartz_R.gdml","3df2a07970bc1ec5d10372a01a6a040a"},
    {"detector/ThinQuartz/DetectorArray/BFSeg5_R5-FF-Quartz_C.gdml","8b6ac081a2bb174df9a99dd4143f538a"},
    {"detector/ThinQuartz/DetectorArray/BFSeg5_R6-BF-Quartz.gdml","f39a11469d5d90092965c948f2879864"},
    {"detector/ThinQuartz/DetectorArray/BFSeg6_R1-BF-Quartz.gdml","e22b4ca8edf363a75e7049a25e435e8f"},
    {"detector/ThinQuartz/DetectorArray/BFSeg6_R2-BF-Quartz.gdml","f9697d95d7ce5fdece7c88ec4e55fdd2"},
    {"detector/ThinQuartz/DetectorArray/BFSeg6_R3-BF-Quartz.gdml","ea5f3301f65015ded59ee89428b13f00"},
    {"detector/ThinQuartz/DetectorArray/BFSeg6_R4-BF-Quartz.gdml","582fef4e866a8f2d5c8b86177214f63b"},
    {"detector/ThinQuartz/DetectorArray/BFSeg6_R5-BF-Quartz_L.gdml","c25f638538c87bef2de340c80472d2bf"},
    {"detector/ThinQuartz/DetectorArray/BFSeg6_R5-BF-Quartz_R.gdml","a12baaff148f4cf70b08c3a4db8f1526"},
    {"detector/ThinQuartz/DetectorArray/BFSeg6_R5-FF-Quartz_C.gdml","542dac86724704075049643a519d134b"},
    {"detector/ThinQuartz/DetectorArray/BFSeg6_R6-BF-Quartz.gdml","27867f2218a44c0bb3df039fe05f8f6e"},
    {"detector/ThinQuartz/DetectorArray/BFSeg7_R1-BF-Quartz.gdml","253e956f434b5700cff9711bb6a6cbad"},
    {"detector/ThinQuartz/DetectorArray/BFSeg7_R2-BF-Quartz.gdml","aabeead025ef810c419bb6bb133ff0b2"},
    {"detector/ThinQuartz/DetectorArray/BFSeg7_R3-BF-Quartz.gdml","d88288a761e3c2f22a39b562756f8fb1"},
    {"detector/ThinQuartz/DetectorArray/BFSeg7_R4-BF-Quartz.gdml","54e2e9855b43ac171cf0aec2b71add15"},
    {"detector/ThinQuartz/DetectorArray/BFSeg7_R5-BF-Quartz_L.gdml","8bf040fd3aa52cce34c47f0b5d37168e"},
    {"detector/ThinQuartz/DetectorArray/BFSeg7_R5-BF-Quartz_R.gdml","e1486544a0e7f4f781594f6f4f1dffc2"},
    {"detector/ThinQuartz/DetectorArray/BFSeg7_R5-FF-Quartz_C.gdml","d1a54cc8582878af988d7ebd844c364b"},
    {"detector/ThinQuartz/DetectorArray/BFSeg7_R6-BF-Quartz.gdml","e5cc48e46486a3a656f8a849d64af646"},
    {"detector/ThinQuartz/DetectorArray/BFSeg8_R1-BF-Quartz.gdml","69e0ebda64267a4301d3b28174fa7369"},
    {"detector/ThinQuartz/DetectorArray/BFSeg8_R2-BF-Quartz.gdml","a7dfe1ed6e1ed8487ac8425008b939ff"},
    {"detector/ThinQuartz/DetectorArray/BFSeg8_R3-BF-Quartz.gdml","2a90c1d3181703812b45e0bde5233807"},
    {"detector/ThinQuartz/DetectorArray/BFSeg8_R4-BF-Quartz.gdml","ab035a000df0dcbd288431294fcc0e35"},
    {"detector/ThinQuartz/DetectorArray/BFSeg8_R5-BF-Quartz_L.gdml","5c8546c6c0f92876d48bc9edff6df78f"},
    {"detector/ThinQuartz/DetectorArray/BFSeg8_R5-BF-Quartz_R.gdml","24c3e2a94b9454de6e8f478416904ff0"},
    {"detector/ThinQuartz/DetectorArray/BFSeg8_R5-FF-Quartz_C.gdml","42cb819fb2e4a43206fc3e1b35f94e3f"},
    {"detector/ThinQuartz/DetectorArray/BFSeg8_R6-BF-Quartz.gdml","3a4432c5778bcb63fa645b3f62462a17"},
    {"detector/ThinQuartz/DetectorArray/BFSeg9_R1-BF-Quartz.gdml","42e052beea546703ef3e562d294868e3"},
    {"detector/ThinQuartz/DetectorArray/BFSeg9_R2-BF-Quartz.gdml","bc524abc271958e217c31fb19cfdc174"},
    {"detector/ThinQuartz/DetectorArray/BFSeg9_R3-BF-Quartz.gdml","5d32cfd2e90a6ebe50db8130b5d5f7de"},
    {"detector/ThinQuartz/DetectorArray/BFSeg9_R4-BF-Quartz.gdml","ae5e6978ac4ff57a0d73d4774df11a85"},
    {"detector/ThinQuartz/DetectorArray/BFSeg9_R5-BF-Quartz_L.gdml","1975a6292197d1310f7de8f9c2e0a8f8"},
    {"detector/ThinQuartz/DetectorArray/BFSeg9_R5-BF-Quartz_R.gdml","26feabc6d6545d656ed696dffaaaf1e1"},
    {"detector/ThinQuartz/DetectorArray/BFSeg9_R5-FF-Quartz_C.gdml","a6a3624a9e1e005128fd36e4d76d9e8a"},
    {"detector/ThinQuartz/DetectorArray/BFSeg9_R6-BF-Quartz.gdml","5ba452603d4f8d3b7866184d27177ba0"},
    {"detector/ThinQuartz/DetectorArray/DetectorArray.gdml","d59f47a32196a28491654b34db0b5742"},
    {"detector/ThinQuartz/DetectorArray/FFSeg10_R1-FF-Quartz.gdml","25d39ee4cfd2eb55dc54acbfa4d0c7af"},
    {"detector/ThinQuartz/DetectorArray/FFSeg10_R2-FF-Quartz.gdml","791fc805d29d5db1acf4de248d17a44a"},
    {"detector/ThinQuartz/DetectorArray/FFSeg10_R3-FF-Quartz.gdml","6a9bec9ed77ab11884cfaf4105eb1838"},
    {"detector/ThinQuartz/DetectorArray/FFSeg10_R4-FF-Quartz.gdml","9c3afdae79b25f408435d24a7994bc61"},
    {"detector/ThinQuartz/DetectorArray/FFSeg10_R5-BF-Quartz_C.gdml","3808cec3261229be5f7aa25aca26530d"},
    {"detector/ThinQuartz/DetectorArray/FFSeg10_R5-FF-Quartz_L.gdml","57e0d595729248a9f1ee4132c19c8d98"},
    {"detector/ThinQuartz/DetectorArray/FFSeg10_R5-FF-Quartz_R.gdml","e7ce54134023e65dd61678b2ade71e24"},
    {"detector/ThinQuartz/DetectorArray/FFSeg10_R6-FF-Quartz.gdml","f6894077ec086cfa56c9660a77f9361e"},
    {"detector/ThinQuartz/DetectorArray/FFSeg11_R1-FF-Quartz.gdml","ed98f670602f6cf7393656e9c8c3ffac"},
    {"detector/ThinQuartz/DetectorArray/FFSeg11_R2-FF-Quartz.gdml","fd7874c098c50dd2012b25a95852451b"},
    {"detector/ThinQuartz/DetectorArray/FFSeg11_R3-FF-Quartz.gdml","cf3cbb396bcd41d565878b8621349b29"},
    {"detector/ThinQuartz/DetectorArray/FFSeg11_R4-FF-Quartz.gdml","89a1316dce946580935fc5026e816d0a"},
    {"detector/ThinQuartz/DetectorArray/FFSeg11_R5-BF-Quartz_C.gdml","ef6df577f693cddbf456a2f85a29f696"},
    {"detector/ThinQuartz/DetectorArray/FFSeg11_R5-FF-Quartz_L.gdml","7d143bc722603fd385bf754ff843bc96"},
    {"detector/ThinQuartz/DetectorArray/FFSeg11_R5-FF-Quartz_R.gdml","2e9662eb074d396b3a8ff67d2b737971"},
    {"detector/ThinQuartz/DetectorArray/FFSeg11_R6-FF-Quartz.gdml","08e848850c371705f71e93ac24321543"},
    {"detector/ThinQuartz/DetectorArray/FFSeg12_R1-FF-Quartz.gdml","7a624fe2408e6e7cb13ef85f201848c9"},
    {"detector/ThinQuartz/DetectorArray/FFSeg12_R2-FF-Quartz.gdml","f2627cc1a971be141802104f569b2b13"},
    {"detector/ThinQuartz/DetectorArray/FFSeg12_R3-FF-Quartz.gdml","47a03b0a80ff30d8f994424e77e5b62a"},
    {"detector/ThinQuartz/DetectorArray/FFSeg12_R4-FF-Quartz.gdml","a1b841527f637ff42f7544a3eb19fa7b"},
    {"detector/ThinQuartz/DetectorArray/FFSeg12_R5-BF-Quartz_C.gdml","c4b04ec8dbb0ae3915f1d59ce69994d8"},
    {"detector/ThinQuartz/DetectorArray/FFSeg12_R5-FF-Quartz_L.gdml","9937864a26bfa1c1a242d4942dc8f6cd"},
    {"detector/ThinQuartz/DetectorArray/FFSeg12_R5-FF-Quartz_R.gdml","4515702a843aaf1f0584968671de0288"},
    {"detector/ThinQuartz/DetectorArray/FFSeg12_R6-FF-Quartz.gdml","0e58e56ecc432cd9984d6e97cf768516"},
    {"detector/ThinQuartz/DetectorArray/FFSeg13_R1-FF-Quartz.gdml","0daf5e74c47faccb2f7861ddcc89c2b0"},
    {"detector/ThinQuartz/DetectorArray/FFSeg13_R2-FF-Quartz.gdml","aa046d257ba2e2011e2a63c79be26002"},
    {"detector/ThinQuartz/DetectorArray/FFSeg13_R3-FF-Quartz.gdml","82f7f15d153c28f785579ee075dc61b1"},
    {"detector/ThinQuartz/DetectorArray/FFSeg13_R4-FF-Quartz.gdml","8cd0edf4ca4df9743555039e857fe16f"},
    {"detector/ThinQuartz/DetectorArray/FFSeg13_R5-BF-Quartz_C.gdml","fcfe1a06eb7e9798b01606967abe3ebe"},
    {"detector/ThinQuartz/DetectorArray/FFSeg13_R5-FF-Quartz_L.gdml","95efea5e08c34483b4b4d4caad686cc7"},
    {"detector/ThinQuartz/DetectorArray/FFSeg13_R5-FF-Quartz_R.gdml","dced241c350143bd3366037923636979"},
    {"detector/ThinQuartz/DetectorArray/FFSeg13_R6-FF-Quartz.gdml","2b28593316c7feea778098dec52319b4"},
    {"detector/ThinQuartz/DetectorArray/FFSeg14_R1-FF-Quartz.gdml","1916ca986834636246bd999997674bd9"},
    {"detector/ThinQuartz/DetectorArray/FFSeg14_R2-FF-Quartz.gdml","be15b3fac272690fbe9ae64500734d25"},
    {"detector/ThinQuartz/DetectorArray/FFSeg14_R3-FF-Quartz.gdml","1ba375efc19fa6f3f22c4bbe665a888b"},
    {"detector/ThinQuartz/DetectorArray/FFSeg14_R4-FF-Quartz.gdml","9e3d5b9b294ad0067fe5789b81b4e777"},
    {"detector/ThinQuartz/DetectorArray/FFSeg14_R5-BF-Quartz_C.gdml","f3a79f05f2e3ce244ce550275e616fa4"},
    {"detector/ThinQuartz/DetectorArray/FFSeg14_R5-FF-Quartz_L.gdml","79a8c02285becdf157e0d8c7f1e6f635"},
    {"detector/ThinQuartz/DetectorArray/FFSeg14_R5-FF-Quartz_R.gdml","2a3bd71f584eee01b1028bef802092fb"},
    {"detector/ThinQuartz/DetectorArray/FFSeg14_R6-FF-Quartz.gdml","e1f0e36d745fc835840a8e10d3844737"},
    {"detector/ThinQuartz/DetectorArray/FFSeg1_R1-FF-Quartz.gdml","fbbb3e410ed75b77d18dee62eda6c734"},
    {"detector/ThinQuartz/DetectorArray/FFSeg1_R2-FF-Quartz.gdml","e87b8d60e77f776f227f545a415f55b2"},
    {"detector/ThinQuartz/DetectorArray/FFSeg1_R3-FF-Quartz.gdml","eefef669599f64f27da6732a88fa921d"},
    {"detector/ThinQuartz/DetectorArray/FFSeg1_R4-FF-Quartz.gdml","c7bca9151f445f52ca755f42821b8566"},
    {"detector/ThinQuartz/DetectorArray/FFSeg1_R5-BF-Quartz_C.gdml","455ac87a8d93acb21a1db581a2d88aee"},
    {"detector/ThinQuartz/DetectorArray/FFSeg1_R5-FF-Quartz_L.gdml","f3482e10bca1aacf439579be406bc4be"},
    {"detector/ThinQuartz/DetectorArray/FFSeg1_R5-FF-Quartz_R.gdml","6941c9d6c14107586bed056f6d99e279"},
    {"detector/ThinQuartz/DetectorArray/FFSeg1_R6-FF-Quartz.gdml","dc430f78285e642f58166ff9afaba952"},
    {"detector/ThinQuartz/DetectorArray/FFSeg2_R1-FF-Quartz.gdml","70a65c8290b9682c4e0b17074ce483d2"},
    {"detector/ThinQuartz/DetectorArray/FFSeg2_R2-FF-Quartz.gdml","35eb671cb273fc22b120a3efbda95f61"},
    {"detector/ThinQuartz/DetectorArray/FFSeg2_R3-FF-Quartz.gdml","ed81ba4d7c9eb22bacd7a850668a1a39"},
    {"detector/ThinQuartz/DetectorArray/FFSeg2_R4-FF-Quartz.gdml","35f099c4893f3c20d10fa2d6109247a2"},
    {"detector/ThinQuartz/DetectorArray/FFSeg2_R5-BF-Quartz_C.gdml","843c8ca642626f5499de00aea858a53f"},
    {"detector/ThinQuartz/DetectorArray/FFSeg2_R5-FF-Quartz_L.gdml","0520b0a8f7a85d0632b6619eb62a364b"},
    {"detector/ThinQuartz/DetectorArray/FFSeg2_R5-FF-Quartz_R.gdml","74d8a62797b126a967cc1aac8bc93d3a"},
    {"detector/ThinQuartz/DetectorArray/FFSeg2_R6-FF-Quartz.gdml","23d47d80d962edf449df52aef8075b59"},
    {"detector/ThinQuartz/DetectorArray/FFSeg3_R1-FF-Quartz.gdml","4010326023cc74e1c9cdb86b34613de9"},
    {"detector/ThinQuartz/DetectorArray/FFSeg3_R2-FF-Quartz.gdml","0d48b75c653a22fc2332eac96ad58b02"},
    {"detector/ThinQuartz/DetectorArray/FFSeg3_R3-FF-Quartz.gdml","4e16dd86a6867b7aca83a44ceb6c0035"},
    {"detector/ThinQuartz/DetectorArray/FFSeg3_R4-FF-Quartz.gdml","8024f9a85ab86e8c5ff7cfbb9291f5a5"},
    {"detector/ThinQuartz/DetectorArray/FFSeg3_R5-BF-Quartz_C.gdml","b1692af44e5cac71f958d5152ad8af30"},
    {"detector/ThinQuartz/DetectorArray/FFSeg3_R5-FF-Quartz_L.gdml","29167674bf3115683155c67ffd96966c"},
    {"detector/ThinQuartz/DetectorArray/FFSeg3_R5-FF-Quartz_R.gdml","bb3389a8a65d71a1a5d80e981b74b94b"},
    {"detector/ThinQuartz/DetectorArray/FFSeg3_R6-FF-Quartz.gdml","df7babdc97e184734ba5c7acd205ec95"},
    {"detector/ThinQuartz/DetectorArray/FFSeg4_R1-FF-Quartz.gdml","1aa1f0aa73dd0a7fc8ae9ec7ff4ed71e"},
    {"detector/ThinQuartz/DetectorArray/FFSeg4_R2-FF-Quartz.gdml","a2c9524ed36cef5b36879f1ce6f3ead9"},
    {"detector/ThinQuartz/DetectorArray/FFSeg4_R3-FF-Quartz.gdml","37ae24f3da6ffeb40745e0d96753233d"},
    {"detector/ThinQuartz/DetectorArray/FFSeg4_R4-FF-Quartz.gdml","50974b47082a948cb5fc1e87fbf11f93"},
    {"detector/ThinQuartz/DetectorArray/FFSeg4_R5-BF-Quartz_C.gdml","8d5bb9a6c5d3aff89a75641aa2482bc9"},
    {"detector/ThinQuartz/DetectorArray/FFSeg4_R5-FF-Quartz_L.gdml","0539498112fdb5037a2a307dccfdbf9a"},
    {"detector/ThinQuartz/DetectorArray/FFSeg4_R5-FF-Quartz_R.gdml","016c36079c236c1a222ac403ee7e1a53"},
    {"detector/ThinQuartz/DetectorArray/FFSeg4_R6-FF-Quartz.gdml","bb46dce81ac0c75ad04be35099f46052"},
    {"detector/ThinQuartz/DetectorArray/FFSeg5_R1-FF-Quartz.gdml","fed13afabe6faf962537d966d73a5cf9"},
    {"detector/ThinQuartz/DetectorArray/FFSeg5_R2-FF-Quartz.gdml","b17866d2323703031b7e1e1c4254e4ba"},
    {"detector/ThinQuartz/DetectorArray/FFSeg5_R3-FF-Quartz.gdml","3611bc6a3dbcb131aa149e02dcccc8e0"},
    {"detector/ThinQuartz/DetectorArray/FFSeg5_R4-FF-Quartz.gdml","6795305814cb1de800e5a5bb9c1832a0"},
    {"detector/ThinQuartz/DetectorArray/FFSeg5_R5-BF-Quartz_C.gdml","14f6a03c014fcd8ad3094d82fcb1e4eb"},
    {"detector/ThinQuartz/DetectorArray/FFSeg5_R5-FF-Quartz_L.gdml","5140afc7d31c6e01d1fcb87a6a0f9d57"},
    {"detector/ThinQuartz/DetectorArray/FFSeg5_R5-FF-Quartz_R.gdml","09d9261e068d773010eaf54b9916d79c"},
    {"detector/ThinQuartz/DetectorArray/FFSeg5_R6-FF-Quartz.gdml","a125a902d8573739379815926712054a"},
    {"detector/ThinQuartz/DetectorArray/FFSeg6_R1-FF-Quartz.gdml","9c7fa441a6f42754203f64b7977d55fc"},
    {"detector/ThinQuartz/DetectorArray/FFSeg6_R2-FF-Quartz.gdml","7a42b647066294ac9ef60020064d3088"},
    {"detector/ThinQuartz/DetectorArray/FFSeg6_R3-FF-Quartz.gdml","5bfa73e95a59f932078fbf90c2b38dc1"},
    {"detector/ThinQuartz/DetectorArray/FFSeg6_R4-FF-Quartz.gdml","cc8fca055c56361718ac885911c3cc39"},
    {"detector/ThinQuartz/DetectorArray/FFSeg6_R5-BF-Quartz_C.gdml","f122ec5a49dfc89f500e7931b403ac79"},
    {"detector/ThinQuartz/DetectorArray/FFSeg6_R5-FF-Quartz_L.gdml","4b776188baa315ec82510c7f864487b9"},
    {"detector/ThinQuartz/DetectorArray/FFSeg6_R5-FF-Quartz_R.gdml","fe01de88a5ccc4ea489cac3148b8714c"},
    {"detector/ThinQuartz/DetectorArray/FFSeg6_R6-FF-Quartz.gdml","c0f66bfec385bcc10abcebabb20bb7cd"},
    {"detector/ThinQuartz/DetectorArray/FFSeg7_R1-FF-Quartz.gdml","f1ab12e4f43a23a842af434bfa10fed2"},
    {"detector/ThinQuartz/DetectorArray/FFSeg7_R2-FF-Quartz.gdml","fd50f53bce50436f13ba9deb81386bba"},
    {"detector/ThinQuartz/DetectorArray/FFSeg7_R3-FF-Quartz.gdml","4826c0f07410e91d480876deb7004129"},
    {"detector/ThinQuartz/DetectorArray/FFSeg7_R4-FF-Quartz.gdml","3529b389ceaa40791c46994b8cbd5b8b"},
    {"detector/ThinQuartz/DetectorArray/FFSeg7_R5-BF-Quartz_C.gdml","4e0491674dde10bb1677516cac733415"},
    {"detector/ThinQuartz/DetectorArray/FFSeg7_R5-FF-Quartz_L.gdml","93e658434f00a93df49701a510ce64d1"},
    {"detector/ThinQuartz/DetectorArray/FFSeg7_R5-FF-Quartz_R.gdml","7e24cf0ec7e64b8f8e9ab85fbc0ff25c"},
    {"detector/ThinQuartz/DetectorArray/FFSeg7_R6-FF-Quartz.gdml","2269c5cb96b5ac3b27b1cc801d899e4a"},
    {"detector/ThinQuartz/DetectorArray/FFSeg8_R1-FF-Quartz.gdml","66d49da96a9a214a4bfb0008afe56efa"},
    {"detector/ThinQuartz/DetectorArray/FFSeg8_R2-FF-Quartz.gdml","8f775c6850c6d224bb7bd3f93a3704ce"},
    {"detector/ThinQuartz/DetectorArray/FFSeg8_R3-FF-Quartz.gdml","e05e69553fcffee30f2b66f80809b510"},
    {"detector/ThinQuartz/DetectorArray/FFSeg8_R4-FF-Quartz.gdml","06bcc30e0347a8aac1bc1c50e6452573"},
    {"detector/ThinQuartz/DetectorArray/FFSeg8_R5-BF-Quartz_C.gdml","0279222d80352ebe911b920a6cbe8c6d"},
    {"detector/ThinQuartz/DetectorArray/FFSeg8_R5-FF-Quartz_L.gdml","e357148d49588ce15fa597514b63dae5"},
    {"detector/ThinQuartz/DetectorArray/FFSeg8_R5-FF-Quartz_R.gdml","77453569b27a10702d5016b111e2f582"},
    {"detector/ThinQuartz/DetectorArray/FFSeg8_R6-FF-Quartz.gdml","4bc9c37d42f6785e45ae1edc80cf8f60"},
    {"detector/ThinQuartz/DetectorArray/FFSeg9_R1-FF-Quartz.gdml","54aeb97754f4ce57c484a92cfba2c6d2"},
    {"detector/ThinQuartz/DetectorArray/FFSeg9_R2-FF-Quartz.gdml","89e8fbafa50fb6250efc1e61ff3d0ea6"},
    {"detector/ThinQuartz/DetectorArray/FFSeg9_R3-FF-Quartz.gdml","5d5a901e385d12b360e47815172d8757"},
    {"detector/ThinQuartz/DetectorArray/FFSeg9_R4-FF-Quartz.gdml","dafaeafe13dd089b8e89fb7533d814aa"},
    {"detector/ThinQuartz/DetectorArray/FFSeg9_R5-BF-Quartz_C.gdml","b8f3dafefe3ad4b46342cfe251b1155b"},
    {"detector/ThinQuartz/DetectorArray/FFSeg9_R5-FF-Quartz_L.gdml","9b867f5567fb340924314b9572abfd97"},
    {"detector/ThinQuartz/DetectorArray/FFSeg9_R5-FF-Quartz_R.gdml","0735b30e0e67b0efb8f41f46fdf06a03"},
    {"detector/ThinQuartz/DetectorArray/FFSeg9_R6-FF-Quartz.gdml","d6766598a700acc2530a1e33e47fec65"},
};}

inline std::vector<fs::path> input_files(const std::vector<std::string>& patterns){
    std::vector<fs::path> files;std::set<fs::path> unique;
    for(const auto& pattern:patterns){
        if(pattern.empty())throw std::runtime_error("Empty input file pattern");
        TChain chain("T");
        if(!chain.Add(pattern.c_str()))throw std::runtime_error("No ROOT files match: "+pattern);
        auto* list=chain.GetListOfFiles();
        if(!list||!list->GetEntries())throw std::runtime_error("No ROOT files match: "+pattern);
        for(int i=0;i<list->GetEntries();i++){
            auto* element=dynamic_cast<TChainElement*>(list->At(i));
            if(!element)throw std::runtime_error("Invalid ROOT file list");
            const auto path=fs::absolute(element->GetTitle()).lexically_normal();
            if(!unique.insert(fs::canonical(path)).second)throw std::runtime_error("Duplicate input ROOT file: "+path.string());
            files.push_back(path);
        }
    }
    if(files.empty())throw std::runtime_error("Set fnms to your input ROOT files");
    return files;
}
struct Tables {fs::path main,showermax;};
inline Tables find_tables(const std::vector<fs::path>& inputs,const char* requested,bool needMain=true,bool needSM=true){
    std::vector<fs::path> roots;
    if(requested&&*requested)roots.emplace_back(requested);
    else if(const char* env=gSystem->Getenv("MOLLER_PE_TABLES")){
        if(*env)roots.emplace_back(env);
    }
    if(roots.empty()){
        const fs::path macro=fs::absolute(__FILE__).parent_path();
        roots={fs::current_path(),macro,macro.parent_path(),
               fs::current_path()/"_local_work/experiments/detector_launcher/data",
               macro.parent_path()/"experiments/detector_launcher/data"};
        for(const auto& input:inputs)roots.push_back(input.parent_path());
    }
    const char* species[]={"e-","e+","gamma","mu-","mu+","pi-","pi+","neutron"};
    const char* relative[]={"","data","tables","table","main_pe","response",
        "data/main_pe","data/response","tables/main_pe","tables/response",
        "response/main_pe","remoll_deconvolution/response"};
    auto main_complete=[](const fs::path& folder){
        for(int r=1;r<=6;r++)for(const char* suffix:{"_NewLookupTable.csv","_EnergyDepTable.csv"})
            if(!fs::is_regular_file(folder/("r"+std::to_string(r)+suffix)))return false;
        return true;
    };
    auto sm_complete=[&](const fs::path& folder){
        for(const char* name:species)if(!fs::is_regular_file(folder/(std::string("fit_param_xy_")+name+"_ifarm.csv")))return false;
        return true;
    };
    Tables found;std::set<fs::path> checked;
    for(const auto& root:roots)for(const char* rel:relative){
        auto folder=(fs::absolute(root)/rel).lexically_normal();
        if(!checked.insert(folder).second)continue;
        if(found.main.empty()&&main_complete(folder))found.main=fs::canonical(folder);
        if(found.showermax.empty()&&sm_complete(folder))found.showermax=fs::canonical(folder);
        if((!needMain||!found.main.empty())&&(!needSM||!found.showermax.empty()))return found;
    }
    std::string missing=needMain&&found.main.empty()?"main-quartz":"";
    if(needSM&&found.showermax.empty())missing+=(missing.empty()?"":" and ")+std::string("ShowerMax");
    std::ostringstream message;
    message<<"PE enabled, but "<<missing<<" tables were not found. Pass tableDirectory or set MOLLER_PE_TABLES. Searched:";
    for(const auto& folder:checked)message<<"\n  "<<folder.string();
    throw std::runtime_error(message.str());
}

// Develop's report layout, with the correlated estimator and separate MC term.
inline void print_report(const Fit& fit,const TMatrixD& sums,const std::vector<int>& ids,
                         double referenceNoise,double days,double polarization,bool usePE,bool useShowerMaxPE=true){
    const char* names[]={"moller","epElastic","epInelasticW1","epInelasticW2","epInelasticW3"};
    const int n=ids.size(),d=n+1;const double seconds=days*86400;
    const auto oldFlags=std::cout.flags();const auto oldPrecision=std::cout.precision();
    std::cout<<std::setprecision(7);
    std::cout<<"\nR S";
    for(int c=0;c<5;c++)std::cout<<"\tA"<<names[c]<<"\tf"<<names[c];
    std::cout<<"\tAm[ppb]\td(Am)[ppb]\trTot["<<(usePE?"PE/s":useShowerMaxPE?"main Hz; SM PE/s":"Hz")<<"]\n";
    for(int ring=0;ring<7;ring++)for(int sector=0;sector<3;sector++){
        std::vector<int> rows;
        for(int i=0;i<n;i++){
            const int det=ids[i],suffix=det%100;
            const int r=det>=73001&&det<=73028?6:suffix>=51&&suffix<=53?4:suffix/10-1;
            if(r==ring&&standalone_detector_selection::region(det)==sector)rows.push_back(i);
        }
        double rate[5]={0},weighted[5]={0},total=0,mixed=0,noise=0;
        for(int c=0;c<5;c++)for(int i:rows){rate[c]+=sums(c,i);weighted[c]+=sums(c,d+i);}
        for(int c=0;c<5;c++){total+=rate[c];mixed+=weighted[c];}
        for(int i:rows)for(int j:rows)noise+=fit.B(i,j);
        std::cout<<ring<<' '<<sector;
        for(int c=0;c<5;c++)std::cout<<'\t'<<(rate[c]>0?weighted[c]/rate[c]:std::numeric_limits<double>::quiet_NaN())<<'\t'<<(total>0?rate[c]/total:std::numeric_limits<double>::quiet_NaN());
        std::cout<<'\t'<<(total>0?mixed/total:std::numeric_limits<double>::quiet_NaN())
                 <<'\t'<<(total>0?1e9*std::sqrt(noise)/(polarization*std::sqrt(seconds)*total):std::numeric_limits<double>::quiet_NaN())<<'\t'<<total<<'\n';
    }
    std::cout<<"R=0..5: Rings 1..6; R=6: ShowerMax. S=0 closed, 1 transition, 2 open. Physical detector membership.\n";
    auto matrix=[](const TMatrixD& value){for(int i=0;i<value.GetNrows();i++){std::cout<<'\t';for(int j=0;j<value.GetNcols();j++)std::cout<<std::setw(15)<<value(i,j);std::cout<<'\n';}};
    std::cout<<"\nBefore Inversion\n\tDet: "<<fit.M.Determinant()<<"\n\tM matrix [ppb^-2]\n";matrix(fit.M);
    std::cout<<"\nAfter Inversion\n\tDet: "<<fit.experimental_covariance.Determinant()<<"\n\tM^-1 matrix [ppb^2]\n";matrix(fit.experimental_covariance);
    std::cout<<"\n\t\t\tAsymmetry extractor (simulation closure)\n";
    for(int c=0;c<5;c++)std::cout<<"Asymmetry "<<names[c]<<'\t'<<fit.extracted_asymmetry_ppb[c]<<'\n';
    std::cout<<"\n\t\t\toverall\nName\tAsymmetry[ppb]\tuncert[ppb]\trelative uncert[%]\tMC[ppb]\tcombined[ppb]\tcombined relative[%]\n";
    for(int c=0;c<5;c++){
        const double stat=std::sqrt(fit.experimental_covariance(c,c)),mc=std::sqrt(std::max(0.,fit.mc_covariance(c,c))),combined=std::sqrt(fit.combined_covariance(c,c));
        std::cout<<names[c]<<'\t'<<fit.reference_asymmetry_ppb[c]<<'\t'<<stat<<'\t'<<100*stat/std::abs(fit.reference_asymmetry_ppb[c])<<'\t'<<mc<<'\t'<<combined<<'\t'<<100*combined/std::abs(fit.reference_asymmetry_ppb[c])<<'\n';
    }
    const double selected=1e9*std::sqrt(referenceNoise)/(polarization*std::sqrt(seconds)*sums(0,n));
    std::cout<<"\n\t\t\tin the selected bin (Ring 5 open, rate reference)\nstat moller\t"<<selected<<" ppb\t"<<100*selected/std::abs(fit.reference_asymmetry_ppb[0])<<"%\n";
    std::cout.flags(oldFlags);std::cout.precision(oldPrecision);
}

inline void run(const std::vector<std::string>& patterns,bool usePE,const char* tableDirectory,
                const char* output,double beamDays,double polarization,double current_uA,double minimumK,bool useShowerMaxPE){
    if(!output||!*output)throw std::runtime_error("Choose an output ROOT filename");
    if(fs::exists(output))throw std::runtime_error("Output already exists; choose another filename");
    if(!(std::isfinite(minimumK)&&minimumK>=0&&std::isfinite(beamDays)&&beamDays>0&&
         std::isfinite(polarization)&&polarization>0&&polarization<=1&&std::isfinite(current_uA)&&current_uA>0))
        throw std::runtime_error("Invalid rate cut, exposure, polarization or current");
    if(usePE&&!useShowerMaxPE)throw std::runtime_error("Main PE mode requires ShowerMax PE to match the existing analysis");
    const auto paths=input_files(patterns);
    const auto ids=channel_ids();
    if(usePE&&std::abs(minimumK-1)>1e-12)throw std::runtime_error("Existing main PE mode requires minimumK=1 MeV");
    const auto thin=usePE?thin_geometry():std::map<std::string,std::string>{};
    std::vector<Input> files;std::set<std::pair<int,long long>> seeds;long long totals[3]={0,0,0};
    for(const auto& path:paths){
        Input input=inspect(path,ids,usePE,thin);
        if(!seeds.insert({input.process,input.seed}).second)throw std::runtime_error("Repeated random seed in one interaction");
        const double energy=beam_energy_MeV(input);
        if(!files.empty()&&std::abs(energy-beam_energy_MeV(files.front()))>1e-10*energy)
            throw std::runtime_error("Physical beam energy differs between inputs");
        if(allowMixedSources&&!trustedMixedTransportMD5.count(input.transport))
            throw std::runtime_error("Mixed-source setup is not in the reviewed cohort. Confirm source, geometry, target, physics and saved commands before adding its transport MD5 to trustedMixedTransportMD5: "+input.path.string()+" MD5 "+input.transport);
        if(!files.empty()&&!allowMixedSources){
            check_transport(files[0],input);
            if(input.commit!=files[0].commit||input.geometry!=files[0].geometry)
                throw std::runtime_error("Simulation source or geometry differs between inputs");
            for(const char* key:{"/remoll/beamene","/remoll/target/mother","/remoll/target/volume","/remoll/physlist/set"})
                if(input.settings[key]!=files[0].settings[key])throw std::runtime_error(std::string("Incompatible setting: ")+key);
        }
        for(const auto& previous:files)if(previous.process==input.process){
            check_generation_settings(previous,input);
            if(!allowMixedSources)check_transport(previous,input);
        }
        totals[input.process]+=input.events;files.push_back(std::move(input));
    }
    for(int p=0;p<3;p++)if(!totals[p])throw std::runtime_error("Provide raw Moller, ep elastic and ep inelastic ROOT files");
    Tables tables;standalone_portable::Fits sm;std::unique_ptr<MainModel> main;
    if(usePE||useShowerMaxPE){
        tables=find_tables(paths,tableDirectory,usePE,useShowerMaxPE);
        if(usePE){std::cout<<"Main PE tables: "<<tables.main<<"\n";main=std::make_unique<MainModel>(tables.main);}
        if(useShowerMaxPE){std::cout<<"ShowerMax PE tables: "<<tables.showermax<<"\n";sm=standalone_portable::load(tables.showermax.string());}
    }
    const int n=ids.size(),d=n+1,dim=2*d,components=5;
    TMatrixD sums(components,dim),mc(components*dim,components*dim),b(n,n);
    std::vector<double> coverage(n);
    std::ostringstream provenance;provenance<<std::setprecision(17);
    provenance<<"Standalone adaptation of gdml-geometry-viewer 5077001171756458963314261943a0dddd9bfa40\nROOT "<<gROOT->GetVersion()<<"\n";
    provenance<<"Mode: "<<(usePE?"PE (main quartz + ShowerMax)":useShowerMaxPE?"rate (main rate + ShowerMax PE)":"pure rate (no PE tables)")<<"\n";
    provenance<<"Main charged K >= "<<minimumK<<" MeV; Ring 5 open reference is rate weighted\n";
    provenance<<(usePE?"Main PE: forward e+/- K > 1.1 MeV, total-energy lookup, polar capped at 4 degrees; ShowerMax forward supported PE\n":
                      useShowerMaxPE?"Rate: main both directions; ShowerMax forward supported expected PE, matching saved rate mode\n":
                      "Pure rate: main both directions; ShowerMax forward charged K >= minimumK; no PE tables used; distinct from saved rate mode\n");
    provenance<<"Recorded entries and re-entries, with event-level sums before moments\n";
    provenance<<"Forecast "<<beamDays<<" days; polarization "<<polarization<<"; current "<<current_uA<<" uA\n";
    provenance<<"Template MC propagated through H at the simulated reference, with C held fixed\n";
    provenance<<"Additional optical/electronics/beam noise and model systematics excluded\n";
    provenance<<"Macro MD5 "<<digest(fs::absolute(__FILE__))<<"\n";
    for(const auto& folder:{usePE?tables.main:fs::path{},useShowerMaxPE?tables.showermax:fs::path{}})
        if(!folder.empty())for(const auto& p:fs::directory_iterator(folder))if(p.path().extension()==".csv")
            provenance<<"PE table "<<p.path().string()<<" MD5 "<<digest(p.path())<<"\n";
    provenance<<"Mixed-source cohort explicitly enabled: "<<(allowMixedSources?"yes":"no")<<"\n";
    long long histories=0,candidates=0,unmatchedMain=0,mainCandidates=0,unsupportedSM=0,smCandidates=0;
    double referenceNoise=0;size_t filesDone=0;
    for(const auto& input:files){
        std::cout<<"Reading ["<<++filesDone<<'/'<<files.size()<<"] "<<input.path.filename()<<" ("<<input.events<<" histories)"<<std::endl;
        const auto moments=scan(input,ids,sm,main.get(),minimumK,useShowerMaxPE);
        const int nc=input.process==2?3:1,first=input.process==2?2:input.process;
        const double factor=double(input.events)/totals[input.process]*current_uA/input.current;
        for(int c=0;c<nc;c++)for(int i=0;i<dim;i++)sums(first+c,i)+=factor*moments.y[c*dim+i];
        moments.add_covariance(mc,first,dim,factor*factor);
        for(int i=0;i<n;i++){
            coverage[i]+=moments.coverage[i];
            for(int j=0;j<n;j++)b(i,j)+=factor*moments.noise[i*n+j];
        }
        histories+=input.events;candidates+=input.generated;referenceNoise+=factor*moments.reference_noise;
        unmatchedMain+=moments.unmatched_main;mainCandidates+=moments.main_pe_candidates;
        unsupportedSM+=moments.unsupported_sm;smCandidates+=moments.sm_pe_candidates;
        provenance<<"Input "<<input.path.string()<<"\nMD5 "<<digest(input.path)<<"\nTransport "<<input.transport;
        provenance<<"\nTransported "<<input.events<<"; generated/candidate counter "<<input.generated<<"; seed "<<input.seed<<"; current "<<input.current;
        provenance<<" uA; factor "<<factor<<"; unsupported SM "<<moments.unsupported_sm<<"/"<<moments.sm_pe_candidates
                  <<"; unmatched main PE "<<moments.unmatched_main<<"/"<<moments.main_pe_candidates<<"\nCommit "<<input.commit;
        provenance<<"\nMacro:\n"<<input.macro<<"\nGeometry fingerprints:\n"<<input.geometry;
        if(!input.samplingPlan.empty())provenance<<"Sampling plan "<<input.samplingPlan.string()<<" MD5 "<<digest(input.samplingPlan)<<"\n";
    }
    for(int i=0;i<n;i++)if(!coverage[i])throw std::runtime_error("No recorded hits in detector "+std::to_string(ids[i]));
    auto fit=solve(sums,mc,b,beamDays,polarization);
    TFile result(output,"CREATE");if(result.IsZombie())throw std::runtime_error("Cannot create output ROOT file");
    write_fit(result,fit,ids,sums,provenance.str());
    mc.Write("rate_asymmetry_mc_covariance");
    TVectorD referenceMoment(1);referenceMoment[0]=referenceNoise;referenceMoment.Write("reference_second_moment");
    TVectorD historyCounts(2);historyCounts[0]=histories;historyCounts[1]=candidates;historyCounts.Write("history_counts");result.Close();
    std::cout<<"\n"<<histories<<" histories; "<<n<<" channels; "<<(usePE?"PE":"rate")<<" weighting\n";
    std::cout<<"PE response diagnostics: unmatched main "<<unmatchedMain<<"/"<<mainCandidates
             <<"; unsupported ShowerMax "<<unsupportedSM<<"/"<<smCandidates<<" eligible entries.\n";
    print_report(fit,sums,ids,referenceNoise,beamDays,polarization,usePE,useShowerMaxPE);
    std::cout<<"Saved "<<output<<"\n";
}
} // namespace standalone_deconv

void deconvolution_correlated(bool usePE,const char* tableDirectory,const char* output,
                   double beamDays,double polarization,double current_uA,double minimumK,bool useShowerMaxPE){
    try{standalone_deconv::run(fnms,usePE,tableDirectory,output,beamDays,polarization,current_uA,minimumK,useShowerMaxPE);}
    catch(const std::exception& error){std::cerr<<"Deconvolution stopped: "<<error.what()<<std::endl;throw;}
}
