// ROOT-only entry point. No Python, CMake, Geant4 installation or remoll library.
// root -l -b -q 'Deconvolution.C+("inputs.txt","deconvolution.root")'
// inputs.txt: one ROOT filename per line, relative to that list or absolute.
// mainPE=true: main quartz PE + ShowerMax PE. false: crossings + ShowerMax PE.
#include "metadata_schema.hh"
#include "LinkDef.h"
#include "root_solver.hh"
#include <TMath.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TMD5.h>
#include <TObjString.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TROOT.h>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <regex>
#include <set>
#include <tuple>
#include "../vendor/portable/scan.C"
#include "pe_response.hh"
namespace root_deconv {
namespace fs=std::filesystem;
const char* component_names[]={"Moller","ep elastic","ep inelastic Delta","ep inelastic resonance","ep inelastic continuum"};
inline std::string text(const fs::path& p){std::ifstream f(p);if(!f)throw std::runtime_error("Cannot open "+p.string());return std::string(std::istreambuf_iterator<char>(f),{});}
inline std::string trim(std::string s){auto a=s.find_first_not_of(" \t\r\n");return a==s.npos?"":s.substr(a,s.find_last_not_of(" \t\r\n")-a+1);}
inline std::string digest(const fs::path& p){std::unique_ptr<TMD5> h(TMD5::FileChecksum(p.string().c_str()));if(!h)throw std::runtime_error("Cannot fingerprint "+p.string());return h->AsString();}
inline std::string digest_text(const std::string& s){TMD5 h;h.Update(reinterpret_cast<const UChar_t*>(s.data()),s.size());h.Final();return h.AsString();}
// Shared provenance checks for both PE maps and the deconvolution.
inline void check_sampling(TFile& file,TTree* tree,const fs::path& path,std::string macro){
    std::transform(macro.begin(),macro.end(),macro.begin(),[](unsigned char c){return std::tolower(c);});
    if(macro.find("neyman")!=std::string::npos||macro.find("fixed_quota")!=std::string::npos||file.Get("sampling_plan")||
       tree->GetBranch("hit_weight")||tree->GetBranch("sampling_weight")||fs::exists(path.parent_path()/(path.stem().string()+".neyman.tsv")))
        throw std::runtime_error("Sampled input needs its stratified normalization adapter. Hit counts remain available.");
}
inline std::map<std::string,std::string> saved_geometry(const remollRunData& data){
    if(data.fGitInfo.empty())throw std::runtime_error("Simulation source identity is missing. Hit counts remain available.");
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
    if(gdml.empty())throw std::runtime_error("Saved geometry is missing. Hit counts remain available.");
    return gdml;
}
inline void check_main_geometry(const std::map<std::string,std::string>& gdml,const std::map<std::string,std::string>& expected){
    bool thin=false;
    for(const auto& item:gdml)if(item.first.rfind("detector/ThinQuartz/",0)==0){
        thin=true;
        if(!expected.count(item.first)||expected.at(item.first)!=item.second)
            throw std::runtime_error("Main quartz geometry does not match the PE model. Hit counts remain available.");
    }
    if(!thin)throw std::runtime_error("Main PE needs saved thin-quartz geometry. Hit counts remain available.");
}
// Compare transport across interactions, not only generator settings within one.
// Preserve ordered commands: multiple field scale/current commands are meaningful.
inline std::string transport_identity(const remollRunData& data,const std::string& macro){
    std::string identity=data.fGitInfo+"\n";
    for(const auto& g:saved_geometry(data))identity+=g.first+"="+g.second+"\n";
    std::vector<std::string> fields;
    for(const auto& field:data.fMagData){
        auto bounded=[](const char* value,size_t size){const auto* end=static_cast<const char*>(std::memchr(value,0,size));
            if(!end||end==value)throw std::runtime_error("Invalid saved magnetic-field identity");return std::string(value,end);};
        auto name=bounded(field.filename,sizeof(field.filename));
        std::replace(name.begin(),name.end(),'\\','/');name=name.substr(name.find_last_of('/')+1);
        fields.push_back(name+"="+bounded(field.hashsum,sizeof(field.hashsum)));
    }
    std::sort(fields.begin(),fields.end());for(const auto& f:fields)identity+="field "+f+"\n";
    std::istringstream lines(macro);std::string line;
    while(std::getline(lines,line)){
        line=line.substr(0,line.find('#'));std::istringstream words(line);std::string command,word,value;words>>command;
        if(command=="/run/beamOn")break;
        if(command=="/remoll/addfield"){
            if(fields.empty())throw std::runtime_error("Saved magnetic-field fingerprints are missing");continue;
        }
        if(command.rfind("/remoll/evgen/",0)==0||command.rfind("/remoll/SD/",0)==0||
           command=="/remoll/seed"||command=="/remoll/beamcurr"||command=="/remoll/filename"||
           command=="/remoll/setgeofile"||command=="/remoll/setgeofileParallel"||
           command=="/run/initialize"||command=="/run/numberOfThreads"||command=="/run/verbose"||command=="/run/printProgress")continue;
        if(command.rfind("/remoll/",0)!=0&&command.rfind("/process/",0)!=0&&command.rfind("/physics/",0)!=0&&command.rfind("/run/",0)!=0)continue;
        if(command.find("verbose")!=std::string::npos||command=="/remoll/print"||command=="/remoll/target/print")continue;
        while(words>>word)value+=(value.empty()?"":" ")+word;
        identity+=command+" "+value+"\n";
    }
    return digest_text(identity);
}
struct Input {
    fs::path path;long long events=0,seed=0;int process=-1;double current=0;
    std::string commit,geometry,macro,transport;
    std::map<std::string,std::string> settings;
};
inline void check_transport(const Input& a,const Input& b){
    if(a.transport!=b.transport)throw std::runtime_error("Incompatible transport: source, geometry, magnetic fields or beam/physics settings differ");
}
inline void check_generation_settings(const Input& a,const Input& b){
    if(a.process!=b.process)return;
    std::map<std::string,std::string> left,right;
    for(const auto& setting:a.settings)if(setting.first.rfind("/remoll/evgen/",0)==0)left.insert(setting);
    for(const auto& setting:b.settings)if(setting.first.rfind("/remoll/evgen/",0)==0)right.insert(setting);
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
    TFile file(path.string().c_str(),"READ");auto* tree=file.Get<TTree>("T");auto* data=file.Get<remollRunData>("run_data");
    if(!tree||!data)throw std::runtime_error("Missing T or run_data in "+path.string());
    Input r;r.path=path;r.events=data->fNthrown;r.seed=data->fSeed;r.commit=data->fGitInfo;r.macro=data->fMacro.fBuffer?data->fMacro.fBuffer:"";
    if(r.events<2||r.events!=tree->GetEntries())throw std::runtime_error("Generated/stored event count mismatch in "+path.string());
    if(r.commit.empty())throw std::runtime_error("Simulation source identity is missing");
    check_sampling(file,tree,path,r.macro);
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
    return r;
}
inline portable_deconv::Moments scan(const Input& input,const std::vector<int>& ids,const portable_deconv::Fits& sm,
                                    MainModel* main,double minimumK){
    TFile file(input.path.string().c_str(),"READ");auto* tree=file.Get<TTree>("T");
    auto has=[&](const char* key){return tree->GetLeaf(key)||tree->GetBranch(key);};
    for(auto key:{"rate","ev.A","ev.W2","hit.det","hit.pid","hit.k","hit.e","hit.x","hit.y","hit.pz"})if(!has(key))throw std::runtime_error(std::string("Missing ROOT field: ")+key);
    if(main)for(auto key:{"hit.z","hit.xl","hit.yl","hit.zl","hit.px","hit.py"})if(!has(key))throw std::runtime_error(std::string("Main PE needs ROOT field: ")+key);
    TTreeReader reader(tree);TTreeReaderValue<double> rate(reader,"rate"),asym(reader,"ev.A"),w2(reader,"ev.W2");
    TTreeReaderArray<int> det(reader,"hit.det"),pid(reader,"hit.pid");
    TTreeReaderArray<double> k(reader,"hit.k"),e(reader,"hit.e"),x(reader,"hit.x"),y(reader,"hit.y"),pz(reader,"hit.pz");
    auto optional=[&](const char* key){return main?std::make_unique<TTreeReaderArray<double>>(reader,key):nullptr;};
    auto z=optional("hit.z"),xl=optional("hit.xl"),yl=optional("hit.yl"),zl=optional("hit.zl"),px=optional("hit.px"),py=optional("hit.py");
    const int n=ids.size(),d=n+1,dim=2*d,nc=input.process==2?3:1;
    portable_deconv::Moments out;out.y.resize(nc*dim);out.cross.resize(nc*dim*dim);out.noise.resize(n*n);out.coverage.resize(n);
    std::map<int,int> index;for(int i=0;i<n;i++)index[ids[i]]=i;
    while(reader.Next()){
        ++out.entries;
        if(!std::isfinite(*rate)||*rate<0||!std::isfinite(*asym))throw std::runtime_error("Invalid event rate/asymmetry");
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
                double pe=portable_deconv::pe(sm,det[j],pid[j],e[j],x[j],y[j]);
                if(pe<0){out.unsupported_sm++;continue;}q[i]+=pe;
            }else{
                // Main quartz cuts: charged particles, kinetic energy >= 1 MeV.
                if(!portable_deconv::charged(pid[j])||k[j]<minimumK)continue;
                int suffix=det[j]%100;double phi=std::fmod(std::atan2(y[j],x[j])+2*M_PI,2*M_PI/7);
                // Reference stays crossing-weighted in both PE and crossing modes.
                if(suffix>=51&&suffix<=53&&phi>=3*M_PI/28&&phi<5*M_PI/28)q[n]+=1;
                if(!main){q[i]+=1;continue;}
                // Upstream main-PE model support restrictions.
                if(std::abs(pid[j])!=11||pz[j]<=0||k[j]<=1.1)continue;
                root_main_pe::remollGenericDetectorHit_t h;
                h.det=det[j];h.x=x[j];h.y=y[j];h.z=(*z)[j];h.xl=(*xl)[j];h.yl=(*yl)[j];h.zl=(*zl)[j];h.px=(*px)[j];h.py=(*py)[j];h.pz=pz[j];
                q[i]+=main->pe(h,k[j]);
            }
        }
        int c=0;
        if(input.process==2){double w=*w2>0?std::sqrt(*w2)*.001:-1;c=w>=1&&w<1.4?0:w>=1.4&&w<2.5?1:w>=2.5&&w<6?2:-1;}
        if(c<0){if(std::any_of(q.begin(),q.end(),[](double v){return v>0;}))throw std::runtime_error("Signal outside inelastic W bins [1,6) GeV");continue;}
        std::vector<int> nz;for(int i=0;i<d;i++)if(q[i])nz.push_back(i);
        for(int i:nz){
            double yi=*rate*q[i];out.y[c*dim+i]+=yi;out.y[c*dim+d+i]+=yi*(*asym);
            for(int j:nz){
                double yy=yi*(*rate)*q[j];int off=c*dim*dim;
                out.cross[off+i*dim+j]+=yy;out.cross[off+i*dim+d+j]+=yy*(*asym);
                out.cross[off+(d+i)*dim+j]+=yy*(*asym);out.cross[off+(d+i)*dim+d+j]+=yy*(*asym)*(*asym);
                if(i<n&&j<n)out.noise[i*n+j]+=*rate*q[i]*q[j];
            }
        }
    }
    if(out.entries!=input.events)throw std::runtime_error("Incomplete ROOT read");return out;
}
inline void write_fit(TFile& file,Fit& fit,const std::vector<int>& ids,const TMatrixD& sums,const std::string& provenance,const std::vector<std::string>& labels={}){
    file.cd();
    for(auto item:std::initializer_list<std::pair<const char*,TMatrixD*>>{
        {"rates",&fit.rates},{"F",&fit.F},{"F_mc_sigma",&fit.F_mc_sigma},{"H",&fit.H},{"B",&fit.B},{"C",&fit.C},{"correlation",&fit.correlation},{"M",&fit.M},
        {"experimental_covariance",&fit.experimental_covariance},{"mc_covariance",&fit.mc_covariance},{"combined_covariance",&fit.combined_covariance}})item.second->Write(item.first);
    fit.reference_asymmetry_ppb.Write("reference_asymmetry_ppb");sums.Write("rate_asymmetry_sums");
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
// Public macro entry point. ROOT's .L Deconvolution.C+ also exposes the functions
// above for interactive inspection, different cuts, and alternative weighting.
void Deconvolution(const char* inputs,const char* output="deconvolution.root",bool mainPE=true,
                   double minimumK=1,double beamDays=344,double polarization=.8,double current_uA=65,
                   const char* mainTables="",const char* showermaxTables="",const char* extraGenerators=""){
    using namespace root_deconv;
    try{
        fs::path here=fs::absolute(fs::path(__FILE__)).parent_path().parent_path();
        fs::path tool=here/"data";
        fs::path mainFolder=tool/"main_pe";
        if(mainTables && *mainTables)mainFolder=fs::absolute(mainTables);
        fs::path smFolder=showermaxTables && *showermaxTables?fs::absolute(showermaxTables):tool/"response";
        fs::path channelFile=tool/"channels.json";
        std::string channelText=text(channelFile);std::regex detectorPattern("\\\"detector\\\"\\s*:\\s*([0-9]+)");
        std::vector<int> ids;for(auto i=std::sregex_iterator(channelText.begin(),channelText.end(),detectorPattern);i!=std::sregex_iterator();++i)ids.push_back(std::stoi((*i)[1]));
        if(ids.size()!=252||std::set<int>(ids.begin(),ids.end()).size()!=252)throw std::runtime_error("Expected 224 quartz and 28 ShowerMax channel IDs");
        if(fs::exists(output))throw std::runtime_error("Output already exists; choose another file");
        if(!(std::isfinite(minimumK)&&minimumK>=0&&std::isfinite(beamDays)&&beamDays>0&&std::isfinite(polarization)&&polarization>0&&polarization<=1&&std::isfinite(current_uA)&&current_uA>0))throw std::runtime_error("Invalid cut or forecast settings");
        std::map<std::string,std::string> thinMD5;
        std::istringstream hashes(text(tool/"geometry.md5"));std::string hash,geometryName;
        while(hashes>>hash>>geometryName)thinMD5[geometryName]=hash;
        fs::path list=fs::absolute(inputs);std::istringstream lines(text(list));std::string line;
        std::vector<std::string> extra,labels(component_names,component_names+5);
        std::istringstream extraWords(extraGenerators?extraGenerators:"");std::string name;
        while(extraWords>>name){if(name=="moller"||name=="elastic"||name=="inelastic"||std::find(extra.begin(),extra.end(),name)!=extra.end())throw std::runtime_error("Duplicate background name");extra.push_back(name);labels.push_back(name);}
        if(extra.size()>8)throw std::runtime_error("At most eight additional backgrounds");
        const int components=labels.size();
        std::vector<Input> files;std::set<fs::path> paths;std::set<std::pair<int,long long>> seeds;std::vector<long long> totals(3+extra.size());
        while(std::getline(lines,line)){
            line=trim(line);if(line.empty()||line[0]=='#')continue;
            fs::path path=fs::canonical(list.parent_path()/fs::path(line));
            if(!paths.insert(path).second)throw std::runtime_error("Duplicate ROOT file");
            Input r=inspect(path,ids,mainPE,thinMD5,extra);
            if(!seeds.insert({r.process,r.seed}).second)throw std::runtime_error("Repeated seed in an interaction");
            if(!files.empty()){
                check_transport(files[0],r);
                if(r.commit!=files[0].commit||r.geometry!=files[0].geometry)throw std::runtime_error("Simulation source or geometry differs between inputs");
                for(auto key:{"/remoll/beamene","/remoll/target/mother","/remoll/target/volume","/remoll/physlist/set"})if(r.settings[key]!=files[0].settings[key])throw std::runtime_error(std::string("Incompatible setting: ")+key);
            }
            for(const auto& previous:files)check_generation_settings(previous,r);
            totals[r.process]+=r.events;files.push_back(std::move(r));
        }
        for(size_t p=0;p<totals.size();p++)if(!totals[p])throw std::runtime_error(p<3?"Provide Moller, ep elastic and ep inelastic ROOT files":"Missing configured background: "+extra[p-3]);
        const int n=ids.size(),d=n+1,dim=2*d;
        TMatrixD sums(components,dim),mc(components*dim,components*dim),b(n,n);std::vector<double> coverage(n);
        auto sm=portable_deconv::load(smFolder.string());std::unique_ptr<MainModel> model;if(mainPE)model=std::make_unique<MainModel>(mainFolder);
        std::ostringstream provenance;provenance<<std::setprecision(17)<<"ROOT "<<gROOT->GetVersion()<<"\nMode: "<<(mainPE?"main PE + ShowerMax PE":"main crossings + ShowerMax PE")<<"\nMain charged K >= "<<minimumK<<" MeV; main PE also forward e+/- K > 1.1, polar capped at 4 degrees\nShowerMax forward supported PE\nEntries and re-entries; exits absent by producer recording\nForecast "<<beamDays<<" days, polarization "<<polarization<<", current "<<current_uA<<" uA\nMC propagated through H at simulated reference; experimental C fixed\nExpected mean PE; no additional optical/electronics/beam noise or model systematic uncertainty\n";
        provenance<<"Deconvolution.C MD5 "<<digest(here/"root/Deconvolution.C")<<"\nroot_solver.hh MD5 "<<digest(here/"root/root_solver.hh")<<"\nmetadata_schema.hh MD5 "<<digest(here/"root/metadata_schema.hh")<<"\n";
        provenance<<"Main PE tables "<<mainFolder<<"\nShowerMax tables "<<smFolder<<"\n";
#ifdef MOLLER_MAIN_LOOKUP_HEADER
        provenance<<"Main evaluator MD5 "<<digest(MOLLER_MAIN_LOOKUP_HEADER)<<"\nCoordinate evaluator MD5 "<<digest(MOLLER_MAIN_COORDINATE_HEADER)<<"\n";
#endif
        for(const auto& folder:{mainFolder,smFolder})for(const auto& p:fs::directory_iterator(folder))if(p.path().extension()==".csv")provenance<<p.path().filename().string()<<" MD5 "<<digest(p.path())<<"\n";
        long long events=0;
        for(const auto& file:files){
            std::cout<<"Scanning "<<file.path.filename()<<" ("<<file.events<<" events)"<<std::endl;
            auto m=scan(file,ids,sm,model.get(),minimumK);
            int nc=file.process==2?3:1,first=file.process==2?2:file.process<2?file.process:file.process+2;
            double factor=double(file.events)/totals[file.process]*current_uA/file.current;
            for(int c=0;c<nc;c++)for(int i=0;i<dim;i++)sums(first+c,i)+=factor*m.y[c*dim+i];
            // MC covariance: centered squared-weight moments, including zero scores.
            for(int c=0;c<nc;c++)for(int e=0;e<nc;e++)for(int i=0;i<dim;i++)for(int j=0;j<dim;j++){
                double value=-m.y[c*dim+i]*m.y[e*dim+j]/(file.events-1);
                if(c==e)value+=m.cross[c*dim*dim+i*dim+j]*file.events/(file.events-1);
                mc((first+c)*dim+i,(first+e)*dim+j)+=factor*factor*value;
            }
            // Physical second moment B: one event rate weight.
            for(int i=0;i<n;i++){coverage[i]+=m.coverage[i];for(int j=0;j<n;j++)b(i,j)+=factor*m.noise[i*n+j];}
            events+=file.events;
            provenance<<"Input "<<file.path.string()<<"\nMD5 "<<digest(file.path)<<"\nTransport "<<file.transport<<"\nGenerated "<<file.events<<"; seed "<<file.seed<<"; current "<<file.current<<" uA; factor "<<factor<<"; unsupported SM "<<m.unsupported_sm<<"\nCommit "<<file.commit<<"\nMacro:\n"<<file.macro<<"\nGeometry fingerprints:\n"<<file.geometry;
        }
        for(int i=0;i<n;i++)if(!coverage[i])throw std::runtime_error("No recorded hits in detector "+std::to_string(ids[i]));
        auto fit=solve(sums,mc,b,beamDays,polarization);
        TFile result(output,"CREATE");if(result.IsZombie())throw std::runtime_error("Cannot create output ROOT file");
        write_fit(result,fit,ids,sums,provenance.str(),labels);result.Close();
        std::cout<<"\n"<<events<<" events; "<<n<<" detector channels. Uncertainties [ppb]:\n"<<"Interaction                  Statistical       MC       Combined\n";
        for(int c=0;c<components;c++)std::cout<<std::left<<std::setw(28)<<labels[c]<<std::right<<std::setw(12)<<std::sqrt(fit.experimental_covariance(c,c))<<std::setw(12)<<std::sqrt(std::max(0.,fit.mc_covariance(c,c)))<<std::setw(12)<<std::sqrt(fit.combined_covariance(c,c))<<"\n";
        std::cout<<"Saved "<<output<<"\n";
    }catch(const std::exception& e){std::cerr<<"Deconvolution stopped: "<<e.what()<<std::endl;throw;}
}
