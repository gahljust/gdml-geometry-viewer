#include "../root/json.hh"
// Read hits event by event; apply a ROOT selection. No implicit cuts or weights.
// ROOT only: root -n -l -b -q 'load_hits.C("input.root","counts.json","hit.det==73003")'
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TTreeFormula.h>
#include <TH1D.h>
#include <TH2D.h>
#include <functional>
#include <limits>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <array>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <sstream>

// Same charge definition as the campaign main-quartz selection (PDG nuclear Z for ions).
bool charged(int pid) {
    const auto code = std::abs(static_cast<long long>(pid));
    if (code >= 1000000000LL) return ((code / 10000) % 1000) > 0;
    if (!pid) return false;
    auto* particle = TDatabasePDG::Instance()->GetParticle(pid);
    if (!particle) throw std::runtime_error("Unknown particle charge for PDG " + std::to_string(pid));
    return particle->Charge() != 0;
}

struct HitFace {int id,det;std::array<double,3> o,u,v,n;std::vector<std::array<double,2>> polygon;};
double dot(const std::array<double,3>& a,const std::array<double,3>& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}

// Compact ROOT histograms, filled inside the same selected-hit loop as the map.
struct HitHistogram {
    std::string key, label, expression; bool eventValue=false, logBins=false; TLeaf* leaf; TH1D hist; bool available; Long64_t invalid=0;
    HitHistogram(const char* name,TLeaf* field,int n,double low,double high,bool log=false,bool present=false)
      : key(name),leaf(field),hist(name,"",n,low,high),available(field||present) {
        hist.SetDirectory(nullptr); logBins=log;
        if(log){std::vector<double> edges(n+1);for(int i=0;i<=n;++i)edges[i]=low*std::pow(high/low,double(i)/n);hist.SetBins(n,edges.data());}
    }
    void fill(double value,double weight=1){if(std::isfinite(value))hist.Fill(value,weight);else ++invalid;}
    void write(std::ostream& out){
        out<<viewer_json::quote(key)<<":{\"label\":"<<viewer_json::quote(label)<<",\"expression\":"<<viewer_json::quote(expression)<<",\"event_value\":"<<(eventValue?"true":"false")<<",\"log_bins\":"<<(logBins?"true":"false")<<",\"available\":"<<(available?"true":"false")<<",\"invalid\":"<<invalid
           <<",\"underflow\":"<<hist.GetBinContent(0)<<",\"overflow\":"<<hist.GetBinContent(hist.GetNbinsX()+1)<<",\"edges\":[";
        for(int i=1;i<=hist.GetNbinsX()+1;++i){if(i>1)out<<',';out<<hist.GetXaxis()->GetBinLowEdge(i);}out<<"],\"counts\":[";
        for(int i=1;i<=hist.GetNbinsX();++i){if(i>1)out<<',';out<<hist.GetBinContent(i);}out<<"]}";
    }
};

struct JointHistogram {
    std::string key, x, y; TH2D hist; bool available; Long64_t invalid=0,outside=0;
    JointHistogram(const char* name,HitHistogram& a,HitHistogram& b)
      :key(name),x(a.key),y(b.key),hist(name,"",a.hist.GetNbinsX(),0,1,b.hist.GetNbinsX(),0,1),available(a.available&&b.available){
        hist.SetDirectory(nullptr);
        std::vector<double> ax,bx;
        for(int i=1;i<=a.hist.GetNbinsX()+1;++i)ax.push_back(a.hist.GetXaxis()->GetBinLowEdge(i));
        for(int i=1;i<=b.hist.GetNbinsX()+1;++i)bx.push_back(b.hist.GetXaxis()->GetBinLowEdge(i));
        hist.SetBins(ax.size()-1,ax.data(),bx.size()-1,bx.data());
    }
    void fill(double a,double b,double weight=1){
        if(!std::isfinite(a)||!std::isfinite(b)){++invalid;return;}
        if(a<hist.GetXaxis()->GetXmin()||a>=hist.GetXaxis()->GetXmax()||b<hist.GetYaxis()->GetXmin()||b>=hist.GetYaxis()->GetXmax())++outside;
        else hist.Fill(a,b,weight);
    }
    void write(std::ostream& out){
        out<<viewer_json::quote(key)<<":{\"x\":"<<viewer_json::quote(x)<<",\"y\":"<<viewer_json::quote(y)<<",\"available\":"<<(available?"true":"false")
           <<",\"invalid\":"<<invalid<<",\"outside\":"<<outside<<",\"cells\":[";
        bool first=true;
        for(int i=1;i<=hist.GetNbinsX();++i)for(int j=1;j<=hist.GetNbinsY();++j)if(auto n=hist.GetBinContent(i,j)){
            if(!first)out<<',';first=false;out<<'['<<i-1<<','<<j-1<<','<<n<<']';
        }out<<"]}";
    }
};

void load_hits(const char* input = nullptr, const char* output = nullptr, const char* cut = nullptr) {
    try {
        if (!input) input = gSystem->Getenv("MOLLER_HIT_INPUT");
        if (!output) output = gSystem->Getenv("MOLLER_HIT_OUTPUT");
        if (!input || !output) throw std::runtime_error("Provide an input ROOT file and output counts file.");
        if (!cut) cut = gSystem->Getenv("MOLLER_HIT_CUT");
        if (!cut || !*cut) cut = "1";
        std::unique_ptr<TFile> file(TFile::Open(input, "READ"));
        if (!file || file->IsZombie()) throw std::runtime_error("Cannot open ROOT file.");
        auto* tree = file->Get<TTree>("T");
        if (!tree) throw std::runtime_error("Missing event tree T.");
        auto* detector = tree->GetLeaf("hit.det");
        if (!detector) throw std::runtime_error("Missing hit.det; a split remoll hit branch is required.");
        if (std::string(detector->GetTypeName()) != "Int_t")
            throw std::runtime_error("hit.det must contain integer detector IDs.");

        const auto events = tree->GetEntries();
        std::string weighting="hits",normalization="",generator="",transport="";double current=0,signalTotal=0;Long64_t seed=0,unsupportedPE=0;
        std::map<int,double> detectorSignal,particleSignal;
#ifdef MOLLER_WEIGHTED_HITS
        weighting="pe";
        auto* rate=tree->GetLeaf("rate");
        auto* run=file->Get<remollRunData>("run_data");
        if(!rate)throw std::runtime_error("PE rate needs ROOT field rate.");
        if(!run||run->fNthrown!=events||events<1)throw std::runtime_error("PE rate/current needs run_data with the generated event count. Hit counts remain available.");
        std::map<std::string,std::string> settings;std::istringstream macro(run->fMacro.fBuffer?run->fMacro.fBuffer:"");std::string line;int runs=0;
        while(std::getline(macro,line)){
            line=line.substr(0,line.find('#'));std::istringstream words(line);std::string command,value;words>>command;
            if(command=="/run/beamOn"){++runs;continue;}if(runs)continue;
            if(command=="/control/execute"||command=="/control/loop"||command=="/control/foreach")settings.clear();
            std::getline(words,value);if(command.rfind("/remoll/",0)==0)settings[command]=root_deconv::trim(value);
        }
        if(runs!=1||settings["/remoll/evgen/set"].empty())throw std::runtime_error("PE rate needs a saved generator and one effective beamOn; no normalization is guessed.");
        std::string macroText=run->fMacro.fBuffer?run->fMacro.fBuffer:"";
        root_deconv::check_sampling(*file,tree,input,macroText);
        const auto gdml=root_deconv::saved_geometry(*run);
        std::map<std::string,std::string> modelGeometry;
        std::istringstream modelHashes(root_deconv::text(std::filesystem::path(__FILE__).parent_path().parent_path()/"data/geometry.md5"));
        std::string hash,name;while(modelHashes>>hash>>name)modelGeometry[name]=hash;
        bool mainGeometryChecked=false;
        std::set<int> verifiedEntries;
        std::istringstream beam(settings["/remoll/beamcurr"]);std::string unit;beam>>current>>unit;
        if(unit=="mA")current*=1000;else if(unit=="ampere")current*=1e6;else if(unit=="nA")current*=.001;else if(unit!="uA"&&unit!="microampere")current=0;
        if(!(current>0)||!std::isfinite(current))throw std::runtime_error("PE rate/current needs the saved beam current.");
        seed=run->fSeed;generator=settings["/remoll/evgen/set"];
        transport=root_deconv::transport_identity(*run,macroText);
        std::string runIdentity=transport;
        for(const auto& g:gdml)runIdentity+=g.first+g.second;
        for(const auto& entry:settings)if(entry.first!="/remoll/beamcurr"&&entry.first!="/remoll/seed")runIdentity+=entry.first+"="+entry.second+"\n";
        normalization=root_deconv::digest_text(runIdentity);
        root_deconv::MainModel mainPE(gSystem->Getenv("MOLLER_MAIN_PE"));
        auto showerPE=portable_deconv::load(gSystem->Getenv("MOLLER_SM_PE"));
        std::map<std::string,TLeaf*> responseFields;
        for(auto name:{"hit.pid","hit.k","hit.e","hit.x","hit.y","hit.z","hit.xl","hit.yl","hit.zl","hit.px","hit.py","hit.pz"})responseFields[name]=tree->GetLeaf(name);
        auto responseValue=[&](const char* name,int index){auto* leaf=responseFields.at(name);if(!leaf)throw std::runtime_error(std::string("PE response needs ROOT field: ")+name);if(leaf->GetNdata()!=detector->GetNdata())throw std::runtime_error(std::string("Mismatched PE field: ")+name);double v=leaf->GetValue(index);if(!std::isfinite(v))throw std::runtime_error(std::string("Nonfinite PE field: ")+name);return v;};
#endif

        // A cut may use hit fields or scalar event fields present in this file.
        std::string expression(cut);
        // TTreeFormula does not accept the editor's line breaks as whitespace.
        for(auto& c : expression) if(c=='\n'||c=='\r'||c=='\t') c=' ';
        TTreeFormula selection("hit_selection", expression.c_str(), tree);
        if (!selection.GetNdim())
            throw std::runtime_error("Invalid ROOT cut or missing field. Check the expression and available fields.");
        std::vector<TLeaf*> cutFields;
        for (int code = 0; code < selection.GetNcodes(); ++code) {
            auto* leaf = selection.GetLeaf(code);
            if (!leaf) continue;
            auto* branch = leaf->GetBranch();
            const std::string name = branch->GetName();
            const std::string owner = branch->GetMother()->GetClassName();
            if(name.rfind("hit.",0)==0)cutFields.push_back(leaf);
            if (name.rfind("hit.", 0) != 0 &&
                (leaf->GetLeafCount() || leaf->GetLenStatic() > 1 || owner.find("vector<") != std::string::npos))
                throw std::runtime_error("Cut uses a different particle array: " + name + ". Use hit fields or scalar event fields.");
            if (branch->GetEntries() != events)
                throw std::runtime_error("Incomplete cut branch: " + name);
        }
        for (auto* object : *tree->GetListOfLeaves()) {
            auto* branch = static_cast<TLeaf*>(object)->GetBranch();
            const std::string name = branch->GetName();
            if ((name == "hit" || name.rfind("hit.", 0) == 0) && branch->GetEntries() != events)
                throw std::runtime_error("Incomplete hit branch: " + name);
        }

        auto* x = tree->GetLeaf("hit.x");
        auto* y = tree->GetLeaf("hit.y");
        auto* z = tree->GetLeaf("hit.z");
        const bool positions = x && y && z;
        // All plots count the same selected hits. Scalar event values repeat per selected hit.
        std::map<std::string,std::unique_ptr<HitHistogram>> plots;
        std::map<std::string,std::function<double(int)>> getters;
        std::map<TLeaf*,bool> plotFields;
        const double nan=std::numeric_limits<double>::quiet_NaN();
        auto add=[&](const char* key,const char* label,const char* expression,std::initializer_list<const char*> fields,
                     int n,double low,double high,bool log,std::function<double(const std::array<double,3>&)> calc,bool scalar=false){
            std::vector<TLeaf*> leaves; bool present=true;
            for(auto name:fields){auto* leaf=tree->GetLeaf(name);present&=leaf!=nullptr;leaves.push_back(leaf);
                if(leaf){plotFields[leaf]=scalar;
                    const std::string owner=leaf->GetBranch()->GetMother()->GetClassName();
                    if(leaf->GetBranch()->GetEntries()!=events || (scalar&&(leaf->GetLeafCount()||leaf->GetLenStatic()!=1||owner.find("vector<")!=std::string::npos)))
                        throw std::runtime_error(std::string("Invalid plot field: ")+name);
                }
            }
            auto h=std::make_unique<HitHistogram>(key,nullptr,n,low,high,log,present);
            h->label=label;h->expression=expression;h->eventValue=scalar;plots[key]=std::move(h);
            getters[key]=[=](int hit){if(!present)return nan;std::array<double,3> v{};size_t i=0;for(auto* leaf:leaves){double a=leaf->GetValue(scalar?0:hit);if(!std::isfinite(a))return nan;v[i++]=a;}return calc(v);};
        };
        auto identity=[](const std::array<double,3>& v){return v[0];};
        auto direct=[&](const char* key,const char* label,const char* field,int n,double low,double high,bool log=false,bool scalar=false){add(key,label,field,{field},n,low,high,log,identity,scalar);};
        direct("energy","Kinetic energy [MeV]","hit.k",64,.001,20000,true);
        direct("total_energy","Total energy [MeV]","hit.e",64,.001,20000,true);
        direct("edep","Deposited energy [MeV]","hit.edep",64,0,20);
        add("radius","Radius [mm]","hypot(hit.x,hit.y)",{"hit.x","hit.y"},80,0,2000,false,[](const auto& v){return std::hypot(v[0],v[1]);});
        direct("vertex","Creation z [mm]","hit.vz",90,-1000,35000);
        add("vertex_r","Creation radius [mm]","hypot(hit.vx,hit.vy)",{"hit.vx","hit.vy"},80,0,2000,false,[](const auto& v){return std::hypot(v[0],v[1]);});
        direct("x","Hit x [mm]","hit.x",80,-2000,2000);
        direct("y","Hit y [mm]","hit.y",80,-2000,2000);
        direct("z","Hit z [mm]","hit.z",90,-1000,35000);
        direct("time","Hit time [ns]","hit.t",80,0,400);
        direct("momentum","Momentum [MeV/c]","hit.p",64,.001,20000,true);
        direct("px","Momentum x [MeV/c]","hit.px",80,-2000,2000);
        direct("py","Momentum y [MeV/c]","hit.py",80,-2000,2000);
        direct("pz","Momentum z [MeV/c]","hit.pz",80,-12000,12000);
        add("phi","Position azimuth [deg]","atan2(hit.y,hit.x)*180/pi",{"hit.x","hit.y"},72,-180,180,false,[=](const auto& v){return (v[0]==0&&v[1]==0)?nan:std::atan2(v[1],v[0])*180/std::acos(-1.);});
        add("theta","Momentum polar angle [deg]","atan2(hypot(hit.px,hit.py),hit.pz)*180/pi",{"hit.px","hit.py","hit.pz"},90,0,180.000001,false,[=](const auto& v){return (v[0]==0&&v[1]==0&&v[2]==0)?nan:std::atan2(std::hypot(v[0],v[1]),v[2])*180/std::acos(-1.);});
        direct("spin_x","Polarization x","hit.sx",50,-1.000001,1.000001);
        direct("spin_y","Polarization y","hit.sy",50,-1.000001,1.000001);
        direct("spin_z","Polarization z","hit.sz",50,-1.000001,1.000001);
        direct("asymmetry","Generator asymmetry [ppb]","ev.A",80,-2000,2000,false,true);
        add("q2","Generator Q² [GeV²]","ev.Q2/1e6",{"ev.Q2"},80,0,.2,false,[](const auto& v){return v[0]/1e6;},true);
        add("w","Generator W [GeV]","sqrt(ev.W2)/1000",{"ev.W2"},80,0,5,false,[=](const auto& v){return v[0]<0?nan:std::sqrt(v[0])/1000;},true);
        direct("rate","Stored event rate value","rate",80,1e-6,1e12,true,true);
        direct("beam_p","Beam momentum at vertex [MeV/c]","ev.beamp",80,0,12000,false,true);
        // Optional per-quantity bins, before filling either the 1D or 2D histograms.
        if(const char* settings=gSystem->Getenv("MOLLER_HIT_BINS")){
            std::istringstream rows(settings);std::string line;
            while(std::getline(rows,line)){
                if(line.empty())continue;
                std::istringstream row(line);std::string key,extra;int n,log;double low,high;
                if(!(row>>key>>n>>low>>high>>log)||(row>>extra)||!plots.count(key)||n<1||n>256||
                   !std::isfinite(low)||!std::isfinite(high)||low>=high||std::abs(low)>1e15||std::abs(high)>1e15||
                   (log!=0&&log!=1)||(log&&(low<=0||!std::isfinite(high/low))))
                    throw std::runtime_error("Invalid histogram bin settings.");
                std::vector<double> edges(n+1);
                for(int i=0;i<=n;++i)edges[i]=log?low*std::pow(high/low,double(i)/n):low+(high-low)*double(i)/n;
                edges.front()=low;edges.back()=high;
                for(int i=1;i<=n;++i)if(!std::isfinite(edges[i])||edges[i]<=edges[i-1])throw std::runtime_error("Histogram bins are too narrow.");
                plots.at(key)->hist.SetBins(n,edges.data());plots.at(key)->logBins=log;
            }
        }
        std::vector<std::unique_ptr<JointHistogram>> joints;
        auto pair=[&](const char* key,const char* a,const char* b){joints.push_back(std::make_unique<JointHistogram>(key,*plots.at(a),*plots.at(b)));};
        pair("xy","x","y");pair("energy_radius","radius","energy");pair("creation_rz","vertex","vertex_r");
        pair("radius_z","z","radius");pair("radius_phi","phi","radius");pair("energy_theta","theta","energy");
        pair("energy_time","time","energy");pair("energy_creation","vertex","energy");pair("momentum_radius","radius","momentum");
        pair("px_py","px","py");pair("edep_energy","energy","edep");pair("q2_w","w","q2");pair("asymmetry_q2","q2","asymmetry");
        auto* pid=tree->GetLeaf("hit.pid");std::map<int,Long64_t> particles;
        if(pid && std::string(pid->GetTypeName())!="Int_t")
            throw std::runtime_error("hit.pid must contain integer particle IDs.");
        Long64_t invalidParticles=0;
        // Same nearest-face / 0.02 mm tolerance and 10 mm local bins as scan_hit_surfaces.cc.
        std::map<int,std::vector<HitFace>> faces;
        const char* geometry=gSystem->Getenv("MOLLER_HIT_FACES");
        if(geometry && *geometry){
            std::ifstream in(geometry);if(!in)throw std::runtime_error("Cannot read detector faces.");
            HitFace f;int count;
            while(in>>f.id>>f.det){
                for(auto* vec:{&f.o,&f.u,&f.v,&f.n})for(auto& value:*vec)in>>value;
                in>>count;if(count<3||count>64)throw std::runtime_error("Invalid detector face.");
                f.polygon.resize(count);for(auto& p:f.polygon)in>>p[0]>>p[1];
                if(!in)throw std::runtime_error("Incomplete detector face.");faces[f.det].push_back(f);
            }
            if(faces.empty())throw std::runtime_error("No detector faces.");
        }
        std::map<std::array<int,3>, double> bins;
        Long64_t invalidPositions = 0, unmatchedPositions = 0, otherDetectors = 0;
        std::map<int, Long64_t> counts;
        Long64_t hits = 0, stored = 0, selectedEvents = 0;
        std::map<std::string,double> pointValues;for(auto& item:plots)pointValues[item.first]=nan;
        for (Long64_t event = 0; event < events; ++event) {
            if (tree->GetEntry(event) <= 0)
                throw std::runtime_error("Cannot read event " + std::to_string(event) + ".");
            const int count = detector->GetNdata();
            stored += count;
            if (positions && (x->GetNdata()!=count || y->GetNdata()!=count || z->GetNdata()!=count))
                throw std::runtime_error("Position fields do not match hit.det at event " + std::to_string(event));
            for(auto* field:{pid,x,y})if(field && field->GetNdata()!=count)
                throw std::runtime_error(std::string("Plot field does not match hit.det: ")+field->GetName());
            for(const auto& item:plotFields)if(item.first->GetNdata()!=(item.second?1:count))
                throw std::runtime_error(std::string("Plot field length mismatch: ")+item.first->GetName());
            for(auto* leaf:cutFields)if(leaf->GetNdata()!=count)
                throw std::runtime_error(std::string("Cut field length mismatch: ")+leaf->GetName());
            const int values = selection.GetNdata();
            const bool scalar = values == 1; // Constants, event scalars, reductions or explicit hit[index].
            if (count && values != (scalar ? 1 : count))
                throw std::runtime_error("Cut does not yield one value per hit at event " + std::to_string(event) + ".");
            const auto before = hits;
            for (int hit = 0; hit < count; ++hit) {
                const double keep = selection.EvalInstance(scalar ? 0 : hit);
                if (!std::isfinite(keep))
                    throw std::runtime_error("Nonfinite cut value at event " + std::to_string(event) + ", hit " + std::to_string(hit) + ".");
                if (keep == 0) continue; // Selection, not a rate weight.
                double weight=1;
#ifdef MOLLER_WEIGHTED_HITS
                const int id=detector->GetValue(hit);double pe=0;const int suffix=id%100;
                const bool main=detector_selection::main_quartz(id);
                if(main&&!mainGeometryChecked){root_deconv::check_main_geometry(gdml,modelGeometry);mainGeometryChecked=true;}
                if(main||(id>=73001&&id<=73028)){
                    if(!verifiedEntries.count(id)){
                        if(!root_deconv::entry_recorded(macroText,id))throw std::runtime_error("PE rate/current requires enabled entry recording for detector "+std::to_string(id)+". Internal steps cannot be PE weighted; hit counts remain available.");
                        verifiedEntries.insert(id);
                    }
                    int code=responseValue("hit.pid",hit);double pzValue=responseValue("hit.pz",hit);
                    if(id>=73001&&id<=73028){if(pzValue>0){pe=portable_deconv::pe(showerPE,id,code,responseValue("hit.e",hit),responseValue("hit.x",hit),responseValue("hit.y",hit));if(pe<0){++unsupportedPE;pe=0;}}}
                    else if(std::abs(code)==11&&pzValue>0&&responseValue("hit.k",hit)>1.1){
                        root_main_pe::remollGenericDetectorHit_t h;h.det=id;
                        h.x=responseValue("hit.x",hit);h.y=responseValue("hit.y",hit);h.z=responseValue("hit.z",hit);
                        h.xl=responseValue("hit.xl",hit);h.yl=responseValue("hit.yl",hit);h.zl=responseValue("hit.zl",hit);
                        h.px=responseValue("hit.px",hit);h.py=responseValue("hit.py",hit);h.pz=pzValue;
                        pe=mainPE.pe(h,responseValue("hit.k",hit));
                    }else ++unsupportedPE;
                }else ++unsupportedPE;
                if(rate->GetNdata()!=1||!std::isfinite(rate->GetValue())||rate->GetValue()<0)throw std::runtime_error("Invalid stored event rate.");
                weight=pe*rate->GetValue()*65/current*events; // Divide by total histories per compatible cohort in Python.
                if(!std::isfinite(weight)||weight<0)throw std::runtime_error("Invalid PE response weight.");
#endif
                signalTotal+=weight;detectorSignal[static_cast<int>(detector->GetValue(hit))]+=weight;
                if(pid&&std::isfinite(pid->GetValue(hit))&&pid->GetValue(hit)>=-2147483648.&&pid->GetValue(hit)<=2147483647.)particleSignal[int(pid->GetValue(hit))]+=weight;
                ++counts[static_cast<int>(detector->GetValue(hit))];
                ++hits;
                for(auto& item:plots){double value=getters.at(item.first)(hit);pointValues[item.first]=value;if(item.second->available)item.second->fill(value,weight);}
                for(auto& joint:joints)if(joint->available)joint->fill(pointValues.at(joint->x),pointValues.at(joint->y),weight);
                if(pid){const double code=pid->GetValue(hit);if(std::isfinite(code)&&code>=-2147483648.&&code<=2147483647.)++particles[int(code)];else ++invalidParticles;}
                if (positions) {
                    const double px=x->GetValue(hit), py=y->GetValue(hit), pz=z->GetValue(hit);
                    if (!std::isfinite(px)||!std::isfinite(py)||!std::isfinite(pz) ||
                        std::abs(px)>1e12||std::abs(py)>1e12||std::abs(pz)>1e12) { ++invalidPositions; continue; }
                    if(faces.empty())continue;
                    const int det=static_cast<int>(detector->GetValue(hit));
                    if(!faces.count(det)){++otherDetectors;continue;}
                    const std::array<double,3> point{px,py,pz};
                    double best=1e99,bu=0,bv=0;int face=-1;
                    for(const auto& f:faces.at(det)){
                        std::array<double,3> d{px-f.o[0],py-f.o[1],pz-f.o[2]};
                        const double distance=std::abs(dot(d,f.n));if(distance>=best)continue;
                        const double u=dot(d,f.u),v=dot(d,f.v);bool pos=false,neg=false;
                        for(size_t j=0;j<f.polygon.size();++j){const auto& p=f.polygon[j];const auto& q=f.polygon[(j+1)%f.polygon.size()];
                            const double c=(q[0]-p[0])*(v-p[1])-(q[1]-p[1])*(u-p[0]);pos|=c>.02;neg|=c<-.02;}
                        if(pos&&neg)continue;best=distance;face=f.id;bu=u;bv=v;
                    }
                    if(face<0||best>.02){++unmatchedPositions;continue;}
                    bins[{face,int(std::floor(bu/10)),int(std::floor(bv/10))}]+=weight;
                    if(bins.size()>200000)throw std::runtime_error("Too many occupied surface bins. Narrow the detector or hit cut.");
                }
            }
            if (hits != before) ++selectedEvents;
        }

        // Keep only a compact receipt; do not copy the ROOT file or dump its hits.
        std::ofstream out(output);
        out<<std::setprecision(17);
        out << "{\"events\":" << events << ",\"hits\":" << hits << ",\"detectors\":{";
        bool first = true;
        for (const auto& item : counts) {
            if (!first) out << ',';
            first = false;
            out << '"' << item.first << "\":" << item.second;
        }
        out << "},\"stored_hits\":" << stored << ",\"selected_events\":" << selectedEvents << ",\"fields\":[";
        first = true;
        for (auto* object : *tree->GetListOfLeaves()) {
            auto* leaf = static_cast<TLeaf*>(object);
            if (!first) out << ',';
            first = false;
            out << viewer_json::quote(leaf->GetName());
        }
        out << "],\"positions_available\":" << (positions?"true":"false")
            << ",\"surfaces_available\":" << (!faces.empty()?"true":"false")
            << ",\"invalid_positions\":" << invalidPositions << ",\"unmatched_positions\":" << unmatchedPositions
            << ",\"other_detectors\":" << otherDetectors << ",\"cell_mm\":10,\"cells\":[";
        first = true;
        for(const auto& item:bins){
            if(!first)out<<',';first=false;
            out<<'['<<item.first[0]<<','<<item.first[1]<<','<<item.first[2]<<','<<item.second<<']';
        }
        out << "],\"histograms\":{";out<<std::setprecision(17);
        first=true;for(auto& item:plots){if(!first)out<<',';first=false;item.second->write(out);}
        out<<"},\"histograms2d\":{";first=true;for(auto& joint:joints){if(!first)out<<',';first=false;joint->write(out);}
        out<<"},\"particles_available\":"<<(pid?"true":"false")<<",\"invalid_particles\":"<<invalidParticles<<",\"particles\":{";
        first=true;for(const auto& item:particles){if(!first)out<<',';first=false;out<<'"'<<item.first<<"\":"<<item.second;}
        out<<"},\"weighting\":"<<viewer_json::quote(weighting)<<",\"transport\":"<<viewer_json::quote(transport)<<",\"normalization\":"<<viewer_json::quote(normalization)<<",\"generator\":"<<viewer_json::quote(generator)<<",\"seed\":"<<seed<<",\"signal_total\":"<<signalTotal<<",\"unsupported_pe\":"<<unsupportedPE<<",\"detector_signal\":{";
        first=true;for(auto item:detectorSignal){if(!first)out<<',';first=false;out<<'"'<<item.first<<"\":"<<item.second;}out<<"},\"particle_signal\":{";
        first=true;for(auto item:particleSignal){if(!first)out<<',';first=false;out<<'"'<<item.first<<"\":"<<item.second;}out<<"}}\n";
        out.close();
        if (!out) throw std::runtime_error("Cannot write hit counts.");
    } catch (const std::exception& error) {
        std::cerr << "HIT_ERROR " << error.what() << std::endl;
        gSystem->Exit(1);
    }
}
