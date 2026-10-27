#include "../root/json.hh"
// Browser adapter only: calculations remain in the public ROOT macro.
#include "../Deconvolution.C"

void deconvolution_bridge(){
    using namespace root_deconv;
    try{
        auto env=[](const char* name){const char* v=gSystem->Getenv(name);if(!v)throw std::runtime_error(std::string("Missing ")+name);return std::string(v);};
        const fs::path base=fs::path(__FILE__).parent_path().parent_path();
        std::vector<int> ids;std::regex re("\"detector\"\\s*:\\s*([0-9]+)");std::string channels=text(base/"data/channels.json");
        for(auto i=std::sregex_iterator(channels.begin(),channels.end(),re);i!=std::sregex_iterator();++i)ids.push_back(std::stoi((*i)[1]));
        std::map<std::string,std::string> geometry;std::istringstream hashes(text(base/"data/geometry.md5"));std::string hash,name;while(hashes>>hash>>name)geometry[name]=hash;
        std::vector<std::string> extra,labels(component_names,component_names+5);std::istringstream names(env("MOLLER_EXTRA"));while(names>>name){extra.push_back(name);labels.push_back(name);}
        std::istringstream input(text(env("MOLLER_INPUT_LIST")));std::vector<Input> records;std::vector<std::string> paths,errors;std::string line;
        std::ofstream out(env("MOLLER_FIT_JSON"));out<<std::setprecision(17);
        if(env("MOLLER_FIT_ACTION")=="check"){
            std::set<int> present;std::set<std::pair<int,long long>> seeds;bool ready=true;long long total=0;
            out<<"{\"files\":[";bool first=true;
            while(std::getline(input,line)){
                if(line.empty())continue;if(!first)out<<',';first=false;
                out<<"{\"path\":"<<viewer_json::quote(line);
                try{
                    auto r=inspect(fs::path(line),ids,true,geometry,extra);
                    TFile file(line.c_str(),"READ");auto* t=file.Get<TTree>("T");
                    for(auto key:{"rate","ev.A","ev.W2","hit.det","hit.pid","hit.k","hit.e","hit.x","hit.y","hit.z","hit.xl","hit.yl","hit.zl","hit.px","hit.py","hit.pz"}){
                        if(!t->GetLeaf(key)&&!t->GetBranch(key))throw std::runtime_error(std::string("Missing ROOT field: ")+key);
                    }
                    if(!seeds.insert({r.process,r.seed}).second)throw std::runtime_error("Repeated seed in one interaction");
                    if(!records.empty()){
                        check_transport(records[0],r);
                        if(r.commit!=records[0].commit||r.geometry!=records[0].geometry)throw std::runtime_error("Source or geometry differs between files");
                        for(auto key:{"/remoll/beamene","/remoll/target/mother","/remoll/target/volume","/remoll/physlist/set"})if(r.settings[key]!=records[0].settings[key])throw std::runtime_error(std::string("Incompatible ")+key);
                    }
                    for(auto& previous:records)check_generation_settings(previous,r);
                    records.push_back(r);present.insert(r.process);total+=r.events;
                    out<<",\"process\":"<<viewer_json::quote(r.settings["/remoll/evgen/set"])<<",\"events\":"<<r.events<<",\"status\":\"Ready\"";
                }catch(const std::exception& e){ready=false;out<<",\"status\":"<<viewer_json::quote(std::string(e.what()));}
                out<<'}';
            }
            out<<"],\"missing\":[";first=true;
            for(size_t i=0;i<3+extra.size();i++)if(!present.count(i)){ready=false;if(!first)out<<',';first=false;out<<viewer_json::quote(i==0?std::string("Møller"):i==1?std::string("ep elastic"):i==2?std::string("ep inelastic"):extra[i-3]);}
            out<<"],\"ready\":"<<(ready&&!records.empty()?"true":"false")<<",\"events\":"<<total<<'}';
        }else{
            const std::string result=env("MOLLER_FIT_ROOT");
            Deconvolution(env("MOLLER_INPUT_LIST").c_str(),result.c_str(),true,1,
                std::stod(env("MOLLER_DAYS")),std::stod(env("MOLLER_POLARIZATION")),std::stod(env("MOLLER_CURRENT")),
                env("MOLLER_MAIN_PE").c_str(),env("MOLLER_SM_PE").c_str(),env("MOLLER_EXTRA").c_str());
            TFile file(result.c_str(),"READ");
            out<<"{\"labels\":[";for(size_t i=0;i<labels.size();i++){if(i)out<<',';out<<viewer_json::quote(labels[i]);}out<<"],\"detectors\":[";
            for(size_t i=0;i<ids.size();i++){if(i)out<<',';out<<ids[i];}out<<"],\"matrices\":{";bool first=true;
            for(auto key:{"F","F_mc_sigma","H","B","C","correlation","M","rates","experimental_covariance","mc_covariance","combined_covariance"}){
                auto* matrix=file.Get<TMatrixD>(key);if(!matrix)throw std::runtime_error(std::string("Missing result ")+key);
                if(!first)out<<',';first=false;out<<viewer_json::quote(key)<<":[";
                for(int r=0;r<matrix->GetNrows();r++){if(r)out<<',';out<<'[';for(int c=0;c<matrix->GetNcols();c++){if(c)out<<',';double v=(*matrix)(r,c);if(!std::isfinite(v))throw std::runtime_error("Nonfinite result");out<<v;}out<<']';}out<<']';
            }
            auto* reference=file.Get<TVectorD>("reference_asymmetry_ppb");out<<"},\"reference\":[";for(int i=0;i<reference->GetNrows();i++){if(i)out<<',';out<<(*reference)[i];}out<<"]}";
        }
        out.close();if(!out)throw std::runtime_error("Cannot write result JSON");
    }catch(const std::exception& e){std::cerr<<"FIT_ERROR "<<e.what()<<std::endl;gSystem->Exit(1);}
}
