#pragma once
#include <filesystem>
// The upstream evaluator needs only these hit fields, not a remoll I/O class.
namespace root_main_pe {
struct remollGenericDetectorHit_t {int det=0;double x=0,y=0,z=0,xl=0,yl=0,zl=0,px=0,py=0,pz=0;};
using std::vector;
#if defined(MOLLER_MAIN_LOOKUP_HEADER) && defined(MOLLER_MAIN_COORDINATE_HEADER)
#include MOLLER_MAIN_LOOKUP_HEADER
#include MOLLER_MAIN_COORDINATE_HEADER
#else
#include "../vendor/main_response/lookuptable.C"
#include "../vendor/main_response/convert_labframe_to_local_thinquartz.C"
#endif
}
namespace root_deconv {
namespace fs=std::filesystem;
struct MainModel {
    std::array<std::vector<root_main_pe::tableEntry>,6> position;
    std::array<std::vector<root_main_pe::tableEntry_energy>,6> energy;
    explicit MainModel(const fs::path& folder){for(int r=0;r<6;r++){
        auto stem=folder/("r"+std::to_string(r+1));
        // Strict parsing; keep the upstream fit functions unchanged.
        const auto positions=portable_deconv::numeric_table<13>(stem.string()+"_NewLookupTable.csv");
        for(const auto& p:positions)if(p[0]>=p[1]||p[2]>=p[3])throw std::runtime_error("Invalid main PE position bounds: "+stem.string());
        const auto energies=portable_deconv::numeric_table<2>(stem.string()+"_EnergyDepTable.csv");
        if(energies.size()!=500)throw std::runtime_error("Main PE energy table needs 500 rows: "+stem.string());
        for(size_t i=0;i<energies.size();i++)if(std::abs(energies[i][0]-(2.511+2*i))>.001||energies[i][1]<0)
            throw std::runtime_error("Invalid main PE energy grid/scaling: "+stem.string());
        for(const auto& p:positions)position[r].push_back({p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12]});
        for(const auto& e:energies)energy[r].push_back({e[0],e[1]});
    }}
    double pe(root_main_pe::remollGenericDetectorHit_t h,double kinetic){
        if(!detector_selection::main_quartz(h.det))throw std::runtime_error("Unknown main quartz detector: "+std::to_string(h.det));
        int suffix=h.det%100,ring=suffix>=51&&suffix<=53?4:suffix/10-1;
        if(ring<0||ring>=6)throw std::runtime_error("Unknown main quartz detector");
        auto angle=root_main_pe::convert_labframe_to_local_thinquartz(h);
        const double offset[]={15,30,30,60,70,50};
        double value=root_main_pe::cal_PE(h.xl,h.yl+offset[ring],std::min(angle.lpolar,4.),angle.lazimuthal,position[ring]);
        return value<0?0:value*root_main_pe::cal_PE_energy(kinetic,energy[ring]);
    }
};
}
