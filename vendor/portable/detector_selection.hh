#ifndef REMOLL_ANALYSIS_DETECTOR_SELECTION_HH
#define REMOLL_ANALYSIS_DETECTOR_SELECTION_HH
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

// Current ThinQuartz and 28 ShowerMax entrance-plane numbering.
// Region order matches campaign products: closed, transition, open.
namespace detector_selection {
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
