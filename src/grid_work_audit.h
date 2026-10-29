#ifndef APRL_GRID_WORK_AUDIT_H
#define APRL_GRID_WORK_AUDIT_H

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>
#include "eb_field.h"

namespace aprl {
// Opt-in, read-only observer. It measures J dot E independently of particle
// work and of an energy residual. Each rank integrates only its owned cells;
// native shared edges receive the ordinary adjacent-cell quadrature weights.
// Neither the cell-centred nor the native-edge product is asserted to be an
// exact discrete energy invariant of the Cowan operator.
class GridWorkAudit {
public:
  GridWorkAudit(const std::string& directory, int rank, std::size_t cx,
                std::size_t cy) : rank_(rank), cx_(cx), cy_(cy) {
    std::ostringstream name;
    name << directory << "/mesh-work-rank-" << std::setw(5)
         << std::setfill('0') << rank << ".csv";
    std::ifstream existing(name.str());
    if (existing.good()) throw std::runtime_error("Grid-work audit file exists");
    output_.open(name.str());
    if (!output_) throw std::runtime_error("Cannot create grid-work audit file");
    output_ << std::setprecision(21)
      << "step,native_old_J,native_mid_J,native_new_J,"
         "center_old_J,center_mid_J,center_new_J,"
         "native_mid_inlet_overwrite_pairing_J,center_mid_inlet_overwrite_pairing_J\n";
  }
  void begin(const EBFieldGrid& f, std::size_t step) {
    old_[0]=f.ex(); old_[1]=f.ey(); old_[2]=f.ez();
    step_=step; inletCaptured_=false;
  }
  static void beforeInlet(void* context, const EBFieldGrid& f) {
    static_cast<GridWorkAudit*>(context)->captureInlet(f);
  }
  void complete(const EBFieldGrid& f) {
    if (rank_==0 && !inletCaptured_)
      throw std::runtime_error("Missing pre-overwrite electric-plane sample");
    const YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};
    const YeeComponent* current[3]={&f.jx(),&f.jy(),&f.jz()};
    const auto& g=f.geometry();
    long double no=0.,nn=0.,co=0.,cn=0.,di=0.,dc=0.;
    // Parallel summation affects diagnostics only; no field or source is written.
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for reduction(+:no,nn,co,cn,di,dc) schedule(static)
#endif
    for (std::ptrdiff_t kk=0;kk<static_cast<std::ptrdiff_t>(g.nz);++kk) {
      const auto k=static_cast<std::size_t>(kk);
      for (std::size_t j=cy_;j<g.ny-cy_;++j)
        for (std::size_t i=cx_;i<g.nx-cx_;++i)
          for (unsigned a=0;a<3;++a) {
            long double js=0.,os=0.,ns=0.,ds=0.;
            for (unsigned p=0;p<2;++p) for (unsigned q=0;q<2;++q) {
              const std::size_t x=i+(a==0?0:p);
              const std::size_t y=j+(a==0?p:(a==1?0:q));
              const std::size_t z=k+(a==2?0:q);
              const long double vj=(*current[a])(x,y,z);
              const long double vo=old_[a](x,y,z);
              const long double actual=(*e[a])(x,y,z);
              const long double vn=(rank_==0 && a<2 && z==0) ?
                inlet_[a][y*e[a]->nx()+x] : actual;
              no+=.25L*vj*vo; nn+=.25L*vj*vn;
              di+=.125L*vj*(actual-vn);
              js+=vj; os+=vo; ns+=vn; ds+=actual-vn;
            }
            co+=js*os/16.L; cn+=js*ns/16.L; dc+=js*ds/32.L;
          }
    }
    const long double scale=static_cast<long double>(g.dt)*g.dx*g.dy*g.dz;
    output_ << step_+1 << ',' << no*scale << ',' << (no+nn)*.5L*scale
      << ',' << nn*scale << ',' << co*scale << ',' << (co+cn)*.5L*scale
      << ',' << cn*scale << ',' << di*scale << ',' << dc*scale << '\n';
    if (!output_) throw std::runtime_error("Cannot write grid-work audit");
  }
  void close() {
    output_.flush();
    if (!output_) throw std::runtime_error("Cannot flush grid-work audit");
    output_.close();
  }
private:
  void captureInlet(const EBFieldGrid& f) {
    if (rank_!=0) return;
    if (inletCaptured_) throw std::runtime_error("Duplicate pre-inlet sample");
    const YeeComponent* e[2]={&f.ex(),&f.ey()};
    for (unsigned a=0;a<2;++a)
      inlet_[a].assign(e[a]->data(),e[a]->data()+e[a]->nx()*e[a]->ny());
    inletCaptured_=true;
  }
  int rank_;
  std::size_t cx_,cy_,step_=0;
  bool inletCaptured_=false;
  YeeComponent old_[3];
  std::vector<Double> inlet_[2];
  std::ofstream output_;
};
}
#endif
