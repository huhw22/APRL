#ifndef APRL_PARTICLE_WORK_GEOMETRY_H
#define APRL_PARTICLE_WORK_GEOMETRY_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <stdexcept>
#include "eb_particles.h"

namespace aprl {

// Independent diagnostic quadrature of the Whitney/trajectory interpolant
// dual to the existing first-order charge-conserving current deposit. Two
// Gauss points integrate its quadratic restriction to each straight cell leg.
// No J, E, particle state, CPML state or ownership data is written.
inline long double trajectoryElectricWork(const EBFieldGrid& f,
    const FieldVector<Double>& origin,const FieldVector<Double>& begin,
    const FieldVector<Double>& end,Double charge) {
  const auto& g=f.geometry();
  const Double spacing[3]={g.dx,g.dy,g.dz};
  Double x[3],d[3],cuts[32]={0.,1.};unsigned count=2;
  for(unsigned a=0;a<3;++a) {
    x[a]=(begin[a]-origin[a])/spacing[a];d[a]=(end[a]-begin[a])/spacing[a];
    if(d[a]!=0.) {
      const Double lo=std::min(x[a],x[a]+d[a]),hi=std::max(x[a],x[a]+d[a]);
      for(long long k=static_cast<long long>(std::floor(lo))+1;k<hi;++k) {
        if(count==32)throw std::runtime_error("Diagnostic trajectory crosses too many cells");
        const Double t=(k-x[a])/d[a];if(t>0. && t<1.)cuts[count++]=t;
      }
    }
  }
  std::sort(cuts,cuts+count);long double result=0.;
  for(unsigned n=1;n<count;++n) {
    if(cuts[n]<=cuts[n-1])continue;
    const Double middle=.5*(cuts[n]+cuts[n-1]),half=.5*(cuts[n]-cuts[n-1]);
    std::size_t cell[3];
    const std::size_t sizes[3]={g.nx,g.ny,g.nz};
    for(unsigned a=0;a<3;++a) {
      const long long k=static_cast<long long>(std::floor(x[a]+middle*d[a]));
      if(k<0 || k>=static_cast<long long>(sizes[a]))
        throw std::out_of_range("Diagnostic trajectory leaves its native field slab");
      cell[a]=static_cast<std::size_t>(k);
    }
    for(int sign:{-1,1}) {
      const Double t=middle+sign*half/std::sqrt(3.);Double w[3][2];
      for(unsigned a=0;a<3;++a) {w[a][1]=std::max(0.,std::min(1.,x[a]+t*d[a]-cell[a]));w[a][0]=1.-w[a][1];}
      long double electric[3]={};
      for(unsigned a=0;a<2;++a)for(unsigned b=0;b<2;++b) {
        electric[0]+=static_cast<long double>(w[1][a])*w[2][b]*f.ex()(cell[0],cell[1]+a,cell[2]+b);
        electric[1]+=static_cast<long double>(w[0][a])*w[2][b]*f.ey()(cell[0]+a,cell[1],cell[2]+b);
        electric[2]+=static_cast<long double>(w[0][a])*w[1][b]*f.ez()(cell[0]+a,cell[1]+b,cell[2]);
      }
      for(unsigned a=0;a<3;++a)result+=static_cast<long double>(charge)*
        (static_cast<long double>(end[a])-begin[a])*electric[a]*half;
    }
  }
  return result;
}

struct ParticleWorkGeometry {
  const EBFieldGrid* fields=NULL;
  FieldVector<Double> origin,first,last,electric;
  std::size_t cx=0,cy=0;
  bool completeActive=false,interior=true;
  unsigned substeps=0;
  long double boris=0.,gatherPath=0.,subpath=0.,chord=0.;
  Double charge=0.,maximumGamma=1.;
  void reset(const EBFieldGrid* f,const FieldVector<Double>& o,
             std::size_t x,std::size_t y,bool active) {
    fields=f;origin=o;cx=x;cy=y;completeActive=active;interior=true;
    substeps=0;boris=gatherPath=subpath=chord=0.;maximumGamma=1.;
  }
  bool contains(const FieldVector<Double>& p) const {
    const auto& g=fields->geometry();
    const Double v[3]={(p[0]-origin[0])/g.dx,(p[1]-origin[1])/g.dy,(p[2]-origin[2])/g.dz};
    return v[0]>cx+1. && v[0]<g.nx-cx-1. &&
      v[1]>cy+1. && v[1]<g.ny-cy-1. && v[2]>1. && v[2]<g.nz-1.;
  }
  void observe(const RelativisticParticleSI& before,const RelativisticParticleSI& after,
               const BorisWork* work) {
    if(!fields)return;
    if(!work)throw std::runtime_error("Geometry audit requires independently recorded Boris work");
    if(substeps==0) {
      first=before.position;charge=before.charge;
      FieldVector<Double> magnetic(0.);
      fields->sampleFieldsPosition(first,origin,electric,magnetic);
    }
    ++substeps;last=after.position;boris+=work->gridWork;
    maximumGamma=std::max(maximumGamma,std::max(
      BoostFrameTransform::gammaFromProperVelocity(before.properVelocity),
      BoostFrameTransform::gammaFromProperVelocity(after.properVelocity)));
    interior=interior && contains(before.position) && contains(after.position);
    for(unsigned a=0;a<3;++a)gatherPath+=static_cast<long double>(charge)*electric[a]*
      (static_cast<long double>(after.position[a])-before.position[a]);
    if(completeActive && interior)subpath+=trajectoryElectricWork(*fields,origin,before.position,after.position,charge);
  }
  void finish() {
    if(substeps && completeActive && interior)
      chord=trajectoryElectricWork(*fields,origin,first,last,charge);
  }
};

// Rank-local output: no added MPI collectives. Only complete active pushes
// whose entire orbit stays one cell inside this rank's physical interior are
// compared. Excluded populations and their absolute work remain explicit.
class ParticleWorkGeometryWriter {
public:
  ParticleWorkGeometryWriter(const std::string& directory,int rank):steps_(0) {
    std::ostringstream name;name<<directory<<"/particle-geometry-rank-"<<std::setw(5)<<std::setfill('0')<<rank<<".csv";
    path_=name.str();std::ifstream prior(path_);if(prior.good())throw std::runtime_error("Geometry audit output exists");
    out_.open(path_);if(!out_)throw std::runtime_error("Cannot open geometry audit output");
    out_<<std::setprecision(21)<<"step,observed_pushes,eligible_pushes,excluded_release_pushes,excluded_boundary_interface_pushes,eligible_boris_grid_J,eligible_gather_displacement_J,eligible_subpath_old_mesh_J,eligible_chord_old_mesh_J,excluded_push_boris_grid_J,excluded_abs_push_boris_grid_J,max_observed_gamma_box\n";
    beginStep();
  }
  void beginStep() {observed_=eligible_=release_=boundary_=0;boris_=gather_=subpath_=chord_=excluded_=excludedAbs_=0.;maxGamma_=1.;}
  void add(const ParticleWorkGeometry& s) {
    if(!s.substeps)return;
    ++observed_;maxGamma_=std::max(maxGamma_,s.maximumGamma);
    if(!s.completeActive || !s.interior) {
      if(!s.completeActive)++release_;else ++boundary_;
      excluded_+=s.boris;excludedAbs_+=std::abs(s.boris);return;
    }
    ++eligible_;boris_+=s.boris;gather_+=s.gatherPath;subpath_+=s.subpath;chord_+=s.chord;
  }
  void finishStep(std::size_t step) {
    steps_=step;
    out_<<step<<','<<observed_<<','<<eligible_<<','<<release_<<','<<boundary_<<','<<boris_<<','<<gather_<<','<<subpath_<<','<<chord_<<','<<excluded_<<','<<excludedAbs_<<','<<maxGamma_<<'\n';
    if(!out_)throw std::runtime_error("Cannot write geometry audit");
  }
  void close() {
    out_.flush();if(!out_)throw std::runtime_error("Cannot flush geometry audit");out_.close();
    std::ofstream marker(path_.substr(0,path_.size()-4)+".complete");marker<<"1 "<<steps_<<'\n';
    if(!marker)throw std::runtime_error("Cannot close geometry audit marker");
  }
private:
  std::ofstream out_;std::string path_;std::size_t steps_;
  unsigned long long observed_,eligible_,release_,boundary_;
  long double boris_,gather_,subpath_,chord_,excluded_,excludedAbs_;Double maxGamma_;
};
}
#endif
