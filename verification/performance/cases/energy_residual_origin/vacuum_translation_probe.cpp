// Particle-free test of the existing Cowan / moving-cut implementation.
// The central aperture contains a 1-D right-going pulse. Its exact discrete
// solution shifts by two cells per complete moving-window step. Transverse
// CPML is far enough away that its numerical stencil cannot reach this
// aperture during these 28 steps. No production operator is substituted.
#include "moving_window_fields.h"
#include "compute_threads.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <mpi.h>

namespace {
using namespace aprl;
const std::size_t lo=72,hi=88,nz=64;
const double amplitude=1.e6;
Double zeroReference(void*,unsigned,std::size_t,std::size_t,std::ptrdiff_t,std::size_t) {return 0.;}
double profile(double cell,int start,int width) {
  const double x=(cell-start)/width;
  return x>0.&&x<1. ? std::pow(std::sin(std::acos(-1.)*x),2) : 0.;
}
long double energy(const EBFieldGrid& f) {
  const auto& g=f.geometry();
  const long double volume=static_cast<long double>(g.dx)*g.dy*g.dz;
  long double u=0.;
  for(std::size_t k=0;k<g.nz;++k)
    for(std::size_t j=lo;j<hi;++j)for(std::size_t i=lo;i<hi;++i)
      u+=static_cast<long double>(f.fieldEnergyDensityCell(i,j,k))*volume;
  return u;
}
long double power(const EBFieldGrid& f) {
  const auto& g=f.geometry();
  const long double xy=static_cast<long double>(g.dx)*g.dy;
  const long double xz=static_cast<long double>(g.dx)*g.dz,yz=static_cast<long double>(g.dy)*g.dz;
  long double p=0.;
  for(std::size_t j=lo;j<hi;++j)for(std::size_t i=lo;i<hi;++i)
    p+=xy*(static_cast<long double>(f.poyntingComponentCell(i,j,g.nz-1,2))-f.poyntingComponentCell(i,j,0,2));
  for(std::size_t k=0;k<g.nz;++k) {
    for(std::size_t j=lo;j<hi;++j)
      p+=yz*(static_cast<long double>(f.poyntingComponentCell(hi-1,j,k,0))-f.poyntingComponentCell(lo,j,k,0));
    for(std::size_t i=lo;i<hi;++i)
      p+=xz*(static_cast<long double>(f.poyntingComponentCell(i,hi-1,k,1))-f.poyntingComponentCell(i,lo,k,1));
  }
  return p;
}
double translationError(const EBFieldGrid& f,int step,int start,int width) {
  double error=0.;
  for(std::size_t k=0;k<=nz;++k) {
    const double expected=amplitude*profile(static_cast<double>(k)-2*step,start,width);
    error=std::max(error,std::abs(f.ex()(80,80,k)-expected)/amplitude);
    error=std::max(error,std::abs(f.ey()(80,80,k))/amplitude);
    if(k<nz) {
      const double expectedB=amplitude/SI::c*profile(static_cast<double>(k)+1-2*step,start,width);
      error=std::max(error,std::abs(f.by()(80,80,k)-expectedB)*SI::c/amplitude);
      error=std::max(error,std::abs(f.bx()(80,80,k))*SI::c/amplitude);
      error=std::max(error,std::abs(f.ez()(80,80,k))/amplitude);
    }
    error=std::max(error,std::abs(f.bz()(80,80,k))*SI::c/amplitude);
  }
  return error;
}
void run(const std::string& directory,int start,int width) {
  const double dz=1.2543262415359034e-6;
  const EBGridGeometry g(160,160,nz,5.e-6,5.e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
  EBFieldGrid f(g);EBCPMLParameters p;p.cells[0]=p.cells[1]=10;p.cells[2]=0;
  MovingWindowReference reference;reference.value=&zeroReference;
  MovingWindowFields moving(g,nz,0,MPI_COMM_SELF,p,reference);moving.installBoundary(f);
  for(std::size_t k=0;k<f.ex().nz();++k)
    for(std::size_t j=0;j<f.ex().ny();++j)for(std::size_t i=0;i<f.ex().nx();++i)
      f.ex()(i,j,k)=amplitude*profile(k,start,width);
  for(std::size_t k=0;k<f.by().nz();++k)
    for(std::size_t j=0;j<f.by().ny();++j)for(std::size_t i=0;i<f.by().nx();++i)
      f.by()(i,j,k)=amplitude/SI::c*profile(k+1,start,width);
  std::ostringstream filename;filename<<directory<<"/pulse-start-"<<start<<"-width-"<<width<<".csv";
  std::ofstream out(filename.str());if(!out)throw std::runtime_error("Cannot create pulse trace");
  out<<std::setprecision(21)<<"step,field_J,outward_J,inlet_J,remap_J,residual_J,max_translation_relative_error\n";
  const long double initial=energy(f);long double q=0.,inlet=0.,remap=0.;double maxError=0.;
  out<<0<<','<<initial<<",0,0,0,0,0\n";
  for(int step=0;step<28;++step) {
    const long double oldPower=power(f);
    moving.beginStep(f,step,-step);
    f.advanceMagnetic(moving.physicalElectricHalos(f,step));
    f.advanceElectric();
    // Ex/Ey at the inlet remain exactly zero for this compact outgoing pulse.
    // The before/after energy of the inlet overwrite is therefore zero.
    q+=.5L*g.dt*(oldPower+power(f));
    const long double beforeShift=energy(f);
    moving.shift(f,step+1,-step-1);
    const long double afterShift=energy(f);remap+=afterShift-beforeShift;
    maxError=std::max(maxError,translationError(f,step+1,start,width));
    out<<step+1<<','<<afterShift<<','<<q<<','<<inlet<<','<<remap<<','<<afterShift-initial+q-inlet-remap<<','<<maxError<<'\n';
  }
  if(maxError>1.e-9)throw std::runtime_error("Central field does not match exact two-cell translation");
  const long double residual=energy(f)-initial+q-inlet-remap;
  std::cout<<std::setprecision(21)<<"{\"type\":\"result\",\"start\":"<<start<<",\"width_cells\":"<<width
    <<",\"initial_energy_J\":"<<initial<<",\"final_energy_J\":"<<energy(f)<<",\"residual_J\":"<<residual
    <<",\"residual_over_initial_energy\":"<<residual/initial<<",\"maximum_translation_relative_error\":"<<maxError
    <<",\"steps\":28,\"particles\":0,\"current\":0,\"pass\":true}\n";
}
}
int main(int argc,char** argv) {
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  try {
    if(argc!=2)throw std::runtime_error("Expected a new output directory");
    configureComputeThreads(provided);
    run(argv[1],16,8);run(argv[1],48,8);run(argv[1],8,16);run(argv[1],8,32);
    MPI_Finalize();return 0;
  }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
}
