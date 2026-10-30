#include "blocked_cowan_reference.h"
#include "eb_cpml.h"
#include "compute_threads.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace aprl;
using aprl::verification::BlockedCowanReference;
struct BoundaryProbe : EBConvolutionalPML {
  std::ofstream& out;BlockedCowanReference& reference;
  std::size_t offset,k;double origin,beta;
  BoundaryProbe(const EBGridGeometry& g,std::size_t nz,std::size_t start,
      const EBCPMLParameters& p,std::ofstream& stream,BlockedCowanReference& ref,double x0,double speed)
    :EBConvolutionalPML(g,nz,start,p,false,false),out(stream),reference(ref),offset(start),k(8),origin(x0),beta(speed) {}
  double translated(unsigned a,std::size_t i,std::size_t j,double z) {
    const auto lo=static_cast<std::size_t>(std::floor(z));const double f=z-lo;
    return (1-f)*reference.field(a,i,j,lo)+f*reference.field(a,i,j,lo+1);
  }
  void sample(const char* stage,const EBFieldGrid& f,bool endElectric,bool endMagnetic) {
    const auto& g=f.geometry();const auto j=g.ny/2;
    for(std::size_t i=0;i<=g.nx;++i) {
      const double by=i<g.nx?f.by()(i,j,k):0.;
      const double bzref=i<g.nx?translated(4,i,j,offset+k-(endMagnetic?beta:0.)):0.;
      out<<stage<<','<<i<<','<<origin+i*g.dx<<','<<f.ez()(i,j,k)<<','
        <<translated(2,i,j,offset+k-(endElectric?beta:0.))<<','<<by<<','<<bzref<<'\n';
    }
  }
  void afterMagneticUpdate(EBFieldGrid& f,const EBElectricHaloView& h) override {
    sample("magnetic_curl",f,false,true);
    EBConvolutionalPML::afterMagneticUpdate(f,h);
    sample("magnetic_cpml",f,false,true);
  }
  void afterElectricUpdate(EBFieldGrid& f) override {
    sample("electric_curl",f,true,true);
    EBConvolutionalPML::afterElectricUpdate(f);
    sample("electric_cpml",f,true,true);
  }
};
int main(int argc,char** argv) {
  int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  try {
    if(argc!=8 && argc!=9)throw std::runtime_error("reference-dir output nx cpml beta reference-nz reference-origin [reference-guard] required");
    configureComputeThreads(provided,8);
    const std::size_t nx=std::stoul(argv[3]),guard=std::stoul(argv[4]),nz=std::stoul(argv[6]);
    const double beta=std::stod(argv[5]),z0=std::stod(argv[7]),dz=1.2543262415359034e-6;
    const EBGridGeometry global(nx,nx,nz,5.e-6,5.e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
    FieldVector<Double> origin(0.);origin[0]=origin[1]=-.5*nx*5.e-6;origin[2]=z0;
    const std::size_t referenceGuard=argc==9?std::stoul(argv[8]):guard;
    const std::size_t guards[3]={referenceGuard,referenceGuard,0},start=nz/2-8;
    auto reference=BlockedCowanReference::openPrepared(global,origin,beta,guards,32,argv[1],MPI_COMM_SELF);
    const EBGridGeometry local(nx,nx,18,global.dx,global.dy,global.dz,global.dt,global.solver);
    EBFieldGrid f(local);reference->fill(f,start);
    EBCPMLParameters p;p.cells[0]=p.cells[1]=guard;p.cells[2]=0;
    p.polynomialOrder=3.;p.targetReflection=1.e-8;p.kappaMax=8.;p.alphaFraction=0.;
    std::ofstream out(argv[2]);if(!out)throw std::runtime_error("Cannot create probe CSV");
    out<<std::setprecision(21)<<"stage,i,x_m,Ez_V_per_m,translated_reference_Ez_V_per_m,By_T,translated_reference_By_T\n";
    auto boundary=std::unique_ptr<BoundaryProbe>(new BoundaryProbe(local,nz,start,p,out,*reference,origin[0],beta));
    auto* observer=boundary.get();f.setBoundary(std::move(boundary));
    observer->sample("initial",f,false,false);f.clearCurrent();f.advanceMagnetic();f.advanceElectric();
    // Only |x|>=150um is used: this region has no source particles/current.
    // The central z sample is 8 cells from an artificial slab end, so one
    // original curl step cannot sample those ends. No trajectory is advanced.
    out.flush();if(!out)throw std::runtime_error("Cannot finish probe CSV");
    std::cout<<"FIRST_STEP_BOUNDARY_PROBE_COMPLETE outside-source-only\n";
  }catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
  MPI_Finalize();return 0;
}
