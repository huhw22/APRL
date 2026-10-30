#include "eb_cpml.h"
#include "compute_threads.h"
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
using namespace aprl;

// A normal-x plane mode, periodic in y/z, using the original 3D Cowan curl
// and CPML. The extra boundary copies enforce this 1D symmetry only; the x
// absorber and x termination are never modified.
struct PlaneBoundary : EBBoundaryOperator {
  std::unique_ptr<EBBoundaryOperator> absorber;
  PlaneBoundary(const EBGridGeometry& g,bool cpml) {
    if(cpml) {
      EBCPMLParameters p;p.cells[0]=10;p.cells[1]=p.cells[2]=0;
      absorber.reset(new EBConvolutionalPML(g,g.nz,0,p,false,false));
    } else absorber.reset(new PerfectElectricConductorBoundary(false,false));
  }
  void afterMagneticUpdate(EBFieldGrid& f,const EBElectricHaloView& h) override {
    absorber->afterMagneticUpdate(f,h);
  }
  std::size_t memoryBytes() const override {return absorber->memoryBytes();}
  void afterElectricUpdate(EBFieldGrid& f) override {
    absorber->afterElectricUpdate(f);
    for(YeeComponent* a:{&f.ex(),&f.ey(),&f.ez()})
      for(std::size_t i=0;i<a->nx();++i) {
        const double v=(*a)(i,a->ny()/2,a->nz()/2);
        for(std::size_t k=0;k<a->nz();++k)for(std::size_t j=0;j<a->ny();++j)(*a)(i,j,k)=v;
      }
  }
};
std::vector<double> run(bool cpml,bool background,bool wave,double& maxStaticError) {
  const double dz=1.2543262415359034e-6;
  EBGridGeometry g(160,2,8,5.e-6,5.e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
  EBFieldGrid f(g);f.setBoundary(std::unique_ptr<EBBoundaryOperator>(new PlaneBoundary(g,cpml)));
  const auto pulse=[](double x){const double r=x-45.;return 1.e6*std::exp(-r*r/(2.*8.*8.))*std::cos(2.*std::acos(-1.)*r/16.);};
  std::vector<double> ex0(f.ex().nx());
  for(std::size_t i=0;i<f.ex().nx();++i) {
    ex0[i]=background?2.e6*std::sin(std::acos(-1.)*(i+.5)/g.nx):0.;
    for(std::size_t k=0;k<f.ex().nz();++k)for(std::size_t j=0;j<f.ex().ny();++j)f.ex()(i,j,k)=ex0[i];
  }
  for(std::size_t i=0;i<f.ey().nx();++i)
    for(std::size_t k=0;k<f.ey().nz();++k)for(std::size_t j=0;j<f.ey().ny();++j)f.ey()(i,j,k)=wave?pulse(i):0.;
  for(std::size_t i=0;i<f.bz().nx();++i)
    for(std::size_t k=0;k<f.bz().nz();++k)for(std::size_t j=0;j<f.bz().ny();++j)
      f.bz()(i,j,k)=wave?pulse(i+.5+.5*SI::c*g.dt/g.dx)/SI::c:0.;
  std::vector<double> trace;
  for(unsigned n=0;n<1300;++n) {
    EBElectricHaloView h;h.lowerEx=h.upperEx=f.ex().data();h.lowerEy=h.upperEy=f.ey().data();h.lowerEz=h.upperEz=f.ez().data();
    f.advanceMagnetic(h);f.advanceElectric();trace.push_back(f.ey()(75,0,4));
    for(std::size_t i=0;i<f.ex().nx();++i)maxStaticError=std::max(maxStaticError,std::abs(f.ex()(i,1,4)-ex0[i]));
  }
  return trace;
}
int main(int argc,char** argv) {
  int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  try {
    if(argc!=2)throw std::runtime_error("output directory required");
    configureComputeThreads(provided,1);
    double staticError=0.;auto pec=run(false,false,true,staticError),pml=run(true,false,true,staticError);
    auto bg=run(true,true,false,staticError),total=run(true,true,true,staticError);
    double incident=0.,reflected=0.,pecIncident=0.,pecReflected=0.,linear=0.;
    std::ofstream csv(std::string(argv[1])+"/wave_trace.csv");csv<<std::setprecision(21)<<"step,PEC_Ey,CPML_Ey,background_Ey,total_Ey\n";
    for(std::size_t n=0;n<pec.size();++n) {
      if(n<300){incident=std::max(incident,std::abs(pml[n]));pecIncident=std::max(pecIncident,std::abs(pec[n]));}
      if(n>550){reflected=std::max(reflected,std::abs(pml[n]));pecReflected=std::max(pecReflected,std::abs(pec[n]));}
      linear=std::max(linear,std::abs(total[n]-bg[n]-pml[n]));
      csv<<n<<','<<pec[n]<<','<<pml[n]<<','<<bg[n]<<','<<total[n]<<'\n';
    }
    const bool pass=incident>5.e5&&pecReflected/pecIncident>.5&&reflected/incident<.05&&linear<1.e-6&&staticError<1.e-6;
    std::ofstream out(std::string(argv[1])+"/WAVE_CHECK.json");out<<std::setprecision(21)
      <<"{\"checks_pass\":"<<(pass?"true":"false")<<",\"CPML_reflected_amplitude_ratio\":"<<reflected/incident
      <<",\"PEC_reflected_amplitude_ratio\":"<<pecReflected/pecIncident<<",\"linear_superposition_error_V_per_m\":"<<linear
      <<",\"static_background_max_error_V_per_m\":"<<staticError<<",\"scope\":\"normal-x plane pulse, production dx/dy/dz/dt, Cowan-z and ten-cell CPML; transverse periodic plane reduction; no oblique-wave or beam-dynamics acceptance\"}\n";
    if(!pass)throw std::runtime_error("Wave/background absorption check failed");
    std::cout<<"WAVE_ABSORPTION_PASS ratio="<<reflected/incident<<'\n';
  }catch(const std::exception& ex){std::cerr<<ex.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
  MPI_Finalize();return 0;
}
