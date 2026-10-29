#include "discrete_field_ledger.h"
#include "moving_window_fields.h"
#include "eb_mpi.h"
#include "compute_threads.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <sys/stat.h>

namespace {
using namespace aprl;
Double zero(void*,unsigned,std::size_t,std::size_t,std::ptrdiff_t,std::size_t){return 0.;}
void require(bool b,const char* s){if(!b)throw std::runtime_error(s);}
double pulse(double k){return k>16.&&k<24.?std::pow(std::sin(std::acos(-1.)*(k-16.)/8.),2):0.;}
}
int main(int argc,char** argv) {
  int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  try {
    require(argc==3,"directory and manufactured/pulse/injected required");
    configureComputeThreads(provided);
    const std::string directory=argv[1],mode=argv[2];
    const bool isPulse=mode=="pulse",bad=mode=="injected";
    require(isPulse||bad||mode=="manufactured","Unknown case");
    int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    if(rank==0)require(::mkdir(directory.c_str(),0777)==0,"Output directory exists");
    MPI_Barrier(MPI_COMM_WORLD); // Test-only setup; external launch is bounded.
    const std::size_t globalNz=isPulse?64:65,nx=isPulse?160:24,ny=isPulse?160:26;
    const std::size_t offset=(globalNz/size)*rank+std::min<std::size_t>(rank,globalNz%size);
    const std::size_t nz=globalNz/size+(static_cast<std::size_t>(rank)<globalNz%size);
    const double dz=1.2543262415359034e-6;
    const EBGridGeometry g(nx,ny,nz,5.e-6,isPulse?5.e-6:6.e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
    EBFieldGrid f(g);EBCPMLParameters p;p.cells[0]=p.cells[1]=isPulse?10:3;p.cells[2]=0;
    MovingWindowReference reference;reference.value=&zero;
    MovingWindowFields moving(g,globalNz,offset,MPI_COMM_WORLD,p,reference);moving.installBoundary(f);
    EBZSlabHaloExchange halo(MPI_COMM_WORLD);
    const std::size_t aperture=isPulse?72:3;
    DiscreteFieldLedger ledger(directory,g,globalNz,offset,aperture,aperture,MPI_COMM_WORLD);
    moving.setDiscreteObserver(&DiscreteFieldLedger::observe,&ledger);
    YeeComponent* fields[6]={&f.ex(),&f.ey(),&f.ez(),&f.bx(),&f.by(),&f.bz()};
    for(unsigned a=0;a<6;++a) {
      auto& c=*fields[a];
      for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)for(std::size_t i=0;i<c.nx();++i) {
        const double z=static_cast<double>(offset+k);
        c(i,j,k)=isPulse ? (a==0?1.e6*pulse(z):(a==4?1.e6/SI::c*pulse(z+1):0.)) :
          (a<3?1.e3:1.e-5)*std::sin(.47*i+.38*j+.21*z+.31*a);
      }
    }
    long double maximumNormalized=0.,sumCowan=0.,maxPrediction=0.,maxTranslation=0.,maxFlux=0.;
    const unsigned steps=isPulse?28:(bad?1:4);
    for(unsigned step=0;step<steps;++step) {
      f.clearCurrent();ledger.begin(f,step);moving.beginStep(f,step,-static_cast<std::ptrdiff_t>(step));
      halo.advanceMagnetic(f,moving.physicalElectricHalos(f,step));
      long double qAnalytic=0.;
      if(isPulse) {
        const long double area=static_cast<long double>(g.dx)*g.dy;
        for(std::size_t j=72;j<88;++j)for(std::size_t i=72;i<88;++i) {
          if(rank==size-1)qAnalytic+=area*g.dt/SI::mu0*f.ex()(i,j,nz)*f.by()(i,j,nz-1);
          if(rank==0)qAnalytic-=area*g.dt/SI::mu0*f.ex()(i,j,0)*f.by()(i,j,0);
        }
      }else {
        YeeComponent* current[3]={&f.jx(),&f.jy(),&f.jz()};
        for(unsigned a=0;a<3;++a) {
          auto& c=*current[a];
          for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)for(std::size_t i=0;i<c.nx();++i) {
            double w=1.;if(a<2&&((k==0&&rank>0)||(k==nz&&rank<size-1)))w=.5;
            c(i,j,k)=w*100.*std::cos(.23*i+.19*j+.17*(static_cast<double>(offset+k)-step)+a);
          }
        }
      }
      halo.advanceElectric(f,EBCurrentFilter::None);
      if(bad&&offset<=32&&32<offset+nz)f.ex()(nx/2,ny/2,32-offset)+=1.;
      ledger.finishElectric(f);moving.shift(f,step+1,-static_cast<std::ptrdiff_t>(step)-1);ledger.afterRemap(f);
      const auto& r=ledger.lastRecord();
      const long double scale=std::abs(r.begin)+std::abs(r.end)+std::abs(r.exchange)+std::abs(r.jMid)+std::abs(r.timeE)+std::abs(r.timeB)+std::abs(r.cowan)+std::abs(r.inlet)+std::abs(r.remap);
      long double values[4]={r.residual,scale,r.cowan,r.exchange-qAnalytic},global[4];
      MPI_Allreduce(values,global,4,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      maximumNormalized=std::max(maximumNormalized,std::abs(global[0])/global[1]);sumCowan+=std::abs(global[2]);
      maxFlux=std::max(maxFlux,std::abs(global[3])/global[1]);
      maxPrediction=std::max(maxPrediction,r.eScale>0?r.eError/r.eScale:0.L);
      if(!bad) {
        require(std::abs(global[0])<512*std::numeric_limits<double>::epsilon()*global[1],"Operator energy identity failed");
        require(r.eError<512*std::numeric_limits<double>::epsilon()*std::max(r.eScale,1.e-100L),"Electric stencil prediction failed");
        require(r.bError<512*std::numeric_limits<double>::epsilon()*std::max(r.bScale,1.e-100L),"Magnetic stencil prediction failed");
      }
      if(isPulse)for(std::size_t k=0;k<=nz;++k) {
        const double z=static_cast<double>(offset+k)-2*(step+1);
        maxTranslation=std::max(maxTranslation,static_cast<long double>(std::abs(f.ex()(80,80,k)/1.e6-pulse(z))));
        if(k<nz)maxTranslation=std::max(maxTranslation,static_cast<long double>(std::abs(f.by()(80,80,k)*SI::c/1.e6-pulse(z+1))));
      }
    }
    long double local[3]={maxPrediction,maxTranslation,sumCowan},global[3];MPI_Allreduce(local,global,3,MPI_LONG_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(bad)require(maximumNormalized>1.e-10&&global[0]>1.e-7,"Injected field error was hidden");
    else if(isPulse)require(global[1]<1.e-9&&maxFlux<1.e-10,"Pulse/analytic boundary-flux check failed");
    else require(global[2]>1.e-30,"Manufactured 3D case did not exercise Cowan operator term");
    ledger.close(true);
    if(rank==0) {
      std::ofstream summary(directory+"/summary.json");summary<<std::setprecision(21)
        <<"{\"pass\":true,\"mode\":\""<<mode<<"\",\"ranks\":"<<size<<",\"steps\":"<<steps
        <<",\"max_normalized_operator_residual\":"<<maximumNormalized<<",\"max_electric_prediction_relative_error\":"<<global[0]
        <<",\"maximum_translation_relative_error\":"<<global[1]<<",\"maximum_analytic_flux_relative_error\":"<<maxFlux
        <<",\"sum_abs_Cowan_term_J\":"<<global[2]<<"}\n";
      std::cout<<"PASS "<<mode<<" ranks="<<size<<" max_normalized_residual="<<maximumNormalized<<'\n';
    }
    MPI_Finalize();return 0;
  }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
}
