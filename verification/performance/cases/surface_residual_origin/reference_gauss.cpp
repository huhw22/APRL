#include "blocked_cowan_reference.h"
#include "compute_threads.h"
#include <mpi.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace aprl;
int main(int argc,char** argv) {
  int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  try {
    if(argc!=7)throw std::runtime_error("reference output solver-beta translation-beta initial-offset steps required");
    configureComputeThreads(provided,8);const double dz=1.2543262415359034e-6;
    EBGridGeometry g(160,160,3328,5.e-6,5.e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
    FieldVector<Double> origin(0.);origin[0]=origin[1]=-.0004;origin[2]=-.0020871988659157434;
    const std::size_t guards[3]={10,10,0};const double beta=std::stod(argv[3]),translation=std::stod(argv[4]);
    auto reference=verification::BlockedCowanReference::openPrepared(g,origin,beta,guards,32,argv[1],MPI_COMM_SELF);
    std::ifstream rho(std::string(argv[1])+"/rho.bin",std::ios::binary);if(!rho)throw std::runtime_error("rho missing");
    std::ofstream out(argv[2]);out<<std::setprecision(21)<<"step,native_k,initial_gauss_signed_C,reference_divE_charge_C,rigid_rho_charge_C,initial_rho_charge_C\n";
    auto charge=[&](std::size_t k){long double sum=0.;std::vector<double> row(49);for(std::size_t j=56;j<=104;++j){rho.seekg(((k-1)*139*139+(j-11)*139+56-11)*sizeof(double));rho.read(reinterpret_cast<char*>(row.data()),row.size()*sizeof(double));if(!rho)throw std::runtime_error("rho short read");for(double x:row)sum+=x;}return sum*g.dx*g.dy*g.dz;};
    const long double volume=static_cast<long double>(g.dx)*g.dy*g.dz;
    const auto steps=std::stoul(argv[6]),offset=std::stoul(argv[5]);
    const auto first=static_cast<std::size_t>(std::floor(offset+1-(1+translation)*(steps-1)))-1,last=offset+1;
    std::vector<long double> divcharge(last-first+1),rhocharge(last-first+1);
    long double previousEz=0.;
    for(std::size_t k=first-1;k<=last;++k) {
      long double ez=0.,transverse=0.;
      for(std::size_t j=56;j<=104;++j)for(std::size_t i=56;i<=104;++i)ez+=reference->field(2,i,j,k);
      if(k>=first) {
        for(std::size_t j=56;j<=104;++j)transverse+=(static_cast<long double>(reference->field(0,104,j,k))-reference->field(0,55,j,k))/g.dx;
        for(std::size_t i=56;i<=104;++i)transverse+=(static_cast<long double>(reference->field(1,i,104,k))-reference->field(1,i,55,k))/g.dy;
        divcharge[k-first]=SI::epsilon0*volume*(transverse+(ez-previousEz)/g.dz);rhocharge[k-first]=charge(k);
      }
      previousEz=ez;
    }
    for(std::size_t n=0;n<steps;++n) {
      const auto k=offset-n+1;const double z=k-translation*n;const auto lo=static_cast<std::size_t>(std::floor(z));const double f=z-lo;
      const auto initial=divcharge[k-first]-rhocharge[k-first];const auto hi=std::min(lo+1,last);
      const auto advected=(1-f)*divcharge[lo-first]+f*divcharge[hi-first];
      const auto qref=(1-f)*rhocharge[lo-first]+f*rhocharge[hi-first];
      out<<n<<','<<k<<','<<initial<<','<<advected<<','<<qref<<','<<rhocharge[k-first]<<'\n';
    }
    out.flush();if(!out)throw std::runtime_error("Gauss reference output failed");std::cout<<"REFERENCE_GAUSS_COMPLETE\n";
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
  MPI_Finalize();return 0;
}
