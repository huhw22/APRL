#include "momentum_z_gather.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
using namespace aprl;using namespace source_coupling_math;
int main(int argc,char** argv){MPI_Init(&argc,&argv);
 EBGridGeometry geometry(16,15,18,5e-6,5e-6,1.2543262415359034e-6,1.2543262415359034e-6/SI::c,EBMaxwellSolver::CowanZ);
 EBFieldGrid f(geometry);NativeHalo halo;halo.exchange(f,MPI_COMM_WORLD,true);const Vec origin(0.);const auto& g=f.geometry();
 for(unsigned a=0;a<3;++a){auto& b=a==0?f.bx():(a==1?f.by():f.bz());for(std::size_t k=0;k<b.nz();++k)for(std::size_t j=0;j<b.ny();++j)for(std::size_t i=0;i<b.nx();++i)b(i,j,k)=.003*(std::sin(.71*i+.39*j+.53*k+.7*a)+.3*std::cos(.11*i-.47*j+.31*k));}
 ChargeConservingCurrentDepositor depositor(f,origin);long double worst=0.,oldWorst=0.,energyWorst=0.;std::vector<std::pair<Vec,Vec>> paths;
 for(unsigned seed=0;seed<256;++seed){Vec a(0.),b(0.);const double spacing[3]={g.dx,g.dy,g.dz};for(unsigned c=0;c<3;++c){a[c]=(7.4+.21*std::sin(.71*seed+.43*c))*spacing[c];b[c]=a[c]+(seed%2?-.83:1.75)*std::cos(.29*seed+.81*c)*spacing[c];}paths.push_back({a,b});f.clearCurrent();const double q=-SI::elementaryCharge;depositor.depositSegment(a,b,q);
   long double mesh=0.,scale=0.;for(std::size_t k=1;k<g.nz;++k)for(std::size_t j=1;j<g.ny;++j)for(std::size_t i=1;i<g.nx;++i){const long double x=f.jx()(i,j,k)*.5L*(f.by()(i,j,k-1)+f.by()(i,j,k)),y=f.jy()(i,j,k)*.5L*(f.bx()(i,j,k-1)+f.bx()(i,j,k));mesh+=x-y;scale+=std::abs(x)+std::abs(y);}
   const long double factor=(long double)g.dt*g.dx*g.dy*g.dz;mesh*=factor;scale*=factor;const auto candidate=momentumOrbitMagnetic(f,origin,a,b,&halo),old=openOrbitMagnetic(f,origin,a,b,&halo);
   auto impulse=[&](const Vec& field){return (long double)q*(((long double)b[0]-a[0])*field[1]-((long double)b[1]-a[1])*field[0]);};worst=std::max(worst,std::abs(impulse(candidate)-mesh)/scale);oldWorst=std::max(oldWorst,std::abs(impulse(old)-mesh)/scale);
   Vec u(0.);u[0]=.4;u[1]=-.17;u[2]=.23;const Vec before(u);Vec e(0.);e[0]=10.;e[1]=-17.;momentum(u,e,candidate,q,SI::electronMass,g.dt);const long double v=SI::electronMass*SI::c*SI::c*deltaGamma(before,u);const double gn=BoostFrameTransform::gammaFromProperVelocity(u),go=BoostFrameTransform::gammaFromProperVelocity(before);Vec delta(before);delta+=u;delta*=g.dt*SI::c/(gn+go);energyWorst=std::max(energyWorst,std::abs(v-chargeWork(q,e,Vec(0.),delta))/(SI::electronMass*SI::c*SI::c));
 }
 if(!(worst<2e-12&&oldWorst>1e-4&&energyWorst<1e-14))throw std::runtime_error("Adjoint momentum/work check failed");
 volatile double sink=0.;double times[2]={};for(unsigned mode=0;mode<2;++mode){auto start=std::chrono::steady_clock::now();for(unsigned repeat=0;repeat<500;++repeat)for(const auto& p:paths){const auto b=mode?momentumOrbitMagnetic(f,origin,p.first,p.second,&halo):openOrbitMagnetic(f,origin,p.first,p.second,&halo);sink+=b[0]+b[1]+b[2];}times[mode]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
 std::cout<<std::setprecision(18)<<"{\"cases\":256,\"adjoint_relative_max\":"<<worst<<",\"old_gather_relative_max\":"<<oldWorst<<",\"magnetic_no_work_normalized_max\":"<<energyWorst<<",\"old_gather_seconds\":"<<times[0]<<",\"candidate_gather_seconds\":"<<times[1]<<",\"scalar_gather_time_ratio\":"<<times[1]/times[0]<<"}\n";MPI_Finalize();
}
