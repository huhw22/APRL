// Independent vacuum experiment using the immutable production Maxwell/CPML
// archive. No particles, reference supply, moving window, current or forcing.
#include "eb_field.h"
#include "eb_cpml.h"
#include "eb_mpi.h"
#include "compute_threads.h"
#include "shutdown_watchdog.h"
#include "boostframe.h"
#include "yee_sbp_observer.h"
#include <mpi.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <array>
using namespace aprl;
int main(int argc,char**argv){int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rank=0,size=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);initializeShutdownWatchdog(rank,true);configureShutdownWatchdog(true,600,60);int status=0;
try{
 if(argc!=5)throw std::runtime_error("Need directory cowan|yee axis|oblique xy-refinement");configureComputeThreads(provided,8);
 const std::string root=argv[1],solver=argv[2],mode=argv[3];if(solver!="yee")throw std::runtime_error("This exact SBP control supports Yee only");const unsigned refine=std::stoul(argv[4]);if(refine!=1&&refine!=2)throw std::runtime_error("Invalid refinement");
 const std::size_t nx=160*refine,ny=160*refine,globalNz=512;if(globalNz%size)throw std::runtime_error("Unequal test slabs");const std::size_t nz=globalNz/size,off=rank*nz;
 const double dx=5e-6/refine,dy=dx,dz=1.2543262415359034e-6;const double dt=solver=="cowan"?dz/SI::c:.98/(SI::c*std::sqrt(1/(dx*dx)+1/(dy*dy)+1/(dz*dz)));const auto kind=solver=="cowan"?EBMaxwellSolver::CowanZ:EBMaxwellSolver::Yee;
 EBGridGeometry g(nx,ny,nz,dx,dy,dz,dt,kind);EBFieldGrid f(g),view(g);EBZSlabHaloExchange halo(MPI_COMM_WORLD);halo.installPhysicalBoundaryMask(f);
 EBCPMLParameters parameters;parameters.cells[0]=parameters.cells[1]=10*refine;parameters.cells[2]=16;parameters.targetReflection=1e-8;parameters.kappaMax=8;parameters.alphaFraction=0.;std::unique_ptr<EBConvolutionalPML> boundary(new EBConvolutionalPML(g,globalNz,off,parameters,rank==0,rank+1==size));auto* cpml=boundary.get();f.setBoundary(std::move(boundary));
 const long double pi=std::acos(-1.L),wave=2*pi/25e-6L,cosine=mode=="axis"?1.L:.23L,sine=std::sqrt(1-cosine*cosine),radius=120e-6L,sigma=40e-6L;const double ox=-.0004,oy=-.0004,oz=-256*dz;
 auto potential=[&](long double x,long double y,long double z){const long double rr=x*x+y*y+z*z;if(rr>=radius*radius)return 0.L;return std::sin(wave*(sine*x+cosine*z))*std::exp(-rr/(2*sigma*sigma)-rr/(radius*radius-rr))/(SI::c*wave);};
 auto Cx=[&](long double x,long double y,long double z){return SI::c*cosine*potential(x,y,z);};auto Cz=[&](long double x,long double y,long double z){return -SI::c*sine*potential(x,y,z);};
 auto nativeE=[&](unsigned a,std::size_t i,std::size_t j,std::size_t k){const long double x=ox+i*dx,y=oy+j*dy,z=oz+(off+k)*dz;if(a==0)return (Cz(x+.5L*dx,y+.5L*dy,z)-Cz(x+.5L*dx,y-.5L*dy,z))/dy;if(a==1)return (Cx(x,y+.5L*dy,z+.5L*dz)-Cx(x,y+.5L*dy,z-.5L*dz))/dz-(Cz(x+.5L*dx,y+.5L*dy,z)-Cz(x-.5L*dx,y+.5L*dy,z))/dx;return -(Cx(x,y+.5L*dy,z+.5L*dz)-Cx(x,y-.5L*dy,z+.5L*dz))/dy;};
 YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};YeeComponent* b[3]={&f.bx(),&f.by(),&f.bz()};
 for(unsigned a=0;a<3;++a)for(std::size_t k=0;k<e[a]->nz();++k)for(std::size_t j=0;j<e[a]->ny();++j)for(std::size_t i=0;i<e[a]->nx();++i)(*e[a])(i,j,k)=(double)nativeE(a,i,j,k);
 // Obtain the exact implemented initial Faraday increment with E0 and B=0;
 // compact support guarantees the CPML histories remain identically zero.
 halo.advanceMagnetic(f);for(std::size_t h=0;h<cpml->checkpointHistoryCount();++h)for(double v:cpml->checkpointHistory(h))if(v!=0.)throw std::runtime_error("Pulse touches initial CPML");
 for(unsigned a=0;a<3;++a)for(std::size_t k=0;k<b[a]->nz();++k)for(std::size_t j=0;j<b[a]->ny();++j)for(std::size_t i=0;i<b[a]->nx();++i){const long double x=ox+i*dx,y=oy+j*dy,z=oz+(off+k)*dz;long double physical=0.;if(a==0)physical=-(potential(x,y+.5L*dy,z+dz)-potential(x,y+.5L*dy,z))/dz;if(a==2)physical=(potential(x+dx,y+.5L*dy,z)-potential(x,y+.5L*dy,z))/dx;(*b[a])(i,j,k)=(double)(physical-.5L*(*b[a])(i,j,k));}
 const std::size_t gx=12*refine,gy=12*refine,gz=66;if(gz%nz==0||(globalNz-gz)%nz==0)throw std::runtime_error("Test face coincides with MPI interface");
 BoostFrameTransform frame;frame.setOriginsFromGamma(929.1305492858543,SI::c,0.,0.,0.);const long double gamma=frame.gamma(),beta=frame.beta(),volume=(long double)dx*dy*dz;
 if(!rank){std::ofstream params(root+"/frame_parameters.json");params<<std::setprecision(21)<<"{\"gamma\":"<<gamma<<",\"beta\":"<<beta<<",\"dt_s\":"<<dt<<"}\n";if(!params)throw std::runtime_error("Cannot write frame parameters");}
 YeeSBPObserver sbp(root,g,off,globalNz,gx,gy,gz,MPI_COMM_WORLD);
 std::array<YeeComponent,3> oldE,oldB;std::ofstream csv;if(!rank){csv.open(root+"/history.csv");csv<<std::setprecision(21)<<"step,time_box_s,Uprime_J,cPzprime_J,lab_charge_J,lab_outward_J,boost_outward_J,lab_balance_pct,boost_balance_pct,divB_fixed_scale\n";}
 // Collocate native components directly at the geometric face. Averaging
 // already-centred cells would apply an extra longitudinal low-pass filter.
 auto facePoint=[&](const EBFieldGrid& q,unsigned axis,std::size_t i,std::size_t j,std::size_t k,std::array<long double,6>& v){
  const YeeComponent* components[6]={&q.ex(),&q.ey(),&q.ez(),&q.bx(),&q.by(),&q.bz()};
  for(unsigned a=0;a<6;++a){std::size_t lo[3]={i,j,k},count[3]={1,1,1};for(unsigned d=0;d<3;++d){const bool targetHalf=d!=axis,sourceHalf=a<3?d==a:d!=a-3;if(targetHalf!=sourceHalf){count[d]=2;if(sourceHalf)--lo[d];}}
   long double sum=0.;for(unsigned z=0;z<count[2];++z)for(unsigned y=0;y<count[1];++y)for(unsigned x=0;x<count[0];++x)sum+=(*components[a])(lo[0]+x,lo[1]+y,lo[2]+z);v[a]=sum/(count[0]*count[1]*count[2]);}
 };
 auto energies=[&](const EBFieldGrid& q){long double u=0.,p=0.,div=0.,bnorm=0.;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for reduction(+:u,p,div,bnorm) schedule(static)
#endif
  for(std::ptrdiff_t kk=0;kk<(std::ptrdiff_t)nz;++kk){const auto k=(std::size_t)kk;if(off+k<gz||off+k>=globalNz-gz)continue;
   const YeeComponent* ee[3]={&q.ex(),&q.ey(),&q.ez()};const YeeComponent* bb[3]={&q.bx(),&q.by(),&q.bz()};
   for(std::size_t j=gy;j<ny-gy;++j)for(std::size_t i=gx;i<nx-gx;++i){
    for(unsigned a=0;a<3;++a){for(unsigned p0=0;p0<2;++p0)for(unsigned q0=0;q0<2;++q0){const long double v=(*ee[a])(i+(a==0?0:p0),j+(a==0?p0:(a==1?0:q0)),k+(a==2?0:q0));u+=volume*SI::epsilon0*v*v/8;}
     for(unsigned p0=0;p0<2;++p0){const long double v=(*bb[a])(i+(a==0?p0:0),j+(a==1?p0:0),k+(a==2?p0:0));u+=volume*v*v/(4*SI::mu0);bnorm+=volume*v*v/2;}}
    const long double d=(q.bx()(i+1,j,k)-q.bx()(i,j,k))/dx+(q.by()(i,j+1,k)-q.by()(i,j,k))/dy+(q.bz()(i,j,k+1)-q.bz()(i,j,k))/dz;div+=volume*d*d;
   }
   for(std::size_t i=gx;i<nx-gx;++i)for(std::size_t j=gy;j<=ny-gy;++j){const long double w=(j==gy||j==ny-gy)?.5L:1.L;p+=volume*SI::epsilon0*SI::c*w*.5L*(q.ex()(i,j,k)+q.ex()(i,j,k+1))*q.by()(i,j,k);}
   for(std::size_t j=gy;j<ny-gy;++j)for(std::size_t i=gx;i<=nx-gx;++i){const long double w=(i==gx||i==nx-gx)?.5L:1.L;p-=volume*SI::epsilon0*SI::c*w*.5L*(q.ey()(i,j,k)+q.ey()(i,j,k+1))*q.bx()(i,j,k);}
  }
  long double local[4]={u,p,div,bnorm},all[4]={};MPI_Allreduce(local,all,4,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);return std::array<long double,4>{{all[0],all[1],all[2],all[3]}};};
 auto flux=[&](const EBFieldGrid& q){long double lab=0.,boost=0.,pressure=0.;auto face=[&](unsigned axis,int sign,std::size_t i,std::size_t j,std::size_t k,long double area){std::array<long double,6> a;facePoint(q,axis,i+(axis==0),j+(axis==1),k+(axis==2),a);long double u=0.;for(unsigned n=0;n<3;++n)u+=.5L*(SI::epsilon0*a[n]*a[n]+a[n+3]*a[n+3]/SI::mu0);const auto x=(axis+1)%3,y=(axis+2)%3;const long double s=(a[x]*a[y+3]-a[y]*a[x+3])/SI::mu0,t=(axis==2?u:0.)-SI::epsilon0*a[2]*a[axis]-a[5]*a[axis+3]/SI::mu0;boost+=sign*area*s;lab+=sign*area*gamma*(s+beta*SI::c*t);pressure+=sign*area*SI::c*t;};
  for(std::size_t k=0;k<nz;++k){if(off+k<gz||off+k>=globalNz-gz)continue;for(std::size_t j=gy;j<ny-gy;++j){face(0,-1,gx-1,j,k,(long double)dy*dz);face(0,1,nx-gx-1,j,k,(long double)dy*dz);}for(std::size_t i=gx;i<nx-gx;++i){face(1,-1,i,gy-1,k,(long double)dx*dz);face(1,1,i,ny-gy-1,k,(long double)dx*dz);}}
  for(unsigned side=0;side<2;++side){const auto k=side?globalNz-gz:gz;if(k>off&&k<off+nz)for(std::size_t j=gy;j<ny-gy;++j)for(std::size_t i=gx;i<nx-gx;++i)face(2,side?1:-1,i,j,k-off-1,(long double)dx*dy);}
  long double local[3]={lab,boost,pressure},all[3]={};MPI_Allreduce(local,all,3,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);return std::array<long double,3>{{all[0],all[1],all[2]}};};
 const double duration=400*dz/SI::c;const unsigned steps=(unsigned)std::ceil(duration/dt);long double firstG=0.,firstU=0.,firstB=0.,outLab=0.,outBoost=0.;std::array<long double,3> previousFlux{{0.,0.,0.}};
 for(unsigned n=0;n<=steps;++n){for(unsigned a=0;a<3;++a){oldE[a]=*e[a];oldB[a]=*b[a];}halo.advanceMagnetic(f);YeeComponent* ve[3]={&view.ex(),&view.ey(),&view.ez()};YeeComponent* vb[3]={&view.bx(),&view.by(),&view.bz()};for(unsigned a=0;a<3;++a){*ve[a]=oldE[a];for(std::size_t k=0;k<vb[a]->size();++k)vb[a]->data()[k]=.5*(oldB[a].data()[k]+b[a]->data()[k]);}
  sbp.state(n,view,oldB,f);const auto energy=energies(view);const auto gLab=gamma*(energy[0]+beta*energy[1]);if(!n){firstG=gLab;firstU=energy[0];firstB=energy[3];if(!(firstG>0.&&firstU>0.&&firstB>0.))throw std::runtime_error("Zero pulse energy");}else{outLab+=dt*previousFlux[0];outBoost+=dt*previousFlux[1];}
  if(!rank){csv<<n<<','<<n*dt<<','<<energy[0]<<','<<energy[1]<<','<<gLab<<','<<outLab<<','<<outBoost<<','<<100*(gLab-firstG+outLab)/firstG<<','<<100*(energy[0]-firstU+outBoost)/firstU<<','<<dz*std::sqrt(energy[2]/firstB)<<'\n';if(!csv)throw std::runtime_error("Pulse CSV failure");if(n%40==0)std::cout<<"step "<<n<<'/'<<steps<<" lab balance % "<<(double)(100*(gLab-firstG+outLab)/firstG)<<std::endl;}
  if(n==steps)break;halo.advanceElectric(f);for(unsigned a=0;a<3;++a){*vb[a]=*b[a];for(std::size_t k=0;k<ve[a]->size();++k)ve[a]->data()[k]=.5*(oldE[a].data()[k]+e[a]->data()[k]);}sbp.transition(n,oldE,f);previousFlux=flux(view);sbp.oldFlux(n,dt*previousFlux[1],dt*previousFlux[2]);
 }
 sbp.close();if(!rank){csv.flush();if(!csv)throw std::runtime_error("Pulse flush failed");csv.close();std::ofstream(root+"/complete")<<"1 "<<steps<<'\n';}
}catch(const std::exception&e){std::cerr<<"rank "<<rank<<": "<<e.what()<<std::endl;armShutdownWatchdog("vacuum-failure");MPI_Abort(MPI_COMM_WORLD,2);status=2;}
armShutdownWatchdogForFinalize("vacuum-finalize");MPI_Finalize();disarmShutdownWatchdog();shutdownShutdownWatchdog();return status;}
