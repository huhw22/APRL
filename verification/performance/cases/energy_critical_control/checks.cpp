#include "longitudinal_shape.h"
#include "source_coupling.h"
#include "source_coupling_math.h"
#include "blocked_cowan_reference.h"
#include "charge_matched_inlet.h"
#include "eb_deposition.h"
#include "eb_mpi.h"
#include "runtime_util.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace aprl;
using Vec=FieldVector<Double>;
using Shape=longitudinal_shape::Shape;
namespace {
const double coeff[7]={1./64.,-6./64.,15./64.,44./64.,15./64.,-6./64.,1./64.};
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
double sign(long long k){return k%2?-1.:1.;}
YeeComponent& current(EBFieldGrid& f,unsigned a){return a==0?f.jx():(a==1?f.jy():f.jz());}
double independentWeight(double z,long long k,bool edge){
 double value=0.;
 for(int shift=-3;shift<=3;++shift){
  const double distance=z-(k-shift)-(edge?.5:0.);
  const double base=edge?(z>=k-shift&&z<k-shift+1?1.:0.):std::max(0.,1.-std::abs(distance));
  value+=coeff[shift+3]*base;
 }
 return value;
}
double actualWeight(double z,long long k,bool edge){
 const auto b=longitudinal_shape::basis(z,edge,Shape::LinearNotch);
 return k<b.first||k>=b.first+b.count?0.:b.v[k-b.first];
}
void charge(YeeComponent& rho,const EBGridGeometry& g,const Vec& origin,const Vec& p,double q){
 const double x=(p[0]-origin[0])/g.dx,y=(p[1]-origin[1])/g.dy,z=(p[2]-origin[2])/g.dz;
 const auto i=(long long)std::floor(x),j=(long long)std::floor(y),k=(long long)std::floor(z);
 for(long long zz=k-4;zz<=k+5;++zz)for(unsigned b=0;b<2;++b)for(unsigned a=0;a<2;++a){
  const double value=q/(g.dx*g.dy*g.dz)*(a?x-i:1.-x+i)*(b?y-j:1.-y+j)*independentWeight(z,zz,false);
  if(value!=0.){require(zz>=0&&zz<(long long)rho.nz(),"Independent charge support outside test");rho(i+a,j+b,zz)+=value;}
 }
}
void fill(EBFieldGrid& f,std::size_t offset,double phase){
 YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};
 for(unsigned a=0;a<3;++a)for(std::size_t k=0;k<e[a]->nz();++k)for(std::size_t j=0;j<e[a]->ny();++j)for(std::size_t i=0;i<e[a]->nx();++i)
  (*e[a])(i,j,k)=1e5*std::sin(.17*i+.23*j+.31*(k+offset)+.37*a+phase);
}
void legacyConvolution(EBFieldGrid& f,const Vec& origin,const Vec& begin,const Vec& end,double q){
 for(int shift=-3;shift<=3;++shift){Vec a(begin),b(end);a[2]+=shift*f.geometry().dz;b[2]+=shift*f.geometry().dz;
  ChargeConservingCurrentDepositor(f,origin).depositSegment(a,b,q*coeff[shift+3]);}
}
double currentError(EBFieldGrid& f,EBFieldGrid& oracle,std::size_t offset){
 double error=0.,scale=0.;
 for(unsigned a=0;a<3;++a){auto& v=current(f,a);auto& r=current(oracle,a);
  for(std::size_t k=0;k<v.nz();++k)for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i){
   error=std::max(error,std::abs(v(i,j,k)-r(i,j,k+offset)));scale=std::max(scale,std::abs(r(i,j,k+offset)));}}
 return error/std::max(scale,1e-300);
}
long double continuity(EBFieldGrid& f,const YeeComponent& before,const YeeComponent& after,unsigned particles){
 const auto& g=f.geometry();long double error=0.;
 for(std::size_t k=1;k<g.nz;++k)for(std::size_t j=1;j<g.ny;++j)for(std::size_t i=1;i<g.nx;++i){
  const long double div=((long double)f.jx()(i,j,k)-f.jx()(i-1,j,k))/g.dx+
   ((long double)f.jy()(i,j,k)-f.jy()(i,j-1,k))/g.dy+((long double)f.jz()(i,j,k)-f.jz()(i,j,k-1))/g.dz;
  error+=std::abs((long double)after(i,j,k)-before(i,j,k)+g.dt*div)*g.dx*g.dy*g.dz/(particles*SI::elementaryCharge);
 }
 return error;
}
double filterAndBasis(){
 double maxError=0.;
 for(unsigned n=0;n<401;++n){const double theta=std::acos(-1.)*n/400.;double symbol=0.;
  for(int r=-3;r<=3;++r)symbol+=coeff[r+3]*std::cos(r*theta);
  maxError=std::max(maxError,std::abs(symbol-(1.-std::pow(std::sin(theta/2.),6))));
 }
 for(unsigned n=0;n<400;++n){const double z=-3.71+.071*n;
  for(bool edge:{false,true}){double sum=0.,nyquist=0.;
   for(long long k=(long long)std::floor(z)-6;k<=(long long)std::floor(z)+7;++k){
    const double a=actualWeight(z,k,edge),b=independentWeight(z,k,edge);
    maxError=std::max(maxError,std::abs(a-b));sum+=a;nyquist+=sign(k)*a;
   }
   require(std::abs(sum-1.)<2e-14,"Notch basis partition of unity");
   require(std::abs(nyquist)<2e-14,"Notch source does not remove exact Nyquist");
  }
 }
 require(maxError<5e-14,"Notch basis or Fourier response differs from independent convolution");return maxError;
}
struct DistributedResult{double current=0.,continuity=0.,work=0.,nyquist=0.;};
DistributedResult distributed(MPI_Comm comm,const std::string& output){
 int rank=0,size=1;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
 const double dz=8.293564128124092e-6;
 EBGridGeometry g(10,11,24,4e-5,4e-5,dz,dz/SI::c,EBMaxwellSolver::CowanZ),gg=g;gg.nz*=size;
 EBFieldGrid f(g),oracle(gg);Vec origin(0.),local(0.);local[2]=rank*g.nz*dz;
 fill(f,rank*g.nz,0.);const auto old=source_coupling_math::copyElectric(f);
 createDirectories(output,comm);ConservativeSourceCoupling coupling(comm,output,false,true);
 coupling.begin(f,local,std::vector<RelativisticParticleSI>());fill(f,rank*g.nz,.31);
 YeeComponent rho0(gg.nx+1,gg.ny+1,gg.nz+1),rho1(gg.nx+1,gg.ny+1,gg.nz+1);
 const unsigned particles=67;
 for(unsigned n=0;n<particles;++n){Vec a(0.);a[0]=(3.13+.017*(n%47))*g.dx;a[1]=(4.29+.019*(n%41))*g.dy;
  a[2]=(5.23+(gg.nz-10.56)*((n*37)%particles)/particles)*dz;Vec b(a);
  b[0]+=(n%2?.37:-.23)*g.dx;b[1]+=(n%3?.29:-.31)*g.dy;b[2]+=(n%2?.43:-.34)*dz;
  const double q=-SI::elementaryCharge;legacyConvolution(oracle,origin,a,b,q);charge(rho0,gg,origin,a,q);charge(rho1,gg,origin,b,q);
  const double za=(a[2]-local[2])/dz,zb=(b[2]-local[2])/dz;double lo=0.,hi=1.;
  if(zb>za){lo=std::max(lo,-za/(zb-za));hi=std::min(hi,(g.nz-za)/(zb-za));}
  else{lo=std::max(lo,(g.nz-za)/(zb-za));hi=std::min(hi,-za/(zb-za));}
  if(hi>lo){Vec first(a),last(a);for(unsigned d=0;d<3;++d){first[d]=a[d]+lo*(b[d]-a[d]);last[d]=a[d]+hi*(b[d]-a[d]);}
   coupling.depositFrozenSourceSupport(first,last,q,true);}
 }
 coupling.completeCurrent();EBZSlabHaloExchange exchange(comm);exchange.sumSharedTangentialCurrent(f);
 DistributedResult result;result.current=currentError(f,oracle,rank*g.nz);result.continuity=(double)continuity(oracle,rho0,rho1,particles);
 // Every global current component's exact alternating sum must vanish.
 double nyquist=0.,currentScale=0.;
 for(unsigned a=0;a<3;++a){auto& v=current(oracle,a);for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i){long double s=0.;
  for(std::size_t k=0;k<v.nz();++k){s+=sign(k)*v(i,j,k);currentScale=std::max(currentScale,std::abs(v(i,j,k)));}
  nyquist=std::max(nyquist,(double)std::abs(s));}}
 result.nyquist=nyquist/std::max(currentScale,1e-300);
 const auto recorded=coupling.finalizeConstraints();long double localWork=0.,absolute=0.;
 const YeeComponent* electric[3]={&f.ex(),&f.ey(),&f.ez()};
 for(unsigned a=0;a<3;++a){auto& j=current(f,a);for(std::size_t k=0;k<j.nz();++k)for(std::size_t y=0;y<j.ny();++y)for(std::size_t x=0;x<j.nx();++x){
  long double weight=1.;if(a<2&&((rank>0&&k==0)||(rank+1<size&&k==g.nz)))weight=.5;
  const long double w=weight*.5L*(old[a](x,y,k)+(*electric[a])(x,y,k))*j(x,y,k)*g.dt*g.dx*g.dy*g.dz;
  localWork+=w;absolute+=std::abs(w);}}
 const long double localValues[3]={localWork,recorded[0],absolute};long double globalValues[3]={};
 MPI_Allreduce(localValues,globalValues,3,MPI_LONG_DOUBLE,MPI_SUM,comm);
 result.work=(double)(std::abs(globalValues[0]+globalValues[1])/std::max(globalValues[2],1e-100L));
 coupling.close(1);return result;
}
double magneticNyquist(const std::string& output,bool notch){
 const double dz=8.293564128124092e-6;EBGridGeometry g(8,9,32,4e-5,4e-5,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
 EBFieldGrid f(g);for(std::size_t k=0;k<f.by().nz();++k)for(std::size_t j=0;j<f.by().ny();++j)for(std::size_t i=0;i<f.by().nx();++i)f.by()(i,j,k)=sign(k)*.1;
 RelativisticParticleSI p;p.id=1;p.sourceId=1;p.mass=SI::electronMass;p.charge=-SI::elementaryCharge;
 p.position=Vec(0.);p.position[0]=3.27*g.dx;p.position[1]=4.31*g.dy;p.position[2]=15.23*dz;
 p.properVelocity=Vec(0.);p.properVelocity[0]=.23;p.properVelocity[2]=.51;const Vec u0(p.properVelocity);
 createDirectories(output,MPI_COMM_SELF);ConservativeSourceCoupling coupling(MPI_COMM_SELF,output,false,notch);
 coupling.begin(f,Vec(0.),std::vector<RelativisticParticleSI>(1,p));SIFieldSourceSet sources;BoostFrameTransform frame;frame.setFromGamma(1.,SI::c,0.);
 coupling.push(p,sources,frame,0.,g.dt,7,nullptr,nullptr,nullptr);Vec difference(p.properVelocity);difference-=u0;
 coupling.close(0);return std::sqrt(difference.norm2());
}
double admissionSupport(const std::string& output){
 const double dz=8.293564128124092e-6,q=-SI::elementaryCharge;
 EBGridGeometry g(8,9,24,4e-5,4e-5,dz,dz/SI::c,EBMaxwellSolver::CowanZ),gg=g;gg.nz+=16;
 EBFieldGrid f(g),oracle(gg),old(g);fill(f,0,0.);Vec origin(0.),padded(0.);padded[2]=-8*dz;
 Vec begin(0.);begin[0]=3.27*g.dx;begin[1]=4.31*g.dy;begin[2]=-.2*dz;
 Vec end(begin);end[2]+=.51479762393803552*dz;Vec entry(begin);entry[2]=0.;
 createDirectories(output,MPI_COMM_SELF);ConservativeSourceCoupling coupling(MPI_COMM_SELF,output,false,true);
 coupling.begin(f,origin,std::vector<RelativisticParticleSI>());
 coupling.depositFrozenSourceSupport(begin,entry,q,true);
 ChargeConservingCurrentDepositor(f,origin,coupling.longitudinalCurrent()).depositSegment(entry,end,q);
 coupling.recordConstraint(entry,end,q,false);coupling.completeCurrent();
 legacyConvolution(oracle,padded,begin,end,q);
 longitudinal_shape::Current clipped(old,MPI_COMM_SELF,Shape::LinearNotch);clipped.clear();clipped.deposit(origin,entry,end,q);
 YeeComponent rho0(gg.nx+1,gg.ny+1,gg.nz+1),rho1(gg.nx+1,gg.ny+1,gg.nz+1);
 charge(rho0,gg,padded,begin,q);charge(rho1,gg,padded,end,q);
 long double fullDefect=0.,oldDefect=0.;
 for(std::size_t k=1;k<g.nz;++k)for(std::size_t j=1;j<g.ny;++j)for(std::size_t i=1;i<g.nx;++i){
  const long double dr=(long double)rho1(i,j,k+8)-rho0(i,j,k+8);
  fullDefect+=std::abs(dr+g.dt*((long double)f.jz()(i,j,k)-f.jz()(i,j,k-1))/dz)*g.dx*g.dy*g.dz/std::abs(q);
  oldDefect+=std::abs(dr+g.dt*((long double)old.jz()(i,j,k)-old.jz()(i,j,k-1))/dz)*g.dx*g.dy*g.dz/std::abs(q);
 }
 require(oldDefect>.01L,"Notch admission negative control did not retain an internal defect");
 const auto work=coupling.finalizeConstraints();long double mesh=0.,absolute=0.;
 for(std::size_t k=0;k<f.jz().nz();++k)for(std::size_t j=0;j<f.jz().ny();++j)for(std::size_t i=0;i<f.jz().nx();++i){
  const long double value=(long double)f.ez()(i,j,k)*f.jz()(i,j,k)*g.dt*g.dx*g.dy*g.dz;mesh+=value;absolute+=std::abs(value);}
 const double error=std::max(currentError(f,oracle,8),std::max((double)fullDefect,(double)(std::abs(mesh+work[0])/std::max(absolute,1e-100L))));
 require(error<2e-10,"Notch admission source or work support is incomplete");coupling.close(1);return error;
}
double referenceAndInlet(const std::string& output){
 const double dz=8.293564128124092e-6;
 EBGridGeometry g(16,16,48,4e-5,4e-5,dz,dz/SI::c,EBMaxwellSolver::CowanZ);Vec origin(0.);
 std::vector<RelativisticParticleSI> p(3);
 for(unsigned n=0;n<p.size();++n){p[n].position=Vec(0.);p[n].position[0]=(6.1+.3*n)*g.dx;
  p[n].position[1]=(7.2+.2*n)*g.dy;p[n].position[2]=(14.2+5.4*n)*g.dz;
  p[n].charge=-SI::elementaryCharge;p[n].mass=SI::electronMass;}
 const std::size_t guards[3]={0,0,0};
 verification::BlockedCowanReference ref(g,origin,p,.51479762393803552,guards,8,output,MPI_COMM_SELF,false,true);
 YeeComponent rho(17,17,49);for(const auto& particle:p)charge(rho,g,origin,particle.position,particle.charge);
 double error=0.,scale=0.;
 for(std::size_t k=1;k<g.nz;++k)for(std::size_t j=1;j<g.ny;++j)for(std::size_t i=1;i<g.nx;++i){
  error=std::max(error,std::abs(rho(i,j,k)-ref.chargeDensity(i,j,k)));scale=std::max(scale,std::abs(rho(i,j,k)));}
 require(error/std::max(scale,1e-300)<2e-13,"Prepared reference does not use matched notch charge");
 auto reader=verification::BlockedCowanReference::openPrepared(g,origin,.51479762393803552,guards,8,output,MPI_COMM_SELF,false,true);
 bool rejected=false;try{auto wrong=verification::BlockedCowanReference::openPrepared(g,origin,.51479762393803552,guards,8,output,MPI_COMM_SELF,false,false);}catch(...){rejected=true;}
 require(rejected,"Prepared reference accepted wrong source kernel");
 EBFieldGrid f(g);ChargeMatchedInlet inlet(g,0,0,false,true);inlet.begin(origin);
 p[0].position[2]=-.2*dz;p[1].position[2]=1.4*dz;
 inlet.particle(p[0]);inlet.particle(p[1]);inlet.complete(f,true);
 double gauss=0.,rhoScale=0.;
 for(std::size_t j=1;j<g.ny;++j)for(std::size_t i=1;i<g.nx;++i){
  long double expected=0.;
  for(unsigned n=0;n<2;++n){const auto& a=p[n];const double x=a.position[0]/g.dx,y=a.position[1]/g.dy;
   expected+=(long double)a.charge/(g.dx*g.dy*g.dz)*std::max(0.,1.-std::abs(x-i))*std::max(0.,1.-std::abs(y-j))*independentWeight(a.position[2]/dz,1,false);}
  gauss=std::max(gauss,(double)std::abs(SI::epsilon0*((long double)f.ez()(i,j,1)-f.ez()(i,j,0))/dz-expected));
  rhoScale=std::max(rhoScale,(double)std::abs(expected));
 }
 require(gauss/std::max(rhoScale,1e-300)<2e-13,"Inlet does not use matched notch charge");
 return std::max(error/std::max(scale,1e-300),gauss/std::max(rhoScale,1e-300));
}
void bulkSourceControl(double& ordinary,double& filtered,double& phaseError){
 double multiplier=0.;for(int r=-3;r<=3;++r)multiplier+=coeff[r+3]*sign(r);
 double a[2]={0.,0.},b[2]={0.,0.};
 for(unsigned n=0;n<1000;++n){const double drive=sign(n),ea=-3*a[0]-2*a[1]+drive,ha=2*a[0]+a[1];
  const double eb=-3*b[0]-2*b[1]+multiplier*drive,hb=2*b[0]+b[1];a[0]=ea;a[1]=ha;b[0]=eb;b[1]=hb;}
 ordinary=a[0]*a[0]+a[1]*a[1];filtered=b[0]*b[0]+b[1]*b[1];
 require(ordinary>1e8&&filtered==0.,"Notch did not suppress a resonant source in the unchanged bulk recurrence");
 const auto c=EBCowanCoefficients::forZDispersion(4e-5,4e-5,8.293564128124092e-6);phaseError=0.;
 const double normalized=c.alpha[2]+2*c.beta[0]+2*c.beta[1]+4*c.deltaXY;
 for(unsigned n=1;n<100;++n){const double theta=std::acos(-1.)*n/100.,q=4*std::pow(std::sin(theta/2.),2)*normalized;
  phaseError=std::max(phaseError,std::abs(std::acos(1-q/2)-theta));}
 require(phaseError<1e-12,"Unchanged axial Cowan dispersion is not exact");
}
}
int main(int argc,char** argv){MPI_Init(&argc,&argv);int rank=0,size=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
 try{require(argc==2,"Need independent output directory");createDirectories(argv[1],MPI_COMM_WORLD);const double basisError=filterAndBasis();
  auto d=distributed(MPI_COMM_WORLD,std::string(argv[1])+"/distributed");double local[4]={d.current,d.continuity,d.work,d.nyquist},global[4]={};
  MPI_Allreduce(local,global,4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  require(global[0]<2e-10&&global[1]<2e-10,"Filtered current/charge does not preserve the expected source");
  require(global[2]<2e-10,"Filtered E gather is not adjoint to actual mesh current");require(global[3]<2e-10,"Filtered mesh source retains a Nyquist component");
  double ordinary=0.,filtered=0.,phaseError=0.;bulkSourceControl(ordinary,filtered,phaseError);
  if(rank==0){const double oldB=magneticNyquist(std::string(argv[1])+"/magnetic-old",false),newB=magneticNyquist(std::string(argv[1])+"/magnetic-notch",true);
   const double referenceError=referenceAndInlet(std::string(argv[1])+"/reference");
   const double admissionError=admissionSupport(std::string(argv[1])+"/admission");
   require(oldB>1e-7&&newB<1e-13,"Notch magnetic gather did not remove grid-scale magnetic force");
   std::cout<<std::setprecision(18)<<"{\"passed\":true,\"physical_acceptance\":false,\"ranks\":"<<size<<",\"basis_error\":"<<basisError
    <<",\"current_relative_error\":"<<global[0]<<",\"continuity_L1_over_abs_charge\":"<<global[1]<<",\"mesh_gather_work_relative_error\":"<<global[2]
    <<",\"nyquist_current_relative_residual\":"<<global[3]<<",\"ordinary_nyquist_magnetic_delta_u\":"<<oldB<<",\"notched_nyquist_magnetic_delta_u\":"<<newB
    <<",\"admission_support_relative_error\":"<<admissionError<<",\"reference_inlet_charge_relative_error\":"<<referenceError<<",\"forced_bulk_native_norm_ordinary\":"<<ordinary<<",\"forced_bulk_native_norm_notch\":"<<filtered<<",\"axial_dispersion_phase_error\":"<<phaseError<<"}\n";
  }
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}MPI_Finalize();return 0;
}
