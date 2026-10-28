#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>
#include "blocked_cowan_reference.h"
#include "boostframe.h"
#include "compute_threads.h"
#include "eb_deposition.h"
#include "eb_initial_field.h"

namespace {
using aprl::EBFieldGrid;using aprl::EBGridGeometry;using aprl::FieldVector;
using aprl::RelativisticParticleSI;using aprl::YeeComponent;
using aprl::verification::BlockedCowanReference;
const double CG_TOLERANCE=1e-12,GAUSS_TOLERANCE=5e-11,CIC_TOLERANCE=2e-14;
void require(bool value,const std::string& message){if(!value)throw std::runtime_error(message);}
void freshDirectory(const std::string& path){require(mkdir(path.c_str(),0755)==0,"Cannot create fresh directory "+path);}
YeeComponent& component(EBFieldGrid& f,unsigned a) {
  if(a==0)return f.ex();if(a==1)return f.ey();if(a==2)return f.ez();
  if(a==3)return f.bx();if(a==4)return f.by();return f.bz();
}
bool unchanged(const std::vector<RelativisticParticleSI>& a,const std::vector<RelativisticParticleSI>& b) {
  if(a.size()!=b.size())return false;
  for(std::size_t n=0;n<a.size();++n) {
    for(unsigned d=0;d<3;++d)if(a[n].position[d]!=b[n].position[d]||a[n].properVelocity[d]!=b[n].properVelocity[d])return false;
    if(a[n].id!=b[n].id||a[n].sourceId!=b[n].sourceId||a[n].charge!=b[n].charge||
       a[n].mass!=b[n].mass||a[n].weight!=b[n].weight||a[n].frozen!=b[n].frozen||
       a[n].frozenBetaZ!=b[n].frozenBetaZ)return false;
  }
  return true;
}
struct Spectrum {double minimum=std::numeric_limits<double>::infinity(),maximum=0.,gradientMaximum=0.;};
Spectrum spectrum(const EBGridGeometry& g,const std::size_t guards[3],double velocity) {
  const auto c=aprl::EBCowanCoefficients::forZDispersion(g.dx,g.dy,g.dz);
  const double delta[3]={c.deltaYZ,c.deltaZX,c.deltaXY},h[3]={g.dx,g.dy,g.dz};
  const std::size_t n[3]={g.nx-2*guards[0],g.ny-2*guards[1],g.nz-2*guards[2]};
  const double eta=velocity*velocity*aprl::SI::epsilon0*aprl::SI::mu0;
  Spectrum result;
  for(std::size_t k=1;k<n[2];++k)for(std::size_t j=1;j<n[1];++j)for(std::size_t i=1;i<n[0];++i) {
    const std::size_t mode[3]={i,j,k};double cosine[3],wave2[3];
    for(unsigned d=0;d<3;++d){const double t=std::acos(-1.)*mode[d]/n[d];cosine[d]=std::cos(t);const double q=2*std::sin(t/2)/h[d];wave2[d]=q*q;}
    double eigenvalue=0.,gradient2=0.;
    for(unsigned d=0;d<3;++d){const unsigned p=(d+1)%3,q=(d+2)%3;
      const double s=c.alpha[d]+2*c.beta[p]*cosine[p]+2*c.beta[q]*cosine[q]+4*delta[d]*cosine[p]*cosine[q]-(d==2?eta:0.);
      eigenvalue+=s*wave2[d];gradient2+=s*s*wave2[d];}
    result.minimum=std::min(result.minimum,eigenvalue);result.maximum=std::max(result.maximum,eigenvalue);
    result.gradientMaximum=std::max(result.gradientMaximum,std::sqrt(gradient2));
  }
  require(result.minimum>0.,"Nonpositive independent spectrum");return result;
}
std::vector<RelativisticParticleSI> source(const EBGridGeometry& g,const std::size_t guard[3],double beta) {
  const double fraction[8]={.013,.071,.183,.207,.397,.414,.797,.983};
  const double h[3]={g.dx,g.dy,g.dz};const std::size_t n[3]={g.nx,g.ny,g.nz};
  std::vector<RelativisticParticleSI> result;
  for(unsigned z=0;z<8;++z)for(unsigned y=0;y<3;++y)for(unsigned x=0;x<4;++x) {
    RelativisticParticleSI p;const unsigned index[3]={x+1,y+2,z};
    for(unsigned d=0;d<3;++d){const double lo=guard[d]+1.1,hi=n[d]-guard[d]-1.1;p.position[d]=(lo+(hi-lo)*fraction[index[d]])*h[d];}
    p.properVelocity=FieldVector<double>(0.);p.properVelocity[2]=beta/std::sqrt(1-beta*beta);
    p.weight=500.+x+3*y+7*z;p.mass=aprl::SI::electronMass*p.weight;p.charge=-aprl::SI::elementaryCharge*p.weight;
    p.id=p.sourceId=result.size();p.frozen=false;p.frozenBetaZ=0.;result.push_back(p);
  }
  return result;
}
struct Comparison {double electricRelative=0,magneticRelative=0,gauss=0,cic=0,seamMaximum=0;};
Comparison run(const std::string& root,double beta,bool zGuard,std::size_t block) {
  const double dz=1.2543262415359034e-6;
  const EBGridGeometry g(34,36,81,5e-6,5e-6,dz,dz/aprl::SI::c,aprl::EBMaxwellSolver::CowanZ);
  const std::size_t guards[3]={10,10,zGuard?2U:0U};
  auto particles=source(g,guards,beta);const auto original=particles;
  EBFieldGrid cg(g),candidate(g);
  const auto report=aprl::EBGaussFieldInitializer::initialize(cg,g,0,FieldVector<double>(0.),particles,
    aprl::InitialSelfFieldModel::CowanRelativisticPoisson,3e-5,guards,CG_TOLERANCE,10000,MPI_COMM_SELF);
  BlockedCowanReference generator(g,FieldVector<double>(0.),particles,report.meanBetaZ,guards,block,root);
  generator.fill(candidate);
  require(unchanged(original,particles),"Initializer changed original particle state");
  require(generator.report().complete,"Incomplete generator report");
  YeeComponent rho(g.nx+1,g.ny+1,g.nz+1);rho.fill(0.);
  aprl::ChargeConservingCurrentDepositor depositor(cg,FieldVector<double>(0.));
  for(const auto& p:particles)depositor.depositCharge(p.position,p.charge,rho);
  long double source2=0.,residual2=0.,cicDifference2=0.,charge2=0.;
  for(std::size_t k=0;k<=g.nz;++k)for(std::size_t j=0;j<=g.ny;++j)for(std::size_t i=0;i<=g.nx;++i) {
    const double q=rho(i,j,k),difference=generator.chargeDensity(i,j,k)-q;
    charge2+=static_cast<long double>(q)*q;cicDifference2+=static_cast<long double>(difference)*difference;
    if(i<=guards[0]||j<=guards[1]||k<=guards[2]||i>=g.nx-guards[0]||j>=g.ny-guards[1]||k>=g.nz-guards[2])continue;
    const double rhs=q/aprl::SI::epsilon0;
    const double residual=(candidate.ex()(i,j,k)-candidate.ex()(i-1,j,k))/g.dx+
      (candidate.ey()(i,j,k)-candidate.ey()(i,j-1,k))/g.dy+
      (candidate.ez()(i,j,k)-candidate.ez()(i,j,k-1))/g.dz-rhs;
    source2+=static_cast<long double>(rhs)*rhs;residual2+=static_cast<long double>(residual)*residual;
  }
  Comparison result;result.gauss=std::sqrt(residual2/source2);result.cic=std::sqrt(cicDifference2/charge2);
  require(result.gauss<GAUSS_TOLERANCE,"Independent native Gauss audit failed");
  require(result.cic<CIC_TOLERANCE,"Blocked charge differs from original CIC");
  const auto spec=spectrum(g,guards,report.meanBetaZ*aprl::SI::c);
  const double roundoff=256*std::numeric_limits<double>::epsilon()*spec.maximum/spec.minimum;
  const double absoluteElectricAllowance=spec.gradientMaximum/spec.minimum*
    (4*CG_TOLERANCE+GAUSS_TOLERANCE)*std::sqrt(source2);
  long double errors[2]={},norms[2]={};
  for(unsigned a=0;a<6;++a) {
    auto& expected=component(cg,a);auto& actual=component(candidate,a);long double error=0.,norm=0.;
    for(std::size_t k=0;k<actual.nz();++k)for(std::size_t j=0;j<actual.ny();++j)for(std::size_t i=0;i<actual.nx();++i) {
      const double difference=actual(i,j,k)-expected(i,j,k);
      require(std::isfinite(difference),"Nonfinite native field difference");
      error+=static_cast<long double>(difference)*difference;norm+=static_cast<long double>(expected(i,j,k))*expected(i,j,k);
      const bool seam=k==0||k+1==actual.nz()||k%block<=1||(k+1)%block==0;
      if(seam)result.seamMaximum=std::max(result.seamMaximum,std::abs(difference));
      if(a==5||(beta==0&&a>=3))require(actual(i,j,k)==0.&&expected(i,j,k)==0.,"Expected identically zero B component");
    }
    const double scale=a<3?1.:std::abs(report.meanBetaZ*aprl::SI::c*aprl::SI::epsilon0*aprl::SI::mu0);
    const double bound=scale*absoluteElectricAllowance+roundoff*std::sqrt(norm);
    require(std::sqrt(error)<=bound,"Native field exceeds independent residual/spectrum bound");
    errors[a>=3]+=error;norms[a>=3]+=norm;
  }
  result.electricRelative=std::sqrt(errors[0]/norms[0]);
  result.magneticRelative=norms[1]>0?std::sqrt(errors[1]/norms[1]):0.;
  // Direct provider queries cross block seams in both directions; fill() is
  // not allowed to hide an incorrect cache or a seam-specific reconstruction.
  for(unsigned a=0;a<6;++a){const auto& values=component(candidate,a);
    for(std::size_t k=values.nz();k-- >0;) {
      if(k!=0&&k+1!=values.nz()&&k%block>1&&(k+1)%block!=0)continue;
      const std::size_t i=values.nx()/2,j=values.ny()/2;
      require(generator.field(a,i,j,k)==values(i,j,k),"Repeated provider query changed a native value");
    }}
  bool invalidComponent=false,invalidPlane=false;
  try{generator.field(6,0,0,0);}catch(const std::exception&){invalidComponent=true;}
  try{generator.field(2,0,0,g.nz);}catch(const std::exception&){invalidPlane=true;}
  require(invalidComponent&&invalidPlane,"Out-of-range provider query was not rejected");
  require(unchanged(original,particles),"Provider changed source particle state");
  generator.writeReport(root+"/report.json");
  std::cout<<"case beta="<<beta<<" z_guard="<<zGuard<<" block="<<block
    <<" E_relative="<<result.electricRelative<<" B_relative="<<result.magneticRelative
    <<" gauss="<<result.gauss<<" CIC="<<result.cic<<" seam_max="<<result.seamMaximum<<'\n';
  return result;
}
void faultChecks(const std::string& root) {
  const double dz=1.2543262415359034e-6,beta=.23;
  const EBGridGeometry g(34,36,81,5e-6,5e-6,dz,dz/aprl::SI::c,aprl::EBMaxwellSolver::CowanZ);
  const std::size_t guards[3]={10,10,2},block=7;
  const auto particles=source(g,guards,beta);
  BlockedCowanReference generator(g,FieldVector<double>(0.),particles,beta,guards,block,root);
  const double originalPotential=generator.potential(guards[0]+2,guards[1]+2,guards[2]+2);
  require(std::isfinite(originalPotential)&&originalPotential!=0.,"Fault fixture is not nonzero");
  bool existsRejected=false;
  try{BlockedCowanReference duplicate(g,FieldVector<double>(0.),particles,beta,guards,block,root);}
  catch(const std::exception&){existsRejected=true;}
  require(existsRejected,"Existing generator directory was not rejected");
  generator.clearCache();
  require(generator.potential(guards[0]+2,guards[1]+2,guards[2]+2)==originalPotential,
    "Rejected overwrite altered an existing potential");
  generator.clearCache();const auto builds=generator.report().fieldCacheBuilds;
  for(unsigned iteration=0;iteration<1000;++iteration) {
    const std::size_t i=11+iteration%10,j=11+(iteration/10)%12;
    require(std::isfinite(generator.field(0,i,j,block-1))&&std::isfinite(generator.field(0,i,j,block)),
      "Nonfinite adjacent-plane cache query");
  }
  const auto repeatedBuilds=generator.report().fieldCacheBuilds-builds;
  require(repeatedBuilds<=1,"Adjacent interpolation planes repeatedly rebuilt the field cache");
  // This is a dedicated failure fixture. Never truncate a normal-case file.
  generator.clearCache();require(truncate((root+"/phi.bin").c_str(),0)==0,"Cannot prepare short-read negative control");
  bool shortReadRejected=false;
  try{generator.potential(guards[0]+1,guards[1]+1,guards[2]+1);}
  catch(const std::exception&){shortReadRejected=true;}
  require(shortReadRejected,"Truncated potential file did not fail closed");
  std::cout<<"fault_checks existing_directory_rejected=1 short_read_rejected=1 repeated_adjacent_cache_builds="<<repeatedBuilds<<'\n';
}
}
int main(int argc,char** argv) {
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  try {
    aprl::configureComputeThreads(provided);int ranks=0;MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    require(ranks==1&&argc==2,"Usage: aprl_blocked_cowan_reference_checks NEW_OUTPUT_DIRECTORY (one rank)");
    freshDirectory(argv[1]);double maxE=0.,maxB=0.,maxGauss=0.,maxCic=0.;unsigned count=0;
    std::cout<<std::setprecision(17);
    for(double beta:{0.,.23})for(bool guard:{false,true})for(std::size_t block:{7U,16U,32U}) {
      const std::string path=std::string(argv[1])+"/case-"+std::to_string(count);
      const auto result=run(path,beta,guard,block);++count;
      maxE=std::max(maxE,result.electricRelative);maxB=std::max(maxB,result.magneticRelative);
      maxGauss=std::max(maxGauss,result.gauss);maxCic=std::max(maxCic,result.cic);
    }
    faultChecks(std::string(argv[1])+"/fault-fixture");
    std::ofstream output(std::string(argv[1])+"/checks.json");output<<std::setprecision(17)
      <<"{\"status\":\"PASS\",\"cases\":"<<count<<",\"cg_tolerance\":"<<CG_TOLERANCE
      <<",\"gauss_tolerance\":"<<GAUSS_TOLERANCE<<",\"cic_tolerance\":"<<CIC_TOLERANCE
      <<",\"maximum_electric_relative_l2\":"<<maxE<<",\"maximum_magnetic_relative_l2\":"<<maxB
      <<",\"maximum_gauss_relative_l2\":"<<maxGauss<<",\"maximum_cic_relative_l2\":"<<maxCic
      <<",\"existing_directory_rejected\":true,\"truncated_potential_rejected\":true,\"adjacent_plane_cache_reuse\":true"
      <<",\"scope\":\"Serial fixed-source native E/B, CIC, finite guards, block seams, endpoints and source immutability; no moving PIC, distributed IO or production speed claim.\"}\n";
    require(output.good(),"Cannot write check result");MPI_Finalize();return 0;
  }catch(const std::exception& error){std::cerr<<error.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
}
