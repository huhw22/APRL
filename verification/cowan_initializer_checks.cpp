#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

#include "cowan_reference_field.h"
#include "eb_deposition.h"
#include "eb_initial_field.h"

namespace {

using aprl::Double;
using aprl::EBFieldGrid;
using aprl::EBGridGeometry;
using aprl::FieldVector;
using aprl::RelativisticParticleSI;
using aprl::YeeComponent;

struct Case {
  const char* name;
  std::size_t nz;
  double dx;
  double dy;
  double beta;
  bool guard;
  bool padding;
};

std::vector<RelativisticParticleSI> particlesFor(
    const Case& test, const EBGridGeometry& geometry) {
  const double fraction[5] = {.071, .263, .499, .727, .931};
  const double lower[3] = {test.guard ? 3.1 : 1.1,
    test.guard ? 3.1 : 1.1, test.guard ? 3.1 : 1.1};
  const double cells[3] = {static_cast<double>(geometry.nx),
    static_cast<double>(geometry.ny), static_cast<double>(geometry.nz)};
  const double spacing[3] = {geometry.dx, geometry.dy, geometry.dz};
  std::vector<RelativisticParticleSI> particles;
  for (unsigned iz=0; iz<5; ++iz)
    for (unsigned iy=0; iy<3; ++iy)
      for (unsigned ix=0; ix<4; ++ix) {
        RelativisticParticleSI particle;
        const unsigned index[3] = {ix, iy+1, iz};
        for (unsigned axis=0; axis<3; ++axis)
          particle.position[axis] = spacing[axis]*(lower[axis]+
            (cells[axis]-2*lower[axis])*fraction[index[axis]]);
        // Exercise charges actually inside the transverse CPML, including
        // a CIC cloud on its outermost original-grid vertex.
        if (test.padding && ix==0) particle.position[0]=.37*geometry.dx;
        if (test.padding && iy==0) particle.position[1]=1.43*geometry.dy;
        particle.properVelocity=FieldVector<double>(0.);
        particle.properVelocity[2]=test.beta/std::sqrt(1-test.beta*test.beta);
        particle.weight=1000.+ix+3*iy+7*iz;
        particle.charge=-aprl::SI::elementaryCharge*particle.weight;
        particle.mass=aprl::SI::electronMass*particle.weight;
        particle.id=particle.sourceId=particles.size();
        particle.frozen=false;
        particles.push_back(particle);
      }
  return particles;
}

bool unchanged(const std::vector<RelativisticParticleSI>& before,
               const std::vector<RelativisticParticleSI>& after) {
  if (before.size()!=after.size()) return false;
  for (std::size_t n=0; n<before.size(); ++n) {
    for (unsigned axis=0; axis<3; ++axis)
      if (before[n].position[axis]!=after[n].position[axis] ||
          before[n].properVelocity[axis]!=after[n].properVelocity[axis])
        return false;
    if (before[n].charge!=after[n].charge || before[n].mass!=after[n].mass ||
        before[n].weight!=after[n].weight || before[n].id!=after[n].id ||
        before[n].sourceId!=after[n].sourceId ||
        before[n].frozen!=after[n].frozen ||
        before[n].frozenBetaZ!=after[n].frozenBetaZ) return false;
  }
  return true;
}

struct Spectrum {
  double minimum=std::numeric_limits<double>::infinity();
  double maximum=0.;
  double gradientMaximum=0.;
};

Spectrum operatorSpectrum(const Case& test, const EBGridGeometry& geometry) {
  const auto coefficient=aprl::EBCowanCoefficients::forZDispersion(
    geometry.dx,geometry.dy,geometry.dz);
  const double delta[3]={coefficient.deltaYZ,coefficient.deltaZX,coefficient.deltaXY};
  const std::size_t cells[3]={
    geometry.nx+(test.padding?4:0)-(test.guard&&!test.padding?4:0),
    geometry.ny+(test.padding?4:0)-(test.guard&&!test.padding?4:0),
    geometry.nz-(test.guard?4:0)};
  const double spacing[3]={geometry.dx,geometry.dy,geometry.dz};
  const double velocity=test.beta*aprl::SI::c;
  const double a=velocity*velocity*aprl::SI::epsilon0*aprl::SI::mu0;
  const double pi=std::acos(-1.);
  Spectrum result;
  for(std::size_t k=1;k<cells[2];++k)
    for(std::size_t j=1;j<cells[1];++j)
      for(std::size_t i=1;i<cells[0];++i) {
        const std::size_t mode[3]={i,j,k};
        double cosine[3],wave2[3];
        for(unsigned axis=0;axis<3;++axis) {
          const double angle=pi*mode[axis]/cells[axis];
          cosine[axis]=std::cos(angle);
          const double wave=2*std::sin(.5*angle)/spacing[axis];
          wave2[axis]=wave*wave;
        }
        double eigenvalue=0.,gradient2=0.;
        for(unsigned axis=0;axis<3;++axis) {
          const unsigned p=(axis+1)%3,q=(axis+2)%3;
          const double symbol=coefficient.alpha[axis]+
            2*coefficient.beta[p]*cosine[p]+2*coefficient.beta[q]*cosine[q]+
            4*delta[axis]*cosine[p]*cosine[q]-(axis==2?a:0.);
          eigenvalue+=symbol*wave2[axis];
          gradient2+=symbol*symbol*wave2[axis];
        }
        result.minimum=std::min(result.minimum,eigenvalue);
        result.maximum=std::max(result.maximum,eigenvalue);
        result.gradientMaximum=std::max(result.gradientMaximum,std::sqrt(gradient2));
      }
  if(!(result.minimum>0.)) throw std::runtime_error("Invalid oracle operator spectrum");
  return result;
}

double sourceNorm(EBFieldGrid& full,
                  const std::vector<RelativisticParticleSI>& particles) {
  const auto& geometry=full.geometry();
  YeeComponent rho(geometry.nx+1,geometry.ny+1,geometry.nz+1);
  rho.fill(0.);
  aprl::ChargeConservingCurrentDepositor depositor(full,FieldVector<double>(0.));
  for(const auto& particle:particles)
    depositor.depositCharge(particle.position,particle.charge,rho);
  long double norm=0.;
  for(std::size_t n=0;n<rho.size();++n) {
    const long double value=rho.data()[n]/aprl::SI::epsilon0;
    norm+=value*value;
  }
  return std::sqrt(norm);
}

struct Difference {
  double l2=0.;
  double maximum=0.;
  double referenceL2=0.;
  double referenceMaximum=0.;
  double sharedMaximum=0.;
};

Difference compareFields(const std::array<const YeeComponent*,3>& local,
                         const std::array<const YeeComponent*,3>& reference,
                         std::size_t offset, std::size_t localCells,
                         int rank, int ranks) {
  long double localSums[2]={};
  double localMaxima[3]={};
  for(unsigned axis=0;axis<3;++axis) {
    const auto& field=*local[axis];
    const auto& oracle=*reference[axis];
    for(std::size_t k=0;k<field.nz();++k)
      for(std::size_t j=0;j<field.ny();++j)
        for(std::size_t i=0;i<field.nx();++i) {
          const double wanted=oracle(i,j,k+offset);
          const double delta=field(i,j,k)-wanted;
          if(!std::isfinite(delta)||!std::isfinite(wanted))
            throw std::runtime_error("Non-finite field in initializer comparison");
          const bool duplicateTop=k==localCells && rank+1<ranks;
          if(!duplicateTop) {
            localSums[0]+=static_cast<long double>(delta)*delta;
            localSums[1]+=static_cast<long double>(wanted)*wanted;
          }
          localMaxima[0]=std::max(localMaxima[0],std::abs(delta));
          localMaxima[1]=std::max(localMaxima[1],std::abs(wanted));
          if((k==0&&rank>0)||duplicateTop)
            localMaxima[2]=std::max(localMaxima[2],std::abs(delta));
        }
  }
  long double sums[2]={};double maxima[3]={};
  MPI_Allreduce(localSums,sums,2,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(localMaxima,maxima,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  Difference result;
  result.l2=std::sqrt(sums[0]);result.referenceL2=std::sqrt(sums[1]);
  result.maximum=maxima[0];result.referenceMaximum=maxima[1];
  result.sharedMaximum=maxima[2];return result;
}

void runCase(const Case& test,double tolerance,int rank,int ranks) {
  const double dz=1e-6;
  const EBGridGeometry global(12,14,test.nz,test.dx*dz,test.dy*dz,dz,
    dz/aprl::SI::c,aprl::EBMaxwellSolver::CowanZ);
  EBGridGeometry slab(global);slab.nz=global.nz/static_cast<std::size_t>(ranks);
  const std::size_t offset=static_cast<std::size_t>(rank)*slab.nz;
  if(slab.nz*static_cast<std::size_t>(ranks)!=global.nz)
    throw std::runtime_error("Test mesh does not divide into requested ranks");
  const auto all=particlesFor(test,global);
  std::vector<RelativisticParticleSI> particles;
  for(const auto& particle:all) {
    const std::size_t cell=static_cast<std::size_t>(std::floor(particle.position[2]/dz));
    if(cell>=offset&&cell<offset+slab.nz)particles.push_back(particle);
  }
  const auto before=particles;
  FieldVector<double> origin(0.);origin[2]=offset*dz;
  const std::size_t guard[3]={test.guard?2U:0U,test.guard?2U:0U,test.guard?2U:0U};
  EBFieldGrid distributed(slab),reference(global);
  const auto report=aprl::EBGaussFieldInitializer::initialize(distributed,global,
    offset,origin,particles,aprl::InitialSelfFieldModel::CowanRelativisticPoisson,
    1e-8,guard,tolerance,5000,MPI_COMM_WORLD,test.padding);
  double referenceGauss=0.;
  if(test.guard) {
    const auto oracle=aprl::EBGaussFieldInitializer::initialize(reference,global,
      0,FieldVector<double>(0.),all,aprl::InitialSelfFieldModel::CowanRelativisticPoisson,
      1e-8,guard,tolerance,5000,MPI_COMM_SELF,test.padding);
    referenceGauss=oracle.relativeResidual;
  } else {
    const auto oracle=aprl::verification::initializeCowanReference(reference,
      FieldVector<double>(0.),all,test.beta,MPI_COMM_SELF);
    referenceGauss=oracle.relativeGaussL2;
  }
  const auto electric=compareFields({{&distributed.ex(),&distributed.ey(),&distributed.ez()}},
    {{&reference.ex(),&reference.ey(),&reference.ez()}},offset,slab.nz,rank,ranks);
  const auto magnetic=compareFields({{&distributed.bx(),&distributed.by(),&distributed.bz()}},
    {{&reference.bx(),&reference.by(),&reference.bz()}},offset,slab.nz,rank,ranks);
  const auto spectrum=operatorSpectrum(test,global);
  const double condition=spectrum.maximum/spectrum.minimum;
  // A e_phi=r and E_error=G e_phi imply ||E_error||_2 <=
  // ||G||_2/lambda_min(A) * ||r||_2. The production post-Gauss gate
  // permits 4*tolerance; two CG solutions need twice that allowance.
  // The 256*epsilon*kappa floor covers small-grid transform/stencil
  // roundoff conservatively; it is specified before collecting results.
  const double residualAllowance=(test.guard?8.:4.)*tolerance;
  const double roundoff=256*std::numeric_limits<double>::epsilon()*condition;
  const double electricBound=spectrum.gradientMaximum/spectrum.minimum*
    residualAllowance*sourceNorm(reference,all)+roundoff*electric.referenceL2;
  const double magneticScale=std::abs(test.beta*aprl::SI::c*
    aprl::SI::epsilon0*aprl::SI::mu0);
  // Convex z interpolation has l2 operator norm <=1, so B inherits this
  // bound with the magnetic scale. Static B must be exactly zero.
  const double magneticBound=magneticScale*electricBound+
    roundoff*magnetic.referenceL2;
  int localOk=unchanged(before,particles)&&electric.l2<=electricBound&&
    electric.maximum<=electricBound&&magnetic.l2<=magneticBound&&
    magnetic.maximum<=magneticBound&&report.relativeResidual<=4*tolerance&&
    std::abs(report.meanBetaZ-test.beta)<=32*std::numeric_limits<double>::epsilon();
  int passed=0;MPI_Allreduce(&localOk,&passed,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if(rank==0)std::cout<<std::setprecision(17)
    <<"{\"case\":\""<<test.name<<"\",\"ranks\":"<<ranks
    <<",\"local_z_cells\":"<<slab.nz
    <<",\"oracle\":\""<<(test.guard?"serial-production-CG":"independent-DST")
    <<"\",\"iterations\":"<<report.iterations
    <<",\"gauss_l2\":"<<report.relativeResidual
    <<",\"oracle_gauss_l2\":"<<referenceGauss
    <<",\"condition\":"<<condition
    <<",\"E_error_l2\":"<<electric.l2<<",\"E_bound_l2\":"<<electricBound
    <<",\"E_relative_l2\":"<<electric.l2/electric.referenceL2
    <<",\"E_shared_error_max\":"<<electric.sharedMaximum
    <<",\"B_error_l2\":"<<magnetic.l2<<",\"B_bound_l2\":"<<magneticBound
    <<",\"B_shared_error_max\":"<<magnetic.sharedMaximum
    <<",\"passed\":"<<(passed?"true":"false")<<"}\n";
  if(!passed)throw std::runtime_error(std::string("Cowan initializer verification failed: ")+test.name);
}

void checkRejections(double tolerance,int rank) {
  const std::size_t guard[3]={0,0,0};
  for(unsigned check=0;check<2;++check) {
    const Case test={"negative",16,1.,1.,check==0?.23:.9,false,false};
    const EBGridGeometry geometry(12,14,16,1e-6,1e-6,1e-6,
      (check==0?.5e-6:1e-6)/aprl::SI::c,
      check==0?aprl::EBMaxwellSolver::Yee:aprl::EBMaxwellSolver::CowanZ);
    EBFieldGrid field(geometry);
    const auto particles=particlesFor(test,geometry);
    const std::string expected=check==0?"requires the cowan-z Maxwell solver":
      "positive-definiteness bound";
    int rejected=0;
    try {
      aprl::EBGaussFieldInitializer::initialize(field,geometry,0,
        FieldVector<double>(0.),particles,
        aprl::InitialSelfFieldModel::CowanRelativisticPoisson,
        1e-8,guard,tolerance,5000,MPI_COMM_SELF);
    } catch(const std::exception& error) {
      rejected=std::string(error.what()).find(expected)!=std::string::npos;
    }
    int passed=0;MPI_Allreduce(&rejected,&passed,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    if(rank==0)std::cout<<"{\"case\":\"reject-"
      <<(check==0?"yee":"non-spd")<<"\",\"passed\":"
      <<(passed?"true":"false")<<"}\n";
    if(!passed)throw std::runtime_error("Missing expected Cowan initializer rejection: "+expected);
  }
}

} // namespace

int main(int argc,char** argv) {
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int rank=0,ranks=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  try {
    double tolerance=1e-11;
    if(argc==3&&std::string(argv[1])=="--tolerance")tolerance=std::stod(argv[2]);
    else if(argc!=1)throw std::invalid_argument("usage: aprl_cowan_initializer_checks [--tolerance 1e-11]");
    if(!(tolerance>0.&&tolerance<1e-5)||!std::isfinite(tolerance))
      throw std::invalid_argument("tolerance must be finite and between zero and 1e-5");
    if(ranks!=1&&ranks!=2&&ranks!=4)
      throw std::invalid_argument("verification requires 1, 2, or 4 MPI ranks");
    const Case cases[]={
      {"isotropic-static",16,1.,1.,0.,false,false},
      {"isotropic-moving",16,1.,1.,.23,false,false},
      {"anisotropic-static",16,4.,2.,0.,false,false},
      {"anisotropic-moving",16,4.,2.,.23,false,false},
      // EBFieldGrid requires at least two cells per local direction.
      // Eight global z cells exercise that minimum with four ranks.
      {"thin-static",8,1.,1.,0.,false,false},
      {"thin-moving",8,4.,2.,.23,false,false},
      {"cpml-guard-moving",16,4.,2.,.23,true,false},
      {"cpml-padding-moving",16,4.,2.,.23,true,true}};
    for(const auto& test:cases)runCase(test,tolerance,rank,ranks);
    checkRejections(tolerance,rank);
    if(rank==0)std::cout<<"Cowan initializer checks PASS: 8 comparisons + 2 rejections, "<<ranks<<" ranks\n";
    MPI_Finalize();return 0;
  } catch(const std::exception& error) {
    std::cerr<<"rank "<<rank<<": "<<error.what()<<'\n';
    MPI_Abort(MPI_COMM_WORLD,1);return 1;
  }
}
