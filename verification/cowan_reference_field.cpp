#include "cowan_reference_field.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <fftw3.h>
#include <mpi.h>

#include "eb_deposition.h"

namespace aprl {
namespace verification {
namespace {

struct Nodes {
  std::size_t cells[3];
  std::vector<double> data;
  explicit Nodes(const EBGridGeometry& g)
      : cells{g.nx, g.ny, g.nz}, data((g.nx+1)*(g.ny+1)*(g.nz+1), 0.) {}
  std::size_t index(std::size_t i, std::size_t j, std::size_t k) const {
    return (k*(cells[1]+1)+j)*(cells[0]+1)+i;
  }
  double& operator()(std::size_t i, std::size_t j, std::size_t k) {
    return data[index(i,j,k)];
  }
  double operator()(std::size_t i, std::size_t j, std::size_t k) const {
    return data[index(i,j,k)];
  }
  // Odd extension of zero Dirichlet potential. Its forward normal derivative
  // has exactly the even/constant first-ghost extension used for normal E;
  // tangential E has the required odd extension at each PEC surface.
  double extended(std::ptrdiff_t i, std::ptrdiff_t j, std::ptrdiff_t k) const {
    std::ptrdiff_t q[3]={i,j,k};
    double sign=1.;
    for(unsigned a=0;a<3;++a) {
      const auto n=static_cast<std::ptrdiff_t>(cells[a]);
      if(q[a]<0) { q[a]=-q[a]; sign=-sign; }
      if(q[a]>n) { q[a]=2*n-q[a]; sign=-sign; }
      if(q[a]<0||q[a]>n)
        throw std::out_of_range("Cowan reference potential extension exceeds one halo");
      if(q[a]==0||q[a]==n) return 0.;
    }
    return sign*(*this)(q[0],q[1],q[2]);
  }
};

double deltaFor(const EBCowanCoefficients& c,unsigned axis) {
  return axis==0?c.deltaYZ:(axis==1?c.deltaZX:c.deltaXY);
}

double smoothed(const Nodes& phi,const EBCowanCoefficients& c,unsigned axis,
                std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k) {
  const unsigned a=(axis+1)%3,b=(axis+2)%3;
  const std::ptrdiff_t q[3]={i,j,k};
  double value=c.alpha[axis]*phi.extended(i,j,k);
  for(int sign:{-1,1}) {
    std::ptrdiff_t qa[3]={i,j,k},qb[3]={i,j,k};
    qa[a]+=sign; qb[b]+=sign;
    value+=c.beta[a]*phi.extended(qa[0],qa[1],qa[2]);
    value+=c.beta[b]*phi.extended(qb[0],qb[1],qb[2]);
  }
  for(int sa:{-1,1}) for(int sb:{-1,1}) {
    std::ptrdiff_t point[3]={q[0],q[1],q[2]};
    point[a]+=sa; point[b]+=sb;
    value+=deltaFor(c,axis)*phi.extended(point[0],point[1],point[2]);
  }
  return value;
}

class SineTransform {
 public:
  SineTransform(const EBGridGeometry& g,std::vector<double>& work)
      : plan_(fftw_plan_r2r_3d(static_cast<int>(g.nz-1),
          static_cast<int>(g.ny-1),static_cast<int>(g.nx-1),
          work.data(),work.data(),FFTW_RODFT00,FFTW_RODFT00,FFTW_RODFT00,
          FFTW_ESTIMATE)) {
    if(!plan_) throw std::runtime_error("Cannot create Cowan-reference DST-I plan");
  }
  ~SineTransform() { fftw_destroy_plan(plan_); }
  void execute() { fftw_execute(plan_); }
 private:
  fftw_plan plan_;
  SineTransform(const SineTransform&)=delete;
  SineTransform& operator=(const SineTransform&)=delete;
};

}  // namespace

CowanMagneticDivergenceReport cowanMagneticDivergence(const EBFieldGrid& fields) {
  const auto& g=fields.geometry();
  const double h[3]={g.dx,g.dy,g.dz};
  const YeeComponent* magnetic[3]={&fields.bx(),&fields.by(),&fields.bz()};
  const auto c=EBCowanCoefficients::forZDispersion(g.dx,g.dy,g.dz);
  const auto smooth=[&c](const YeeComponent& b,unsigned axis,
      std::size_t i,std::size_t j,std::size_t k) {
    const unsigned a=(axis+1)%3,d=(axis+2)%3;
    const std::ptrdiff_t q[3]={static_cast<std::ptrdiff_t>(i),
      static_cast<std::ptrdiff_t>(j),static_cast<std::ptrdiff_t>(k)};
    double value=c.alpha[axis]*b(i,j,k);
    for(int sign:{-1,1}) {
      std::ptrdiff_t qa[3]={q[0],q[1],q[2]},qd[3]={q[0],q[1],q[2]};
      qa[a]+=sign;qd[d]+=sign;
      value+=c.beta[a]*b(qa[0],qa[1],qa[2]);
      value+=c.beta[d]*b(qd[0],qd[1],qd[2]);
    }
    for(int sa:{-1,1}) for(int sd:{-1,1}) {
      std::ptrdiff_t point[3]={q[0],q[1],q[2]};
      point[a]+=sa;point[d]+=sd;
      value+=deltaFor(c,axis)*b(point[0],point[1],point[2]);
    }
    return value;
  };
  CowanMagneticDivergenceReport report;
  long double ordinary2=0,modified2=0,scale2=0;
  for(std::size_t k=2;k+2<g.nz;++k) for(std::size_t j=2;j+2<g.ny;++j)
    for(std::size_t i=2;i+2<g.nx;++i) {
      double ordinary=0,modified=0;
      long double localScale2=0;
      for(unsigned axis=0;axis<3;++axis) {
        const std::size_t upper[3]={i+(axis==0),j+(axis==1),k+(axis==2)};
        const auto& b=*magnetic[axis];
        const double lowerValue=b(i,j,k),upperValue=b(upper[0],upper[1],upper[2]);
        ordinary+=(upperValue-lowerValue)/h[axis];
        modified+=(smooth(b,axis,upper[0],upper[1],upper[2])-smooth(b,axis,i,j,k))/h[axis];
        localScale2+=(static_cast<long double>(lowerValue)*lowerValue+
          static_cast<long double>(upperValue)*upperValue)/(2.L*h[axis]*h[axis]);
      }
      ++report.samples;
      ordinary2+=static_cast<long double>(ordinary)*ordinary;
      modified2+=static_cast<long double>(modified)*modified;
      scale2+=localScale2;
      report.ordinaryMaximum=std::max(report.ordinaryMaximum,std::abs(ordinary));
      report.modifiedMaximum=std::max(report.modifiedMaximum,std::abs(modified));
      report.normalizationMaximum=std::max(report.normalizationMaximum,
        static_cast<double>(std::sqrt(localScale2)));
    }
  if(report.samples)report.normalizationRms=std::sqrt(scale2/report.samples);
  const auto ratio=[](long double numerator,long double denominator) {
    if(denominator>0)return static_cast<double>(std::sqrt(numerator/denominator));
    return numerator==0?0.:std::numeric_limits<double>::infinity();
  };
  report.ordinaryRelativeL2=ratio(ordinary2,scale2);
  report.modifiedRelativeL2=ratio(modified2,scale2);
  const double maxScale=report.normalizationMaximum;
  report.ordinaryRelativeMaximum=maxScale>0?report.ordinaryMaximum/maxScale:
    (report.ordinaryMaximum==0?0.:std::numeric_limits<double>::infinity());
  report.modifiedRelativeMaximum=maxScale>0?report.modifiedMaximum/maxScale:
    (report.modifiedMaximum==0?0.:std::numeric_limits<double>::infinity());
  return report;
}

CowanReferenceReport initializeCowanReference(EBFieldGrid& fields,
    const FieldVector<double>& origin,
    const std::vector<RelativisticParticleSI>& particles,double beta,
    MPI_Comm communicator) {
  const auto start=std::chrono::steady_clock::now();
  int ranks=0;
  if(communicator==MPI_COMM_NULL)
    throw std::invalid_argument("Cowan reference communicator cannot be null");
  MPI_Comm_size(communicator,&ranks);
  const auto& g=fields.geometry();
  if(ranks!=1||g.solver!=EBMaxwellSolver::CowanZ)
    throw std::invalid_argument("Cowan reference requires one MPI rank and CowanZ");
  if(!std::isfinite(beta)||std::abs(beta)>=1.)
    throw std::invalid_argument("Cowan reference beta must be finite and subluminal");
  for(auto n:{g.nx,g.ny,g.nz})
    if(n<3||n>static_cast<std::size_t>(INT_MAX))
      throw std::invalid_argument("Cowan reference requires 3 <= cells <= INT_MAX");
  const std::size_t maximum=std::numeric_limits<std::size_t>::max()/sizeof(double);
  if((g.nx+1)>maximum/(g.ny+1)||(g.nx+1)*(g.ny+1)>maximum/(g.nz+1))
    throw std::overflow_error("Cowan reference allocation overflows size_t");
  const double h[3]={g.dx,g.dy,g.dz};
  const std::size_t cells[3]={g.nx,g.ny,g.nz};
  for(const auto& p:particles) for(unsigned axis=0;axis<3;++axis)
    if(!(p.position[axis]>origin[axis]+h[axis]&&
         p.position[axis]<origin[axis]+(cells[axis]-1)*h[axis]))
      throw std::invalid_argument("Cowan reference requires particles one cell inside PEC");

  fields.clearFields();
  fields.clearCurrent();
  YeeComponent rho(g.nx+1,g.ny+1,g.nz+1);
  rho.fill(0.);
  ChargeConservingCurrentDepositor depositor(fields,origin);
  CowanReferenceReport report;
  long double charge=0;
  for(const auto& p:particles) {
    depositor.depositCharge(p.position,p.charge,rho);
    charge+=p.charge;
  }
  report.sourceCharge=static_cast<double>(charge);
  const std::size_t nx=g.nx-1,ny=g.ny-1,nz=g.nz-1;
  const auto index=[nx,ny](std::size_t i,std::size_t j,std::size_t k) {
    return (k*ny+j)*nx+i;
  };
  std::vector<double> work(nx*ny*nz);
  for(std::size_t k=1;k<g.nz;++k) for(std::size_t j=1;j<g.ny;++j)
    for(std::size_t i=1;i<g.nx;++i)
      work[index(i-1,j-1,k-1)]=rho(i,j,k)/SI::epsilon0;
  SineTransform transform(g,work);
  transform.execute();

  const auto coefficient=EBCowanCoefficients::forZDispersion(g.dx,g.dy,g.dz);
  const double v=beta*SI::c;
  const double a=v*v*SI::epsilon0*SI::mu0;
  const double pi=std::acos(-1.);
  report.minimumDenominator=std::numeric_limits<double>::infinity();
  report.minimumSmoothingSymbol=std::numeric_limits<double>::infinity();
  report.minimumLongitudinalSymbol=std::numeric_limits<double>::infinity();
  for(std::size_t k=1;k<g.nz;++k) for(std::size_t j=1;j<g.ny;++j)
    for(std::size_t i=1;i<g.nx;++i) {
      const double angle[3]={pi*i/g.nx,pi*j/g.ny,pi*k/g.nz};
      double wave2[3],cosine[3],symbol[3];
      for(unsigned d=0;d<3;++d) {
        const double wave=2*std::sin(.5*angle[d])/h[d];
        wave2[d]=wave*wave;
        cosine[d]=std::cos(angle[d]);
      }
      for(unsigned d=0;d<3;++d) {
        const unsigned p=(d+1)%3,q=(d+2)%3;
        symbol[d]=coefficient.alpha[d]+2*coefficient.beta[p]*cosine[p]+
          2*coefficient.beta[q]*cosine[q]+4*deltaFor(coefficient,d)*cosine[p]*cosine[q];
        report.minimumSmoothingSymbol=std::min(report.minimumSmoothingSymbol,symbol[d]);
      }
      const double longitudinal=symbol[2]-a;
      const double denominator=symbol[0]*wave2[0]+symbol[1]*wave2[1]+longitudinal*wave2[2];
      report.minimumDenominator=std::min(report.minimumDenominator,denominator);
      report.minimumLongitudinalSymbol=std::min(report.minimumLongitudinalSymbol,longitudinal);
      const double scale=wave2[0]+wave2[1]+wave2[2];
      if(!(denominator>64*std::numeric_limits<double>::epsilon()*scale))
        throw std::runtime_error("Cowan reference elliptic denominator is nonpositive or nearly singular");
      work[index(i-1,j-1,k-1)]/=denominator;
    }
  transform.execute();
  report.transforms=2;
  const double normalization=8.*static_cast<double>(g.nx)*g.ny*g.nz;
  Nodes phi(g);
  for(std::size_t k=1;k<g.nz;++k) for(std::size_t j=1;j<g.ny;++j)
    for(std::size_t i=1;i<g.nx;++i)
      phi(i,j,k)=work[index(i-1,j-1,k-1)]/normalization;

  YeeComponent* electric[3]={&fields.ex(),&fields.ey(),&fields.ez()};
  for(unsigned d=0;d<3;++d) {
    Nodes filtered(g);
    for(std::size_t k=0;k<=g.nz;++k) for(std::size_t j=0;j<=g.ny;++j)
      for(std::size_t i=0;i<=g.nx;++i)
        filtered(i,j,k)=smoothed(phi,coefficient,d,i,j,k)-(d==2?a*phi(i,j,k):0.);
    auto& e=*electric[d];
    for(std::size_t k=0;k<e.nz();++k) for(std::size_t j=0;j<e.ny();++j)
      for(std::size_t i=0;i<e.nx();++i) {
        const std::size_t next[3]={i+(d==0),j+(d==1),k+(d==2)};
        e(i,j,k)=-(filtered(next[0],next[1],next[2])-filtered(i,j,k))/h[d];
      }
  }

  // The translating reference is evaluated at B's stored time -dt/2, using
  // the existing initializer's spatial half-cell plus time-shift interpolation.
  const double upperWeight=.5+.5*v*g.dt/g.dz;
  if(upperWeight<0.||upperWeight>1.)
    throw std::invalid_argument("Cowan reference B half-step shift exceeds one cell");
  const double lowerWeight=1.-upperWeight,magneticScale=v*SI::epsilon0*SI::mu0;
  for(std::size_t k=0;k<g.nz;++k) for(std::size_t j=0;j<g.ny;++j)
    for(std::size_t i=0;i<=g.nx;++i)
      fields.bx()(i,j,k)=-magneticScale*(lowerWeight*fields.ey()(i,j,k)+upperWeight*fields.ey()(i,j,k+1));
  for(std::size_t k=0;k<g.nz;++k) for(std::size_t j=0;j<=g.ny;++j)
    for(std::size_t i=0;i<g.nx;++i)
      fields.by()(i,j,k)=magneticScale*(lowerWeight*fields.ex()(i,j,k)+upperWeight*fields.ex()(i,j,k+1));
  fields.bz().fill(0.);

  long double error2=0,source2=0;
  for(std::size_t k=1;k<g.nz;++k) for(std::size_t j=1;j<g.ny;++j)
    for(std::size_t i=1;i<g.nx;++i) {
      const double source=rho(i,j,k)/SI::epsilon0;
      const double residual=(fields.ex()(i,j,k)-fields.ex()(i-1,j,k))/g.dx+
        (fields.ey()(i,j,k)-fields.ey()(i,j-1,k))/g.dy+
        (fields.ez()(i,j,k)-fields.ez()(i,j,k-1))/g.dz-source;
      report.maximumGaussResidual=std::max(report.maximumGaussResidual,std::abs(residual));
      error2+=static_cast<long double>(residual)*residual;
      source2+=static_cast<long double>(source)*source;
    }
  report.relativeGaussL2=source2>0?std::sqrt(error2/source2):0.;
  report.initializationSeconds=std::chrono::duration<double>(
      std::chrono::steady_clock::now()-start).count();
  return report;
}

}  // namespace verification
}  // namespace aprl
