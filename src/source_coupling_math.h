#ifndef APRL_SOURCE_COUPLING_MATH_H
#define APRL_SOURCE_COUPLING_MATH_H
// Local equations validated by the independent source-coupling prototype.
#include "eb_particles.h"
#include "eb_deposition.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <mpi.h>
namespace source_coupling_math {
using namespace aprl;
using Vec=FieldVector<Double>;
using Electric=std::array<YeeComponent,3>;
inline Double scalar(const Vec& a,const Vec& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline long double gammaLD(const Vec& u) {return std::sqrt(1.L+(long double)u[0]*u[0]+(long double)u[1]*u[1]+(long double)u[2]*u[2]);}
inline long double deltaGamma(const Vec& u,const Vec& v) {
 long double s=0;for(unsigned a=0;a<3;++a)s+=((long double)v[a]-u[a])*((long double)v[a]+u[a]);
 return s/(gammaLD(u)+gammaLD(v));
}
inline long double chargeWork(Double q,const Vec& e,const Vec& x,const Vec& y) {
 long double s=0;for(unsigned a=0;a<3;++a)s+=(long double)q*e[a]*((long double)y[a]-x[a]);return s;
}
inline Electric copyElectric(const EBFieldGrid& f) {return {{f.ex(),f.ey(),f.ez()}};}
class NativeHalo {
 bool lower_=false,upper_=false;
 std::vector<Double> down_[3],up_[3];
public:
 void exchange(const EBFieldGrid& f,MPI_Comm comm,bool magnetic=false) {
  int rank=0,size=1;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);lower_=rank>0;upper_=rank+1<size;
  if(size==1)return;
  const YeeComponent* v[3]={magnetic?&f.bx():&f.ex(),magnetic?&f.by():&f.ey(),magnetic?&f.bz():&f.ez()};
  for(unsigned a=0;a<3;++a){const auto& q=*v[a];const std::size_t count=q.nx()*q.ny();if(count>std::size_t(std::numeric_limits<int>::max()))throw std::overflow_error("Coupling halo too large");
   down_[a].resize(count);up_[a].resize(count);const bool node=magnetic?a==2:a<2;
   MPI_Sendrecv(q.data()+(q.nz()-(node?2:1))*count,int(count),MPI_DOUBLE,upper_?rank+1:MPI_PROC_NULL,11000+2*a,down_[a].data(),int(count),MPI_DOUBLE,lower_?rank-1:MPI_PROC_NULL,11000+2*a,comm,MPI_STATUS_IGNORE);
   MPI_Sendrecv(q.data()+(node?1:0)*count,int(count),MPI_DOUBLE,lower_?rank-1:MPI_PROC_NULL,11001+2*a,up_[a].data(),int(count),MPI_DOUBLE,upper_?rank+1:MPI_PROC_NULL,11001+2*a,comm,MPI_STATUS_IGNORE);
  }
 }
 bool contains(long long k,std::size_t nz)const{return(k==-1&&lower_)||(k==static_cast<long long>(nz)&&upper_);}
 Double value(const YeeComponent& v,unsigned a,long long i,long long j,long long k)const {
  if(i<0||j<0||i>=static_cast<long long>(v.nx())||j>=static_cast<long long>(v.ny()))throw std::out_of_range("Coupling transverse native sample outside grid");
  if(k>=0&&k<static_cast<long long>(v.nz()))return v(i,j,k);
  if(!contains(k,v.nz()))throw std::out_of_range("Coupling sample outside one native MPI collar");
  return (k<0?down_[a]:up_[a])[j*v.nx()+i];
 }
};
struct MagneticSample {Vec value=Vec(0.),maximum=Vec(0.),point=Vec(0.);int cell[2][3];};
inline MagneticSample magneticSample(const EBFieldGrid& f,const Vec& origin,const Vec& position,bool bounds=true,const NativeHalo* halo=nullptr) {
 MagneticSample result;result.point=position;const auto& g=f.geometry();const Double ds[3]={g.dx,g.dy,g.dz};Double fraction[2][3];
 for(unsigned half=0;half<2;++half)for(unsigned d=0;d<3;++d){const Double x=(position[d]-origin[d])/ds[d]-.5*half;result.cell[half][d]=(int)std::floor(x);fraction[half][d]=x-result.cell[half][d];}
 const YeeComponent* fields[3]={&f.bx(),&f.by(),&f.bz()};
 for(unsigned axis=0;axis<3;++axis) {
  const unsigned hx=axis==0?0:1,hy=axis==1?0:1,hz=axis==2?0:1;
  for(unsigned k=0;k<2;++k)for(unsigned j=0;j<2;++j)for(unsigned i=0;i<2;++i) {
   const Double x=i?fraction[hx][0]:1-fraction[hx][0],y=j?fraction[hy][1]:1-fraction[hy][1],z=k?fraction[hz][2]:1-fraction[hz][2];
   const Double b=halo?halo->value(*fields[axis],axis,result.cell[hx][0]+int(i),result.cell[hy][1]+int(j),result.cell[hz][2]+int(k)):(*fields[axis])(result.cell[hx][0]+i,result.cell[hy][1]+j,result.cell[hz][2]+k);
   result.value[axis]+=x*y*z*b;if(bounds)result.maximum[axis]=std::max(result.maximum[axis],std::abs(b));
  }
 }
 return result;
}
// A rigorous trilinear Lipschitz bound while both points use the same cell.
// Otherwise evaluate the actual new sample. B nodes do not change in this loop.
inline Double magneticDifferenceBound(const EBFieldGrid& f,const Vec& origin,const Vec& point,const MagneticSample& old) {
 const auto& g=f.geometry();const Double ds[3]={g.dx,g.dy,g.dz};bool same=true;
 for(unsigned half=0;half<2;++half)for(unsigned d=0;d<3;++d)if((int)std::floor((point[d]-origin[d])/ds[d]-.5*half)!=old.cell[half][d])same=false;
 if(!same){const auto fresh=magneticSample(f,origin,point);Vec difference=fresh.value;difference-=old.value;return std::sqrt(difference.norm2())+64*std::numeric_limits<Double>::epsilon()*(std::sqrt(old.maximum.norm2())+std::sqrt(fresh.maximum.norm2()));}
 Double sum=0;for(unsigned d=0;d<3;++d)sum+=std::abs(point[d]-old.point[d])/ds[d];return (2*sum+64*std::numeric_limits<Double>::epsilon())*std::sqrt(old.maximum.norm2());
}

// Exact analytic orbit average of the first-order depositor-adjoint electric
// interpolant. Inside each cell the transverse product is quadratic along a
// straight chord: midpoint product plus signed d_eta*d_zeta/12.
inline Vec orbitElectric(const EBFieldGrid& f,const Electric* old,
 const Vec& origin,const Vec& begin,const Vec& end,const NativeHalo* halo=nullptr,const NativeHalo* oldHalo=nullptr) {
 const auto& g=f.geometry();const Double spacing[3]={g.dx,g.dy,g.dz};
 Double x[3],d[3],cuts[32]={0.,1.};unsigned count=2;
 for(unsigned a=0;a<3;++a) {
  x[a]=(begin[a]-origin[a])/spacing[a];d[a]=(end[a]-begin[a])/spacing[a];
  if(d[a]!=0.) {
   const Double lo=std::min(x[a],x[a]+d[a]),hi=std::max(x[a],x[a]+d[a]);
   for(long long k=(long long)std::floor(lo)+1;k<hi;++k) {
    if(count==32)throw std::runtime_error("Prototype orbit exceeds cell-crossing budget");
    const Double t=(k-x[a])/d[a];if(t>0.&&t<1.)cuts[count++]=t;
   }
  }
 }
 std::sort(cuts,cuts+count);Vec result(0.);const YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};
 for(unsigned n=1;n<count;++n) {
  const Double width=cuts[n]-cuts[n-1];if(width<=0.)continue;
  const Double middle=.5*(cuts[n]+cuts[n-1]);long long cell[3];Double w[3][2];
  for(unsigned a=0;a<3;++a) {
   cell[a]=(long long)std::floor(x[a]+middle*d[a]);
   if((cell[a]<0||cell[a]>=(long long)(a==0?g.nx:(a==1?g.ny:g.nz))) && !(a==2&&halo&&halo->contains(cell[a],g.nz)))throw std::out_of_range("Coupling orbit outside grid collar");
   w[a][1]=x[a]+middle*d[a]-cell[a];w[a][0]=1.-w[a][1];
  }
  for(unsigned a=0;a<3;++a) {
   const unsigned b=(a+1)%3,c=(a+2)%3;
   for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j) {
    long long node[3]={cell[0],cell[1],cell[2]};node[b]+=i;node[c]+=j;
    Double value=halo?halo->value(*e[a],a,node[0],node[1],node[2]):(*e[a])(node[0],node[1],node[2]);
    if(old)value=.5*(value+(oldHalo?oldHalo->value((*old)[a],a,node[0],node[1],node[2]):(*old)[a](node[0],node[1],node[2])));
    const Double correction=(i==j?1.:-1.)*d[b]*d[c]*width*width/12.;
    result[a]+=width*(w[b][i]*w[c][j]+correction)*value;
   }
  }
 }
 return result;
}

// Lapenta--Markidis mean velocity, with the local scalar equation solved to
// roundoff. It is the momentum equation in RelSIM Appendix A, not its approximate
// global response. The magnetic force and position share the energy velocity.
inline unsigned momentum(Vec& u,const Vec& e,const Vec& magnetic,Double q,Double m,Double h,Double knownGamma=0.) {
 const Vec before(u);const Double g0=knownGamma>0.?knownGamma:BoostFrameTransform::gammaFromProperVelocity(before);
 Vec kick(e);kick*=q*h/(2*m*SI::c);Vec b(magnetic);b*=q*h/(2*m);
 Vec a(before);a+=kick;const Vec axb=cross(a,b);
 const Double b2=scalar(b,b),ka=scalar(kick,a),kb=scalar(kick,b),ab=scalar(a,b),kcross=scalar(kick,axb);
 const Double first=scalar(kick,before)/g0;
 Double G=g0+first+(scalar(kick,kick)+kcross/g0-first*first)/g0;
 G=std::max(G,.5*(g0+1.));unsigned iteration=0;
 for(;iteration<12;++iteration) {
  const Double den=G*G+b2,numerator=G*G*ka+G*kcross+ab*kb;
  const Double f=G*(G-g0)-numerator/den;
  const Double derivative=2*G-g0-((2*G*ka+kcross)*den-2*G*numerator)/(den*den);
  const Double correction=f/derivative;const Double next=G-correction;
  if(!std::isfinite(next)||next<.5*(g0+1.)-1e-13)throw std::runtime_error("Nonphysical local mean-gamma root");
  G=std::max(next,.5*(g0+1.));
  if(std::abs(correction)<=8*std::numeric_limits<Double>::epsilon()*std::max(1.,G))break;
 }
 if(iteration==12)throw std::runtime_error("Local mean-gamma root did not converge");
 const Double den=G*G+b2;
 for(unsigned k=0;k<3;++k)u[k]=2*(G*G*a[k]+ab*b[k]+G*axb[k])/den-before[k];
 return iteration+1;
}
struct PushWork {
 long double grid=0.,prescribed=0.,absoluteGrid=0.;unsigned rootIterations=0;
};
inline PushWork push(RelativisticParticleSI& particle,const Vec& gridE,const Vec& gridB,
 const SIFieldSourceSet& sources,const BoostFrameTransform& frame,Double time,Double dt,unsigned substeps,bool prescribedWork=true,bool gridWork=true) {
 PushWork work;const Double h=dt/substeps;const Vec firstPosition(particle.position);Double absoluteGrid=0.;
 for(unsigned s=0;s<substeps;++s) {
  const Vec beforeU(particle.properVelocity),beforeX(particle.position);
  const Double g0=BoostFrameTransform::gammaFromProperVelocity(beforeU);
  Vec knownMidpoint(beforeX);knownMidpoint.pmv(.5*h*SI::c/g0,beforeU);
  Vec externalE(0.),externalB(0.);
  // Explicit second-order geometry for a prescribed analytic device only.
  // No future self field is predicted; the self field is solved by the outer
  // current/particle fixed point and its residual is checked independently.
  sources.addPrescribedBox(knownMidpoint,time+(s+.5)*h,frame,externalE,externalB);
  Vec totalE(gridE);totalE+=externalE;Vec totalB(gridB);totalB+=externalB;
  work.rootIterations=std::max(work.rootIterations,momentum(particle.properVelocity,totalE,totalB,particle.charge,particle.mass,h,g0));
  const Double g1=BoostFrameTransform::gammaFromProperVelocity(particle.properVelocity);
  Vec sum(beforeU);sum+=particle.properVelocity;particle.position.pmv(h*SI::c/(g0+g1),sum);
  // An absolute leaf-work bound is needed for stopping, not its extended
  // precision decomposition. The actual net held-field work telescopes.
  if(gridWork){Double subwork=0.;for(unsigned a=0;a<3;++a)subwork+=particle.charge*gridE[a]*(particle.position[a]-beforeX[a]);absoluteGrid+=std::abs(subwork);}
  if(prescribedWork)work.prescribed+=chargeWork(particle.charge,externalE,beforeX,particle.position);
 }
 if(gridWork)work.grid=chargeWork(particle.charge,gridE,firstPosition,particle.position);
 work.absoluteGrid=absoluteGrid;
 return work;
}
}
#endif
