#ifndef APRL_EXPERIMENTAL_SOURCE_COUPLING_H
#define APRL_EXPERIMENTAL_SOURCE_COUPLING_H
// Independent prototype only. Nothing in the production Simulation includes it.
#include "eb_particles.h"
#include "eb_deposition.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace source_trial {
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

// Exact analytic orbit average of the first-order depositor-adjoint electric
// interpolant. Inside each cell the transverse product is quadratic along a
// straight chord: midpoint product plus signed d_eta*d_zeta/12.
inline Vec orbitElectric(const EBFieldGrid& f,const Electric* old,
 const Vec& origin,const Vec& begin,const Vec& end) {
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
   if(cell[a]<0||cell[a]>=(long long)(a==0?g.nx:(a==1?g.ny:g.nz)))throw std::out_of_range("Prototype orbit outside grid");
   w[a][1]=x[a]+middle*d[a]-cell[a];w[a][0]=1.-w[a][1];
  }
  for(unsigned a=0;a<3;++a) {
   const unsigned b=(a+1)%3,c=(a+2)%3;
   for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j) {
    std::size_t node[3]={(std::size_t)cell[0],(std::size_t)cell[1],(std::size_t)cell[2]};node[b]+=i;node[c]+=j;
    Double value=(*e[a])(node[0],node[1],node[2]);
    if(old)value=.5*(value+(*old)[a](node[0],node[1],node[2]));
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
inline unsigned momentum(Vec& u,const Vec& e,const Vec& magnetic,Double q,Double m,Double h) {
 const Vec before(u);const Double g0=BoostFrameTransform::gammaFromProperVelocity(before);
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
 const SIFieldSourceSet& sources,const BoostFrameTransform& frame,Double time,Double dt,unsigned substeps) {
 PushWork work;const Double h=dt/substeps;
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
  work.rootIterations=std::max(work.rootIterations,momentum(particle.properVelocity,totalE,totalB,particle.charge,particle.mass,h));
  const Double g1=BoostFrameTransform::gammaFromProperVelocity(particle.properVelocity);
  Vec sum(beforeU);sum+=particle.properVelocity;particle.position.pmv(h*SI::c/(g0+g1),sum);
  const auto grid=chargeWork(particle.charge,gridE,beforeX,particle.position);
  work.grid+=grid;work.absoluteGrid+=std::abs(grid);
  work.prescribed+=chargeWork(particle.charge,externalE,beforeX,particle.position);
 }
 return work;
}
}
#endif
