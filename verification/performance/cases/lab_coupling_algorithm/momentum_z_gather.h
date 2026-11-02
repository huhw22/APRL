#ifndef APRL_TEST_MOMENTUM_Z_GATHER_H
#define APRL_TEST_MOMENTUM_Z_GATHER_H
#include "source_coupling_math.h"
namespace source_coupling_math {
// Test-only alternative: make the longitudinal magnetic impulse adjoint to
// the existing whole-step chord current and native Yee cPz functional.
// The ordinary three-vector Lorentz force is retained, so magnetic work is
// still zero. No force/residual subtraction or particle-momentum projection.
inline Vec momentumOrbitMagnetic(const EBFieldGrid& f,const Vec& origin,
 const Vec& begin,const Vec& end,const NativeHalo* halo) {
 const auto clipped=clipOpenOrbit(f.geometry(),origin,begin,end,halo);
 // Preserve the established open-face/handoff treatment in this first test.
 if(!clipped.valid||clipped.first!=0.||clipped.last!=1.)
   return openOrbitMagnetic(f,origin,begin,end,halo);
 const auto& g=f.geometry();const Double ds[3]={g.dx,g.dy,g.dz};
 Double x[3],d[3],cuts[32]={0.,1.};unsigned count=2;
 for(unsigned a=0;a<3;++a){x[a]=(begin[a]-origin[a])/ds[a];d[a]=(end[a]-begin[a])/ds[a];
   if(d[a]!=0.){const Double lo=std::min(x[a],x[a]+d[a]),hi=std::max(x[a],x[a]+d[a]);
     for(long long k=(long long)std::floor(lo)+1;k<hi;++k){
       if(count==32)throw std::runtime_error("Magnetic adjoint chord crossing budget");
       const Double t=(k-x[a])/d[a];if(t>0.&&t<1.)cuts[count++]=t;
     }
   }
 }
 std::sort(cuts,cuts+count);Vec result(0.);
 for(unsigned n=1;n<count;++n){const Double width=cuts[n]-cuts[n-1];if(width<=0.)continue;
   const Double middle=.5*(cuts[n]+cuts[n-1]);long long cell[3];Double w[3][2];
   for(unsigned a=0;a<3;++a){cell[a]=(long long)std::floor(x[a]+middle*d[a]);w[a][1]=x[a]+middle*d[a]-cell[a];w[a][0]=1.-w[a][1];}
   for(unsigned axis=0;axis<2;++axis){const unsigned b=(axis+1)%3,c=(axis+2)%3,magnetic=1-axis;const auto& v=magnetic==0?f.bx():f.by();
     for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j){long long p[3]={cell[0],cell[1],cell[2]};p[b]+=i;p[c]+=j;
       auto sample=[&](long long k){return halo?halo->openMagneticValue(v,magnetic,p[0],p[1],k):v(p[0],p[1],k);};
       const Double value=.5*(sample(p[2]-1)+sample(p[2]));
       const Double correction=(i==j?1.:-1.)*d[b]*d[c]*width*width/12.;
       result[magnetic]+=width*(w[b][i]*w[c][j]+correction)*value;
     }
   }
 }
 // Keep the old Bz interpolation exactly; only Bx/By enter longitudinal force.
 Vec middle(begin);middle+=end;middle*=.5;long long cell[3];Double w[3][2];
 for(unsigned a=0;a<3;++a){const Double v=(middle[a]-origin[a])/ds[a]-(a<2?.5:0.);cell[a]=(long long)std::floor(v);w[a][1]=v-cell[a];w[a][0]=1.-w[a][1];}
 for(unsigned k=0;k<2;++k)for(unsigned j=0;j<2;++j)for(unsigned i=0;i<2;++i){
   const Double v=halo?halo->openMagneticValue(f.bz(),2,cell[0]+i,cell[1]+j,cell[2]+k):f.bz()(cell[0]+i,cell[1]+j,cell[2]+k);
   result[2]+=w[0][i]*w[1][j]*w[2][k]*v;
 }
 return result;
}
}
#endif
