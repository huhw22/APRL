#ifndef APRL_TEST_COWAN_MATCHED_GATHER_H
#define APRL_TEST_COWAN_MATCHED_GATHER_H
#include "source_coupling_math.h"
#include <map>
namespace source_coupling_math {
struct CowanMetric {
  struct Term {int x,y,z;Double value;};
  using Key=std::array<int,3>;
  using Polynomial=std::map<Key,Double>;
  std::array<std::vector<Term>,3> electric,magnetic;
  Double dx,dy,dz;
  static Polynomial smooth(const EBCowanCoefficients& c,unsigned axis){
    Polynomial p;p[{{0,0,0}}]=c.alpha[axis];const unsigned j=(axis+1)%3,k=(axis+2)%3;const Double delta=axis==0?c.deltaYZ:(axis==1?c.deltaZX:c.deltaXY);
    for(int sign:{-1,1}){Key a{{0,0,0}},b{{0,0,0}};a[j]=sign;b[k]=sign;p[a]+=c.beta[j];p[b]+=c.beta[k];}
    for(int u:{-1,1})for(int v:{-1,1}){Key a{{0,0,0}};a[j]=u;a[k]=v;p[a]+=delta;}return p;
  }
  static std::vector<Term> terms(const Polynomial& p){std::vector<Term> v;for(const auto& q:p)if(q.second!=0.)v.push_back({q.first[0],q.first[1],q.first[2],q.second});return v;}
  explicit CowanMetric(const EBGridGeometry& g):dx(g.dx),dy(g.dy),dz(g.dz){
    const auto c=g.solver==EBMaxwellSolver::CowanZ?EBCowanCoefficients::forZDispersion(dx,dy,dz):EBCowanCoefficients();Polynomial s[3];for(unsigned a=0;a<3;++a){s[a]=smooth(c,a);magnetic[a]=terms(s[a]);}
    for(unsigned a=0;a<3;++a){Polynomial p;for(const auto& u:s[(a+1)%3])for(const auto& v:s[(a+2)%3]){Key q{{u.first[0]+v.first[0],u.first[1]+v.first[1],u.first[2]+v.first[2]}};p[q]+=u.second*v.second;}electric[a]=terms(p);}
  }
  template<class Getter>Double apply(const std::vector<Term>& v,Getter f,unsigned a,long long i,long long j,long long k)const{Double s=0.;for(const auto& q:v)s+=q.value*f(a,i+q.x,j+q.y,k+q.z);return s;}
};
inline const CowanMetric& matchedMetric(const EBGridGeometry& g){
  static thread_local std::unique_ptr<CowanMetric> metric;
  if(!metric||metric->dx!=g.dx||metric->dy!=g.dy||metric->dz!=g.dz)metric.reset(new CowanMetric(g));return *metric;
}
// Same real-global electric continuation as EBFieldGrid::electricValue.
// MPI extensions always read the true donor, never clamp at a rank interface.
inline Double matchedElectricValue(const YeeComponent& v,unsigned a,long long i,long long j,long long k,const NativeHalo* halo){
  long long p[3]={i,j,k},n[3]={(long long)v.nx(),(long long)v.ny(),(long long)v.nz()};Double sign=1.;
  if(i>=0&&j>=0&&k>=0&&i<n[0]&&j<n[1]&&k<n[2])return v(i,j,k);
  for(unsigned d=0;d<3;++d)if(p[d]<0||p[d]>=n[d]){
    if(d==2&&halo&&((p[d]<0&&halo->hasLower())||(p[d]>=n[d]&&halo->hasUpper())))continue;
    if(a==d)p[d]=std::max(0LL,std::min(p[d],n[d]-1));
    else{p[d]=p[d]<0?-p[d]:2*(n[d]-1)-p[d];sign=-sign;}
  }
  return sign*(halo?halo->value(v,a,p[0],p[1],p[2]):v(p[0],p[1],p[2]));
}
template<class Getter>inline Vec matchedOrbit(const EBFieldGrid& f,const Vec& origin,const Vec& begin,const Vec& end,Getter get,unsigned axes=3){
  const auto& g=f.geometry();const Double ds[3]={g.dx,g.dy,g.dz};Double x[3],d[3],cuts[32]={0.,1.};unsigned count=2;
  for(unsigned a=0;a<3;++a){x[a]=(begin[a]-origin[a])/ds[a];d[a]=(end[a]-begin[a])/ds[a];if(d[a]!=0.){const Double lo=std::min(x[a],x[a]+d[a]),hi=std::max(x[a],x[a]+d[a]);for(long long k=(long long)std::floor(lo)+1;k<hi;++k){if(count==32)throw std::runtime_error("Matched orbit crossing budget");const Double t=(k-x[a])/d[a];if(t>0.&&t<1.)cuts[count++]=t;}}}
  std::sort(cuts,cuts+count);Vec out(0.);
  for(unsigned n=1;n<count;++n){const Double width=cuts[n]-cuts[n-1];if(width<=0.)continue;const Double mid=.5*(cuts[n]+cuts[n-1]);long long cell[3];Double w[3][2];for(unsigned a=0;a<3;++a){cell[a]=(long long)std::floor(x[a]+mid*d[a]);w[a][1]=x[a]+mid*d[a]-cell[a];w[a][0]=1.-w[a][1];}
    for(unsigned a=0;a<axes;++a){const unsigned b=(a+1)%3,c=(a+2)%3;for(unsigned i=0;i<2;++i)for(unsigned j=0;j<2;++j){long long p[3]={cell[0],cell[1],cell[2]};p[b]+=i;p[c]+=j;const Double correction=(i==j?1.:-1.)*d[b]*d[c]*width*width/12.;out[a]+=width*(w[b][i]*w[c][j]+correction)*get(a,p[0],p[1],p[2]);}}
  }return out;
}
inline Vec matchedOrbitElectric(const EBFieldGrid& f,const Electric* old,const Vec& origin,const Vec& begin,const Vec& end,const NativeHalo* halo=nullptr,const NativeHalo* oldHalo=nullptr){
  const auto& metric=matchedMetric(f.geometry());const YeeComponent* fields[3]={&f.ex(),&f.ey(),&f.ez()};
  auto raw=[&](unsigned a,long long i,long long j,long long k){const Double v=matchedElectricValue(*fields[a],a,i,j,k,halo);return old?.5*(v+matchedElectricValue((*old)[a],a,i,j,k,oldHalo)):v;};
  auto value=[&](unsigned a,long long i,long long j,long long k){return metric.apply(metric.electric[a],raw,a,i,j,k);};return matchedOrbit(f,origin,begin,end,value);
}
inline Vec matchedOpenOrbitElectric(const EBFieldGrid& f,const Electric* old,const Vec& origin,const Vec& begin,const Vec& end,const NativeHalo* halo,const NativeHalo* oldHalo){
  const auto r=clipOpenOrbit(f.geometry(),origin,begin,end,halo);if(!r.valid)throw std::out_of_range("No in-domain matched orbit");auto v=matchedOrbitElectric(f,old,origin,r.begin,r.end,halo,oldHalo);v*=r.last-r.first;return v;
}
inline Vec matchedOrbitMagnetic(const EBFieldGrid& f,const Vec& origin,const Vec& begin,const Vec& end,const NativeHalo* halo){
  const auto r=clipOpenOrbit(f.geometry(),origin,begin,end,halo);if(!r.valid)throw std::out_of_range("No matched magnetic orbit");const auto& metric=matchedMetric(f.geometry());const YeeComponent* fields[3]={&f.bx(),&f.by(),&f.bz()};
  auto raw=[&](unsigned a,long long i,long long j,long long k){return halo?halo->openMagneticValue(*fields[a],a,i,j,k):(*fields[a])(i,j,k);};
  auto value=[&](unsigned a,long long i,long long j,long long k){const unsigned b=1-a;auto sample=[&](long long z){return metric.apply(metric.magnetic[b],raw,b,i,j,z);};return .5*(sample(k-1)+sample(k));};
  const auto transverse=matchedOrbit(f,origin,r.begin,r.end,value,2);Vec out(0.);out[0]=transverse[1]*(r.last-r.first);out[1]=transverse[0]*(r.last-r.first);
  Vec middle(r.begin);middle+=r.end;middle*=.5;out[2]=magneticSample(f,origin,middle,false,halo,true).value[2];return out;
}
}
#endif
