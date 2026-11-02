#ifndef APRL_VERIFICATION_YEE_SBP_GEOMETRY_H
#define APRL_VERIFICATION_YEE_SBP_GEOMETRY_H
// Diagnostic geometry only. No field updates. All ports are evaluated from
// native values, without using measured energy/momentum residuals.
#include "paired_momentum_geometry.h"
namespace aprl {
struct YeeSBPGeometry:PairedMomentumGeometry {
  YeeSBPGeometry(const EBGridGeometry& grid,std::size_t x,std::size_t y,std::size_t off,std::size_t nz):PairedMomentumGeometry(grid,x,y,off,nz){}
  template<class E,class B> long double electricInner(E e,B v)const {
    long double out=0.;const long double dv=(long double)g.dx*g.dy*g.dz;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for reduction(+:out) schedule(static)
#endif
    for(std::ptrdiff_t k=0;k<(std::ptrdiff_t)g.nz;++k)for(std::size_t j=cy;j<g.ny-cy;++j)for(std::size_t i=cx;i<g.nx-cx;++i)for(unsigned a=0;a<3;++a)for(unsigned p=0;p<2;++p)for(unsigned q=0;q<2;++q){const auto x=i+(a==0?0:p),y=j+(a==0?p:(a==1?0:q));const auto z=k+(a==2?0:q);out+=e(a,x,y,z)*v(a,x,y,z);}
    return out*dv/4;
  }
  template<class B,class C> long double magneticInner(B b,C c)const {
    long double out=0.;const long double dv=(long double)g.dx*g.dy*g.dz;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for reduction(+:out) schedule(static)
#endif
    for(std::ptrdiff_t k=0;k<(std::ptrdiff_t)g.nz;++k)for(std::size_t j=cy;j<g.ny-cy;++j)for(std::size_t i=cx;i<g.nx-cx;++i)for(unsigned a=0;a<3;++a)for(unsigned p=0;p<2;++p)out+=b(a,i+(a==0?p:0),j+(a==1?p:0),k+(a==2?p:0))*c(a,i+(a==0?p:0),j+(a==1?p:0),k+(a==2?p:0));
    return out*dv/2;
  }
  template<class E> long double faraday(E e,unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k)const {
    if(a==0)return (long double)g.dt*((e(1,i,j,k+1)-e(1,i,j,k))/g.dz-(e(2,i,j+1,k)-e(2,i,j,k))/g.dy);
    if(a==1)return (long double)g.dt*((e(2,i+1,j,k)-e(2,i,j,k))/g.dx-(e(0,i,j,k+1)-e(0,i,j,k))/g.dz);
    return (long double)g.dt*((e(0,i,j+1,k)-e(0,i,j,k))/g.dy-(e(1,i+1,j,k)-e(1,i,j,k))/g.dx);
  }
  template<class B> long double ampere(B b,unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k)const {
    const long double s=(long double)g.dt/(SI::mu0*SI::epsilon0);
    if(a==0)return s*((b(2,i,j,k)-b(2,i,j-1,k))/g.dy-(b(1,i,j,k)-b(1,i,j,k-1))/g.dz);
    if(a==1)return s*((b(0,i,j,k)-b(0,i,j,k-1))/g.dz-(b(2,i,j,k)-b(2,i-1,j,k))/g.dx);
    return s*((b(1,i,j,k)-b(1,i-1,j,k))/g.dx-(b(0,i,j,k)-b(0,i,j-1,k))/g.dy);
  }
  template<class E,class B> std::array<long double,6> energyPort(E e,B b)const {
    std::array<long double,6> q{};
    for(unsigned side=0;side<2;++side){const auto i=side?g.nx-cx:cx;const long double sign=side?1:-1;
      for(std::ptrdiff_t k=0;k<=(std::ptrdiff_t)g.nz;++k)for(std::size_t j=cy;j<g.ny-cy;++j)q[side]+=sign*weight(k,0,g.nz)*e(1,i,j,k)*.5L*(b(2,i-1,j,k)+b(2,i,j,k))*g.dy*g.dz;
      for(std::ptrdiff_t k=0;k<(std::ptrdiff_t)g.nz;++k)for(std::size_t j=cy;j<=g.ny-cy;++j)q[side]-=sign*weight(j,cy,g.ny-cy)*e(2,i,j,k)*.5L*(b(1,i-1,j,k)+b(1,i,j,k))*g.dy*g.dz;
      const auto j=side?g.ny-cy:cy;
      for(std::ptrdiff_t k=0;k<(std::ptrdiff_t)g.nz;++k)for(std::size_t x=cx;x<=g.nx-cx;++x)q[2+side]+=sign*weight(x,cx,g.nx-cx)*e(2,x,j,k)*.5L*(b(0,x,j-1,k)+b(0,x,j,k))*g.dx*g.dz;
      for(std::ptrdiff_t k=0;k<=(std::ptrdiff_t)g.nz;++k)for(std::size_t x=cx;x<g.nx-cx;++x)q[2+side]-=sign*weight(k,0,g.nz)*e(0,x,j,k)*.5L*(b(2,x,j-1,k)+b(2,x,j,k))*g.dx*g.dz;
      if((side==0&&offset==0)||(side==1&&offset+g.nz==globalNz)){const std::ptrdiff_t k=side?g.nz:0;
        for(std::size_t y=cy;y<=g.ny-cy;++y)for(std::size_t x=cx;x<g.nx-cx;++x)q[4+side]+=sign*weight(y,cy,g.ny-cy)*e(0,x,y,k)*.5L*(b(1,x,y,k-1)+b(1,x,y,k))*g.dx*g.dy;
        for(std::size_t y=cy;y<g.ny-cy;++y)for(std::size_t x=cx;x<=g.nx-cx;++x)q[4+side]-=sign*weight(x,cx,g.nx-cx)*e(1,x,y,k)*.5L*(b(0,x,y,k-1)+b(0,x,y,k))*g.dx*g.dy;
      }
    }for(auto& x:q)x*=g.dt/SI::mu0;return q;
  }
  template<class E,class B> long double curlEnergyWork(E e,B b)const {
    auto de=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return ampere(b,a,i,j,k);};
    auto db=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return faraday(e,a,i,j,k);};
    return SI::epsilon0*electricInner(e,de)+magneticInner(b,db)/SI::mu0;
  }
  static long double sum(const std::array<long double,6>& a){long double out=0.;for(auto x:a)out+=x;return out;}
};
}
#endif
