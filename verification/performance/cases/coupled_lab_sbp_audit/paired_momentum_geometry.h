#ifndef APRL_PAIRED_MOMENTUM_GEOMETRY_H
#define APRL_PAIRED_MOMENTUM_GEOMETRY_H
// Independent summation-by-parts momentum geometry. All products are
// colocated on native Yee sites BEFORE transverse trapezoid averaging.
// E/B getters include one longitudinal halo; global B ghosts may be zero,
// provided the unapplied end-node curl is reported as a held-node operation.
#include "eb_field.h"
#include <array>
namespace aprl {
struct PairedMomentumGeometry {
  EBGridGeometry g;std::size_t cx,cy,offset,globalNz;
  PairedMomentumGeometry(const EBGridGeometry& grid,std::size_t x,std::size_t y,std::size_t off,std::size_t nz):g(grid),cx(x),cy(y),offset(off),globalNz(nz){}
  static long double weight(std::size_t i,std::size_t lo,std::size_t hi){return i==lo||i==hi?.5L:1.L;}
  template<class E,class B> long double momentum(E e,B b)const {
    long double p=0.;const long double dv=(long double)g.dx*g.dy*g.dz;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for reduction(+:p) schedule(static)
#endif
    for(std::ptrdiff_t k=0;k<(std::ptrdiff_t)g.nz;++k){
      for(std::size_t j=cy;j<=g.ny-cy;++j)for(std::size_t i=cx;i<g.nx-cx;++i)p+=weight(j,cy,g.ny-cy)*.5L*(e(0,i,j,k)+e(0,i,j,k+1))*b(1,i,j,k);
      for(std::size_t j=cy;j<g.ny-cy;++j)for(std::size_t i=cx;i<=g.nx-cx;++i)p-=weight(i,cx,g.nx-cx)*.5L*(e(1,i,j,k)+e(1,i,j,k+1))*b(0,i,j,k);
    }return p*dv*SI::epsilon0*SI::c;
  }
  template<class E> long double gaussElectricImpulse(E e)const {
    auto rho=[&](std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return SI::epsilon0*((e(0,i,j,k)-e(0,i-1,j,k))/g.dx+(e(1,i,j,k)-e(1,i,j-1,k))/g.dy+(e(2,i,j,k)-e(2,i,j,k-1))/g.dz);};
    long double f=0.;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for reduction(+:f) schedule(static)
#endif
    for(std::ptrdiff_t k=0;k<(std::ptrdiff_t)g.nz;++k)for(std::size_t j=cy;j<=g.ny-cy;++j)for(std::size_t i=cx;i<=g.nx-cx;++i)f+=weight(i,cx,g.nx-cx)*weight(j,cy,g.ny-cy)*e(2,i,j,k)*.5L*(rho(i,j,k)+rho(i,j,k+1));
    return f*SI::c*g.dt*(long double)g.dx*g.dy*g.dz;
  }
  template<class B> long double gaussMagneticImpulse(B b)const {
    long double f=0.;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for reduction(+:f) schedule(static)
#endif
    for(std::ptrdiff_t k=0;k<(std::ptrdiff_t)g.nz;++k)for(std::size_t j=cy;j<g.ny-cy;++j)for(std::size_t i=cx;i<g.nx-cx;++i){
      const long double div=(b(0,i+1,j,k)-b(0,i,j,k))/g.dx+(b(1,i,j+1,k)-b(1,i,j,k))/g.dy+(b(2,i,j,k+1)-b(2,i,j,k))/g.dz;
      f+=.5L*(b(2,i,j,k)+b(2,i,j,k+1))*div;
    }return f*SI::c/SI::mu0*g.dt*(long double)g.dx*g.dy*g.dz;
  }
  template<class E> long double gaussInteriorElectricImpulse(E e)const {
    auto rho=[&](std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return SI::epsilon0*((e(0,i,j,k)-e(0,i-1,j,k))/g.dx+(e(1,i,j,k)-e(1,i,j-1,k))/g.dy+(e(2,i,j,k)-e(2,i,j,k-1))/g.dz);};
    long double f=0.;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for reduction(+:f) schedule(static)
#endif
    for(std::ptrdiff_t k=0;k<=(std::ptrdiff_t)g.nz;++k){if((k==0&&offset==0)||(k==(std::ptrdiff_t)g.nz&&offset+g.nz==globalNz))continue;const long double wz=k==0||k==(std::ptrdiff_t)g.nz?.5L:1.L;for(std::size_t j=cy;j<=g.ny-cy;++j)for(std::size_t i=cx;i<=g.nx-cx;++i)f+=wz*weight(i,cx,g.nx-cx)*weight(j,cy,g.ny-cy)*rho(i,j,k)*.5L*(e(2,i,j,k-1)+e(2,i,j,k));}
    return f*SI::c*g.dt*(long double)g.dx*g.dy*g.dz;
  }
  template<class E> std::array<long double,6> electricPortStress(E e)const {
    // This is evaluated directly from INNER normal E and face tangential E.
    // No exterior Ez ghost or inferred residual enters this port exchange.
    auto out=electricStress(e);
    for(unsigned side=0;side<2;++side)if((side==0&&offset==0)||(side==1&&offset+g.nz==globalNz)){
      const std::ptrdiff_t k=side?g.nz:0,inside=side?g.nz-1:0;const long double sign=side?1:-1;long double tangent=0.,normal=0.,surface=0.;
      for(std::size_t j=cy;j<=g.ny-cy;++j)for(std::size_t i=cx;i<g.nx-cx;++i)tangent+=weight(j,cy,g.ny-cy)*e(0,i,j,k)*e(0,i,j,k);
      for(std::size_t j=cy;j<g.ny-cy;++j)for(std::size_t i=cx;i<=g.nx-cx;++i)tangent+=weight(i,cx,g.nx-cx)*e(1,i,j,k)*e(1,i,j,k);
      for(std::size_t j=cy;j<=g.ny-cy;++j)for(std::size_t i=cx;i<=g.nx-cx;++i){const long double w=weight(i,cx,g.nx-cx)*weight(j,cy,g.ny-cy),ez=e(2,i,j,inside);normal+=w*ez*ez;surface+=w*ez*((e(0,i,j,k)-e(0,i-1,j,k))/g.dx+(e(1,i,j,k)-e(1,i,j-1,k))/g.dy);}
      out[4+side]=SI::epsilon0*SI::c*g.dt*.5L*g.dx*g.dy*(sign*(tangent-normal)+g.dz*surface);
    }return out;
  }
  template<class E> std::array<long double,6> electricStress(E e)const {
    std::array<long double,6> out{};
    for(unsigned side=0;side<2;++side){const auto i=side?g.nx-cx:cx;const long double sign=side?1:-1;
      for(std::ptrdiff_t k=0;k<(std::ptrdiff_t)g.nz;++k)for(std::size_t j=cy;j<=g.ny-cy;++j)out[side]-=sign*weight(j,cy,g.ny-cy)*e(2,i,j,k)*.25L*(e(0,i-1,j,k)+e(0,i,j,k)+e(0,i-1,j,k+1)+e(0,i,j,k+1))*g.dy*g.dz;
      const auto j=side?g.ny-cy:cy;
      for(std::ptrdiff_t k=0;k<(std::ptrdiff_t)g.nz;++k)for(std::size_t x=cx;x<=g.nx-cx;++x)out[2+side]-=sign*weight(x,cx,g.nx-cx)*e(2,x,j,k)*.25L*(e(1,x,j-1,k)+e(1,x,j,k)+e(1,x,j-1,k+1)+e(1,x,j,k+1))*g.dx*g.dz;
      if((side==0&&offset==0)||(side==1&&offset+g.nz==globalNz)){const std::ptrdiff_t k=side?g.nz:0;
        long double t=0.;for(std::size_t y=cy;y<=g.ny-cy;++y)for(std::size_t x=cx;x<g.nx-cx;++x)t+=weight(y,cy,g.ny-cy)*e(0,x,y,k)*e(0,x,y,k);
        for(std::size_t y=cy;y<g.ny-cy;++y)for(std::size_t x=cx;x<=g.nx-cx;++x)t+=weight(x,cx,g.nx-cx)*e(1,x,y,k)*e(1,x,y,k);
        for(std::size_t y=cy;y<=g.ny-cy;++y)for(std::size_t x=cx;x<=g.nx-cx;++x)t-=weight(x,cx,g.nx-cx)*weight(y,cy,g.ny-cy)*e(2,x,y,k-1)*e(2,x,y,k);
        out[4+side]+=.5L*sign*t*g.dx*g.dy;
      }
    }for(auto& v:out)v*=SI::epsilon0*SI::c*g.dt;return out;
  }
  template<class B> std::array<long double,6> magneticStress(B b)const {
    std::array<long double,6> out{};
    for(unsigned side=0;side<2;++side){const auto i=side?g.nx-cx:cx;const long double sign=side?1:-1;
      for(std::ptrdiff_t k=0;k<(std::ptrdiff_t)g.nz;++k)for(std::size_t j=cy;j<g.ny-cy;++j)out[side]-=sign*b(0,i,j,k)*.25L*(b(2,i-1,j,k)+b(2,i,j,k)+b(2,i-1,j,k+1)+b(2,i,j,k+1))*g.dy*g.dz;
      const auto j=side?g.ny-cy:cy;
      for(std::ptrdiff_t k=0;k<(std::ptrdiff_t)g.nz;++k)for(std::size_t x=cx;x<g.nx-cx;++x)out[2+side]-=sign*b(1,x,j,k)*.25L*(b(2,x,j-1,k)+b(2,x,j,k)+b(2,x,j-1,k+1)+b(2,x,j,k+1))*g.dx*g.dz;
      if((side==0&&offset==0)||(side==1&&offset+g.nz==globalNz)){const std::ptrdiff_t k=side?g.nz:0;
        long double t=0.;for(std::size_t y=cy;y<g.ny-cy;++y)for(std::size_t x=cx;x<=g.nx-cx;++x)t+=weight(x,cx,g.nx-cx)*b(0,x,y,k-1)*b(0,x,y,k);
        for(std::size_t y=cy;y<=g.ny-cy;++y)for(std::size_t x=cx;x<g.nx-cx;++x)t+=weight(y,cy,g.ny-cy)*b(1,x,y,k-1)*b(1,x,y,k);
        for(std::size_t y=cy;y<g.ny-cy;++y)for(std::size_t x=cx;x<g.nx-cx;++x)t-=b(2,x,y,k)*b(2,x,y,k);
        out[4+side]+=.5L*sign*t*g.dx*g.dy;
      }
    }for(auto& v:out)v*=SI::c/SI::mu0*g.dt;return out;
  }
};
}
#endif
