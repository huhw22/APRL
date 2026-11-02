#ifndef APRL_VERIFICATION_COWAN_SBP_GEOMETRY_H
#define APRL_VERIFICATION_COWAN_SBP_GEOMETRY_H
#include "yee_sbp_geometry.h"
namespace aprl {
// Independent evaluation of the existing Cowan stencil, not a field update.
struct CowanSBPGeometry:YeeSBPGeometry {
  EBCowanCoefficients coefficients;
  CowanSBPGeometry(const EBGridGeometry& grid,std::size_t x,std::size_t y,std::size_t off,std::size_t nz):YeeSBPGeometry(grid,x,y,off,nz),coefficients(grid.solver==EBMaxwellSolver::CowanZ?EBCowanCoefficients::forZDispersion(grid.dx,grid.dy,grid.dz):EBCowanCoefficients()){}
  template<class E> long double smooth(E e,unsigned a,unsigned derivative,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k)const {
    const unsigned first=(derivative+1)%3,second=(derivative+2)%3;
    const long double delta=derivative==0?coefficients.deltaYZ:(derivative==1?coefficients.deltaZX:coefficients.deltaXY);
    long double value=coefficients.alpha[derivative]*e(a,i,j,k);const std::ptrdiff_t base[3]={i,j,k};
    for(int sign=-1;sign<=1;sign+=2){std::ptrdiff_t p[3]={base[0],base[1],base[2]},q[3]={base[0],base[1],base[2]};p[first]+=sign;q[second]+=sign;value+=coefficients.beta[first]*e(a,p[0],p[1],p[2]);value+=coefficients.beta[second]*e(a,q[0],q[1],q[2]);}
    for(int p=-1;p<=1;p+=2)for(int q=-1;q<=1;q+=2){std::ptrdiff_t point[3]={base[0],base[1],base[2]};point[first]+=p;point[second]+=q;value+=delta*e(a,point[0],point[1],point[2]);}return value;
  }
  template<class E> long double cowan(E e,unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k)const {
    auto v=[&](unsigned c,unsigned d,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return smooth(e,c,d,x,y,z);};
    if(a==0)return (long double)g.dt*((v(1,2,i,j,k+1)-v(1,2,i,j,k))/g.dz-(v(2,1,i,j+1,k)-v(2,1,i,j,k))/g.dy);
    if(a==1)return (long double)g.dt*((v(2,0,i+1,j,k)-v(2,0,i,j,k))/g.dx-(v(0,2,i,j,k+1)-v(0,2,i,j,k))/g.dz);
    return (long double)g.dt*((v(0,1,i,j+1,k)-v(0,1,i,j,k))/g.dy-(v(1,0,i+1,j,k)-v(1,0,i,j,k))/g.dx);
  }
  template<class E> long double extra(E e,unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k)const{return cowan(e,a,i,j,k)-faraday(e,a,i,j,k);}
};
}
#endif
