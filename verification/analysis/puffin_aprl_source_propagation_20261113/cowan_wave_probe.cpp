// Exact staggered plane wave against the actual production Cowan/Yee bodies.
#include "eb_field.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

int main(int argc,char** argv) {
  if(argc!=4)return 2;
  const double ky=std::strtod(argv[1],nullptr),q=std::strtod(argv[2],nullptr),phi=std::strtod(argv[3],nullptr);
  const double dx=4e-5,dy=4e-5,dz=8.293564128124092e-6,dt=dz/aprl::SI::c;
  aprl::EBGridGeometry geometry(16,16,20,dx,dy,dz,dt,aprl::EBMaxwellSolver::CowanZ);
  aprl::EBFieldGrid f(geometry);const auto a=f.cowanCoefficients();
  const double cy=std::cos(ky*dy),cz=std::cos(q),
    Cy=a.alpha[1]+2*a.beta[2]*cz+2*a.beta[0]+4*a.deltaZX*cz,
    Cz=a.alpha[2]+2*a.beta[0]+2*a.beta[1]*cy+4*a.deltaXY*cy,
    S=2*std::sin(phi/2)/dt,By=2*std::sin(q/2)/dz*Cz/S,
    Bz=-2*std::sin(ky*dy/2)/dy*Cy/S;
  for(std::size_t k=0;k<f.ex().nz();++k)
    for(std::size_t j=0;j<f.ex().ny();++j)
      for(std::size_t i=0;i<f.ex().nx();++i)f.ex()(i,j,k)=std::cos(ky*j*dy+q*k);
  for(std::size_t k=0;k<f.by().nz();++k)
    for(std::size_t j=0;j<f.by().ny();++j)
      for(std::size_t i=0;i<f.by().nx();++i)f.by()(i,j,k)=By*std::cos(ky*j*dy+q*(k+.5)+phi/2);
  for(std::size_t k=0;k<f.bz().nz();++k)
    for(std::size_t j=0;j<f.bz().ny();++j)
      for(std::size_t i=0;i<f.bz().nx();++i)f.bz()(i,j,k)=Bz*std::cos(ky*(j+.5)*dy+q*k+phi/2);
  f.advanceMagnetic();f.advanceElectric();double e=0,b=0;
  for(std::size_t k=4;k<15;++k)for(std::size_t j=4;j<11;++j)for(std::size_t i=4;i<11;++i) {
    e=std::max(e,std::abs(f.ex()(i,j,k)-std::cos(ky*j*dy+q*k-phi)));
    b=std::max(b,aprl::SI::c*std::abs(f.by()(i,j,k)-By*std::cos(ky*j*dy+q*(k+.5)-phi/2)));
    b=std::max(b,aprl::SI::c*std::abs(f.bz()(i,j,k)-Bz*std::cos(ky*(j+.5)*dy+q*k-phi/2)));
  }
  std::cout<<std::setprecision(17)<<ky<<','<<q<<','<<phi<<','<<e<<','<<b<<'\n';
  return e<2e-12&&b<2e-12?0:1;
}
