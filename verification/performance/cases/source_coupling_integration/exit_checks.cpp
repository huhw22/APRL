#include "source_coupling_math.h"
#include <iostream>
#include <iomanip>
using namespace source_coupling_math;
int main() {
 try {
  const double dz=1.2543262415359034e-6;
  EBGridGeometry g(24,24,24,5e-6,5e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
  EBFieldGrid f(g),old(g);Vec origin(0.);origin[2]=-.0031;
  YeeComponent* ef[3]={&f.ex(),&f.ey(),&f.ez()};
  YeeComponent* eo[3]={&old.ex(),&old.ey(),&old.ez()};
  YeeComponent* bf[3]={&f.bx(),&f.by(),&f.bz()};
  for(unsigned a=0;a<3;++a) {
   for(std::size_t k=0;k<ef[a]->nz();++k)for(std::size_t j=0;j<ef[a]->ny();++j)for(std::size_t i=0;i<ef[a]->nx();++i){(*ef[a])(i,j,k)=1e5*(.2*i+.03*j*k+a+1);(*eo[a])(i,j,k)=2e5*(.1*i-.02*j*k+a+2);}
   for(std::size_t k=0;k<bf[a]->nz();++k)for(std::size_t j=0;j<bf[a]->ny();++j)for(std::size_t i=0;i<bf[a]->nx();++i)(*bf[a])(i,j,k)=1e-3*(.2*i+.03*j*k+a+1);
  }
  NativeHalo halo;const auto saved=copyElectric(old);const double ds[3]={g.dx,g.dy,g.dz};
  long double worst=0.;unsigned samples=0,strictRejected=0;
  for(unsigned axis=0;axis<3;++axis)for(unsigned side=0;side<2;++side)for(unsigned n=0;n<16;++n) {
   Vec start(origin);for(unsigned a=0;a<3;++a)start[a]+=(10.1+.007*n)*ds[a];
   start[axis]=origin[axis]+(side?23.8:0.2)*ds[axis];Vec end(start);
   end[axis]+=(side?1.:-1.)*.4*ds[axis];end[(axis+1)%3]+=.13*ds[(axis+1)%3];
   try{orbitElectric(f,&saved,origin,start,end,&halo,&halo);}catch(const std::out_of_range&){++strictRejected;}
   const auto clipped=clipOpenOrbit(g,origin,start,end,&halo);if(!clipped.valid)throw std::runtime_error("Lost valid physical prefix");
   const auto avg=openOrbitElectric(f,&saved,origin,start,end,&halo,&halo);
   const auto b=openOrbitMagnetic(f,origin,start,end,&halo);Vec middle(clipped.begin);middle+=clipped.end;middle*=.5;Vec unused(0.),reference(0.);f.sampleFieldsPosition(middle,origin,unused,reference);
   for(unsigned a=0;a<3;++a)if(std::abs(b[a]-reference[a])>1e-13)throw std::runtime_error("Global-face B continuation differs from original sampler");
   f.clearCurrent();ChargeConservingCurrentDepositor(f,origin).depositSegment(clipped.begin,clipped.end,-SI::elementaryCharge);
   const YeeComponent* j[3]={&f.jx(),&f.jy(),&f.jz()};long double mesh=0.;
   for(unsigned a=0;a<3;++a)for(std::size_t i=0;i<ef[a]->size();++i)mesh+=(long double)j[a]->data()[i]*(ef[a]->data()[i]+(long double)eo[a]->data()[i])/2;
   mesh*=(long double)g.dt*g.dx*g.dy*g.dz;
   const auto particle=chargeWork(-SI::elementaryCharge,avg,start,end);
   worst=std::max(worst,std::abs(mesh-particle)/std::max(std::abs(mesh),1e-22L));++samples;
  }
  if(strictRejected!=96||worst>2e-10)throw std::runtime_error("Open-face orbit/deposition check failed");
  EBGridGeometry longSlab(24,24,83,g.dx,g.dy,dz,g.dt,EBMaxwellSolver::CowanZ);
  EBFieldGrid edge(longSlab);YeeComponent* uniform[3]={&edge.ex(),&edge.ey(),&edge.ez()};
  for(unsigned a=0;a<3;++a)std::fill(uniform[a]->data(),uniform[a]->data()+uniform[a]->size(),1e5);
  const auto oldEdge=copyElectric(edge);unsigned translated=0;long double translatedError=0.;
  for(unsigned n=0;n<128;++n)for(double fraction:{.07,.18,.29,.499,.73}) {
   Vec shifted(0.);shifted[2]=.0038723942-n*dz;Vec start(shifted);start[0]+=10.3*g.dx;start[1]+=10.4*g.dy;start[2]+=(83-fraction)*dz;Vec end(start);end[2]+=.99*dz;
   const auto avg=openOrbitElectric(edge,&oldEdge,shifted,start,end,&halo,&halo);
   const auto actual=chargeWork(-SI::elementaryCharge,avg,start,end);
   const long double expected=-SI::elementaryCharge*1e5L*((shifted[2]+83*dz)-start[2]);
   translatedError=std::max(translatedError,std::abs(actual-expected)/std::max(std::abs(expected),1e-25L));++translated;
  }
  if(translatedError>2e-10)throw std::runtime_error("Translated global face loses its in-domain orbit");
  std::cout<<std::setprecision(21)<<"{\"passed\":true,\"six_face_crossings\":"<<samples<<",\"old_strict_failure_reproduced\":"<<strictRejected<<",\"max_scaled_physical_prefix_work_error\":"<<worst<<",\"translated_83_cell_face_samples\":"<<translated<<",\"max_scaled_translated_face_work_error\":"<<translatedError<<"}\n";
  return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
