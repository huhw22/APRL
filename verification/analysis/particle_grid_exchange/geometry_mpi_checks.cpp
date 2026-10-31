#include "particle_work_geometry.h"
#include "eb_deposition.h"
#include <iostream>
using namespace aprl;
using Vec=FieldVector<Double>;
void fill(EBFieldGrid& f,std::size_t offset) {
  YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};
  for(unsigned a=0;a<3;++a)for(std::size_t k=0;k<e[a]->nz();++k)for(std::size_t j=0;j<e[a]->ny();++j)for(std::size_t i=0;i<e[a]->nx();++i)
    (*e[a])(i,j,k)=1e5*(std::sin(.21*i+.17*j+.37*(k+offset)+a)+.03*i*j-.002*(k+offset)*j);
}
long double work(const EBFieldGrid& f) {
  const YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};const YeeComponent* js[3]={&f.jx(),&f.jy(),&f.jz()};
  long double s=0;for(unsigned a=0;a<3;++a)for(std::size_t k=0;k<e[a]->nz();++k)for(std::size_t j=0;j<e[a]->ny();++j)for(std::size_t i=0;i<e[a]->nx();++i)s+=(long double)(*e[a])(i,j,k)*(*js[a])(i,j,k);
  const auto& g=f.geometry();return s*g.dt*g.dx*g.dy*g.dz;
}
int main(int argc,char** argv) {
 MPI_Init(&argc,&argv);int rank=0,size=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
 try {
  const Double dz=1.2543262415359034e-6;
  EBGridGeometry localGeometry(16,16,16,5e-6,5e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
  EBGridGeometry globalGeometry(16,16,16*size,5e-6,5e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
  EBFieldGrid local(localGeometry),global(globalGeometry);fill(local,16*rank);fill(global,0);
  GeometryElectricHalo halo;halo.exchange(local,MPI_COMM_WORLD);
  Vec origin(0.),zero(0.);origin[2]=16*rank*dz;
  ChargeConservingCurrentDepositor deposit(global,zero);int checks=0;long double maxError=0;
  for(int direction:{-1,1}) {
    if((direction<0 && rank==0)||(direction>0 && rank+1==size))continue;
    Vec x(0.),m(0.),y(0.);x[0]=8.21*5e-6;x[1]=7.35*5e-6;m=x;y=x;
    x[2]=(16*rank+(direction<0?.3:15.8))*dz;
    m[2]=x[2]+direction*.35*dz;y[2]=x[2]+direction*.8*dz;
    m[0]+=.11*5e-6;m[1]-=.13*5e-6;y[0]-=.07*5e-6;y[1]+=.08*5e-6;
    for(bool curved:{false,true}) {
      global.clearCurrent();long double calculated=0;
      if(curved) {
        deposit.depositSegment(x,m,-SI::elementaryCharge);deposit.depositSegment(m,y,-SI::elementaryCharge);
        calculated=trajectoryElectricWork(local,origin,x,m,-SI::elementaryCharge,&halo)+trajectoryElectricWork(local,origin,m,y,-SI::elementaryCharge,&halo);
      }else {deposit.depositSegment(x,y,-SI::elementaryCharge);calculated=trajectoryElectricWork(local,origin,x,y,-SI::elementaryCharge,&halo);}
      const long double exact=work(global),error=std::abs(calculated-exact)/std::max(std::abs(exact),1e-20L);
      if(error>2e-12L)throw std::runtime_error("Neighbor geometry differs from independently deposited global grid");
      maxError=std::max(maxError,error);++checks;
    }
  }
  int allChecks=0;long double allError=0;MPI_Reduce(&checks,&allChecks,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);MPI_Reduce(&maxError,&allError,1,MPI_LONG_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(rank==0)std::cout<<std::setprecision(21)<<"{\"passed\":true,\"ranks\":"<<size<<",\"cross_interface_straight_and_curved_checks\":"<<allChecks<<",\"max_scaled_work_difference\":"<<allError<<"}\n";
  MPI_Finalize();return 0;
 }catch(const std::exception& e){std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
}
