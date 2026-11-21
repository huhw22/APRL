#include "longitudinal_shape.h"
#include "source_coupling_math.h"
#include "eb_mpi.h"
#include <iomanip>
#include <iostream>
#include <random>
using namespace aprl;
using namespace source_coupling_math;

const auto f3=longitudinal_shape::Shape::LinearBinomial3;

void checkKernel() {
  EBGridGeometry g(22,22,22,1.,1.,1.,1./SI::c,EBMaxwellSolver::CowanZ);
  Vec origin(0.); EBFieldGrid f(g); NativeHalo halo;
  for(auto* b:{&f.bx(),&f.by(),&f.bz()})for(std::size_t n=0;n<b->size();++n)b->data()[n]=std::sin(.245765*n);
  const auto n3=std::size_t(22*22*22);
  auto idx=[](int i,int j,int k){return (k*22+j)*22+i;};
  std::mt19937 gen(20031);std::uniform_real_distribution<double> loc(6.,15.),d(-2.1,2.1);
  double nativeError=0.,continuity=0.,magError=0.,adjoint=0.,oldError=0.;
  for(int path=0;path<151;++path){
    Vec a,b;for(unsigned k=0;k<3;++k){a[k]=loc(gen);b[k]=a[k]+d(gen);}
    if(path==0){a=Vec(8.);b=Vec(11.);}
    if(path==1){a=Vec(10.);b=a;}
    if(path==2){a=Vec(9.5);b=Vec(7.5);}
    std::array<std::vector<double>,3> predicted,explicitJ,defaultJ,oldJ;
    for(auto* p:{&predicted,&explicitJ,&defaultJ,&oldJ})for(auto& arr:*p)arr.assign(n3,0.);
    longitudinal_shape::weights(g,origin,a,b,[&](unsigned axis,int i,int j,int k,double w){
      oldJ[axis][idx(i,j,k)]+=w;
      for(int dx=-1;dx<=1;++dx)predicted[axis][idx(i+dx,j,k)]+=w*(dx==0?.5:.25);
    },f3,false);
    longitudinal_shape::weights(g,origin,a,b,[&](unsigned axis,int i,int j,int k,double w){explicitJ[axis][idx(i,j,k)]+=w;},f3,true);
    longitudinal_shape::weights(g,origin,a,b,[&](unsigned axis,int i,int j,int k,double w){defaultJ[axis][idx(i,j,k)]+=w;},f3);
    std::vector<double> rho(n3,0.);
    for(int end=0;end<2;++end) {
      const auto& x=end?b:a;int i=std::floor(x[0]),j=std::floor(x[1]);
      auto z=longitudinal_shape::basis(x[2],false,f3);auto wx=transverse_source_shape::nodeWeights(x[0]-i);
      double fy=x[1]-j;
      for(unsigned u=0;u<4;++u)for(unsigned v=0;v<2;++v)for(unsigned t=0;t<z.count;++t)
        rho[idx(i-1+u,j+v,z.first+t)]+=(end?1.:-1.)*wx[u]*(v?fy:1.-fy)*z.v[t];
    }
    long double mesh=0.,source=0.;
    for(unsigned k=1;k<22;++k)for(unsigned j=1;j<22;++j)for(unsigned i=1;i<22;++i){
      auto n=idx(i,j,k);double rr=rho[n];int stride[3]={1,22,22*22};
      for(unsigned u=0;u<3;++u) {
        const auto& v=explicitJ[u];auto old=predicted[u][n];double E=std::cos(.21*n+.44*u);
        rr+=(b[u]-a[u])*(v[n]-v[n-stride[u]]);
        nativeError=std::max(nativeError,std::abs(old-v[n]));
        oldError=std::max(oldError,std::abs(defaultJ[u][n]-oldJ[u][n]));
        mesh+=v[n]*(b[u]-a[u])*E;
        source+=oldJ[u][n]*(b[u]-a[u])*(.25*std::cos(.21*(n-1)+.44*u)+.5*E+.25*std::cos(.21*(n+1)+.44*u));
      }continuity=std::max(continuity,std::abs(rr));
    }
    adjoint=std::max(adjoint,(double)std::abs(mesh-source));
    Vec middle=a;middle+=b;middle*=.5;
    if(path==3){middle[0]=.25;middle[1]=.25;middle[2]=.25;}
    if(path==4){middle[0]=21.95;middle[1]=21.2;middle[2]=21.1;}
    Vec ref(0.);
    for(int k=-1;k<=1;++k)for(int i=-1;i<=1;++i){
      Vec x=middle;x[0]+=i;x[2]+=k;
      ref.pmv((k==0?.5:.25)*(i==0?.5:.25),magneticSample(f,origin,x,false,&halo,true).value);
    }
    auto joint=binomialXZMagneticSample(f,origin,middle,halo);for(unsigned u=0;u<3;++u)magError=std::max(magError,std::abs(ref[u]-joint[u]));
  }
  std::cout<<std::setprecision(17)<<"kernel max errors: "<<nativeError<<","<<continuity<<","<<adjoint<<","<<magError<<","<<oldError<<"\n";
  if(std::max({nativeError,continuity,magError,adjoint,oldError})>5.e-13)throw std::runtime_error("Native XZ shape verification failed");
}

int main(int argc,char**argv){
  MPI_Init(&argc,&argv);int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  if(size!=2){MPI_Finalize();return 1;}
  checkKernel();
  double globalError[4]={};
  constexpr int slab=8;
  {
    const auto shape=longitudinal_shape::Shape::LinearBinomial3;
    EBGridGeometry g(18,18,slab,1.,1.,1.,1./SI::c,EBMaxwellSolver::CowanZ);
    EBGridGeometry gr(18,18,slab*size,1.,1.,1.,1./SI::c,EBMaxwellSolver::CowanZ);
    EBFieldGrid f(g),ref(gr);Vec zero(0.),lo(0.);lo[2]=slab*rank;
    longitudinal_shape::Current j(f,MPI_COMM_WORLD,shape,3,true),jr(ref,MPI_COMM_SELF,shape,3,false);
    std::mt19937 rng(20031);
    std::uniform_real_distribution<double> xy(4.,13.),off(-.21,.21),v(-.8,.8);
    auto fill=[](YeeComponent& b,size_t offset,unsigned a){
      for(size_t k=0;k<b.nz();++k)for(size_t y=0;y<b.ny();++y)for(size_t i=0;i<b.nx();++i)b(i,y,k)=sin(.51*i+.317*y+.127*(k+offset)+.39*a);
    };
    YeeComponent* bb[3]={&f.bx(),&f.by(),&f.bz()};YeeComponent* br[3]={&ref.bx(),&ref.by(),&ref.bz()};
    for(unsigned a=0;a<3;++a){fill(*bb[a],slab*rank,a);fill(*br[a],0,a);}
    NativeHalo halo,unpartitioned;halo.exchange(f,MPI_COMM_WORLD,true,3);
    double error[4]={};
    unsigned crossed=0;
    for(int face=1;face<size;++face)for(int n=0;n<128;++n){
      Vec a,b;a[0]=xy(rng);a[1]=xy(rng);a[2]=slab*face+off(rng);
      for(unsigned d=0;d<3;++d)b[d]=a[d]+v(rng);
      const int owner=static_cast<int>(a[2])/slab;
      if(static_cast<int>(b[2])/slab!=owner)++crossed;
      const double q=(n%3?1:-1)/SI::c;
      if(owner==rank){
        j.deposit(lo,a,b,q);
        Vec probe=a;probe+=b;probe*=.5;
        const auto localB=binomialXZMagneticSample(f,lo,probe,halo);
        const auto target=binomialXZMagneticSample(ref,zero,probe,unpartitioned);
        for(int d=0;d<3;++d)error[3]=std::max(error[3],std::abs(target[d]-localB[d]));
      }
      jr.deposit(zero,a,b,q);
    }
    j.complete();
    EBZSlabHaloExchange exchange(MPI_COMM_WORLD);exchange.sumSharedTangentialCurrent(f);
    // Independent serial reference: convolve the completed original source
    // on its native grid, then select this rank's overlapping Yee slab.
    for(unsigned a=0;a<3;++a){
      const auto& local=j.component(a);const auto& all=jr.component(a);
      for(size_t k=0;k<local.nz();++k)for(size_t y=0;y<local.ny();++y)for(size_t x=0;x<local.nx();++x){
        const size_t K=k+slab*rank;
        double target=.5*all(x,y,K);
        if(x>0)target+=.25*all(x-1,y,K);
        if(x+1<all.nx())target+=.25*all(x+1,y,K);
        error[a]=std::max(error[a],std::abs(local(x,y,k)-target));
      }
    }
    MPI_Allreduce(error,globalError,4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(!rank){
      std::cout<<std::setprecision(17)<<"{\"ranks\":"<<size<<",\"interface_paths\":"<<(size-1)*128<<",\"crossing_paths\":"<<crossed<<",\"current_max_error_by_component\":["<<globalError[0]<<","<<globalError[1]<<","<<globalError[2]<<"],\"magnetic_collar_max_error\":"<<globalError[3]<<"}\n";
    }
  }
  MPI_Finalize();return *std::max_element(globalError,globalError+4)>5e-13 ? 2:0;
}
