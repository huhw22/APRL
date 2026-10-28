// Independent native-index migration and MPI/serial propagation checks.
// No particle, production-length or radiation-accuracy acceptance is implied.
#include "moving_window_fields.h"
#include "eb_mpi.h"
#include "compute_threads.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using namespace aprl;
constexpr std::size_t NX=24,NY=26,NZ=48;
constexpr Double DZ=1.2543262415359034e-6;
std::uint64_t checked=0,causalChecked=0,poisonDifferences=0;
double maximumFieldRelative=0.,maximumHistoryRelative=0.;
void require(bool good,const char* message) {if(!good)throw std::runtime_error(message);}
bool equal(Double a,Double b) {return std::memcmp(&a,&b,sizeof(Double))==0;}
YeeComponent& component(EBFieldGrid& f,unsigned a) {
  switch(a){case 0:return f.ex();case 1:return f.ey();case 2:return f.ez();
    case 3:return f.bx();case 4:return f.by();case 5:return f.bz();
    case 6:return f.jx();case 7:return f.jy();case 8:return f.jz();}
  throw std::out_of_range("Check component");
}
EBGridGeometry geometry(std::size_t nz) {
  return EBGridGeometry(NX,NY,nz,5e-6,5e-6,DZ,DZ/SI::c,EBMaxwellSolver::CowanZ);
}
EBCPMLParameters parameters() {
  EBCPMLParameters p;p.cells[0]=p.cells[1]=10;p.cells[2]=0;
  p.polynomialOrder=3.;p.kappaMax=8.;p.alphaFraction=0.;p.targetReflection=1e-8;return p;
}
struct Shape {std::size_t row,layers,slots;};
Shape historyShape(unsigned h,std::size_t nz) {
  // Explicit independent list in the documented magnetic/electric order.
  // Longitudinal derivative terms are empty with transverse-only CPML.
  switch(h) {
    case 0:return {NX+1,nz,20};case 1:return {0,nz,0};case 2:return {0,nz,0};
    case 3:return {NY+1,nz,20};case 4:return {NY,nz+1,20};case 5:return {NX,nz+1,20};
    case 6:return {NX,nz+1,20};case 7:return {0,nz+1,0};case 8:return {0,nz+1,0};
    case 9:return {NY,nz+1,20};case 10:return {NY+1,nz,20};case 11:return {NX+1,nz,20};
  }
  throw std::out_of_range("Check history");
}
Double marker(void*,unsigned a,std::size_t i,std::size_t j,std::ptrdiff_t k,std::size_t n) {
  return std::ldexp(Double(a+1),40)+std::ldexp(Double(i),25)+
    std::ldexp(Double(j),15)+Double(k+1024)*16+Double(n);
}
Double historyMarker(unsigned h,std::size_t slot,std::size_t q,std::ptrdiff_t k,std::size_t n) {
  return marker(nullptr,h,slot,q,k,n);
}
Double wave(void*,unsigned a,std::size_t i,std::size_t j,std::ptrdiff_t k,std::size_t n) {
  const double time=static_cast<double>(n)-(a>=3?.5:0.);
  return std::sin(.17*i+.11*j+.23*k-.19*time+.31*a)*(a<3?1.:1./SI::c);
}
void fill(EBFieldGrid& f,std::size_t local,std::ptrdiff_t offset,std::size_t n,
          MovingWindowReference::Function fn,bool histories) {
  for(unsigned a=0;a<9;++a) {
    auto& c=component(f,a);
    for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
      for(std::size_t i=0;i<c.nx();++i)c(i,j,k)=fn(nullptr,a,i,j,offset+static_cast<std::ptrdiff_t>(local+k),n);
  }
  if(histories) {
    auto& p=dynamic_cast<EBConvolutionalPML&>(f.boundary());
    for(unsigned h=0;h<12;++h) {
      auto& v=p.checkpointHistory(h);const auto s=historyShape(h,f.geometry().nz);
      require(v.size()==s.row*s.layers*s.slots,"Independent history shape mismatch");
      for(std::size_t slot=0;slot<s.slots;++slot)for(std::size_t k=0;k<s.layers;++k)
        for(std::size_t q=0;q<s.row;++q)v[(slot*s.layers+k)*s.row+q]=historyMarker(h,slot,q,offset+static_cast<std::ptrdiff_t>(local+k),n);
    }
  }
}
void stateRemap(int rank,int size,bool reset) {
  const std::size_t local=NZ/size*rank;EBFieldGrid f(geometry(NZ/size));
  MovingWindowReference ref;ref.value=marker;
  MovingWindowFields window(f.geometry(),NZ,local,MPI_COMM_WORLD,parameters(),ref);
  window.installBoundary(f);EBZSlabHaloExchange mpi(MPI_COMM_WORLD);
  bool rejected=false;try {window.shift(f,1,99,reset);}catch(const std::logic_error&){rejected=true;}
  require(rejected,"Shift before complete field step was accepted");
  std::ptrdiff_t offset=100;fill(f,local,offset,0,marker,true);
  for(std::size_t n=0;n<3;++n) {
    window.beginStep(f,n,offset);
    rejected=false;try {window.beginStep(f,n,offset);}catch(const std::logic_error&){rejected=true;}
    require(rejected,"Repeated beginStep was accepted");
    f.clearCurrent();mpi.advanceMagnetic(f,window.physicalElectricHalos(f,n));mpi.advanceElectric(f);
    // Replace the completed state with exactly representable index labels.
    // Only the remap, not the advance arithmetic, is tested bitwise here.
    fill(f,local,offset,n+1,marker,true);
    --offset;const auto report=window.shift(f,n+1,offset,reset);
    require(report.fieldMovedBytes>0&&report.historyMovedBytes>0,"Empty remap");
    if(size==1)require(report.sentBytes==0&&report.receivedBytes==0,"Single rank communicated");
    for(unsigned a=0;a<9;++a) {
      auto& c=component(f,a);
      for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
        for(std::size_t i=0;i<c.nx();++i) {
          const Double want=(a>=6&&(reset||(rank==0&&k==0)))?0.:
            marker(nullptr,a,i,j,offset+static_cast<std::ptrdiff_t>(local+k),n+1);
          ++checked;require(equal(c(i,j,k),want),"Native component migration differs bitwise");
        }
    }
    auto& p=dynamic_cast<EBConvolutionalPML&>(f.boundary());
    for(unsigned h=0;h<12;++h) {
      const auto& v=p.checkpointHistory(h);const auto s=historyShape(h,f.geometry().nz);
      for(std::size_t slot=0;slot<s.slots;++slot)for(std::size_t k=0;k<s.layers;++k)
        for(std::size_t q=0;q<s.row;++q) {
          const Double want=rank==0&&k==0?0.:historyMarker(h,slot,q,offset+static_cast<std::ptrdiff_t>(local+k),n+1);
          ++checked;require(equal(v[(slot*s.layers+k)*s.row+q],want),"CPML migration differs bitwise");
        }
    }
  }
}
void current(EBFieldGrid& f,std::size_t local,std::ptrdiff_t offset,std::size_t n,int rank,int size) {
  f.clearCurrent();
  for(unsigned a=6;a<9;++a) {
    auto& c=component(f,a);
    for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
      for(std::size_t i=0;i<c.nx();++i) {
        Double value=1e-5*std::sin(.1*i+.2*j+.3*(offset+static_cast<std::ptrdiff_t>(local+k))+.07*n+a);
        if(a<8&&((k==0&&rank>0)||(k==f.geometry().nz&&rank+1<size)))value*=.5;
        c(i,j,k)=value;
      }
  }
}
void propagation(int rank,int size) {
  const auto local=NZ/size*rank;EBFieldGrid f(geometry(NZ/size)),serial(geometry(NZ));
  MovingWindowReference ref;ref.value=wave;
  MovingWindowFields w(f.geometry(),NZ,local,MPI_COMM_WORLD,parameters(),ref);
  MovingWindowFields s(serial.geometry(),NZ,0,MPI_COMM_SELF,parameters(),ref);
  w.installBoundary(f);s.installBoundary(serial);
  EBZSlabHaloExchange mpi(MPI_COMM_WORLD),self(MPI_COMM_SELF);
  std::ptrdiff_t offset=100;fill(f,local,offset,0,wave,false);fill(serial,0,offset,0,wave,false);
  std::vector<Double> poison[3];for(unsigned a=0;a<3;++a)poison[a].resize(component(f,a).nx()*component(f,a).ny(),1e50);
  for(std::size_t n=0;n<8;++n) {
    w.beginStep(f,n,offset);s.beginStep(serial,n,offset);
    auto halo=w.physicalElectricHalos(f,n);
    // Supplying invalid physical values for INTERNAL ends must not replace
    // genuine neighbour data in the MPI magnetic operator.
    if(rank>0){halo.lowerEx=poison[0].data();halo.lowerEy=poison[1].data();halo.lowerEz=poison[2].data();}
    if(rank+1<size){halo.upperEx=poison[0].data();halo.upperEy=poison[1].data();halo.upperEz=poison[2].data();}
    current(f,local,offset,n,rank,size);current(serial,0,offset,n,0,1);
    mpi.advanceMagnetic(f,halo);self.advanceMagnetic(serial,s.physicalElectricHalos(serial,n));
    mpi.advanceElectric(f);self.advanceElectric(serial);
    --offset;w.shift(f,n+1,offset);s.shift(serial,n+1,offset);
    for(unsigned a=0;a<6;++a) {
      auto& c=component(f,a);auto& d=component(serial,a);long double sum=0.,delta=0.;
      for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
        for(std::size_t i=0;i<c.nx();++i) {
          const long double value=d(i,j,local+k),error=static_cast<long double>(c(i,j,k))-value;
          require(std::isfinite(error),"Nonfinite distributed field");sum+=value*value;delta+=error*error;
        }
      const double error=std::sqrt(delta/std::max(sum,1e-90L));
      maximumFieldRelative=std::max(maximumFieldRelative,error);
      require(error<=2e-12,"Distributed evolution exceeds roundoff-scale field threshold");
    }
    auto& a=dynamic_cast<EBConvolutionalPML&>(f.boundary());
    auto& b=dynamic_cast<EBConvolutionalPML&>(serial.boundary());
    for(unsigned h=0;h<12;++h) {
      const auto x=historyShape(h,f.geometry().nz),y=historyShape(h,NZ);
      const auto& av=a.checkpointHistory(h);const auto& bv=b.checkpointHistory(h);
      long double sum=0.,delta=0.;
      for(std::size_t slot=0;slot<x.slots;++slot)for(std::size_t k=0;k<x.layers;++k)
        for(std::size_t q=0;q<x.row;++q) {
          const long double value=bv[(slot*y.layers+local+k)*y.row+q];
          const long double error=static_cast<long double>(av[(slot*x.layers+k)*x.row+q])-value;
          require(std::isfinite(error),"Nonfinite distributed history");sum+=value*value;delta+=error*error;
        }
      const double error=std::sqrt(delta/std::max(sum,1e-90L));
      maximumHistoryRelative=std::max(maximumHistoryRelative,error);
      require(error<=2e-12,"Distributed evolution exceeds roundoff-scale history threshold");
    }
  }
}
void causality(int rank,int size) {
  const std::size_t local=NZ/size*rank;EBFieldGrid a(geometry(NZ/size)),b(geometry(NZ/size));
  MovingWindowReference ref;ref.value=wave;
  MovingWindowFields x(a.geometry(),NZ,local,MPI_COMM_WORLD,parameters(),ref);
  MovingWindowFields y(b.geometry(),NZ,local,MPI_COMM_WORLD,parameters(),ref);
  x.installBoundary(a);y.installBoundary(b);EBZSlabHaloExchange am(MPI_COMM_WORLD),bm(MPI_COMM_WORLD);
  std::ptrdiff_t offset=100;fill(a,local,offset,0,wave,false);fill(b,local,offset,0,wave,false);
  std::vector<Double> poison[3];for(unsigned c=0;c<3;++c)poison[c].resize(component(b,c).nx()*component(b,c).ny(),1e9);
  for(std::size_t n=0;n<12;++n) {
    if(rank+1==size)for(unsigned c:{0u,1u,5u}) {
      auto& v=component(b,c);
      for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i)
        v(i,j,b.geometry().nz)=(n+1)*(c<3?1e7:1e-1)+i+3*j;
    }
    a.clearCurrent();b.clearCurrent();x.beginStep(a,n,offset);y.beginStep(b,n,offset);
    auto halo=y.physicalElectricHalos(b,n);
    if(rank+1==size){halo.upperEx=poison[0].data();halo.upperEy=poison[1].data();halo.upperEz=poison[2].data();}
    am.advanceMagnetic(a,x.physicalElectricHalos(a,n));bm.advanceMagnetic(b,halo);
    am.advanceElectric(a);bm.advanceElectric(b);
    if(rank+1==size)for(unsigned c=0;c<6;++c) {
      auto& av=component(a,c);auto& bv=component(b,c);
      for(std::size_t q=0;q<av.size();++q)if(!equal(av.data()[q],bv.data()[q]))++poisonDifferences;
    }
    --offset;x.shift(a,n+1,offset);y.shift(b,n+1,offset);
    for(unsigned c=0;c<6;++c) {
      auto& av=component(a,c);auto& bv=component(b,c);
      for(std::size_t k=0;k<av.nz()&&local+k<NZ;++k)
        for(std::size_t j=0;j<av.ny();++j)for(std::size_t i=0;i<av.nx();++i) {
          ++causalChecked;require(equal(av(i,j,k),bv(i,j,k)),"Right-end perturbation entered retained field");
        }
    }
    auto& ap=dynamic_cast<EBConvolutionalPML&>(a.boundary());auto& bp=dynamic_cast<EBConvolutionalPML&>(b.boundary());
    for(unsigned h=0;h<12;++h) {
      const auto s=historyShape(h,a.geometry().nz);const auto& av=ap.checkpointHistory(h);const auto& bv=bp.checkpointHistory(h);
      for(std::size_t slot=0;slot<s.slots;++slot)for(std::size_t k=0;k<s.layers&&local+k<NZ;++k)
        for(std::size_t q=0;q<s.row;++q) {
          const auto index=(slot*s.layers+k)*s.row+q;++causalChecked;
          require(equal(av[index],bv[index]),"Right-end perturbation entered retained CPML history");
        }
    }
  }
}
}

int main(int argc,char** argv) {
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int rank=0,size=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  try {
    configureComputeThreads(provided);require(size==1||size==2||size==4,"Use 1, 2 or 4 MPI ranks");
    stateRemap(rank,size,false);stateRemap(rank,size,true);propagation(rank,size);causality(rank,size);
    unsigned long long local[3]={checked,causalChecked,poisonDifferences},global[3]={};
    MPI_Allreduce(local,global,3,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    double errors[2]={maximumFieldRelative,maximumHistoryRelative},allErrors[2]={};
    MPI_Allreduce(errors,allErrors,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    require(global[2]>0,"Right-end perturbation did not affect the sacrifice collar");
    if(rank==0) {
      std::ofstream file;if(argc>1)file.open(argv[1]);
      if(argc>1)require(file.good(),"Cannot open check output");
      std::ostream& out=argc>1?static_cast<std::ostream&>(file):std::cout;
      out<<std::setprecision(17)<<"{\"status\":\"PASS\",\"ranks\":"<<size
        <<",\"remap_values_bitwise_checked\":"<<global[0]<<",\"causal_values_bitwise_checked\":"<<global[1]
        <<",\"poisoned_collar_differences\":"<<global[2]<<",\"maximum_field_relative_l2\":"<<allErrors[0]
        <<",\"maximum_history_relative_l2\":"<<allErrors[1]
        <<",\"evolution_relative_tolerance\":2e-12,\"scope\":\"native fields and transverse CPML; fixed slabs; no particle/production/radiation acceptance\"}\n";
    }
    MPI_Finalize();return 0;
  }catch(const std::exception& e){std::cerr<<"rank "<<rank<<": "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;}
}
