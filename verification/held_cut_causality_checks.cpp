// Independent finite-domain causal-envelope checks for the unchanged Cowan
// and transverse CPML kernels. No reference-fed or moving boundary is used.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <mpi.h>
#include "eb_cpml.h"
#include "moving_window_remap.h"

namespace {
using namespace aprl;
using aprl::verification::MovingWindowRemap;
constexpr std::size_t NX=24,NY=26,NZ=112;
constexpr unsigned int STEPS=40;
constexpr Double DZ=1.2543262415359034e-6;
std::uint64_t assertions=0,fieldValues=0,historyValues=0,affectedValues=0;
std::size_t lowerFirstFront=0,upperFirstDepth=0;
std::vector<std::size_t> lowerFront,upperDepth;

void check(bool pass,const char* message) {
  ++assertions;if(!pass)throw std::runtime_error(message);
}
bool equal(Double a,Double b) {return std::memcmp(&a,&b,sizeof(Double))==0;}
YeeComponent& field(EBFieldGrid& f,unsigned int a) {
  switch(a) {
    case 0:return f.ex();case 1:return f.ey();case 2:return f.ez();
    case 3:return f.bx();case 4:return f.by();case 5:return f.bz();
  }
  throw std::out_of_range("Invalid component");
}
EBGridGeometry geometry() {
  return EBGridGeometry(NX,NY,NZ,5e-6,5e-6,DZ,DZ/SI::c,EBMaxwellSolver::CowanZ);
}
EBCPMLParameters parameters() {
  EBCPMLParameters p;p.cells[0]=p.cells[1]=10;p.cells[2]=0;
  p.polynomialOrder=3.;p.kappaMax=8.;p.alphaFraction=0.;p.targetReflection=1e-8;
  return p;
}
void seed(EBFieldGrid& f) {
  for(unsigned int a=0;a<6;++a) {
    auto& c=field(f,a);
    for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
      for(std::size_t i=0;i<c.nx();++i) {
        const int tag=static_cast<int>((i*17+j*31+k*43+a*53)%127)-63;
        c(i,j,k)=std::ldexp(static_cast<Double>(tag),a<3?-4:-24);
      }
  }
  PerfectElectricConductorBoundary transverse(false,false);
  transverse.afterElectricUpdate(f);
  auto& p=dynamic_cast<EBConvolutionalPML&>(f.boundary());
  for(std::size_t h=0;h<12;++h) {
    auto& v=p.checkpointHistory(h);
    for(std::size_t i=0;i<v.size();++i)
      v[i]=std::ldexp(static_cast<Double>(static_cast<int>((i+7*h)%31)-15),h<6?-7:-27);
  }
}
struct Halos {
  std::vector<Double> lower[3],upper[3];
  Halos(EBFieldGrid& f) {
    for(unsigned int a=0;a<3;++a) {
      const auto& c=field(f,a);lower[a].resize(c.nx()*c.ny());upper[a].resize(c.nx()*c.ny());
    }
  }
  EBElectricHaloView view() const {
    EBElectricHaloView h;
    h.lowerEx=lower[0].data();h.lowerEy=lower[1].data();h.lowerEz=lower[2].data();
    h.upperEx=upper[0].data();h.upperEy=upper[1].data();h.upperEz=upper[2].data();return h;
  }
  void poison(bool left,unsigned int n) {
    for(unsigned int a=0;a<3;++a) {
      auto& v=left?lower[a]:upper[a];
      for(std::size_t i=0;i<v.size();++i)v[i]=(n+1)*1e7+(a+1)*1e6+(i%71)*1e3;
    }
  }
};
void poison(EBFieldGrid& f,const MovingWindowRemap& layout,bool left,unsigned n) {
  // Alter the entire first/last native plane, including all half-grid arrays.
  // This is stronger than unknown exterior halos or integer face data alone.
  for(unsigned int a=0;a<6;++a) {
    auto& c=field(f,a);const std::size_t k=left?0:c.nz()-1;
    for(std::size_t j=0;j<c.ny();++j)for(std::size_t i=0;i<c.nx();++i)
      c(i,j,k)=(a<3?1e7:1e-2)*(1+n+a)+(1+i+3*j)*(a<3?1e3:1e-6);
  }
  auto& p=dynamic_cast<EBConvolutionalPML&>(f.boundary());
  for(std::size_t h=0;h<12;++h) {
    const auto& l=layout.historyLayout(h);if(l.derivativeAxis==2)continue;
    auto& v=p.checkpointHistory(h);const std::size_t width=l.dimensions[l.derivativeAxis==0?1:0];
    const std::size_t k=left?0:l.dimensions[2]-1;
    for(std::size_t slot=0;slot<l.activeAxisIndices.size();++slot)
      for(std::size_t ij=0;ij<width;++ij)
        v[(slot*l.dimensions[2]+k)*width+ij]=(h<6?1e6:1e-3)*(1+n+h+slot)+ij;
  }
}

void run(bool left) {
  const auto g=geometry();const auto p=parameters();
  EBFieldGrid control(g),altered(g);
  control.setBoundary(std::unique_ptr<EBBoundaryOperator>(new EBConvolutionalPML(g,NZ,0,p,false,false)));
  altered.setBoundary(std::unique_ptr<EBBoundaryOperator>(new EBConvolutionalPML(g,NZ,0,p,false,false)));
  MovingWindowRemap layout(g,p,{-1});seed(control);seed(altered);
  Halos controlHalo(control),alteredHalo(altered);
  for(unsigned int n=1;n<=STEPS;++n) {
    poison(altered,layout,left,n);alteredHalo.poison(left,n);
    control.advanceMagnetic(controlHalo.view());control.advanceElectric();
    altered.advanceMagnetic(alteredHalo.view());altered.advanceElectric();
    std::size_t maxK=0,minK=NZ;std::uint64_t alteredCount=0;
    const auto inspect=[&](Double a,Double b,std::size_t k,bool history) {
      check(std::isfinite(a)&&std::isfinite(b),"Non-finite poisoned propagation");
      if(history)++historyValues;else ++fieldValues;
      if(!equal(a,b)) {
        ++affectedValues;++alteredCount;maxK=std::max(maxK,k);minK=std::min(minK,k);
        // Lower k=0 can first reach Ex/Ey at k=2 through Ez -> B_t -> E_t.
        // Thereafter its maximum forward front advances by one layer/step.
        // Upper initial k=NZ-1 retreats by at most one layer/step.
        check(left?k<=n+1:k>=NZ-1-n,"Boundary influence exceeded native causal envelope");
      }
    };
    for(unsigned int a=0;a<6;++a) {
      const auto& x=field(control,a);const auto& y=field(altered,a);
      for(std::size_t k=0;k<x.nz();++k)for(std::size_t j=0;j<x.ny();++j)
        for(std::size_t i=0;i<x.nx();++i)inspect(x(i,j,k),y(i,j,k),k,false);
    }
    const auto& x=dynamic_cast<const EBConvolutionalPML&>(control.boundary());
    const auto& y=dynamic_cast<const EBConvolutionalPML&>(altered.boundary());
    for(std::size_t h=0;h<12;++h) {
      const auto& l=layout.historyLayout(h);
      if(l.derivativeAxis==2) {
        check(x.checkpointHistory(h).empty()&&y.checkpointHistory(h).empty(),"z CPML is enabled");continue;
      }
      const std::size_t width=l.dimensions[l.derivativeAxis==0?1:0];
      for(std::size_t slot=0;slot<l.activeAxisIndices.size();++slot)
        for(std::size_t k=0;k<l.dimensions[2];++k)for(std::size_t ij=0;ij<width;++ij) {
          const std::size_t index=(slot*l.dimensions[2]+k)*width+ij;
          inspect(x.checkpointHistory(h)[index],y.checkpointHistory(h)[index],k,true);
        }
    }
    check(alteredCount>0,"Poison test did not perturb the propagator");
    if(left) {
      lowerFront.push_back(maxK);if(n==1)lowerFirstFront=maxK;
    }else {
      upperDepth.push_back(NZ-1-minK);if(n==1)upperFirstDepth=NZ-1-minK;
    }
  }
}
void outputArray(std::ostream& out,const std::vector<std::size_t>& a) {
  out<<'[';for(std::size_t i=0;i<a.size();++i){if(i)out<<',';out<<a[i];}out<<']';
}
}
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);int size=0;MPI_Comm_size(MPI_COMM_WORLD,&size);
  try {
    check(size==1,"Causal checks require one rank");run(true);run(false);
    check(lowerFirstFront==2,"Fixture did not exercise two-stage lower first-step reach");
    check(upperFirstDepth==1,"Fixture did not exercise upper one-layer reach");
    std::ofstream file;std::ostream* out=&std::cout;
    if(argc==2){file.open(argv[1]);check(file.good(),"Cannot create report");out=&file;}
    *out<<"{\"format\":\"aprl-held-cut-causality-checks-v1\",\"pass\":true,"
        <<"\"steps_per_side\":"<<STEPS<<",\"assertions\":"<<assertions
        <<",\"field_values_checked\":"<<fieldValues<<",\"history_values_checked\":"<<historyValues
        <<",\"affected_values\":"<<affectedValues<<",\"lower_maximum_native_k\":";
    outputArray(*out,lowerFront);*out<<",\"upper_maximum_native_depth\":";outputArray(*out,upperDepth);
    *out<<",\"scientific_geometry\":false,\"scope\":\"Original Cowan and xy CPML field+history finite propagation only; no particle, initialization-boundary or production physical acceptance\"}\n";
    MPI_Finalize();return 0;
  }catch(const std::exception& e){std::cerr<<"causality checks FAILED: "<<e.what()<<'\n';MPI_Finalize();return 1;}
}
