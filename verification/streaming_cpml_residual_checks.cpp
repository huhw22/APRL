// Verification-only operator checks. No particle or production kernel changes.
#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <mpi.h>

#include "eb_cpml.h"
#include "moving_window_remap.h"
#include "streaming_cpml_boundary.h"
#include "streaming_reference_boundary.h"

namespace {
using aprl::Double;
using aprl::EBFieldGrid;
using aprl::EBGridGeometry;
using aprl::YeeComponent;
using aprl::EBConvolutionalPML;
using aprl::EBCPMLParameters;
using aprl::EBElectricHaloView;
using aprl::verification::ReferenceFedCPML;
using aprl::verification::StreamingReferenceCache;
using aprl::verification::StreamingReferenceBoundary;
using aprl::verification::StreamingBoundaryMode;
using aprl::verification::MovingWindowRemap;
using aprl::verification::MovingWindowHistoryPolicy;

constexpr std::size_t NX=8,NY=6,NZ=104,CACHE_NZ=256,PML=20;
constexpr Double DZ=1e-6;
std::uint64_t checks=0,failures=0,outsideValuesChecked=0;

void check(bool value,const std::string& message) {
  ++checks;
  if(!value){++failures;throw std::runtime_error(message);}
}
bool sameDouble(Double a,Double b) {
  return std::memcmp(&a,&b,sizeof(Double))==0;
}
EBGridGeometry geometry(std::size_t nz=NZ) {
  return EBGridGeometry(NX,NY,nz,2*DZ,2*DZ,DZ,DZ/aprl::SI::c,
                       aprl::EBMaxwellSolver::CowanZ);
}
EBCPMLParameters parameters() {
  EBCPMLParameters result;
  result.cells[0]=result.cells[1]=0;result.cells[2]=PML;
  return result;
}
YeeComponent& component(EBFieldGrid& f,unsigned a) {
  switch(a) {
    case 0:return f.ex();case 1:return f.ey();case 2:return f.ez();
    case 3:return f.bx();case 4:return f.by();case 5:return f.bz();
    case 6:return f.jx();case 7:return f.jy();case 8:return f.jz();
  }
  throw std::logic_error("Bad component");
}
const YeeComponent& component(const EBFieldGrid& f,unsigned a) {
  return component(const_cast<EBFieldGrid&>(f),a);
}
using FieldCopy=std::array<std::vector<Double>,9>;
FieldCopy snapshot(const EBFieldGrid& f) {
  FieldCopy result;
  for(unsigned a=0;a<9;++a) {
    const auto& c=component(f,a);result[a].assign(c.data(),c.data()+c.size());
  }
  return result;
}
bool equalArray(const YeeComponent& c,const std::vector<Double>& values) {
  return c.size()==values.size()&&
    std::memcmp(c.data(),values.data(),c.size()*sizeof(Double))==0;
}
void compareFields(const EBFieldGrid& a,const EBFieldGrid& b,const std::string& where) {
  for(unsigned c=0;c<9;++c) {
    const auto& x=component(a,c);const auto& y=component(b,c);
    check(x.size()==y.size()&&std::memcmp(x.data(),y.data(),x.size()*sizeof(Double))==0,
          where+": field component "+std::to_string(c)+" is not bitwise equal");
  }
}
void compareHistory(const EBConvolutionalPML& a,const EBConvolutionalPML& b,
                    const std::string& where) {
  check(a.checkpointHistoryCount()==b.checkpointHistoryCount(),where+": history count");
  for(std::size_t h=0;h<a.checkpointHistoryCount();++h) {
    const auto& x=a.checkpointHistory(h);const auto& y=b.checkpointHistory(h);
    check(x.size()==y.size()&&(x.empty()||
      std::memcmp(x.data(),y.data(),x.size()*sizeof(Double))==0),
      where+": history "+std::to_string(h)+" is not bitwise equal");
  }
}
void seedHistory(EBConvolutionalPML& pml) {
  bool nonzero=false;
  for(std::size_t h=0;h<pml.checkpointHistoryCount();++h) {
    auto& values=pml.checkpointHistory(h);
    for(std::size_t i=0;i<values.size();++i) {
      const int tag=static_cast<int>((i+7*h)%43)-21;
      values[i]=std::ldexp(static_cast<Double>(tag),-32);
      nonzero=nonzero||values[i]!=0.;
    }
  }
  check(nonzero,"Nonzero-history fixture did not exercise memory");
}
void fillPattern(EBFieldGrid& f,unsigned tag,bool transversePec=false) {
  for(unsigned a=0;a<9;++a) {
    auto& c=component(f,a);
    for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
      for(std::size_t i=0;i<c.nx();++i) {
        const int code=static_cast<int>((37*(a+1)+11*(i+2)+19*(j+3)+
                                        23*(k+5)+29*tag)%257)-128;
        c(i,j,k)=std::ldexp(static_cast<Double>(code),a<3?-4:(a<6?-20:-10));
      }
  }
  if(transversePec) {
    aprl::PerfectElectricConductorBoundary transverse(false,false);
    transverse.afterElectricUpdate(f);
  }
}
void assertSourcesRestored(const EBFieldGrid& f,const FieldCopy& before,
                          bool magnetic,const std::string& where) {
  for(unsigned a=0;a<9;++a)if((magnetic?a<3:(a>=3&&a<6))||a>=6)
    check(equalArray(component(f,a),before[a]),where+": source/J modified, component "+std::to_string(a));
}

// Independent physical-coordinate target mask; no implementation edge-mask API.
void assertPmlOnlyTargets(const EBFieldGrid& f,const FieldCopy& before,
                         bool magnetic,const std::string& where) {
  for(unsigned a=magnetic?3:0;a<(magnetic?6:3);++a) {
    const auto& c=component(f,a);bool unchanged=true;
    for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
      for(std::size_t i=0;i<c.nx();++i) {
        bool allowed=false;
        if(magnetic&&a<5) {
          const Double z=static_cast<Double>(k)+.5;
          allowed=z<PML||z>NZ-PML;
        }
        if(!magnetic) {
          if(a<2)allowed=k<PML||k>NZ-PML;
          // Existing PEC faces/reference-fed z faces may change separately.
          if(a==0)allowed=allowed||j==0||j==NY||k==0||k==NZ;
          if(a==1)allowed=allowed||i==0||i==NX||k==0||k==NZ;
          if(a==2)allowed=i==0||i==NX||j==0||j==NY;
        }
        if(!allowed) {
          ++outsideValuesChecked;
          const std::size_t index=(k*c.ny()+j)*c.nx()+i;
          unchanged=unchanged&&sameDouble(c(i,j,k),before[a][index]);
        }
      }
    check(unchanged,where+": CPML changed a nonabsorbing target, component "+std::to_string(a));
  }
}

// Independent oracle: replace ALL native source-array values, invoke the
// unchanged core CPML once, then restore the saved source bits. This deliberately
// does not share the candidate's sparse edge ranges or backup implementation.
class FullDomainOracle : public EBConvolutionalPML {
public:
  FullDomainOracle(const EBGridGeometry& g,const EBCPMLParameters& p,
                   StreamingReferenceCache& reference,const std::size_t& offset,
                   bool wrongBTime=false)
    :EBConvolutionalPML(g,g.nz,0,p,true,true),reference_(reference),offset_(offset),
     feed_(reference,offset,StreamingBoundaryMode::ReferenceClamp),wrongBTime_(wrongBTime) {}
  void beginStep(const EBFieldGrid& f,unsigned step) {step_=step;feed_.beginStep(f,step);}
  EBElectricHaloView electricHalos(const EBFieldGrid& f,unsigned step) {
    return feed_.electricHalos(f,step);
  }
  void afterMagneticUpdate(EBFieldGrid& f,const EBElectricHaloView& halo) override {
    const auto saved=snapshot(f);
    subtractSources(f,0,3,step_);
    std::array<std::vector<Double>,6> storage;
    const Double* input[6]={halo.lowerEx,halo.upperEx,halo.lowerEy,halo.upperEy,
                            halo.lowerEz,halo.upperEz};
    for(unsigned a=0;a<3;++a) {
      const auto& c=component(f,a);
      for(unsigned side=0;side<2;++side) {
        const unsigned h=2*a+side;storage[h].resize(c.nx()*c.ny());
        check(input[h]!=nullptr,"Oracle requires explicit matched E halos");
        const Double globalK=static_cast<Double>(offset_)+(side?static_cast<Double>(c.nz()):-1.);
        for(std::size_t j=0;j<c.ny();++j)for(std::size_t i=0;i<c.nx();++i) {
          const std::size_t index=j*c.nx()+i;
          storage[h][index]=input[h][index]-reference_.value(a,i,j,globalK,step_);
        }
      }
    }
    EBElectricHaloView residual;
    residual.lowerEx=storage[0].data();residual.upperEx=storage[1].data();
    residual.lowerEy=storage[2].data();residual.upperEy=storage[3].data();
    residual.lowerEz=storage[4].data();residual.upperEz=storage[5].data();
    EBConvolutionalPML::afterMagneticUpdate(f,residual);
    restore(f,saved,0,3);feed_.afterMagneticUpdate(f,halo);
  }
  void afterElectricUpdate(EBFieldGrid& f) override {
    const auto saved=snapshot(f);
    subtractSources(f,3,6,step_+(wrongBTime_?.5:1.));
    EBConvolutionalPML::afterElectricUpdate(f);
    restore(f,saved,3,6);feed_.afterElectricUpdate(f);
  }
private:
  void subtractSources(EBFieldGrid& f,unsigned first,unsigned last,Double stateStep) {
    for(unsigned a=first;a<last;++a) {
      auto& c=component(f,a);
      for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
        for(std::size_t i=0;i<c.nx();++i)c(i,j,k)-=reference_.value(a,i,j,offset_+k,stateStep);
    }
  }
  static void restore(EBFieldGrid& f,const FieldCopy& saved,unsigned first,unsigned last) {
    for(unsigned a=first;a<last;++a)
      std::copy(saved[a].begin(),saved[a].end(),component(f,a).data());
  }
  StreamingReferenceCache& reference_;const std::size_t& offset_;
  StreamingReferenceBoundary feed_;unsigned step_=0;bool wrongBTime_;
};

void zeroReferenceDegeneracy() {
  EBFieldGrid cacheFields(geometry(CACHE_NZ));cacheFields.clearFields();
  StreamingReferenceCache reference(cacheFields,.23);std::size_t offset=64;
  EBFieldGrid candidate(geometry()),original(geometry());
  fillPattern(candidate,1);fillPattern(original,1);
  auto* adapted=new ReferenceFedCPML(geometry(),parameters(),reference,offset,false,true);
  auto* native=new EBConvolutionalPML(geometry(),NZ,0,parameters(),true,true);
  candidate.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(adapted));
  original.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(native));
  seedHistory(*adapted);seedHistory(*native);
  for(unsigned n=2;n<6;++n) {
    adapted->beginStep(candidate,n);const auto halo=adapted->electricHalos(candidate,n);
    candidate.advanceMagnetic(halo);original.advanceMagnetic(halo);
    compareFields(candidate,original,"Zero reference B");compareHistory(*adapted,*native,"Zero reference B");
    candidate.advanceElectric();original.advanceElectric();
    compareFields(candidate,original,"Zero reference E");compareHistory(*adapted,*native,"Zero reference E");
  }
}

void equalReferenceNoNewResponse(Double beta) {
  EBFieldGrid cacheFields(geometry(CACHE_NZ));fillPattern(cacheFields,3,true);
  StreamingReferenceCache reference(cacheFields,beta);std::size_t offset=64;
  EBFieldGrid candidate(geometry());
  auto* adapted=new ReferenceFedCPML(geometry(),parameters(),reference,offset,false,true);
  candidate.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(adapted));
  const unsigned n=7;
  for(unsigned a=0;a<6;++a) {
    auto& c=component(candidate,a);
    for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
      for(std::size_t i=0;i<c.nx();++i)c(i,j,k)=reference.value(a,i,j,offset+k,a<3?n:n+1);
  }
  adapted->beginStep(candidate,n);const auto halo=adapted->electricHalos(candidate,n);
  auto saved=snapshot(candidate);adapted->afterMagneticUpdate(candidate,halo);
  for(unsigned a=0;a<9;++a)check(equalArray(component(candidate,a),saved[a]),"Equal reference: magnetic response");
  // The electric-stage target is E(n+1); prepare its prescribed outer faces
  // consistently so the separate reference feed also has no change to make.
  for(unsigned a=0;a<3;++a) {
    auto& c=component(candidate,a);
    for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
      for(std::size_t i=0;i<c.nx();++i)c(i,j,k)=reference.value(a,i,j,offset+k,n+1);
  }
  saved=snapshot(candidate);adapted->afterElectricUpdate(candidate);
  for(unsigned a=0;a<9;++a)check(equalArray(component(candidate,a),saved[a]),"Equal reference: electric response");
  for(std::size_t h=0;h<adapted->checkpointHistoryCount();++h)
    for(Double v:adapted->checkpointHistory(h))check(v==0.,"Equal reference created CPML history");
}

// Independent slow world-position remap: full snapshots and direct coordinate
// lookup, including histories. It does not consume MovingWindowRemap plans.
void naiveMinusOne(EBFieldGrid& f,EBConvolutionalPML& cpml,
                   StreamingReferenceCache& reference,std::size_t newOffset,unsigned stateStep) {
  const auto old=snapshot(f);
  for(unsigned a=0;a<9;++a) {
    auto& c=component(f,a);const std::size_t plane=c.nx()*c.ny();
    for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
      for(std::size_t i=0;i<c.nx();++i) {
        const std::size_t transverse=j*c.nx()+i;
        c(i,j,k)=k?old[a][(k-1)*plane+transverse]:
          (a<6?reference.value(a,i,j,newOffset,stateStep):0.);
      }
  }
  const unsigned targets[6]={0,0,1,1,2,2};
  const unsigned axes[6]={1,2,2,0,0,1};
  for(std::size_t h=0;h<cpml.checkpointHistoryCount();++h) {
    auto& current=cpml.checkpointHistory(h);
    if(axes[h%6]!=2) {check(current.empty(),"Disabled transverse history is nonempty");continue;}
    const bool electric=h>=6;const unsigned a=targets[h%6]+(electric?0:3);
    const auto& c=component(f,a);const std::size_t plane=c.nx()*c.ny();
    std::vector<std::size_t> active;
    for(std::size_t k=0;k<c.nz();++k) {
      const Double coordinate=static_cast<Double>(k)+(electric?0.:.5);
      if(coordinate<PML||coordinate>NZ-PML)active.push_back(k);
    }
    check(current.size()==active.size()*plane,"Independent history layout size mismatch");
    const auto saved=current;
    for(std::size_t slot=0;slot<active.size();++slot) {
      const auto k=active[slot];
      const auto found=k?std::lower_bound(active.begin(),active.end(),k-1):active.end();
      const bool keep=k&&found!=active.end()&&*found==k-1;
      for(std::size_t t=0;t<plane;++t)current[slot*plane+t]=keep?
        saved[static_cast<std::size_t>(found-active.begin())*plane+t]:0.;
    }
  }
}

void independentOracleAndShift() {
  EBFieldGrid cacheFields(geometry(CACHE_NZ));fillPattern(cacheFields,9,true);
  StreamingReferenceCache reference(cacheFields,.23);
  std::size_t candidateOffset=64,oracleOffset=64;
  EBFieldGrid candidate(geometry()),oracle(geometry());fillPattern(candidate,11);fillPattern(oracle,11);
  auto* adapted=new ReferenceFedCPML(geometry(),parameters(),reference,candidateOffset,false,true);
  auto* independent=new FullDomainOracle(geometry(),parameters(),reference,oracleOffset);
  candidate.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(adapted));
  oracle.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(independent));
  seedHistory(*adapted);seedHistory(*independent);
  MovingWindowRemap mover(geometry(),parameters(),{-1},MovingWindowHistoryPolicy::WorldPosition);
  for(unsigned n=4;n<8;++n) {
    adapted->beginStep(candidate,n);independent->beginStep(oracle,n);
    const auto candidateHalo=adapted->electricHalos(candidate,n);
    const auto oracleHalo=independent->electricHalos(oracle,n);
    auto before=snapshot(candidate);
    adapted->afterMagneticUpdate(candidate,candidateHalo);independent->afterMagneticUpdate(oracle,oracleHalo);
    assertSourcesRestored(candidate,before,true,"Sparse B callback");
    assertPmlOnlyTargets(candidate,before,true,"Sparse B callback");
    compareFields(candidate,oracle,"Independent B oracle");compareHistory(*adapted,*independent,"Independent B oracle");
    before=snapshot(candidate);
    adapted->afterElectricUpdate(candidate);independent->afterElectricUpdate(oracle);
    assertSourcesRestored(candidate,before,false,"Sparse E callback");
    assertPmlOnlyTargets(candidate,before,false,"Sparse E callback");
    compareFields(candidate,oracle,"Independent E oracle");compareHistory(*adapted,*independent,"Independent E oracle");
    if(n==5) {
      struct Context {StreamingReferenceCache* reference;std::size_t offset;unsigned state;};
      Context context{&reference,candidateOffset-1,n+1};
      aprl::verification::MovingWindowIncoming incoming;incoming.context=&context;
      incoming.value=[](void* raw,unsigned a,std::size_t i,std::size_t j,std::size_t k,Double)->Double {
        const auto& c=*static_cast<Context*>(raw);
        return a<6?c.reference->value(a,i,j,c.offset+k,c.state):0.;
      };
      const auto report=mover.shift(candidate,-1,(candidateOffset-1)*DZ,incoming);
      check(report.historyMovedBytes>0&&report.historyZeroedBytes>0&&report.historyRetainedBytes==0,
            "World-position residual history migration did not exercise overlap/reset");
      naiveMinusOne(oracle,*independent,reference,oracleOffset-1,n+1);
      --candidateOffset;--oracleOffset;
      compareFields(candidate,oracle,"Independent world remap");compareHistory(*adapted,*independent,"Independent world remap");
    }
  }
}

void remappedClampRefresh(bool residualInput) {
  // Deliberately non-PEC reference data exercises the transverse PEC override
  // where a prescribed z face intersects an x/y face.
  EBFieldGrid cacheFields(geometry(CACHE_NZ));fillPattern(cacheFields,23);
  StreamingReferenceCache reference(cacheFields,.23);
  std::size_t offset=64;
  EBFieldGrid candidate(geometry()),oracle(geometry());
  fillPattern(candidate,29,true);fillPattern(oracle,29,true);
  auto* adapted=new ReferenceFedCPML(geometry(),parameters(),reference,offset,false,residualInput);
  auto* native=new EBConvolutionalPML(geometry(),NZ,0,parameters(),true,true);
  candidate.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(adapted));
  oracle.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(native));
  seedHistory(*adapted);seedHistory(*native);
  const unsigned state=8;
  for(EBFieldGrid* fields:{&candidate,&oracle}) {
    // A completed-step fixture: old faces satisfy the clamp, while an old
    // interior point has a nonzero residual that must not become a new face.
    for(unsigned a=0;a<2;++a) {
      auto& c=component(*fields,a);
      for(unsigned side=0;side<2;++side) {
        const std::size_t k=side?c.nz()-1:0;
        for(std::size_t j=0;j<c.ny();++j)for(std::size_t i=0;i<c.nx();++i)
          c(i,j,k)=reference.value(a,i,j,offset+k,state);
      }
      c(2,2,NZ-1)=reference.value(a,2,2,offset+NZ-1,state)+32.;
    }
    aprl::PerfectElectricConductorBoundary transverse(false,false);
    transverse.afterElectricUpdate(*fields);
  }
  struct Context {StreamingReferenceCache* reference;std::size_t offset;unsigned state;};
  Context context{&reference,offset-1,state};
  aprl::verification::MovingWindowIncoming incoming;incoming.context=&context;
  incoming.value=[](void* raw,unsigned a,std::size_t i,std::size_t j,std::size_t k,Double)->Double {
    const auto& c=*static_cast<Context*>(raw);
    return a<6?c.reference->value(a,i,j,c.offset+k,c.state):0.;
  };
  MovingWindowRemap mover(geometry(),parameters(),{-1},MovingWindowHistoryPolicy::WorldPosition);
  mover.shift(candidate,-1,(offset-1)*DZ,incoming);
  naiveMinusOne(oracle,*native,reference,offset-1,state);--offset;
  compareFields(candidate,oracle,"Refresh fixture before clamp");
  compareHistory(*adapted,*native,"Refresh fixture before clamp");
  for(unsigned a=0;a<2;++a) {
    const auto& c=component(candidate,a);
    check(!sameDouble(c(2,2,NZ),reference.value(a,2,2,offset+NZ,state)),
          "Remap fixture did not expose an unclamped new upper face");
    check(sameDouble(c(2,2,0),reference.value(a,2,2,offset,state)),
          "Incoming lower face has the wrong reference time or coordinate");
    check(!sameDouble(reference.value(a,2,2,offset+NZ,state),
                      reference.value(a,2,2,offset+NZ,state-1)),
          "Refresh fixture cannot detect a stale reference time");
  }
  const auto before=snapshot(candidate);
  const auto boundaryValues=adapted->boundaryValues();
  const auto residualValues=adapted->residualInputValues();
  adapted->refreshAfterRemap(candidate,state);
  compareHistory(*adapted,*native,"Refresh changed CPML history");
  assertSourcesRestored(candidate,before,false,"Refresh");
  check(adapted->residualInputValues()==residualValues,"Refresh invoked the residual CPML hook");
  check(adapted->boundaryValues()-boundaryValues==
        2*(candidate.ex().nx()*candidate.ex().ny()+candidate.ey().nx()*candidate.ey().ny()),
        "Refresh did not account for exactly two Ex/Ey faces");
  for(unsigned a=0;a<3;++a) {
    const auto& c=component(candidate,a);bool correct=true;
    for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
      for(std::size_t i=0;i<c.nx();++i) {
        const std::size_t index=(k*c.ny()+j)*c.nx()+i;
        Double expected=before[a][index];
        const bool transverse=(a==0&&(j==0||j==NY))||
          (a==1&&(i==0||i==NX))||(a==2&&(i==0||i==NX||j==0||j==NY));
        if(transverse)expected=0.;
        else if(a<2&&(k==0||k+1==c.nz()))
          expected=reference.value(a,i,j,offset+k,state);
        else ++outsideValuesChecked;
        correct=correct&&sameDouble(c(i,j,k),expected);
      }
    check(correct,"Refresh changed E interior or violated current-time clamp/transverse PEC");
  }
  // A mid-step refresh must fail without changing any state.
  const auto after=snapshot(candidate);adapted->beginStep(candidate,state);
  bool rejected=false;
  try {adapted->refreshAfterRemap(candidate,state);}
  catch(const std::logic_error&) {rejected=true;}
  check(rejected,"Refresh accepted an unfinished field step");
  for(unsigned a=0;a<9;++a)
    check(equalArray(component(candidate,a),after[a]),"Rejected refresh modified a field");
  compareHistory(*adapted,*native,"Rejected refresh modified CPML history");
}

void diagnosticCallbackRoundTrip() {
  EBFieldGrid cacheFields(geometry(CACHE_NZ));fillPattern(cacheFields,31);
  StreamingReferenceCache reference(cacheFields,.23);
  std::array<Double,6> original;
  for(unsigned a=0;a<6;++a)original[a]=reference.value(a,2,3,80.25,7.5);
  const auto initialCount=reference.sampleCount();
  const auto memory=reference.memoryBytes();
  struct Context {
    unsigned calls=0,component=0;std::size_t i=0,j=0;
    Double k=0.,step=0.,result=0.;bool reject=false;
  } context;
  const auto callback=[](void* raw,unsigned a,std::size_t i,std::size_t j,
                         Double k,Double step)->Double {
    auto& c=*static_cast<Context*>(raw);
    ++c.calls;c.component=a;c.i=i;c.j=j;c.k=k;c.step=step;
    if(c.reject)throw std::out_of_range("Intentional diagnostic replay rejection");
    return c.result;
  };
  reference.setDiagnosticQuery(callback,&context);
  for(unsigned a=0;a<6;++a) {
    context.result=a==0?-0.:std::ldexp(static_cast<Double>(a+17),-4);
    // The callback, rather than frozen-cache translation/range checking,
    // owns these coordinates. Record them to catch hidden time/space changes.
    const Double actual=reference.value(a,2,3,-6.25,13.5);
    check(sameDouble(actual,context.result),"Diagnostic callback changed returned bits");
    check(context.component==a&&context.i==2&&context.j==3&&
          sameDouble(context.k,-6.25)&&sameDouble(context.step,13.5),
          "Diagnostic callback received translated coordinates or wrong component");
    check(reference.sampleCount()==initialCount+a+1,"Diagnostic success count mismatch");
  }
  context.reject=true;bool rejected=false;const auto beforeFailure=reference.sampleCount();
  try {reference.value(0,2,3,-6.25,13.5);}
  catch(const std::out_of_range&) {rejected=true;}
  check(rejected&&reference.sampleCount()==beforeFailure,
        "Diagnostic exception was swallowed or counted as a successful query");
  reference.setDiagnosticQuery(nullptr);
  for(unsigned a=0;a<6;++a)
    check(sameDouble(reference.value(a,2,3,80.25,7.5),original[a]),
          "Disabling diagnostic replay did not restore the original cached result");
  check(context.calls==7,"Disabled callback was invoked again");
  check(reference.sampleCount()==initialCount+12,"Restored frozen query count mismatch");
  check(reference.memoryBytes()==memory,"Callback registration changed owned cache storage");
}

void wrongHalfTimeNegativeControl() {
  EBFieldGrid cacheFields(geometry(CACHE_NZ));fillPattern(cacheFields,17,true);
  StreamingReferenceCache reference(cacheFields,.23);std::size_t offset=64;
  EBFieldGrid candidate(geometry()),wrong(geometry());fillPattern(candidate,19);fillPattern(wrong,19);
  ReferenceFedCPML adapted(geometry(),parameters(),reference,offset,false,true);
  FullDomainOracle shiftedHalfStep(geometry(),parameters(),reference,offset,true);
  const unsigned n=7;adapted.beginStep(candidate,n);shiftedHalfStep.beginStep(wrong,n);
  adapted.afterMagneticUpdate(candidate,adapted.electricHalos(candidate,n));
  shiftedHalfStep.afterMagneticUpdate(wrong,shiftedHalfStep.electricHalos(wrong,n));
  compareFields(candidate,wrong,"Wrong-time control before B-reference use");
  adapted.afterElectricUpdate(candidate);shiftedHalfStep.afterElectricUpdate(wrong);
  bool difference=false;
  for(unsigned a=0;a<3;++a) {
    const auto& x=component(candidate,a);const auto& y=component(wrong,a);
    difference=difference||std::memcmp(x.data(),y.data(),x.size()*sizeof(Double))!=0;
  }
  check(difference,"Dynamic beta=.23 fixture cannot detect erroneous B stateStep=n+.5");
}

void rejectedConfigurations() {
  EBFieldGrid cacheFields(geometry(CACHE_NZ));cacheFields.clearFields();
  StreamingReferenceCache reference(cacheFields,0.);std::size_t offset=64;
  for(unsigned axis=0;axis<2;++axis) {
    auto p=parameters();p.cells[axis]=2;bool rejected=false;
    try {ReferenceFedCPML unsupported(geometry(),p,reference,offset,false,true);}
    catch(const std::invalid_argument&) {rejected=true;}
    check(rejected,"Residual mode accepted transverse CPML");
  }
  bool rejected=false;
  try {ReferenceFedCPML unsupported(geometry(),parameters(),reference,offset,true,true);}
  catch(const std::invalid_argument&) {rejected=true;}
  check(rejected,"Residual mode accepted contradictory zeroFeed");
}

void writeMetrics(const std::string& directory,bool complete) {
  std::ofstream out(directory+"/metrics.json");
  if(!out.good())throw std::runtime_error("Cannot write residual checks metrics");
  out<<"{\"format\":\"aprl-streaming-cpml-residual-checks-v1\",\"completed\":"
     <<(complete?"true":"false")<<",\"checks\":"<<checks<<",\"bitwise_failures\":"<<failures
     <<",\"outside_target_values_checked\":"<<outsideValuesChecked
     <<",\"zero_reference_full_steps\":4,\"equal_reference_beta_cases\":2"
       ",\"independent_oracle_callback_steps\":4,\"world_position_shifts\":1"
       ",\"refresh_world_position_shifts\":2,\"refresh_total_and_residual_modes\":2"
       ",\"diagnostic_callback_success_queries\":6,\"diagnostic_callback_exception_checks\":1"
       ",\"diagnostic_callback_restore_queries\":6"
       ",\"scope\":\"isolated CPML operators and zero-reference full Cowan steps; no physical approximation acceptance\"}\n";
  out.close();if(!out.good())throw std::runtime_error("Cannot close residual checks metrics");
}
} // namespace

int main(int argc,char** argv) {
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int size=0;MPI_Comm_size(MPI_COMM_WORLD,&size);
  std::string directory;
  try {
    check(size==1,"Residual checks require one MPI rank");
    check(argc==2,"Usage: aprl_streaming_cpml_residual_checks OUT_DIRECTORY");
    directory=argv[1];
    check(mkdir(directory.c_str(),0755)==0||errno==EEXIST,"Cannot create checks directory");
    check(!std::ifstream(directory+"/metrics.json").good(),"Refusing to overwrite checks metrics");
    rejectedConfigurations();zeroReferenceDegeneracy();
    equalReferenceNoNewResponse(0.);equalReferenceNoNewResponse(.23);
    independentOracleAndShift();wrongHalfTimeNegativeControl();
    remappedClampRefresh(false);remappedClampRefresh(true);diagnosticCallbackRoundTrip();
    writeMetrics(directory,true);
    std::cout<<"Residual CPML checks PASS: "<<checks<<" checks\n";
    MPI_Finalize();return 0;
  } catch(const std::exception& error) {
    std::cerr<<"Residual CPML operator checks FAILED: "<<error.what()<<'\n';
    if(!directory.empty()&&!std::ifstream(directory+"/metrics.json").good())
      try {writeMetrics(directory,false);} catch(...) {}
    MPI_Finalize();return 1;
  }
}
