#include "moving_window_fields.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace aprl {
namespace {
using Clock = std::chrono::steady_clock;
double seconds(const Clock::time_point& start) {
  return std::chrono::duration<double>(Clock::now()-start).count();
}
std::size_t multiply(std::size_t a,std::size_t b) {
  if(b && a>std::numeric_limits<std::size_t>::max()/b)
    throw std::overflow_error("Moving-window allocation overflow");
  return a*b;
}
YeeComponent& component(EBFieldGrid& f,unsigned a) {
  switch(a) {
    case 0:return f.ex();case 1:return f.ey();case 2:return f.ez();
    case 3:return f.bx();case 4:return f.by();case 5:return f.bz();
    case 6:return f.jx();case 7:return f.jy();case 8:return f.jz();
  }
  throw std::out_of_range("Moving-window component");
}
const YeeComponent& component(const EBFieldGrid& f,unsigned a) {
  return component(const_cast<EBFieldGrid&>(f),a);
}
std::array<std::size_t,3> shape(const EBGridGeometry& g,unsigned a) {
  const bool electric=a<3||a>=6;const unsigned axis=a%3;
  const std::size_t cells[3]={g.nx,g.ny,g.nz};
  std::array<std::size_t,3> result;
  for(unsigned d=0;d<3;++d)result[d]=cells[d]+(electric?d!=axis:d==axis);
  return result;
}
void sameGeometry(const EBGridGeometry& a,const EBGridGeometry& b) {
  if(a.nx!=b.nx||a.ny!=b.ny||a.nz!=b.nz||a.dx!=b.dx||a.dy!=b.dy||
     a.dz!=b.dz||a.dt!=b.dt||a.solver!=b.solver)
    throw std::invalid_argument("Moving-window geometry changed");
}
}

struct MovingWindowFields::Impl {
  struct Boundary : EBConvolutionalPML {
    enum class Phase {Idle, Magnetic, Electric, Complete};
    MovingWindowFields::Impl& owner;
    Phase phase=Phase::Idle;
    bool established=false;
    std::size_t step=0;
    std::ptrdiff_t offset=0;
    std::vector<Double> lower[3],upper[3];
    std::array<long double,6> startPower{};
    std::array<long double,6> overwritePowerDelta{};
    bool energyStepFinished=false;
    std::uint64_t queries=0;

    explicit Boundary(MovingWindowFields::Impl& p):EBConvolutionalPML(p.geometry,p.globalNz,
        p.localOffset,p.parameters,false,false),owner(p) {
      for(unsigned a=0;a<3;++a) {
        const auto s=shape(p.geometry,a);const auto count=multiply(s[0],s[1]);
        if(p.rank==0)lower[a].resize(count);
        if(p.rank+1==p.size)upper[a].resize(count,0.);
      }
    }
    Double reference(unsigned a,std::size_t i,std::size_t j,
                     std::ptrdiff_t k,std::size_t n) {
      const Double result=owner.reference.value(owner.reference.context,a,i,j,k,n);
      if(!std::isfinite(result))throw std::runtime_error("Non-finite moving-window inlet");
      ++queries;return result;
    }
    void begin(EBFieldGrid& f,std::size_t n,std::ptrdiff_t k) {
      if(phase!=Phase::Idle || (established&&(n!=step||k!=offset)))
        throw std::logic_error("Moving-window beginStep phase/origin mismatch");
      if(n==std::numeric_limits<std::size_t>::max() ||
         k==std::numeric_limits<std::ptrdiff_t>::min())
        throw std::overflow_error("Moving-window step/origin overflow");
      step=n;offset=k;established=true;
      if(owner.energy.enabled) {
        const auto diagnosticStart=Clock::now();
        startPower=owner.boundaryPower(f);
        overwritePowerDelta.fill(0.);energyStepFinished=false;
        owner.energy.diagnosticSeconds+=seconds(diagnosticStart);
      }
      if(owner.rank==0)for(unsigned a=0;a<3;++a) {
        const auto& c=component(f,a);
        for(std::size_t j=0;j<c.ny();++j)for(std::size_t i=0;i<c.nx();++i)
          lower[a][j*c.nx()+i]=reference(a,i,j,k-1,n);
      }
      phase=Phase::Magnetic;
    }
    EBElectricHaloView halos(std::size_t n) const {
      if(phase!=Phase::Magnetic||n!=step)
        throw std::logic_error("Moving-window halo phase mismatch");
      EBElectricHaloView h;
      if(owner.rank==0){h.lowerEx=lower[0].data();h.lowerEy=lower[1].data();h.lowerEz=lower[2].data();}
      if(owner.rank+1==owner.size){h.upperEx=upper[0].data();h.upperEy=upper[1].data();h.upperEz=upper[2].data();}
      return h;
    }
    void afterMagneticUpdate(EBFieldGrid& f,const EBElectricHaloView& h) override {
      if(phase!=Phase::Magnetic)throw std::logic_error("Moving-window magnetic phase mismatch");
      if(owner.stageObserver)owner.stageObserver(owner.stageContext,"magnetic_curl",f);
      EBConvolutionalPML::afterMagneticUpdate(f,h);
      if(owner.stageObserver)owner.stageObserver(owner.stageContext,"magnetic_cpml",f);
      phase=Phase::Electric;
    }
    void afterElectricUpdate(EBFieldGrid& f) override {
      if(phase!=Phase::Electric)throw std::logic_error("Moving-window electric phase mismatch");
      if(owner.stageObserver)owner.stageObserver(owner.stageContext,"electric_curl_local",f);
      EBConvolutionalPML::afterElectricUpdate(f);
      if(owner.stageObserver)owner.stageObserver(owner.stageContext,"electric_cpml",f);
      if(owner.beforeInletObserver)
        owner.beforeInletObserver(owner.beforeInletContext,f);
      // Only physical cell k=0 uses the overwritten native Ex/Ey k=0 plane.
      // Measure its collocated energy, including the unchanged k=1 samples;
      // squaring the replaced native plane alone would miss the cross terms.
      long double energyBefore=0.;
      std::array<long double,6> powerBefore{};
      const bool measure=owner.energy.enabled&&owner.rank==0;
      if(measure) {
        const auto diagnosticStart=Clock::now();
        energyBefore=owner.planeEnergy(f,0);
        powerBefore=owner.boundaryPower(f,true);
        owner.energy.diagnosticSeconds+=seconds(diagnosticStart);
      }
      if(owner.rank==0)for(unsigned a=0;a<2;++a) {
        auto& c=component(f,a);
        for(std::size_t j=0;j<c.ny();++j)for(std::size_t i=0;i<c.nx();++i)
          c(i,j,0)=reference(a,i,j,offset,step+1);
      }
      if(measure) {
        const auto diagnosticStart=Clock::now();
        owner.energy.inletOverwriteEnergyJ+=owner.planeEnergy(f,0)-energyBefore;
        const auto powerAfter=owner.boundaryPower(f,true);
        for(unsigned face=0;face<6;++face)
          overwritePowerDelta[face]=powerAfter[face]-powerBefore[face];
        ++owner.energy.overwriteEvents;
        owner.energy.diagnosticSeconds+=seconds(diagnosticStart);
      }
      if(owner.stageObserver)owner.stageObserver(owner.stageContext,"inlet_overwrite",f);
      phase=Phase::Complete;
    }
    std::size_t memoryBytes() const override {
      std::size_t n=EBConvolutionalPML::memoryBytes();
      for(unsigned a=0;a<3;++a)n+=(lower[a].capacity()+upper[a].capacity())*sizeof(Double);
      return n;
    }
  };
  struct History {std::size_t layers=0,row=0,slots=0,values=0;};
  struct Packet {
    unsigned fields=0;
    std::array<std::size_t,9> fieldOffsets;
    std::array<std::size_t,12> historyOffsets;
    std::size_t values=0;
  };
  EBGridGeometry geometry;
  std::size_t globalNz,localOffset;
  MPI_Comm communicator;
  EBCPMLParameters parameters;
  MovingWindowReference reference;
  int rank=0,size=1,lowerRank=MPI_PROC_NULL,upperRank=MPI_PROC_NULL;
  Boundary* boundary=nullptr; // Owned by EBFieldGrid, never deleted here.
  std::array<std::array<std::size_t,3>,9> shapes;
  std::array<History,12> histories;
  Packet resetPacket,preservePacket;
  std::vector<Double> send,receive;
  MovingWindowFieldEnergyLedger energy;
  BeforeInletObserver beforeInletObserver = nullptr;
  void* beforeInletContext = nullptr;
  StageObserver stageObserver = nullptr;
  void* stageContext = nullptr;

  Impl(const EBGridGeometry& g,std::size_t nz,std::size_t local,
       MPI_Comm comm,const EBCPMLParameters& p,MovingWindowReference r)
    :geometry(g),globalNz(nz),localOffset(local),communicator(comm),parameters(p),reference(r) {
    MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
    const auto ranks=static_cast<std::size_t>(size);
    const auto base=nz/ranks,remainder=nz%ranks;
    if(g.solver!=EBMaxwellSolver::CowanZ||g.nz<2||g.nx<2||g.ny<2||
       g.nz!=base+(static_cast<std::size_t>(rank)<remainder)||
       local!=base*rank+std::min(static_cast<std::size_t>(rank),remainder)||
       p.cells[2]!=0||p.cells[0]<2||p.cells[1]<2||
       p.cells[0]>=g.nx||p.cells[0]>=g.nx-p.cells[0]||
       p.cells[1]>=g.ny||p.cells[1]>=g.ny-p.cells[1]||
       (rank==0&&!r.value))
      throw std::invalid_argument("Moving-window requires fixed Cowan z slabs, xy-only CPML and an inlet reference");
    if(!std::isfinite(g.dt)||!std::isfinite(g.dz)||!(g.dz>0)||
       std::abs(SI::c*g.dt-g.dz)>64*std::numeric_limits<Double>::epsilon()*g.dz)
      throw std::invalid_argument("Moving-window requires c*dt=dz");
    lowerRank=rank?rank-1:MPI_PROC_NULL;upperRank=rank+1<size?rank+1:MPI_PROC_NULL;
    for(unsigned a=0;a<9;++a)shapes[a]=shape(g,a);
    const unsigned pairs[6][2]={{0,1},{0,2},{1,2},{1,0},{2,0},{2,1}};
    for(unsigned h=0;h<12;++h) {
      const bool electric=h>=6;const unsigned axis=pairs[h%6][1];
      const auto& s=shapes[pairs[h%6][0]+(electric?0:3)];
      auto& v=histories[h];v.layers=s[2];
      if(axis!=2){v.row=s[axis==0?1:0];v.slots=multiply(2,p.cells[axis]);}
      v.values=multiply(multiply(v.row,v.slots),v.layers);
    }
    resetPacket=packet(6);preservePacket=packet(9);
    if(preservePacket.values>static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::overflow_error("Moving-window neighbour plane exceeds MPI count");
    send.resize(preservePacket.values);receive.resize(preservePacket.values);
  }
  Packet packet(unsigned count) const {
    Packet p;p.fields=count;p.fieldOffsets.fill(0);p.historyOffsets.fill(0);
    const auto append=[&](std::size_t n) {
      if(n>std::numeric_limits<std::size_t>::max()-p.values)
        throw std::overflow_error("Moving-window packet overflow");
      p.values+=n;
    };
    for(unsigned a=0;a<count;++a){p.fieldOffsets[a]=p.values;append(multiply(shapes[a][0],shapes[a][1]));}
    for(unsigned h=0;h<12;++h){p.historyOffsets[h]=p.values;append(multiply(histories[h].row,histories[h].slots));}
    return p;
  }
  void validate(const EBFieldGrid& f) const {
    sameGeometry(geometry,f.geometry());
    if(!boundary||&f.boundary()!=boundary)
      throw std::logic_error("Moving-window boundary is not installed on this grid");
    for(unsigned a=0;a<9;++a) {
      const auto& c=component(f,a);const auto& s=shapes[a];
      if(c.nx()!=s[0]||c.ny()!=s[1]||c.nz()!=s[2])
        throw std::logic_error("Moving-window native field shape changed");
    }
    for(unsigned h=0;h<12;++h)
      if(boundary->checkpointHistory(h).size()!=histories[h].values)
        throw std::logic_error("Moving-window transverse history layout mismatch");
  }
  long double planeEnergy(const EBFieldGrid& f,std::size_t k) {
    const long double volume=static_cast<long double>(geometry.dx)*
      static_cast<long double>(geometry.dy)*static_cast<long double>(geometry.dz);
    long double value=0.;
    for(std::size_t j=parameters.cells[1];j<geometry.ny-parameters.cells[1];++j)
      for(std::size_t i=parameters.cells[0];i<geometry.nx-parameters.cells[0];++i) {
        value+=static_cast<long double>(f.fieldEnergyDensityCell(i,j,k))*volume;
        ++energy.diagnosticCellSamples;
      }
    if(!std::isfinite(value))throw std::runtime_error("Non-finite moving-window energy diagnostic");
    return value;
  }
  std::array<long double,6> boundaryPower(const EBFieldGrid& f,bool inletOnly=false) {
    std::array<long double,6> power{};
    const std::size_t i0=parameters.cells[0],i1=geometry.nx-i0;
    const std::size_t j0=parameters.cells[1],j1=geometry.ny-j0;
    const long double yz=static_cast<long double>(geometry.dy)*geometry.dz;
    const long double xz=static_cast<long double>(geometry.dx)*geometry.dz;
    const long double xy=static_cast<long double>(geometry.dx)*geometry.dy;
    for(std::size_t k=0;k<(inletOnly?1:geometry.nz);++k) {
      for(std::size_t j=j0;j<j1;++j) {
        power[0]-=static_cast<long double>(f.poyntingComponentCell(i0,j,k,0))*yz;
        power[1]+=static_cast<long double>(f.poyntingComponentCell(i1-1,j,k,0))*yz;
      }
      for(std::size_t i=i0;i<i1;++i) {
        power[2]-=static_cast<long double>(f.poyntingComponentCell(i,j0,k,1))*xz;
        power[3]+=static_cast<long double>(f.poyntingComponentCell(i,j1-1,k,1))*xz;
      }
    }
    if(rank==0||rank+1==size)
      for(std::size_t j=j0;j<j1;++j)for(std::size_t i=i0;i<i1;++i) {
        if(rank==0)power[4]-=static_cast<long double>(f.poyntingComponentCell(i,j,0,2))*xy;
        if(!inletOnly&&rank+1==size)power[5]+=static_cast<long double>(f.poyntingComponentCell(i,j,geometry.nz-1,2))*xy;
      }
    for(const auto value:power)
      if(!std::isfinite(value))throw std::runtime_error("Non-finite moving-window boundary power");
    return power;
  }
  void finishEnergyLedgerStep(const EBFieldGrid& f) {
    if(!energy.enabled)return;
    validate(f);
    if(boundary->phase!=Boundary::Phase::Complete||boundary->energyStepFinished)
      throw std::logic_error("Moving-window energy flux requires one completed MPI field step");
    const auto diagnosticStart=Clock::now();
    const auto endPower=boundaryPower(f);
    const long double halfDt=.5L*static_cast<long double>(geometry.dt);
    for(unsigned face=0;face<6;++face)
      energy.outwardEnergyFace[face]+=halfDt*(boundary->startPower[face]+
        endPower[face]-boundary->overwritePowerDelta[face]);
    ++energy.fluxSteps;boundary->energyStepFinished=true;
    energy.diagnosticSeconds+=seconds(diagnosticStart);
  }
  MovingWindowFieldShiftReport shift(EBFieldGrid& f,std::size_t n,std::ptrdiff_t k,bool reset) {
    const auto started=Clock::now();validate(f);
    if(boundary->phase!=Boundary::Phase::Complete||n!=boundary->step+1||k!=boundary->offset-1)
      throw std::logic_error("Moving-window requires one left shift after each complete step");
    if(energy.enabled&&!boundary->energyStepFinished)
      throw std::logic_error("Moving-window energy flux must finish after MPI and before shift");
    long double removedLocal=0.;
    if(energy.enabled) {
      const auto diagnosticStart=Clock::now();
      removedLocal=planeEnergy(f,geometry.nz-1);
      energy.diagnosticSeconds+=seconds(diagnosticStart);
    }
    const auto& p=reset?resetPacket:preservePacket;
    const std::size_t oldK=geometry.nz-1;
    for(unsigned a=0;a<p.fields;++a) {
      const auto& c=component(f,a);const auto plane=multiply(c.nx(),c.ny());
      std::memcpy(send.data()+p.fieldOffsets[a],c.data()+oldK*plane,plane*sizeof(Double));
    }
    for(unsigned h=0;h<12;++h) {
      const auto& v=histories[h];const auto& values=boundary->checkpointHistory(h);
      for(std::size_t slot=0;slot<v.slots;++slot)
        std::memcpy(send.data()+p.historyOffsets[h]+slot*v.row,
          values.data()+(slot*v.layers+oldK)*v.row,v.row*sizeof(Double));
    }
    // Validate the entire prescribed input before modifying any field/history.
    if(rank==0) {
      std::fill(receive.begin(),receive.begin()+p.values,0.);
      for(unsigned a=0;a<std::min(p.fields,6u);++a) {
        const auto& c=component(f,a);
        for(std::size_t j=0;j<c.ny();++j)for(std::size_t i=0;i<c.nx();++i)
          receive[p.fieldOffsets[a]+j*c.nx()+i]=boundary->reference(a,i,j,k,n);
      }
    }
    MovingWindowFieldShiftReport report;
    if(size>1) {
      const auto commStart=Clock::now();
      const int status=MPI_Sendrecv(send.data(),static_cast<int>(p.values),MPI_DOUBLE,upperRank,18331,
        receive.data(),static_cast<int>(p.values),MPI_DOUBLE,lowerRank,18331,communicator,MPI_STATUS_IGNORE);
      report.communicationSeconds=seconds(commStart);
      if(status!=MPI_SUCCESS){MPI_Abort(communicator,status);throw std::runtime_error("Moving-window neighbour transfer failed");}
      if(upperRank!=MPI_PROC_NULL)report.sentBytes=p.values*sizeof(Double);
      if(lowerRank!=MPI_PROC_NULL)report.receivedBytes=p.values*sizeof(Double);
    }
    for(unsigned a=0;a<p.fields;++a) {
      auto& c=component(f,a);const auto plane=multiply(c.nx(),c.ny());
      const auto bytes=(c.nz()-1)*plane*sizeof(Double);
      std::memmove(c.data()+plane,c.data(),bytes);
      std::memcpy(c.data(),receive.data()+p.fieldOffsets[a],plane*sizeof(Double));
      report.fieldMovedBytes+=bytes;
    }
    if(reset)for(unsigned a=6;a<9;++a){auto& c=component(f,a);c.fill(0.);report.currentClearedBytes+=c.size()*sizeof(Double);}
    for(unsigned h=0;h<12;++h) {
      const auto& v=histories[h];auto& values=boundary->checkpointHistory(h);
      for(std::size_t slot=0;slot<v.slots;++slot) {
        Double* base=values.data()+slot*v.layers*v.row;const auto bytes=(v.layers-1)*v.row*sizeof(Double);
        std::memmove(base+v.row,base,bytes);
        std::memcpy(base,receive.data()+p.historyOffsets[h]+slot*v.row,v.row*sizeof(Double));
        report.historyMovedBytes+=bytes;
      }
    }
    if(energy.enabled) {
      const auto diagnosticStart=Clock::now();
      const long double insertedLocal=planeEnergy(f,0);
      // Every retained local cell k>=1 is exactly old cell k-1, including
      // its nz+1 Yee nodes. Thus this two-plane difference measures the full
      // local remap energy change without traversing the volume. Across ranks
      // the entering/leaving cell energies cancel only if shared nodes match;
      // keep that interface contribution instead of assuming cancellation.
      energy.remapEnergyJ+=insertedLocal-removedLocal;
      if(rank==0)energy.insertedEnergyJ+=insertedLocal;
      if(rank+1==size)energy.croppedEnergyJ+=removedLocal;
      energy.remapInterfaceCorrectionJ+=(rank>0?insertedLocal:0.L)-
        (rank+1<size?removedLocal:0.L);
      ++energy.remapSteps;
      energy.diagnosticSeconds+=seconds(diagnosticStart);
    }
    boundary->step=n;boundary->offset=k;boundary->phase=Boundary::Phase::Idle;
    report.totalSeconds=seconds(started);return report;
  }
};

MovingWindowFields::MovingWindowFields(const EBGridGeometry& g,std::size_t nz,
    std::size_t local,MPI_Comm comm,const EBCPMLParameters& p,MovingWindowReference r)
  :impl_(new Impl(g,nz,local,comm,p,r)) {}
MovingWindowFields::~MovingWindowFields() {}
void MovingWindowFields::installBoundary(EBFieldGrid& f) {
  sameGeometry(impl_->geometry,f.geometry());
  if(impl_->boundary)throw std::logic_error("Moving-window boundary already installed");
  std::unique_ptr<Impl::Boundary> p(new Impl::Boundary(*impl_));
  impl_->boundary=p.get();f.setBoundary(std::move(p));impl_->validate(f);
}
void MovingWindowFields::enableEnergyLedger(bool enabled) {
  if(impl_->boundary&&impl_->boundary->established)
    throw std::logic_error("Moving-window energy ledger must be configured before the first step");
  impl_->energy.enabled=enabled;
}
void MovingWindowFields::setBeforeInletObserver(
    BeforeInletObserver observer,void* context) {
  impl_->beforeInletObserver=observer;
  impl_->beforeInletContext=context;
}
void MovingWindowFields::setStageObserver(StageObserver observer,void* context) {
  impl_->stageObserver=observer;impl_->stageContext=context;
}
const MovingWindowFieldEnergyLedger& MovingWindowFields::energyLedger() const {
  return impl_->energy;
}
void MovingWindowFields::finishEnergyLedgerStep(const EBFieldGrid& f) {
  impl_->finishEnergyLedgerStep(f);
}
void MovingWindowFields::beginStep(EBFieldGrid& f,std::size_t n,std::ptrdiff_t k) {
  impl_->validate(f);impl_->boundary->begin(f,n,k);
}
EBElectricHaloView MovingWindowFields::physicalElectricHalos(const EBFieldGrid& f,std::size_t n) const {
  impl_->validate(f);return impl_->boundary->halos(n);
}
MovingWindowFieldShiftReport MovingWindowFields::shift(EBFieldGrid& f,std::size_t n,std::ptrdiff_t k,bool reset) {
  return impl_->shift(f,n,k,reset);
}
std::size_t MovingWindowFields::memoryBytes() const {
  return sizeof(Impl)+(impl_->send.capacity()+impl_->receive.capacity())*sizeof(Double);
}
std::uint64_t MovingWindowFields::referenceValues() const {
  return impl_->boundary?impl_->boundary->queries:0;
}
}
