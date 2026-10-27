#include "held_cut_boundary.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace aprl { namespace verification {
namespace {
YeeComponent& component(EBFieldGrid& f,unsigned int a) {
  switch(a) {
    case 0:return f.ex();case 1:return f.ey();case 2:return f.ez();
    case 3:return f.bx();case 4:return f.by();case 5:return f.bz();
  }
  throw std::out_of_range("Held cut field component");
}
const YeeComponent& component(const EBFieldGrid& f,unsigned int a) {
  return component(const_cast<EBFieldGrid&>(f),a);
}
void compatible(const EBGridGeometry& a,const EBGridGeometry& b) {
  if(a.nx!=b.nx||a.ny!=b.ny||a.dx!=b.dx||a.dy!=b.dy||a.dz!=b.dz||
     a.dt!=b.dt||a.solver!=b.solver)
    throw std::invalid_argument("Held cut reference/grid mismatch");
}
std::size_t endPlane(std::size_t nz,std::size_t first,std::size_t count) {
  if(first>nz||count>nz-first)
    throw std::out_of_range("Held cut slab outside target");
  return count?first+count:nz;
}
}

HeldCutBoundary::HeldCutBoundary(const EBGridGeometry& g,
    const EBCPMLParameters& p,const std::size_t& offset,HeldCutReference reference,
    Double rightValue)
  :EBConvolutionalPML(g,g.nz,0,p,false,false),geometry_(g),offset_(offset),
   stepOffset_(offset),reference_(reference),rightValue_(rightValue) {
  if(g.solver!=EBMaxwellSolver::CowanZ||p.cells[2]!=0||
     p.cells[0]==0||p.cells[1]==0||!reference.field||!std::isfinite(rightValue))
    throw std::invalid_argument("Held cut requires Cowan, two transverse CPML axes and a native reference");
  if(offset>static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max())-g.nz)
    throw std::invalid_argument("Held cut global indices overflow");
  for(unsigned int a=0;a<3;++a) {
    const std::size_t nx=g.nx+(a!=0),ny=g.ny+(a!=1);
    lower_[a].resize(nx*ny);upper_[a].resize(nx*ny);
  }
}

void HeldCutBoundary::validate(const EBFieldGrid& f) const {
  compatible(geometry_,f.geometry());
  if(geometry_.nz!=f.geometry().nz)
    throw std::invalid_argument("Held cut window width changed");
}
void HeldCutBoundary::requirePhase(Phase p) const {
  if(phase_!=p||stepOffset_!=offset_)
    throw std::logic_error("Held cut requires ordered hooks on an unshifted lattice");
}
void HeldCutBoundary::beginStep(const EBFieldGrid& fields,unsigned int n) {
  if(phase_!=Phase::Idle)
    throw std::logic_error("Held cut previous step was not completed");
  validate(fields);step_=n;stepOffset_=offset_;
  for(unsigned int a=0;a<3;++a) {
    const auto& c=component(fields,a);
    for(std::size_t j=0;j<c.ny();++j)for(std::size_t i=0;i<c.nx();++i) {
      const std::size_t ij=j*c.nx()+i;
      lower_[a][ij]=reference_.field(reference_.context,a,i,j,
          static_cast<std::ptrdiff_t>(offset_)-1,n);
      if(!std::isfinite(lower_[a][ij]))
        throw std::runtime_error("Non-finite held inlet halo");
      upper_[a][ij]=rightValue_;
      ++queries_;
    }
  }
  phase_=Phase::Magnetic;
}
EBElectricHaloView HeldCutBoundary::electricHalos(
    const EBFieldGrid& fields,unsigned int n) const {
  validate(fields);requirePhase(Phase::Magnetic);
  if(n!=step_)throw std::logic_error("Held cut halo has wrong native time");
  EBElectricHaloView h;
  h.lowerEx=lower_[0].data();h.lowerEy=lower_[1].data();h.lowerEz=lower_[2].data();
  h.upperEx=upper_[0].data();h.upperEy=upper_[1].data();h.upperEz=upper_[2].data();
  return h;
}
void HeldCutBoundary::afterMagneticUpdate(
    EBFieldGrid& f,const EBElectricHaloView& h) {
  validate(f);requirePhase(Phase::Magnetic);
  EBConvolutionalPML::afterMagneticUpdate(f,h);
  phase_=Phase::Electric;
}
void HeldCutBoundary::afterElectricUpdate(EBFieldGrid& f) {
  validate(f);requirePhase(Phase::Electric);
  EBConvolutionalPML::afterElectricUpdate(f);
  // The native Ampere loop intentionally excludes Ex/Ey z endpoints. Supply
  // only the lower prescribed endpoint; the upper is a sacrificial face.
  for(unsigned int a=0;a<2;++a) {
    auto& c=component(f,a);
    for(std::size_t j=0;j<c.ny();++j)for(std::size_t i=0;i<c.nx();++i) {
      const Double value=reference_.field(reference_.context,a,i,j,
          static_cast<std::ptrdiff_t>(offset_),step_+1);
      if(!std::isfinite(value))throw std::runtime_error("Non-finite held inlet face");
      c(i,j,0)=value;++queries_;
    }
  }
  phase_=Phase::Idle;
}
std::size_t HeldCutBoundary::memoryBytes() const {
  std::size_t bytes=EBConvolutionalPML::memoryBytes();
  for(unsigned int a=0;a<3;++a)
    bytes+=(lower_[a].capacity()+upper_[a].capacity())*sizeof(Double);
  return bytes;
}
std::uint64_t HeldCutBoundary::referenceValues() const {return queries_;}
void HeldCutBoundary::setRightHaloValue(Double value) {
  if(phase_!=Phase::Idle||!std::isfinite(value))
    throw std::invalid_argument("Change finite right halo only between complete steps");
  rightValue_=value;
}

void copyHeldCutFields(EBFieldGrid& target,const EBFieldGrid& source,
    std::size_t offset,std::size_t first,std::size_t count) {
  compatible(target.geometry(),source.geometry());
  // Validate the complete operation before changing a target component.
  for(unsigned int a=0;a<6;++a) {
    const auto& t=component(target,a);const auto& s=component(source,a);
    const std::size_t end=endPlane(t.nz(),first,count);
    if(offset>s.nz()||end>s.nz()-offset)
      throw std::out_of_range("Held cut field slab outside donor");
  }
  for(unsigned int a=0;a<6;++a) {
    auto& t=component(target,a);const auto& s=component(source,a);
    const std::size_t end=endPlane(t.nz(),first,count),plane=t.nx()*t.ny();
    std::memcpy(t.data()+first*plane,s.data()+(offset+first)*plane,
                (end-first)*plane*sizeof(Double));
  }
}

void copyHeldCutHistories(EBFieldGrid& target,const MovingWindowRemap& tLayout,
    const EBFieldGrid& source,const MovingWindowRemap& sLayout,
    std::size_t offset,std::size_t first,std::size_t count) {
  compatible(target.geometry(),source.geometry());
  auto* t=dynamic_cast<EBConvolutionalPML*>(&target.boundary());
  const auto* s=dynamic_cast<const EBConvolutionalPML*>(&source.boundary());
  if(!t||!s)throw std::invalid_argument("Held cut histories require native CPML");
  // Validate every history before copying any of them.
  for(std::size_t h=0;h<12;++h) {
    const auto& a=tLayout.historyLayout(h);const auto& b=sLayout.historyLayout(h);
    if(a.fieldComponent!=b.fieldComponent||a.derivativeAxis!=b.derivativeAxis||
       a.activeAxisIndices!=b.activeAxisIndices||a.dimensions[0]!=b.dimensions[0]||
       a.dimensions[1]!=b.dimensions[1]||a.values!=t->checkpointHistory(h).size()||
       b.values!=s->checkpointHistory(h).size())
      throw std::invalid_argument("Held cut transverse history layouts differ");
    if(a.derivativeAxis==2) {
      if(a.values||b.values)throw std::invalid_argument("Held cut cannot import z CPML histories");
      continue;
    }
    const std::size_t end=endPlane(a.dimensions[2],first,count);
    if(offset>b.dimensions[2]||end>b.dimensions[2]-offset)
      throw std::out_of_range("Held cut history slab outside donor");
  }
  for(std::size_t h=0;h<12;++h) {
    const auto& a=tLayout.historyLayout(h);const auto& b=sLayout.historyLayout(h);
    if(a.derivativeAxis==2)continue;
    const std::size_t q=a.dimensions[a.derivativeAxis==0?1:0];
    const std::size_t end=endPlane(a.dimensions[2],first,count);
    auto& tv=t->checkpointHistory(h);const auto& sv=s->checkpointHistory(h);
    for(std::size_t slot=0;slot<a.activeAxisIndices.size();++slot) {
      const std::size_t ti=(slot*a.dimensions[2]+first)*q;
      const std::size_t si=(slot*b.dimensions[2]+offset+first)*q;
      std::memcpy(tv.data()+ti,sv.data()+si,(end-first)*q*sizeof(Double));
    }
  }
}

}} // namespace aprl::verification
