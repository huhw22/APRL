#include "streaming_reference_boundary.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace aprl {
namespace verification {
namespace {

const YeeComponent& fieldComponent(const EBFieldGrid& fields,unsigned int a) {
  switch(a) {
    case 0:return fields.ex();case 1:return fields.ey();case 2:return fields.ez();
    case 3:return fields.bx();case 4:return fields.by();case 5:return fields.bz();
  }
  throw std::out_of_range("Invalid streaming reference component");
}

YeeComponent& electricComponent(EBFieldGrid& fields,unsigned int a) {
  if(a==0)return fields.ex();
  if(a==1)return fields.ey();
  return fields.ez();
}

}  // namespace

StreamingReferenceCache::StreamingReferenceCache(
    const EBFieldGrid& initial,Double beta)
  :geometry_(initial.geometry()),beta_(beta),samples_(0),
   diagnosticQuery_(nullptr),diagnosticContext_(nullptr) {
  if(geometry_.solver!=EBMaxwellSolver::CowanZ ||
     !std::isfinite(beta_) || std::abs(beta_)>=1.)
    throw std::invalid_argument("Streaming reference requires Cowan and |beta|<1");
  for(unsigned int a=0;a<6;++a) {
    const YeeComponent& source=fieldComponent(initial,a);
    Component& target=components_[a];
    target.nx=source.nx();target.ny=source.ny();target.nz=source.nz();
    target.values.assign(source.data(),source.data()+source.size());
    for(Double value:target.values)
      if(!std::isfinite(value))
        throw std::invalid_argument("Non-finite initial streaming reference field");
  }
}

Double StreamingReferenceCache::value(unsigned int a,std::size_t i,
    std::size_t j,Double globalK,Double stateStep) const {
  if(diagnosticQuery_) {
    const Double result=diagnosticQuery_(diagnosticContext_,a,i,j,globalK,stateStep);
    ++samples_;
    return result;
  }
  if(a>=6 || !std::isfinite(globalK) || !std::isfinite(stateStep))
    throw std::out_of_range("Invalid streaming reference query");
  const Component& component=components_[a];
  if(i>=component.nx || j>=component.ny)
    throw std::out_of_range("Streaming reference query outside transverse mesh");
  Double coordinate=globalK-beta_*SI::c*geometry_.dt/geometry_.dz*stateStep;
  const Double maximum=static_cast<Double>(component.nz-1);
  const Double tolerance=64*std::numeric_limits<Double>::epsilon()*
    std::max(1.,std::max(std::abs(coordinate),maximum));
  if(coordinate < -tolerance || coordinate > maximum+tolerance) {
    std::ostringstream message;
    message.precision(17);
    message<<"Streaming reference query outside frozen cache: component="<<a
      <<", requested_k="<<globalK<<", state_step="<<stateStep
      <<", translated_k="<<coordinate<<", last_index="<<maximum;
    throw std::out_of_range(message.str());
  }
  coordinate=std::max(0.,std::min(maximum,coordinate));
  const std::size_t lower=static_cast<std::size_t>(std::floor(coordinate));
  const std::size_t upper=std::min(lower+1,component.nz-1);
  const Double fraction=coordinate-static_cast<Double>(lower);
  const std::size_t plane=component.nx*component.ny;
  const std::size_t transverse=j*component.nx+i;
  ++samples_;
  // Preserve exact cached values for lattice-aligned queries.
  if(fraction==0.)return component.values[lower*plane+transverse];
  return (1.-fraction)*component.values[lower*plane+transverse]+
    fraction*component.values[upper*plane+transverse];
}

void StreamingReferenceCache::setDiagnosticQuery(
    DiagnosticQuery query,void* context) {
  diagnosticQuery_=query;
  diagnosticContext_=query?context:nullptr;
}

const EBGridGeometry& StreamingReferenceCache::geometry() const {return geometry_;}
Double StreamingReferenceCache::beta() const {return beta_;}
std::size_t StreamingReferenceCache::memoryBytes() const {
  std::size_t bytes=0;
  for(const auto& component:components_)bytes+=component.values.capacity()*sizeof(Double);
  return bytes;
}
std::uint64_t StreamingReferenceCache::sampleCount() const {return samples_;}

StreamingReferenceBoundary::StreamingReferenceBoundary(
    const StreamingReferenceCache& reference,const std::size_t& offset,
    StreamingBoundaryMode mode)
  :reference_(reference),offset_(offset),mode_(mode),transversePec_(false,false),
   step_(0),stepOffset_(offset),prepared_(false),boundaryValues_(0),haloValues_(0) {}

void StreamingReferenceBoundary::validateGeometry(const EBFieldGrid& fields) const {
  const auto& actual=fields.geometry();const auto& expected=reference_.geometry();
  if(actual.nx!=expected.nx || actual.ny!=expected.ny ||
     actual.dx!=expected.dx || actual.dy!=expected.dy ||
     actual.dz!=expected.dz || actual.dt!=expected.dt ||
     actual.solver!=expected.solver || actual.nz<2)
    throw std::invalid_argument("Streaming window/reference geometry mismatch");
}

void StreamingReferenceBoundary::requirePrepared() const {
  if(!prepared_ || offset_!=stepOffset_)
    throw std::logic_error("Streaming boundary needs beginStep and an unshifted step lattice");
}

Double StreamingReferenceBoundary::referenceValue(unsigned int a,std::size_t i,
    std::size_t j,Double k,Double step) const {
  return mode_==StreamingBoundaryMode::ZeroFeed?0.:reference_.value(a,i,j,k,step);
}

void StreamingReferenceBoundary::beginStep(const EBFieldGrid& fields,unsigned int step) {
  validateGeometry(fields);step_=step;stepOffset_=offset_;prepared_=true;
  if(mode_!=StreamingBoundaryMode::ResidualMur)return;
  for(unsigned int a=0;a<2;++a) {
    const YeeComponent& component=fieldComponent(fields,a);
    const std::size_t plane=component.nx()*component.ny();
    for(unsigned int side=0;side<2;++side) {
      oldEnd_[a][side].resize(plane);oldAdjacent_[a][side].resize(plane);
      const std::size_t end=side?component.nz()-1:0;
      const std::size_t adjacent=side?end-1:1;
      for(std::size_t j=0;j<component.ny();++j)
        for(std::size_t i=0;i<component.nx();++i) {
          const std::size_t index=j*component.nx()+i;
          oldEnd_[a][side][index]=component(i,j,end)-
            reference_.value(a,i,j,static_cast<Double>(offset_+end),step_);
          oldAdjacent_[a][side][index]=component(i,j,adjacent)-
            reference_.value(a,i,j,static_cast<Double>(offset_+adjacent),step_);
        }
    }
  }
}

EBElectricHaloView StreamingReferenceBoundary::electricHalos(
    const EBFieldGrid& fields,unsigned int step) {
  requirePrepared();validateGeometry(fields);
  if(step!=step_)throw std::logic_error("Streaming halo time differs from beginStep");
  for(unsigned int a=0;a<3;++a) {
    const YeeComponent& component=fieldComponent(fields,a);
    const std::size_t plane=component.nx()*component.ny();
    lower_[a].resize(plane);upper_[a].resize(plane);
    for(unsigned int side=0;side<2;++side) {
      std::vector<Double>& halo=side?upper_[a]:lower_[a];
      const std::size_t end=side?component.nz()-1:0;
      const std::size_t adjacent=side?end-1:1;
      const Double k=static_cast<Double>(offset_)+(side?
        static_cast<Double>(component.nz()):-1.);
      for(std::size_t j=0;j<component.ny();++j)
        for(std::size_t i=0;i<component.nx();++i) {
          Double value=referenceValue(a,i,j,k,step_);
          if(mode_==StreamingBoundaryMode::ResidualMur) {
            const Double residualEnd=component(i,j,end)-
              reference_.value(a,i,j,static_cast<Double>(offset_+end),step_);
            const Double residualAdjacent=component(i,j,adjacent)-
              reference_.value(a,i,j,static_cast<Double>(offset_+adjacent),step_);
            value+=2.*residualEnd-residualAdjacent;
          }
          halo[j*component.nx()+i]=value;++haloValues_;
        }
    }
  }
  EBElectricHaloView halo;
  halo.lowerEx=lower_[0].data();halo.lowerEy=lower_[1].data();halo.lowerEz=lower_[2].data();
  halo.upperEx=upper_[0].data();halo.upperEy=upper_[1].data();halo.upperEz=upper_[2].data();
  return halo;
}

void StreamingReferenceBoundary::afterMagneticUpdate(
    EBFieldGrid& fields,const EBElectricHaloView& halo) {
  requirePrepared();transversePec_.afterMagneticUpdate(fields,halo);
}

void StreamingReferenceBoundary::afterElectricUpdate(EBFieldGrid& fields) {
  requirePrepared();
  const auto& geometry=fields.geometry();
  const Double velocity=1./std::sqrt(SI::mu0*SI::epsilon0);
  const Double ratio=(velocity*geometry.dt-geometry.dz)/
    (velocity*geometry.dt+geometry.dz);
  for(unsigned int a=0;a<2;++a) {
    YeeComponent& component=electricComponent(fields,a);
    for(unsigned int side=0;side<2;++side) {
      const std::size_t end=side?component.nz()-1:0;
      const std::size_t adjacent=side?end-1:1;
      for(std::size_t j=0;j<component.ny();++j)
        for(std::size_t i=0;i<component.nx();++i) {
          Double value=referenceValue(a,i,j,static_cast<Double>(offset_+end),step_+1.);
          if(mode_==StreamingBoundaryMode::ResidualMur) {
            const std::size_t index=j*component.nx()+i;
            const Double nextAdjacent=component(i,j,adjacent)-
              reference_.value(a,i,j,static_cast<Double>(offset_+adjacent),step_+1.);
            value+=oldAdjacent_[a][side][index]+
              ratio*(nextAdjacent-oldEnd_[a][side][index]);
          }
          component(i,j,end)=value;++boundaryValues_;
        }
    }
  }
  transversePec_.afterElectricUpdate(fields);
  prepared_=false;
}

std::size_t StreamingReferenceBoundary::memoryBytes() const {
  std::size_t bytes=0;
  for(unsigned int a=0;a<2;++a)for(unsigned int side=0;side<2;++side)
    bytes+=(oldEnd_[a][side].capacity()+oldAdjacent_[a][side].capacity())*sizeof(Double);
  for(unsigned int a=0;a<3;++a)
    bytes+=(lower_[a].capacity()+upper_[a].capacity())*sizeof(Double);
  return bytes;
}

void StreamingReferenceBoundary::refreshAfterRemap(
    EBFieldGrid& fields, unsigned int stateStep) {
  validateGeometry(fields);
  if(prepared_ || mode_==StreamingBoundaryMode::ResidualMur)
    throw std::logic_error("Remapped face refresh requires a completed clamp step");
  for(unsigned int a=0;a<2;++a) {
    YeeComponent& c=electricComponent(fields,a);
    for(unsigned side=0;side<2;++side) {
      const std::size_t k=side?c.nz()-1:0;
      for(std::size_t j=0;j<c.ny();++j)
        for(std::size_t i=0;i<c.nx();++i) {
          c(i,j,k)=referenceValue(a,i,j,static_cast<Double>(offset_+k),stateStep);
          ++boundaryValues_;
        }
    }
  }
  transversePec_.afterElectricUpdate(fields);
}
std::uint64_t StreamingReferenceBoundary::boundaryValues() const {return boundaryValues_;}
std::uint64_t StreamingReferenceBoundary::haloValues() const {return haloValues_;}

}  // namespace verification
}  // namespace aprl
