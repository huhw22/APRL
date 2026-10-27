#include "streaming_cpml_boundary.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace aprl {
namespace verification {

struct ReferenceFedCPML::ResidualInputWork {
  using Clock = std::chrono::steady_clock;
  struct Patch {
    unsigned component;
    std::size_t firstK, planes, planeValues, backupOffset;
  };
  enum class Phase { Idle, Magnetic, Electric };

  const StreamingReferenceCache& reference;
  const std::size_t& offset;
  std::size_t stepOffset;
  unsigned int step = 0;
  Phase phase = Phase::Idle;
  std::array<Patch, 4> electric, magnetic;
  std::vector<Double> backup;
  std::uint64_t inputValues = 0;
  double inputSeconds = 0.;

  static std::size_t product(std::size_t a, std::size_t b) {
    if (b && a > std::numeric_limits<std::size_t>::max()/b)
      throw std::overflow_error("Residual CPML backup size overflow");
    return a*b;
  }

  ResidualInputWork(const EBGridGeometry& g, const EBCPMLParameters& p,
                    const StreamingReferenceCache& r, const std::size_t& o,
                    bool zeroFeed)
    : reference(r), offset(o), stepOffset(o) {
    if (zeroFeed || p.cells[0] || p.cells[1] || !p.cells[2] ||
        p.cells[2] >= g.nz || p.cells[2] >= g.nz-p.cells[2])
      throw std::invalid_argument(
        "Residual CPML input requires nonoverlapping z-only CPML and reference feed");
    if (g.solver != EBMaxwellSolver::CowanZ)
      throw std::invalid_argument("Residual CPML input requires the existing Cowan probe");
    const std::size_t pz = p.cells[2];
    std::size_t eValues = 0, bValues = 0;
    for (unsigned component = 0; component < 2; ++component) {
      const std::size_t ePlane = product(component == 0 ? g.nx : g.nx+1,
                                         component == 0 ? g.ny+1 : g.ny);
      const std::size_t bPlane = product(component == 0 ? g.nx+1 : g.nx,
                                         component == 0 ? g.ny : g.ny+1);
      for (unsigned side = 0; side < 2; ++side) {
        const unsigned patch = 2*component+side;
        electric[patch] = {component, side ? g.nz-pz : 0,
                           pz+1, ePlane, eValues};
        magnetic[patch] = {component+3, side ? g.nz-pz : 0,
                           pz, bPlane, bValues};
        const std::size_t eCount = product(pz+1, ePlane);
        const std::size_t bCount = product(pz, bPlane);
        if (eValues > std::numeric_limits<std::size_t>::max()-eCount ||
            bValues > std::numeric_limits<std::size_t>::max()-bCount)
          throw std::overflow_error("Residual CPML backup sum overflow");
        eValues += eCount; bValues += bCount;
      }
    }
    const std::size_t values = std::max(eValues, bValues);
    product(values, sizeof(Double));
    backup.resize(values); // One allocation, reused by both native time layers.
  }

  static YeeComponent& source(EBFieldGrid& fields, unsigned component) {
    if (component == 0) return fields.ex();
    if (component == 1) return fields.ey();
    if (component == 3) return fields.bx();
    return fields.by();
  }

  void requireIdle() const {
    if (phase != Phase::Idle)
      throw std::logic_error("Residual CPML step was begun before its prior electric hook");
  }

  void begin(unsigned int nativeStep) {
    step = nativeStep; stepOffset = offset; phase = Phase::Magnetic;
  }

  void requirePhase(Phase expected) const {
    if (phase != expected || offset != stepOffset)
      throw std::logic_error("Residual CPML requires ordered B/E hooks on an unshifted step lattice");
  }

  void restore(EBFieldGrid& fields, const std::array<Patch, 4>& patches) noexcept {
    const auto start = Clock::now();
    for (const Patch& p : patches)
      std::memcpy(source(fields, p.component).data()+p.firstK*p.planeValues,
                  backup.data()+p.backupOffset, p.planes*p.planeValues*sizeof(Double));
    inputSeconds += std::chrono::duration<double>(Clock::now()-start).count();
  }

  template<typename Hook>
  void withResidualSource(EBFieldGrid& fields, bool magneticSource, Hook hook) {
    requirePhase(magneticSource ? Phase::Electric : Phase::Magnetic);
    const auto& patches = magneticSource ? magnetic : electric;
    // Copy the exact stored values first. Restoration uses memcpy, never
    // addition of a large reference, so source arrays recover all bits.
    auto start = Clock::now();
    for (const Patch& p : patches) {
      const YeeComponent& c = source(fields, p.component);
      if (c.nx()*c.ny() != p.planeValues || p.firstK+p.planes > c.nz())
        throw std::logic_error("Residual CPML source geometry differs from construction");
      std::memcpy(backup.data()+p.backupOffset, c.data()+p.firstK*p.planeValues,
                  p.planes*p.planeValues*sizeof(Double));
    }
    // The backup is complete before any source mutation. Even a failed
    // cache query or native hook restores sources. Targets/history are not
    // rolled back on failure; such a failed simulation must terminate.
    try {
      const Double nativeState = static_cast<Double>(step)+(magneticSource ? 1. : 0.);
      for (const Patch& p : patches) {
        YeeComponent& c = source(fields, p.component);
        for (std::size_t k = p.firstK; k < p.firstK+p.planes; ++k)
          for (std::size_t j = 0; j < c.ny(); ++j)
            for (std::size_t i = 0; i < c.nx(); ++i) {
              c(i,j,k) -= reference.value(p.component, i, j,
                  static_cast<Double>(offset)+static_cast<Double>(k), nativeState);
              ++inputValues;
            }
      }
    } catch (...) {
      inputSeconds += std::chrono::duration<double>(Clock::now()-start).count();
      restore(fields, patches);
      throw;
    }
    inputSeconds += std::chrono::duration<double>(Clock::now()-start).count();
    try {
      hook(); // Original CPML arithmetic; deliberately outside input timing.
    } catch (...) {
      restore(fields, patches);
      throw;
    }
    restore(fields, patches);
  }

  std::size_t memoryBytes() const {
    return sizeof(*this)+backup.capacity()*sizeof(Double);
  }
};

ReferenceFedCPML::ReferenceFedCPML(
    const EBGridGeometry& geometry, const EBCPMLParameters& parameters,
    StreamingReferenceCache& reference, const std::size_t& offset,
    bool zeroFeed, bool residualInput)
  : EBConvolutionalPML(geometry, geometry.nz, 0, parameters, true, true),
    feed_(reference, offset, zeroFeed ? StreamingBoundaryMode::ZeroFeed :
                                     StreamingBoundaryMode::ReferenceClamp),
    residual_(residualInput ? new ResidualInputWork(
        geometry, parameters, reference, offset, zeroFeed) : nullptr) {}

ReferenceFedCPML::~ReferenceFedCPML() = default;

void ReferenceFedCPML::beginStep(const EBFieldGrid& fields, unsigned int step) {
  if (residual_) residual_->requireIdle();
  feed_.beginStep(fields, step);
  if (residual_) residual_->begin(step);
}

EBElectricHaloView ReferenceFedCPML::electricHalos(
    const EBFieldGrid& fields, unsigned int step) {
  return feed_.electricHalos(fields, step);
}

void ReferenceFedCPML::afterMagneticUpdate(
    EBFieldGrid& fields, const EBElectricHaloView& halo) {
  if (residual_) {
    // Only z derivative profiles are enabled. Cowan's z derivative smooths
    // in x/y, so CPML reads Ex/Ey at k=0..P and W-P..W, never a z halo.
    // The preceding total-field Maxwell update still uses the caller's
    // original reference halo; no halo or field time convention is changed.
    residual_->withResidualSource(fields, false, [&]() {
      EBConvolutionalPML::afterMagneticUpdate(fields, halo);
    });
  } else EBConvolutionalPML::afterMagneticUpdate(fields, halo);
  feed_.afterMagneticUpdate(fields, halo);
  if (residual_) residual_->phase = ResidualInputWork::Phase::Electric;
}

void ReferenceFedCPML::afterElectricUpdate(EBFieldGrid& fields) {
  if (residual_) {
    // Native B is now at n+1/2: query cache stateStep=n+1. Original Ampere
    // already used total B and true J; only its CPML correction reads this
    // temporary Bx/By residual. PEC termination is still the original hook.
    residual_->withResidualSource(fields, true, [&]() {
      EBConvolutionalPML::afterElectricUpdate(fields);
    });
  } else EBConvolutionalPML::afterElectricUpdate(fields);
  feed_.afterElectricUpdate(fields);
  if (residual_) residual_->phase = ResidualInputWork::Phase::Idle;
}

std::size_t ReferenceFedCPML::memoryBytes() const {
  return EBConvolutionalPML::memoryBytes() + feed_.memoryBytes()+residualWorkBytes();
}

std::uint64_t ReferenceFedCPML::boundaryValues() const {
  return feed_.boundaryValues();
}

std::uint64_t ReferenceFedCPML::haloValues() const {
  return feed_.haloValues();
}

bool ReferenceFedCPML::residualInput() const { return static_cast<bool>(residual_); }
std::size_t ReferenceFedCPML::residualWorkBytes() const {
  return residual_ ? residual_->memoryBytes() : 0;
}
std::uint64_t ReferenceFedCPML::residualInputValues() const {
  return residual_ ? residual_->inputValues : 0;
}
double ReferenceFedCPML::residualInputSeconds() const {
  return residual_ ? residual_->inputSeconds : 0.;
}

} // namespace verification
} // namespace aprl
