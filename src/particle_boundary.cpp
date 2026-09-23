#include "particle_boundary.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fel
{
  namespace
  {
    Double coordinateTolerance(Double lower, Double upper)
    {
      return 64.0 * std::numeric_limits<Double>::epsilon() *
        std::max(1.0, std::max(std::abs(lower), std::abs(upper)));
    }

    void cellCoordinate(Double coordinate, std::size_t cells,
                        std::size_t& cell, Double& fraction)
    {
      const Double upper = static_cast<Double>(cells);
      const Double clamped = std::max(0.0, std::min(upper, coordinate));
      if (clamped >= upper)
        {
          cell = cells - 1;
          fraction = 1.0;
        }
      else
        {
          cell = static_cast<std::size_t>(std::floor(clamped));
          fraction = clamped - static_cast<Double>(cell);
        }
    }
  }

  const char* particleBoundaryFaceName(ParticleBoundaryFace face)
  {
    switch (face)
      {
      case ParticleBoundaryFace::LowerX: return "x-";
      case ParticleBoundaryFace::UpperX: return "x+";
      case ParticleBoundaryFace::LowerY: return "y-";
      case ParticleBoundaryFace::UpperY: return "y+";
      case ParticleBoundaryFace::LowerZ: return "z-";
      case ParticleBoundaryFace::UpperZ: return "z+";
      }
    throw std::invalid_argument("Unknown particle boundary face");
  }

  ParticleBoundaryHit::ParticleBoundaryHit()
    : position(0.0), fraction(1.0),
      face(ParticleBoundaryFace::LowerX)
  {}

  ParticleOpenBoundary::FaceFlux::FaceFlux()
    : firstNodes(0), secondNodes(0), values()
  {}

  ParticleOpenBoundary::ParticleOpenBoundary(
      const EBGridGeometry& localGeometry,
      const FieldVector<Double>& localOriginSI,
      const EBGridGeometry& globalGeometry,
      const FieldVector<Double>& globalOriginSI,
      bool retainSpatialFlux)
    : localGeometry_(localGeometry), localOriginSI_(localOriginSI),
      globalGeometry_(globalGeometry), globalOriginSI_(globalOriginSI),
      localLowerZSI_(localOriginSI[2]),
      localUpperZSI_(localOriginSI[2] +
        static_cast<Double>(localGeometry.nz) * localGeometry.dz),
      retainSpatialFlux_(retainSpatialFlux)
  {
    const Double spacing[3] = {
      globalGeometry_.dx, globalGeometry_.dy, globalGeometry_.dz
    };
    const std::size_t cells[3] = {
      globalGeometry_.nx, globalGeometry_.ny, globalGeometry_.nz
    };
    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        globalLowerSI_[axis] = globalOriginSI_[axis];
        globalUpperSI_[axis] = globalOriginSI_[axis] +
          static_cast<Double>(cells[axis]) * spacing[axis];
        globalToleranceSI_[axis] = coordinateTolerance(
          globalLowerSI_[axis], globalUpperSI_[axis]);
      }
    for (std::size_t face = 0; face < 6; ++face)
      {
        stepCount_[face] = 0;
        cumulativeCount_[face] = 0;
        stepCharge_[face] = 0.0;
        cumulativeCharge_[face] = 0.0;
      }
    validateGeometry();
  }

  void ParticleOpenBoundary::beginStep()
  {
    for (std::size_t face = 0; face < 6; ++face)
      {
        stepCount_[face] = 0;
        stepCharge_[face] = 0.0;
        if (retainSpatialFlux_ && !faceFlux_[face].values.empty())
          std::fill(faceFlux_[face].values.begin(),
                    faceFlux_[face].values.end(), 0.0);
      }
  }

  bool ParticleOpenBoundary::firstExit(
      const FieldVector<Double>& startSI,
      const FieldVector<Double>& endSI,
      ParticleBoundaryHit& hit) const
  {
    Double bestFraction = std::numeric_limits<Double>::infinity();
    ParticleBoundaryFace bestFace = ParticleBoundaryFace::LowerX;
    bool escaped = false;

    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        if (!std::isfinite(startSI[axis]) ||
            !std::isfinite(endSI[axis]))
          throw std::runtime_error(
            "Particle trajectory contains a non-finite coordinate");
        if (startSI[axis] < globalLowerSI_[axis] -
                              globalToleranceSI_[axis] ||
            startSI[axis] > globalUpperSI_[axis] +
                              globalToleranceSI_[axis])
          throw std::out_of_range(
            "Particle trajectory starts outside the global E/B domain");

        Double candidate = 0.0;
        ParticleBoundaryFace candidateFace =
          ParticleBoundaryFace::LowerX;
        bool axisEscaped = false;
        if (endSI[axis] < globalLowerSI_[axis])
          {
            candidate = (globalLowerSI_[axis] - startSI[axis]) /
              (endSI[axis] - startSI[axis]);
            candidateFace = static_cast<ParticleBoundaryFace>(2 * axis);
            axisEscaped = true;
          }
        else if (endSI[axis] > globalUpperSI_[axis])
          {
            candidate = (globalUpperSI_[axis] - startSI[axis]) /
              (endSI[axis] - startSI[axis]);
            candidateFace = static_cast<ParticleBoundaryFace>(2 * axis + 1);
            axisEscaped = true;
          }

        if (axisEscaped && candidate < bestFraction)
          {
            bestFraction = candidate;
            bestFace = candidateFace;
            escaped = true;
          }
      }

    if (!escaped) return false;
    const Double fractionTolerance = 64.0 *
      std::numeric_limits<Double>::epsilon();
    if (bestFraction < -fractionTolerance ||
        bestFraction > 1.0 + fractionTolerance ||
        !std::isfinite(bestFraction))
      throw std::runtime_error(
        "Cannot locate the first particle-domain exit");
    hit.fraction = std::max(0.0, std::min(1.0, bestFraction));
    hit.face = bestFace;
    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        hit.position[axis] = startSI[axis] + hit.fraction *
          (endSI[axis] - startSI[axis]);
        hit.position[axis] = std::max(globalLowerSI_[axis],
          std::min(globalUpperSI_[axis], hit.position[axis]));
      }
    const unsigned int normal = faceAxis(hit.face);
    hit.position[normal] = upperFace(hit.face) ?
      globalUpperSI_[normal] : globalLowerSI_[normal];
    return true;
  }

  void ParticleOpenBoundary::depositOutgoingFlux(
      const ParticleBoundaryHit& hit, Double chargeCoulomb)
  {
    if (!std::isfinite(chargeCoulomb))
      throw std::invalid_argument(
        "Escaped particle charge must be finite");
    validateHit(hit);
    const std::size_t face = faceIndex(hit.face);
    ++stepCount_[face];
    ++cumulativeCount_[face];
    stepCharge_[face] += chargeCoulomb;
    cumulativeCharge_[face] += chargeCoulomb;
    if (retainSpatialFlux_)
      depositSpatialFlux(hit.face, hit.position, chargeCoulomb);
  }

  unsigned long long ParticleOpenBoundary::stepCount(
      ParticleBoundaryFace face) const
  {
    return stepCount_[faceIndex(face)];
  }

  unsigned long long ParticleOpenBoundary::cumulativeCount(
      ParticleBoundaryFace face) const
  {
    return cumulativeCount_[faceIndex(face)];
  }

  Double ParticleOpenBoundary::stepCharge(ParticleBoundaryFace face) const
  {
    return stepCharge_[faceIndex(face)];
  }

  Double ParticleOpenBoundary::cumulativeCharge(
      ParticleBoundaryFace face) const
  {
    return cumulativeCharge_[faceIndex(face)];
  }

  std::size_t ParticleOpenBoundary::memoryBytes() const
  {
    std::size_t bytes = 0;
    for (std::size_t face = 0; face < 6; ++face)
      bytes += faceFlux_[face].values.capacity() * sizeof(Double);
    return bytes;
  }

  bool ParticleOpenBoundary::retainsSpatialFlux() const
  {
    return retainSpatialFlux_;
  }

  Double ParticleOpenBoundary::outwardCurrentDensity(
      ParticleBoundaryFace face, std::size_t first,
      std::size_t second) const
  {
    if (!retainSpatialFlux_)
      throw std::logic_error(
        "Spatial particle-boundary flux retention is disabled");
    return faceValue(face, first, second);
  }

  Double ParticleOpenBoundary::outwardDivergence(
      std::size_t i, std::size_t j, std::size_t k) const
  {
    if (!retainSpatialFlux_)
      throw std::logic_error(
        "Spatial particle-boundary flux retention is disabled");
    if (i > localGeometry_.nx || j > localGeometry_.ny ||
        k > localGeometry_.nz)
      throw std::out_of_range(
        "Particle-boundary divergence vertex is outside the local grid");

    Double divergence = 0.0;
    if (i == 0)
      divergence += faceValue(ParticleBoundaryFace::LowerX, j, k) /
                    localGeometry_.dx;
    if (i == localGeometry_.nx)
      divergence += faceValue(ParticleBoundaryFace::UpperX, j, k) /
                    localGeometry_.dx;
    if (j == 0)
      divergence += faceValue(ParticleBoundaryFace::LowerY, i, k) /
                    localGeometry_.dy;
    if (j == localGeometry_.ny)
      divergence += faceValue(ParticleBoundaryFace::UpperY, i, k) /
                    localGeometry_.dy;
    if (k == 0 &&
        std::abs(localLowerZSI_ - globalLowerSI_[2]) <=
          globalToleranceSI_[2])
      divergence += faceValue(ParticleBoundaryFace::LowerZ, i, j) /
                    localGeometry_.dz;
    if (k == localGeometry_.nz &&
        std::abs(localUpperZSI_ - globalUpperSI_[2]) <=
          globalToleranceSI_[2])
      divergence += faceValue(ParticleBoundaryFace::UpperZ, i, j) /
                    localGeometry_.dz;
    return divergence;
  }

  Double ParticleOpenBoundary::maxContinuityResidual(
      const EBFieldGrid& fields,
      const YeeComponent& rhoBefore,
      const YeeComponent& rhoAfter) const
  {
    if (!retainSpatialFlux_)
      throw std::logic_error(
        "Continuity diagnostics require spatial boundary flux retention");
    const EBGridGeometry& g = fields.geometry();
    if (g.nx != localGeometry_.nx || g.ny != localGeometry_.ny ||
        g.nz != localGeometry_.nz ||
        rhoBefore.nx() != g.nx + 1 || rhoBefore.ny() != g.ny + 1 ||
        rhoBefore.nz() != g.nz + 1 ||
        rhoAfter.nx() != g.nx + 1 || rhoAfter.ny() != g.ny + 1 ||
        rhoAfter.nz() != g.nz + 1)
      throw std::invalid_argument(
        "Open-boundary continuity diagnostic got incompatible grids");

    Double maximum = 0.0;
    for (std::size_t k = 0; k <= g.nz; ++k)
      for (std::size_t j = 0; j <= g.ny; ++j)
        for (std::size_t i = 0; i <= g.nx; ++i)
          {
            const Double rightX =
              i < g.nx ? fields.jx()(i, j, k) : 0.0;
            const Double leftX =
              i > 0 ? fields.jx()(i - 1, j, k) : 0.0;
            const Double rightY =
              j < g.ny ? fields.jy()(i, j, k) : 0.0;
            const Double leftY =
              j > 0 ? fields.jy()(i, j - 1, k) : 0.0;
            const Double rightZ =
              k < g.nz ? fields.jz()(i, j, k) : 0.0;
            const Double leftZ =
              k > 0 ? fields.jz()(i, j, k - 1) : 0.0;
            const Double residual =
              (rhoAfter(i, j, k) - rhoBefore(i, j, k)) / g.dt +
              (rightX - leftX) / g.dx +
              (rightY - leftY) / g.dy +
              (rightZ - leftZ) / g.dz +
              outwardDivergence(i, j, k);
            maximum = std::max(maximum, std::abs(residual));
          }
    return maximum;
  }

  std::size_t ParticleOpenBoundary::faceIndex(
      ParticleBoundaryFace face)
  {
    const int value = static_cast<int>(face);
    if (value < 0 || value >= 6)
      throw std::invalid_argument("Unknown particle boundary face");
    return static_cast<std::size_t>(value);
  }

  unsigned int ParticleOpenBoundary::faceAxis(
      ParticleBoundaryFace face)
  {
    return static_cast<unsigned int>(faceIndex(face) / 2);
  }

  bool ParticleOpenBoundary::upperFace(ParticleBoundaryFace face)
  {
    return (faceIndex(face) % 2) != 0;
  }

  void ParticleOpenBoundary::validateGeometry() const
  {
    if (localGeometry_.nx == 0 || localGeometry_.ny == 0 ||
        localGeometry_.nz == 0 || globalGeometry_.nx == 0 ||
        globalGeometry_.ny == 0 || globalGeometry_.nz == 0 ||
        !(localGeometry_.dx > 0.0) || !(localGeometry_.dy > 0.0) ||
        !(localGeometry_.dz > 0.0) || !(localGeometry_.dt > 0.0) ||
        !std::isfinite(localGeometry_.dx) ||
        !std::isfinite(localGeometry_.dy) ||
        !std::isfinite(localGeometry_.dz) ||
        !std::isfinite(localGeometry_.dt) ||
        localGeometry_.nx != globalGeometry_.nx ||
        localGeometry_.ny != globalGeometry_.ny ||
        localGeometry_.dx != globalGeometry_.dx ||
        localGeometry_.dy != globalGeometry_.dy ||
        localGeometry_.dz != globalGeometry_.dz ||
        localGeometry_.dt != globalGeometry_.dt)
      throw std::invalid_argument(
        "Particle open boundary requires a compatible z-slab geometry");
    for (unsigned int axis = 0; axis < 3; ++axis)
      if (!std::isfinite(globalOriginSI_[axis]) ||
          !std::isfinite(localOriginSI_[axis]) ||
          !std::isfinite(globalUpperSI_[axis]))
        throw std::invalid_argument(
          "Particle open-boundary origin or extent is not finite");
    if (localLowerZSI_ < globalLowerSI_[2] - globalToleranceSI_[2] ||
        localUpperZSI_ > globalUpperSI_[2] + globalToleranceSI_[2])
      throw std::invalid_argument(
        "Particle open-boundary slab lies outside the global domain");
  }

  void ParticleOpenBoundary::validateHit(
      const ParticleBoundaryHit& hit) const
  {
    if (!std::isfinite(hit.fraction) || hit.fraction < 0.0 ||
        hit.fraction > 1.0)
      throw std::invalid_argument(
        "Particle boundary hit has an invalid trajectory fraction");
    const unsigned int normal = faceAxis(hit.face);
    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        if (!std::isfinite(hit.position[axis]) ||
            hit.position[axis] < globalLowerSI_[axis] -
                                   globalToleranceSI_[axis] ||
            hit.position[axis] > globalUpperSI_[axis] +
                                   globalToleranceSI_[axis])
          throw std::out_of_range(
            "Particle boundary hit lies outside the global face");
        if (axis == normal)
          {
            const Double expected = upperFace(hit.face) ?
              globalUpperSI_[axis] : globalLowerSI_[axis];
            if (std::abs(hit.position[axis] - expected) >
                globalToleranceSI_[axis])
              throw std::invalid_argument(
                "Particle boundary hit is not on its selected face");
          }
      }
    const Double zTolerance = coordinateTolerance(
      localLowerZSI_, localUpperZSI_);
    if (hit.position[2] < localLowerZSI_ - zTolerance ||
        hit.position[2] > localUpperZSI_ + zTolerance)
      throw std::out_of_range(
        "Particle boundary hit does not belong to this MPI slab");
  }

  void ParticleOpenBoundary::ensureFaceStorage(
      ParticleBoundaryFace face)
  {
    FaceFlux& storage = faceFlux_[faceIndex(face)];
    if (!storage.values.empty()) return;
    const unsigned int axis = faceAxis(face);
    if (axis == 0)
      {
        storage.firstNodes = localGeometry_.ny + 1;
        storage.secondNodes = localGeometry_.nz + 1;
      }
    else if (axis == 1)
      {
        storage.firstNodes = localGeometry_.nx + 1;
        storage.secondNodes = localGeometry_.nz + 1;
      }
    else
      {
        storage.firstNodes = localGeometry_.nx + 1;
        storage.secondNodes = localGeometry_.ny + 1;
      }
    if (storage.secondNodes >
        std::numeric_limits<std::size_t>::max() / storage.firstNodes)
      throw std::overflow_error(
        "Particle-boundary diagnostic face size overflows size_t");
    storage.values.assign(storage.firstNodes * storage.secondNodes, 0.0);
  }

  void ParticleOpenBoundary::depositSpatialFlux(
      ParticleBoundaryFace face,
      const FieldVector<Double>& positionSI,
      Double chargeCoulomb)
  {
    ensureFaceStorage(face);
    const unsigned int axis = faceAxis(face);
    const unsigned int firstAxis = axis == 0 ? 1 : 0;
    const unsigned int secondAxis = axis == 2 ? 1 : 2;
    const Double spacing[3] = {
      localGeometry_.dx, localGeometry_.dy, localGeometry_.dz
    };
    const std::size_t cells[3] = {
      localGeometry_.nx, localGeometry_.ny, localGeometry_.nz
    };
    std::size_t firstCell = 0;
    std::size_t secondCell = 0;
    Double firstFraction = 0.0;
    Double secondFraction = 0.0;
    cellCoordinate((positionSI[firstAxis] - localOriginSI_[firstAxis]) /
                     spacing[firstAxis],
                   cells[firstAxis], firstCell, firstFraction);
    cellCoordinate((positionSI[secondAxis] - localOriginSI_[secondAxis]) /
                     spacing[secondAxis],
                   cells[secondAxis], secondCell, secondFraction);
    const Double scale = chargeCoulomb /
      (spacing[firstAxis] * spacing[secondAxis] * localGeometry_.dt);
    FaceFlux& storage = faceFlux_[faceIndex(face)];
    for (unsigned int second = 0; second < 2; ++second)
      for (unsigned int first = 0; first < 2; ++first)
        {
          const Double firstWeight = first ? firstFraction :
            1.0 - firstFraction;
          const Double secondWeight = second ? secondFraction :
            1.0 - secondFraction;
          const std::size_t index =
            (secondCell + second) * storage.firstNodes +
            firstCell + first;
          storage.values[index] += scale * firstWeight * secondWeight;
        }
  }

  Double ParticleOpenBoundary::faceValue(
      ParticleBoundaryFace face, std::size_t first,
      std::size_t second) const
  {
    const FaceFlux& storage = faceFlux_[faceIndex(face)];
    if (storage.values.empty()) return 0.0;
    if (first >= storage.firstNodes || second >= storage.secondNodes)
      throw std::out_of_range(
        "Particle-boundary face index is outside the diagnostic lattice");
    return storage.values[second * storage.firstNodes + first];
  }
}
