#ifndef APRL_VERIFICATION_BLOCKED_COWAN_REFERENCE_H
#define APRL_VERIFICATION_BLOCKED_COWAN_REFERENCE_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "eb_field.h"
#include "eb_particles.h"

namespace aprl { namespace verification {

struct BlockedCowanReferenceReport {
  bool complete = false;
  std::size_t particles = 0;
  std::size_t interiorNodes[3] = {0,0,0};
  std::size_t blockLayers = 0;
  std::size_t blocks = 0;
  double beta = 0.;
  double sourceCharge = 0.;
  double depositedCharge = 0.;
  double minimumPivot = 0.;
  double minimumModeEigenvalue = 0.;
  double preparationSeconds = 0.;
  double sourceSeconds = 0.;
  double forwardSeconds = 0.;
  double reconstructionSeconds = 0.;
  double querySeconds = 0.;
  std::uint64_t sourceFileBytes = 0;
  std::uint64_t checkpointFileBytes = 0;
  std::uint64_t potentialFileBytes = 0;
  std::uint64_t bytesRead = 0;
  std::uint64_t bytesWritten = 0;
  std::uint64_t fieldQueries = 0;
  std::uint64_t potentialCacheHits = 0;
  std::uint64_t potentialCacheMisses = 0;
  std::uint64_t fieldCacheBuilds = 0;
  std::size_t sourceIndexBytes = 0;
  std::size_t sourceSortingAllowanceBytes = 0;
  std::size_t solverWorkspaceBytes = 0;
  std::size_t potentialCacheCapacityBytes = 0;
  std::size_t sourceCacheCapacityBytes = 0;
  std::size_t fieldCacheCapacityBytes = 0;
  std::size_t queryWorkspaceBytes = 0;
  std::size_t residentArrayBytes = 0;
};

// Verification-only, one MPI rank, fixed common-beta finite Dirichlet source.
// No CPML-source padding or moving/time-dependent source is implemented.
// The caller supplies the same mass-weighted beta used by the CG initializer.
// Construction refuses an existing work directory and leaves rho.bin,
// checkpoints.bin and phi.bin on disk, including after destruction.
// Input particles are never modified or retained; O(P) temporary contribution
// indices preserve original CIC deposition order and are reported.
//
// field() returns native E(0),B(-dt/2), component order Ex,Ey,Ez,Bx,By,Bz.
// i,j,k are GLOBAL original-grid native indices; invalid component/grid
// indices throw. Finite-domain guard/boundary values match the production
// Cowan initializer, including zero filtered potentials outside its interior.
// A maximum of three potential blocks and one native E/B block (with one
// right-hand E/B plane for cross-block interpolation) are cached. Scalar
// charge diagnostics additionally cache one source block, only on demand;
// calls are serial, not thread safe. Internal bulk construction may parallelize.
class BlockedCowanReference {
public:
  BlockedCowanReference(const EBGridGeometry& globalGeometry,
      const FieldVector<Double>& globalOrigin,
      const std::vector<RelativisticParticleSI>& particles,
      Double beta, const std::size_t guards[3], std::size_t blockLayers,
      const std::string& directory);
  ~BlockedCowanReference();
  BlockedCowanReference(const BlockedCowanReference&) = delete;
  BlockedCowanReference& operator=(const BlockedCowanReference&) = delete;

  Double field(unsigned component, std::size_t i, std::size_t j,
               std::size_t globalK);
  // Writes only E/B. J, particles and boundary auxiliary state remain owned by
  // the caller. x/y grid must match; local z slices may be filled by offset.
  void fill(EBFieldGrid& fields, std::size_t globalZOffset = 0);
  // Read-only raw scalar diagnostics, zero on/outside static Dirichlet faces.
  Double chargeDensity(std::size_t i, std::size_t j, std::size_t globalK);
  Double potential(std::size_t i, std::size_t j, std::size_t globalK);
  void clearCache();
  const BlockedCowanReferenceReport& report() const;
  void writeReport(const std::string& path) const; // Refuses overwrite.

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}} // namespace aprl::verification
#endif
