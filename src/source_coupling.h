#ifndef APRL_CONSERVATIVE_SOURCE_COUPLING_H
#define APRL_CONSERVATIVE_SOURCE_COUPLING_H
#include <array>
#include <memory>
#include <string>
#include <vector>
#include <mpi.h>
#include "eb_particles.h"
namespace aprl {
class ConservativeSourceCoupling {
public:
  ConservativeSourceCoupling(MPI_Comm,const std::string& directory);
  ~ConservativeSourceCoupling();
  void begin(EBFieldGrid&,const FieldVector<Double>& origin,const std::vector<RelativisticParticleSI>&);
  void beginTrial();
  void setReplay(bool enabled);
  void exchangeCandidate();
  void applyCurrentSource(Double sign);
  std::array<Double,4> localResidual(); // force, absolute work defect, work scale, active calls
  void swapElectric();
  void push(RelativisticParticleSI&,const SIFieldSourceSet&,const BoostFrameTransform&,Double time,Double dt,unsigned substeps,Double* prescribed,BorisWork*,ParticleSubstepObserver*);
  void recordConstraint(const FieldVector<Double>& begin,const FieldVector<Double>& end,Double charge,bool carrier);
  std::array<long double,6> finalizeConstraints(); // frozen work/abs/count, carrier work/abs/count
  void finish(std::size_t step,unsigned iterations,Double force,Double work);
  void close(std::size_t step);
private:
  struct Impl;std::unique_ptr<Impl> impl_;
};
}
#endif
