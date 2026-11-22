#include "simulation.h"
#include <algorithm>
#include <array>
#include <stdexcept>
#include <sstream>
namespace aprl {
bool Simulation::conservativeLongitudinalHandoff(const ParticleBoundaryHit& hit,
    const FieldVector<Double>& end) const {
 if(!sourceCoupling_||!movingWindow_||static_cast<unsigned>(hit.face)<4)return false;
 // Preserve transverse CPML damping and its original first-entry policy.
 // A pure longitudinal handoff completes the already performed field step.
 for(unsigned a=0;a<2;++a)if(end[a]<particleCPML_->physicalLower(a)||end[a]>particleCPML_->physicalUpper(a))return false;
 return true;
}
void Simulation::solveConservativeSourceStep() {
 auto& coupling=*sourceCoupling_;const Double cacheStart=performanceNow();coupling.begin(*fields_,localOriginBox_,particles_,config_.boundary.cpml.cells[0],config_.boundary.cpml.cells[1]);performanceSeconds_[CouplingCacheStage]+=performanceNow()-cacheStart;
 // Only the zero-current electric curl is evaluated. The reversible boundary
 // preview keeps its phase and CPML histories for the one later real commit.
 movingWindow_->particleHalo->clear(*fields_);
 const Double previewStart=performanceNow();movingWindow_->fields->previewElectricBoundary(true);
 halo_->advanceElectric(*fields_,EBCurrentFilter::None);
 movingWindow_->fields->previewElectricBoundary(false);
 performanceSeconds_[FieldStage]+=performanceNow()-previewStart;
 movingWindow_->particleHalo->exchangeAndInstall(*fields_);
 const auto originalParticles=particles_;const auto originalCarriers=pmlCarriers_;
 const auto lost=localLostParticles_,ballistic=ballisticBoundaryEntries_,activations=frozenPlaneActivations_,crossings=crossRankActivationPushes_;
 const auto peak=peakPmlCarriers_;const auto exits=movingWindow_->exits;const auto exitCharge=movingWindow_->exitCharge;
 unsigned long long boundaryCounts[6];Double boundaryCharges[6];particleBoundary_->checkpointCumulative(boundaryCounts,boundaryCharges);
 std::array<unsigned long long,6> cpml,direct,carrier;
 std::array<Double,6> cpmlCharge,directCharge,carrierCharge;
 std::copy(cpmlEntryCount_,cpmlEntryCount_+6,cpml.begin());std::copy(directOuterCount_,directOuterCount_+6,direct.begin());std::copy(carrierOuterCount_,carrierOuterCount_+6,carrier.begin());
 std::copy(cpmlEntryCharge_,cpmlEntryCharge_+6,cpmlCharge.begin());std::copy(directOuterCharge_,directOuterCharge_+6,directCharge.begin());std::copy(carrierOuterCharge_,carrierOuterCharge_+6,carrierCharge.begin());
 auto restore=[&](){particles_=originalParticles;pmlCarriers_=originalCarriers;localLostParticles_=lost;ballisticBoundaryEntries_=ballistic;frozenPlaneActivations_=activations;crossRankActivationPushes_=crossings;peakPmlCarriers_=peak;movingWindow_->exits=exits;movingWindow_->exitCharge=exitCharge;
  std::copy(cpml.begin(),cpml.end(),cpmlEntryCount_);std::copy(direct.begin(),direct.end(),directOuterCount_);std::copy(carrier.begin(),carrier.end(),carrierOuterCount_);std::copy(cpmlCharge.begin(),cpmlCharge.end(),cpmlEntryCharge_);std::copy(directCharge.begin(),directCharge.end(),directOuterCharge_);std::copy(carrierCharge.begin(),carrierCharge.end(),carrierOuterCharge_);particleBoundary_->restoreCumulative(boundaryCounts,boundaryCharges);particleBoundary_->beginStep();};
 const bool ledgerEnabled=config_.energyLedger.enabled,trajectoryEnabled=config_.trajectory.enabled;
 auto savedScreen=std::move(radiationScreen_);
 auto savedDetectors=std::move(detectors_);config_.energyLedger.enabled=false;config_.trajectory.enabled=false;
 bool converged=false;unsigned iterations=0;Double force=0.,work=0.;
 try {
  for(unsigned n=0;n<8;++n){if(n)coupling.applyCurrentSource(1.);fields_->clearCurrent();restore();coupling.beginTrial();coupling.exchangeCandidate();pushDepositAndMigrate();coupling.completeCurrent();
   if(config_.mesh.movingWindow.referenceSurfaceCurrent){applyMovingTransverseReferenceCurrent();finishMovingTransverseReferenceCurrent();}
   halo_->sumSharedTangentialCurrent(*fields_);coupling.applyCurrentSource(-1.);coupling.exchangeCandidate();
   const Double residualStart=performanceNow();const auto r=coupling.localResidual();performanceSeconds_[CouplingResidualStage]+=performanceNow()-residualStart;const Double local[3]={r[1],r[2],r[3]};Double total[3]={};
   MPI_Allreduce(&r[0],&force,1,MPI_DOUBLE,MPI_MAX,communicator_);MPI_Allreduce(local,total,3,MPI_DOUBLE,MPI_SUM,communicator_);
   work=total[2]==0.?0.:(total[1]>0.?total[0]/total[1]:total[0]);iterations=n+1;
   if(iterations>=(config_.mesh.oneWaySourceDiagnostic?1u:(total[2]>0.?2u:1u))&&force<1e-13&&work<1e-8){converged=true;break;}
  }
 }catch(...){config_.energyLedger.enabled=ledgerEnabled;config_.trajectory.enabled=trajectoryEnabled;detectors_=std::move(savedDetectors);radiationScreen_=std::move(savedScreen);throw;}
 config_.energyLedger.enabled=ledgerEnabled;config_.trajectory.enabled=trajectoryEnabled;detectors_=std::move(savedDetectors);radiationScreen_=std::move(savedScreen);
 if(!converged){std::ostringstream message;message<<"Conservative source coupling failed at step "<<step_+1<<": force="<<force<<" work="<<work<<" iterations="<<iterations;throw std::runtime_error(message.str());}
 // Remove the trial current, retain the zero-current wave by swapping it into
 // the old-E buffer, and replay the converged substeps against the real old E.
 coupling.applyCurrentSource(1.);fields_->clearCurrent();restore();movingWindow_->particleHalo->clear(*fields_);coupling.swapElectric();movingWindow_->particleHalo->exchangeAndInstall(*fields_);
 coupling.setReplay(true);pushDepositAndMigrate();coupling.completeCurrent();coupling.setReplay(false);coupling.finish(step_+1,iterations,force,work);
}
void Simulation::commitConservativeElectric() {
 sourceCoupling_->swapElectric();movingWindow_->fields->undoPreviewElectricBoundary(*fields_);
 halo_->sumSharedTangentialCurrent(*fields_);sourceCoupling_->applyCurrentSource(-1.);
 movingWindow_->fields->commitPreparedElectricBoundary(*fields_);
}
}
