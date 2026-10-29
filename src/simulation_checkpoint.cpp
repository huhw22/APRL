#include "simulation.h"
#include "checkpoint_io.h"
#include "runtime_util.h"
#include "shutdown_watchdog.h"

#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unistd.h>

namespace aprl
{
  namespace
  {
    // Small restart bookkeeping uses native arithmetic representations so
    // long-double energy sums retain every bit. The explicit ABI contract is
    // checked before decoding. Large arrays have portable HDF5 numeric types.
    class StateBuffer
    {
    public:
      StateBuffer(std::vector<unsigned char>& bytes, bool reading)
        : bytes_(bytes), reading_(reading), offset_(0)
      { if (!reading_) bytes_.clear(); }
      template<class T> void value(T& item)
      {
        static_assert(std::is_arithmetic<T>::value,
                      "checkpoint scalar must be arithmetic");
        if (reading_)
          {
            if (offset_ > bytes_.size() || sizeof(T) > bytes_.size()-offset_)
              throw std::runtime_error("Truncated checkpoint scalar state");
            std::memcpy(&item, bytes_.data()+offset_, sizeof(T));
          }
        else
          {
            const unsigned char* p =
              reinterpret_cast<const unsigned char*>(&item);
            bytes_.insert(bytes_.end(), p, p+sizeof(T));
          }
        offset_ += sizeof(T);
      }
      void value(bool& item)
      {
        unsigned char v = item ? 1 : 0;
        value(v);
        if (v > 1) throw std::runtime_error("Invalid checkpoint boolean");
        if (reading_) item = v != 0;
      }
      template<class T, std::size_t N> void value(T (&items)[N])
      { for (std::size_t i=0; i<N; ++i) value(items[i]); }
      void finish() const
      {
        if (reading_ && offset_ != bytes_.size())
          throw std::runtime_error("Unexpected checkpoint scalar state tail");
      }
    private:
      std::vector<unsigned char>& bytes_;
      bool reading_;
      std::size_t offset_;
    };

    std::string nativeAbi()
    {
      const unsigned int endian = 1;
      std::ostringstream out;
      out << "aprl-step-checkpoint-v2;endian="
          << static_cast<unsigned int>(
               *reinterpret_cast<const unsigned char*>(&endian))
          << ";size=" << sizeof(std::size_t)
          << ";uint=" << sizeof(unsigned int)
          << ";ull=" << sizeof(unsigned long long)
          << ";double=" << sizeof(Double) << ':'
          << std::numeric_limits<Double>::digits << ':'
          << std::numeric_limits<Double>::max_exponent
          << ";long-double=" << sizeof(long double) << ':'
          << std::numeric_limits<long double>::digits << ':'
          << std::numeric_limits<long double>::max_exponent;
      return out.str();
    }

    std::string canonical(const YAML::Node& node)
    {
      if (!node || node.IsNull()) return "null";
      if (node.IsScalar()) return "s" + YAML::Dump(node);
      std::ostringstream out;
      if (node.IsSequence())
        {
          out << '[';
          for (std::size_t i=0; i<node.size(); ++i)
            { const std::string v=canonical(node[i]); out<<v.size()<<':'<<v; }
          out << ']';
        }
      else
        {
          std::vector<std::string> keys;
          for (YAML::const_iterator i=node.begin(); i!=node.end(); ++i)
            keys.push_back(i->first.as<std::string>());
          std::sort(keys.begin(), keys.end());
          out << '{';
          for (const std::string& key: keys)
            { const std::string v=canonical(node[key]);
              out<<key.size()<<':'<<key<<v.size()<<':'<<v; }
          out << '}';
        }
      return out.str();
    }

    void removeOutputTuning(YAML::Node node)
    {
      if (!node || !node.IsMap()) return;
      for (const char* key: {"directory", "filename", "basename",
                            "buffer_records", "buffer_samples", "compression",
                            "flush_every_samples", "warning_relative_tolerance"})
        node.remove(key);
    }

    std::string continuationContract(const std::string& text,
                                     bool retirementChange, bool resetDiagnostics)
    {
      YAML::Node node = YAML::Load(text);
      for (const char* key: {"checkpoint", "output", "runtime", "stop",
                            "radiation_resolution"}) node.remove(key);
      node["mesh"].remove("duration");
      if (retirementChange) node.remove("particle_retirement");
      if (resetDiagnostics)
        { node.remove("detectors"); node.remove("trajectory"); }
      else
        {
          removeOutputTuning(node["trajectory"]);
          YAML::Node detectors = node["detectors"];
          removeOutputTuning(detectors);
          if (detectors)
            for (const char* list: {"field_planes", "particle_planes"})
              {
                YAML::Node planes=detectors[list];
                for (std::size_t i=0; planes && i<planes.size(); ++i)
                  { removeOutputTuning(planes[i]);
                    removeOutputTuning(planes[i]["particle_background"]); }
              }
        }
      removeOutputTuning(node["energy_ledger"]);
      return canonical(node);
    }

    std::string readText(CheckpointArchive& file, const std::string& name)
    {
      const std::vector<unsigned char> bytes=file.readLocalBytes(name);
      return std::string(bytes.begin(),bytes.end());
    }
    void writeText(CheckpointArchive& file, const std::string& name,
                   const std::string& text)
    { file.writeBytes(name,text.data(),text.size()); }

    std::string absolutePath(const std::string& path)
    {
      if (!path.empty() && path[0]=='/') return path;
      char* cwd=::getcwd(NULL,0);
      if (!cwd) throw std::runtime_error("Cannot resolve checkpoint output path");
      const std::string result=std::string(cwd)+"/"+path;
      std::free(cwd);
      return result;
    }

    void requireCheckpointWatchdog(unsigned int seconds)
    {
      if (!armCheckpointWatchdog(seconds))
        throw std::runtime_error(
          "Checkpoint requires the process-local I/O watchdog; worker unavailable");
    }
  }

  void Simulation::transferCheckpointScalars(
      std::vector<unsigned char>& bytes, bool reading, bool geometry)
  {
    StateBuffer state(bytes,reading);
    unsigned int version=2;
    state.value(version);
    if (version!=2) throw std::runtime_error("Unsupported simulation checkpoint schema");
    if (geometry)
      {
        state.value(config_.mesh.cells);
        for (unsigned int a=0;a<3;++a) state.value(config_.mesh.extent[a]);
        Double origins[3]={frame_.labTimeOrigin(),frame_.labZOrigin(),frame_.boxZOrigin()};
        state.value(origins);
        if (reading) frame_.setOriginsFromGamma(config_.mesh.boostGamma,SI::c,
                                                origins[0],origins[1],origins[2]);
        state.value(referenceCenterBoxZ_);
        state.value(particleSubsteps_);
        state.value(initialMaximumLabGamma_);
        state.value(initialGlobalParticleCount_);
        state.value(hasLongitudinalLayoutReport_);
        LongitudinalLayout& l=longitudinalLayoutReport_.layout;
        state.value(l.totalCells); state.value(l.physicalCells);
        state.value(l.minimumPhysicalCells); state.value(l.alignmentExtraCells);
        state.value(l.extent); state.value(l.origin); state.value(l.outerLower);
        state.value(l.outerUpper); state.value(l.physicalLower); state.value(l.physicalUpper);
        state.value(l.translation); state.value(l.initialMinimum); state.value(l.initialMaximum);
        state.value(l.initialMean); state.value(l.ballisticMinimum); state.value(l.ballisticMaximum);
        state.value(l.radiationUpper); state.value(l.selfFieldLower); state.value(l.selfFieldUpper);
        state.value(l.contentLower); state.value(l.contentUpper); state.value(l.requiredLower);
        state.value(l.requiredUpper); state.value(l.rightRoundingPadding);
        state.value(longitudinalLayoutReport_.observationTargetLab);
        state.value(longitudinalLayoutReport_.initialHeadLab);
        state.value(longitudinalLayoutReport_.ballisticStopTime);
        state.value(longitudinalLayoutReport_.plannedBallisticTime);
        state.value(longitudinalLayoutReport_.radiationAdvance);
      }
    else
      {
        state.value(timeBoxSI_); state.value(step_);
        state.value(nextTrajectorySampleTime_); state.value(trajectorySamplesSinceFlush_);
        state.value(energyLedgerPreviousPower_); state.value(energyLedgerOutwardEnergy_);
        state.value(energyLedgerPrescribedWork_); state.value(energyLedgerBoundarySourceWork_);
        state.value(energyLedgerBoundarySourceAbsWork_); state.value(energyLedgerRemovedKinetic_);
        state.value(energyLedgerRemovedTotal_); state.value(energyLedgerInitialKinetic_);
        state.value(energyLedgerBorisWork_.gridWork);
        state.value(energyLedgerBorisWork_.prescribedWork);
        state.value(energyLedgerBorisWork_.magneticRotationDeltaK);
        state.value(energyLedgerBorisWork_.totalDeltaK);
        state.value(energyLedgerBorisWork_.roundoffResidual);
        for (ParticleEnergyMoment* moment : {&energyLedgerRemovedMoment_,
             &energyLedgerReleasedMoment_, &energyLedgerInitialMoment_}) {
          state.value(moment->kineticBox); state.value(moment->totalBox);
          state.value(moment->cPzBox); state.value(moment->kineticLab);
          state.value(moment->totalLab); state.value(moment->cPzLab);
        }
        state.value(energyLedgerReleasedCharge_); state.value(energyLedgerPartialWorkSubsteps_);
        state.value(energyLedgerInitialField_); state.value(energyLedgerRemovedMacroparticles_);
        state.value(energyLedgerLastSampleStep_); state.value(energyLedgerHasReference_);
        state.value(localLostParticles_); state.value(lostParticles_);
        state.value(ballisticBoundaryEntries_); state.value(frozenPlaneActivations_);
        state.value(crossRankActivationPushes_); state.value(cpmlEntryCount_);
        state.value(cpmlEntryCharge_); state.value(directOuterCount_);
        state.value(directOuterCharge_); state.value(carrierOuterCount_);
        state.value(carrierOuterCharge_); state.value(peakPmlCarriers_);
        state.value(retirementEntryCount_); state.value(retirementEntryCharge_);
        state.value(retirementExitCount_); state.value(retirementExitResidualCharge_);
        state.value(peakRetirementCarriers_);
        unsigned long long counts[6]={}; Double charges[6]={};
        if (!reading) particleBoundary_->checkpointCumulative(counts,charges);
        state.value(counts); state.value(charges);
        if (reading) particleBoundary_->restoreCumulative(counts,charges);
        EnergyLedgerRecord& r=lastEnergyLedgerRecord_;
        state.value(r.step); state.value(r.activeMacroparticles); state.value(r.removedMacroparticles);
        state.value(r.boundaryInjectionEntries); state.value(r.closureValid); state.value(r.timeBox);
        state.value(r.activeRepresentedElectrons); state.value(r.particleKineticActive);
        state.value(r.particleTotalActive); state.value(r.particleKineticRemoved);
        state.value(r.particleTotalRemoved); state.value(r.fieldEnergyInterior);
        state.value(r.outwardFieldEnergyFace); state.value(r.outwardFieldEnergy);
        state.value(r.prescribedWork); state.value(r.boundarySourceWork);
        state.value(r.boundarySourceAbsWork); state.value(r.balanceResidualRaw);
        state.value(r.balanceResidual); state.value(r.relativeBalanceExchange);
        state.value(r.relativeBalanceInitial); state.value(r.fieldFractionKinetic);
        state.value(r.fieldFractionIncludingRest); state.value(r.meanGammaLab);
        state.value(r.sigmaGammaLab); state.value(r.relativeEnergySpreadLab);
        state.value(r.minimumGammaLab); state.value(r.maximumGammaLab);
        state.value(r.linearChirpGammaPerM); state.value(r.uncorrelatedSigmaGammaLab);
        state.value(r.meanGammaBox); state.value(r.sigmaGammaBox);
        state.value(r.movingInsertedField);
        state.value(r.movingCroppedField);
        state.value(r.movingRemapFieldChange);
        state.value(r.movingRemapInterfaceCorrection);
        state.value(r.inletOverwriteWork);
        state.value(r.longitudinalBoundaryOutward);
        state.value(r.movingEnergyDiagnosticRankSeconds);
        state.value(r.frozenConstraintWork);
        state.value(r.carrierConstraintWork);
        state.value(r.constraintEnergyDiagnosticRankSeconds);
        state.value(r.particleKineticBoxRemovedEvent);
        state.value(r.particleRemovalInterpolationDelta);
        state.value(r.particleGridWork);
        state.value(r.prescribedBorisWork);
        state.value(r.magneticRotationDeltaK);
        state.value(r.borisPushDeltaK);
        state.value(r.borisRoundoffResidual);
        state.value(r.particleWorkBalanceResidual);
        state.value(r.fieldExchangeResidual);
        state.value(r.balanceResidualBoris);
        state.value(r.particleKineticLabActive);
        state.value(r.particleKineticLabRemoved);
        state.value(r.particleCPzBoxActive);
        state.value(r.particleCPzBoxRemoved);
        state.value(r.particleCPzLabActive);
        state.value(r.particleCPzLabRemoved);
        state.value(r.deltaParticleKineticBox);
        state.value(r.deltaParticleKineticLab);
        state.value(r.deltaParticleCPzBox);
        state.value(r.particleLorentzDeltaResidual);
        state.value(r.deltaFieldEnergy);
        state.value(r.releasedKineticBox);
        state.value(r.releasedKineticLab);
        state.value(r.releasedCPzBox);
        state.value(r.releasedCPzLab);
        state.value(r.releasedCharge);
        state.value(r.partialWorkSubsteps);
        state.value(r.boostGamma);
        state.value(r.boostBeta);

      }
    state.finish();
  }

  void Simulation::loadCheckpointGeometry()
  {
    requireCheckpointWatchdog(config_.checkpoint.ioTimeoutSeconds);
    config_.checkpoint.restartFrom=absolutePath(config_.checkpoint.restartFrom);
    if (config_.output.overwrite)
      throw std::runtime_error("Restart requires fresh output paths and output.overwrite: false");
    checkpointReader_.reset(new CheckpointArchive(
      communicator_,config_.checkpoint.restartFrom,false));
    if (readText(*checkpointReader_,"metadata/native_abi")!=nativeAbi())
      throw std::runtime_error("Checkpoint native scalar ABI differs from this executable");
    const std::string parent=readText(*checkpointReader_,"metadata/configuration");
    if (continuationContract(parent,config_.checkpoint.allowRetirementChange,
          config_.checkpoint.resetDiagnostics)!=
        continuationContract(config_.configurationText,config_.checkpoint.allowRetirementChange,
          config_.checkpoint.resetDiagnostics))
      throw std::runtime_error(
        "Restart physics/diagnostic contract differs; preserve mesh, boost, beam, sources, "
        "boundary, feedback, initial-field and energy settings. Use explicit branch options "
        "only for future retirement or new diagnostic segments.");
    retirementChangedOnRestart_=canonical(YAML::Load(parent)["particle_retirement"])!=
      canonical(YAML::Load(config_.configurationText)["particle_retirement"]);
    std::vector<unsigned char> bytes=checkpointReader_->readLocalBytes("state/geometry");
    const std::size_t requestedCells[3]={config_.mesh.cells[0],config_.mesh.cells[1],config_.mesh.cells[2]};
    transferCheckpointScalars(bytes,true,true);
    for (unsigned int a=0;a<3;++a)
      if ((a!=2 || !config_.mesh.longitudinal.automatic()) &&
          config_.mesh.cells[a]!=requestedCells[a])
        throw std::runtime_error("Checkpoint resolved mesh does not match requested mesh");
    if (!particleSubsteps_ || particleSubsteps_>config_.mesh.maximumParticleSubsteps ||
        !std::isfinite(initialMaximumLabGamma_) || initialMaximumLabGamma_<1)
      throw std::runtime_error("Invalid checkpoint initial beam/substep state");
    initializeGeometry();
    particleCPML_.reset(new ParticleCPMLRegion(globalGeometry_,globalOriginBox_,
      config_.boundary.type==EBBoundaryType::Cpml ? config_.boundary.cpml : EBCPMLParameters()));
    initializeFieldDetectorRegions();
    validateBeamlineExclusionRules();
    checkpointReader_->readParticles("particles",particles_);
    checkpointReader_->readCarriers("carriers/cpml",pmlCarriers_);
    checkpointReader_->readCarriers("carriers/retirement",retirementCarriers_);
    // Scalars include time and cumulative open-boundary state. Construct the
    // latter early; initialize() will reconstruct it before the second restore.
    particleBoundary_.reset(new ParticleOpenBoundary(localGeometry_,localOriginBox_,
      globalGeometry_,globalOriginBox_,localZOffset_));
    bytes=checkpointReader_->readLocalBytes("state/scalars");
    transferCheckpointScalars(bytes,true,false);
    if (!std::isfinite(timeBoxSI_) || timeBoxSI_<0 || !(totalTimeBoxSI_>timeBoxSI_))
      throw std::runtime_error("Restart mesh.duration must exceed the saved absolute box time");
    const Double stepsReal=std::ceil(totalTimeBoxSI_/globalGeometry_.dt);
    if (!std::isfinite(stepsReal) || stepsReal<=0 ||
        stepsReal>=static_cast<Double>(std::numeric_limits<std::size_t>::max()))
      throw std::runtime_error("Restart duration produces an invalid field-step count");
    estimatedMaximumSteps_=static_cast<std::size_t>(stepsReal);
    config_.checkpoint.restartStep=step_;
    config_.checkpoint.restartTimeBox=timeBoxSI_;
    validateCheckpointBranch();
    if (configuredStopReached())
      throw std::runtime_error(
        "Checkpoint already satisfies the requested physical stop; extend stop explicitly to continue");
    nextCheckpointStep_=std::upper_bound(config_.checkpoint.steps.begin(),
      config_.checkpoint.steps.end(),step_)-config_.checkpoint.steps.begin();
    nextCheckpointTime_=std::upper_bound(config_.checkpoint.times.begin(),
      config_.checkpoint.times.end(),timeBoxSI_)-config_.checkpoint.times.begin();
    logRoot(communicator_,"Checkpoint restart: input reconstruction and initial Poisson solve skipped; "
      "saved geometry, particle order and absolute boost time retained.");
  }

  void Simulation::restoreCheckpointFields()
  {
    YeeComponent* arrays[6]={&fields_->ex(),&fields_->ey(),&fields_->ez(),
                            &fields_->bx(),&fields_->by(),&fields_->bz()};
    const char* names[6]={"ex","ey","ez","bx","by","bz"};
    for (unsigned int a=0;a<6;++a)
      checkpointReader_->readDoubles(std::string("fields/")+names[a],arrays[a]->data(),arrays[a]->size());
    if (config_.boundary.type==EBBoundaryType::Cpml)
      {
        EBConvolutionalPML& pml=dynamic_cast<EBConvolutionalPML&>(fields_->boundary());
        for (std::size_t i=0;i<pml.checkpointHistoryCount();++i)
          {
            std::vector<Double>& h=pml.checkpointHistory(i);
            checkpointReader_->readDoubles("cpml/history_"+std::to_string(i),h.data(),h.size());
          }
      }
    std::vector<unsigned char> bytes=checkpointReader_->readLocalBytes("state/scalars");
    transferCheckpointScalars(bytes,true,false);
  }

  void Simulation::validateCheckpointBranch() const
  {
    int invalid=0;
    if (retirementChangedOnRestart_)
      {
        if (retirementEntryCount_ || !retirementCarriers_.empty()) invalid=1;
        if (config_.particleRetirement.enabled)
          for (const RelativisticParticleSI& p:particles_)
            if (frame_.labZFromBoxZT(p.position[2],timeBoxSI_)>=
                config_.particleRetirement.entranceZ) invalid=1;
      }
    int anyInvalid=0;
    MPI_Allreduce(&invalid,&anyInvalid,1,MPI_INT,MPI_MAX,communicator_);
    if (anyInvalid)
      throw std::runtime_error(
        "Retirement branch must start before any original retirement and before any "
        "particle reaches the new entrance; choose an earlier checkpoint. "
        "Particles may still be inside magnets or held.");
    if (config_.checkpoint.resetDiagnostics)
      {
        invalid=0;
        for (const FieldDetectorPlaneConfig& plane:detectorConfig_.fieldPlanes)
          if (plane.particleBackgroundReference)
            {
              if (localLostParticles_ || retirementEntryCount_ || !pmlCarriers_.empty()) invalid=1;
              for (const RelativisticParticleSI& p:particles_)
                if (frame_.labZFromBoxZT(p.position[2],timeBoxSI_)>=plane.referenceEntranceZ)
                  invalid=1;
            }
        MPI_Allreduce(&invalid,&anyInvalid,1,MPI_INT,MPI_MAX,communicator_);
        if (anyInvalid)
          throw std::runtime_error(
            "New ballistic-reference diagnostics require a checkpoint before particle "
            "losses and before every particle crosses the reference entrance; missing "
            "historical crossings cannot be reconstructed from a later snapshot");
      }
  }

  void Simulation::restoreCheckpointDiagnostics()
  {
    if (detectors_ && !config_.checkpoint.resetDiagnostics)
      detectors_->restoreCheckpointState(checkpointReader_->readLocalBytes("state/detectors"));
    if (config_.checkpoint.resetDiagnostics && config_.trajectory.enabled)
      { nextTrajectorySampleTime_=timeBoxSI_;
        trajectorySamplesSinceFlush_=0; sampleTrajectory(); }
    if (config_.checkpoint.resetDiagnostics && detectors_)
      detectors_->sampleFieldPlanes(*fields_,timeBoxSI_);
    checkpointReader_->close();
    checkpointReader_.reset();
    disarmShutdownWatchdog();
    std::ostringstream message;
    message<<std::setprecision(17)<<"[checkpoint] phase=restore step="<<step_
           <<" time_box_s="<<timeBoxSI_<<" path="<<config_.checkpoint.restartFrom
           <<" diagnostics="<<(config_.checkpoint.resetDiagnostics?"new-segment":"continued")
           <<" retirement_changed="<<(retirementChangedOnRestart_?1:0);
    logRoot(communicator_,message.str());
  }

  void Simulation::maybeWriteCheckpoint()
  {
    bool due=false;
    while (nextCheckpointStep_<config_.checkpoint.steps.size() &&
           config_.checkpoint.steps[nextCheckpointStep_]<=step_)
      { ++nextCheckpointStep_; due=true; }
    while (nextCheckpointTime_<config_.checkpoint.times.size() &&
           config_.checkpoint.times[nextCheckpointTime_]<=timeBoxSI_)
      { ++nextCheckpointTime_; due=true; }
    if (!due) return;
    std::ostringstream filename;
    filename<<"step-"<<std::setw(12)<<std::setfill('0')<<step_<<".h5";
    writeCheckpoint(joinPath(config_.checkpoint.directory,filename.str()));
  }

  void Simulation::writeCheckpoint(const std::string& path)
  {
    const double start=MPI_Wtime();
    requireCheckpointWatchdog(config_.checkpoint.ioTimeoutSeconds);
    if (!boundaryInjectionWorkSamples_.empty())
      throw std::logic_error("Checkpoint requested before boundary work was finalized");
    createDirectories(config_.checkpoint.directory,communicator_);
    if (config_.trajectory.enabled) trajectoryWriter_.flushForCheckpoint();
    if (detectors_) detectors_->flushForCheckpoint();
    if (config_.energyLedger.enabled && rank_==0) energyLedgerWriter_.flushForCheckpoint();
    CheckpointArchive file(communicator_,path,true);
    writeText(file,"metadata/native_abi",nativeAbi());
    writeText(file,"metadata/configuration",config_.configurationText);
    writeText(file,"metadata/source_manifest",absolutePath(runMetadata_.manifestPath));
    writeText(file,"metadata/source_run_id",runMetadata_.runId);
    writeText(file,"metadata/source_revision",runMetadata_.sourceRevision);
    writeText(file,"metadata/parent_checkpoint",config_.checkpoint.restartFrom);
    std::vector<unsigned char> bytes;
    transferCheckpointScalars(bytes,false,true);
    file.writeBytes("state/geometry",bytes.data(),bytes.size());
    transferCheckpointScalars(bytes,false,false);
    file.writeBytes("state/scalars",bytes.data(),bytes.size());
    YeeComponent* arrays[6]={&fields_->ex(),&fields_->ey(),&fields_->ez(),
                            &fields_->bx(),&fields_->by(),&fields_->bz()};
    const char* names[6]={"ex","ey","ez","bx","by","bz"};
    for (unsigned int a=0;a<6;++a)
      file.writeDoubles(std::string("fields/")+names[a],arrays[a]->data(),arrays[a]->size());
    if (config_.boundary.type==EBBoundaryType::Cpml)
      {
        const EBConvolutionalPML& pml=dynamic_cast<const EBConvolutionalPML&>(fields_->boundary());
        for (std::size_t i=0;i<pml.checkpointHistoryCount();++i)
          { const std::vector<Double>& h=pml.checkpointHistory(i);
            file.writeDoubles("cpml/history_"+std::to_string(i),h.data(),h.size()); }
      }
    file.writeParticles("particles",particles_);
    file.writeCarriers("carriers/cpml",pmlCarriers_);
    file.writeCarriers("carriers/retirement",retirementCarriers_);
    if (detectors_) bytes=detectors_->checkpointState(); else bytes.clear();
    file.writeBytes("state/detectors",bytes.data(),bytes.size());
    YAML::Emitter prefix;
    prefix<<YAML::BeginMap<<YAML::Key<<"step"<<YAML::Value<<step_
          <<YAML::Key<<"time_box_s"<<YAML::Value<<YAML::DoublePrecision(17)<<timeBoxSI_
          <<YAML::Key<<"parent_checkpoint"<<YAML::Value<<config_.checkpoint.restartFrom
          <<YAML::Key<<"outputs"<<YAML::Value<<YAML::BeginSeq;
    if (detectors_)
      for (const DetectorCheckpointPrefix& p:detectors_->checkpointPrefixes())
        prefix<<YAML::BeginMap<<YAML::Key<<"file"<<YAML::Value<<absolutePath(p.filename)
              <<YAML::Key<<"group"<<YAML::Value<<p.group
              <<YAML::Key<<"records"<<YAML::Value<<p.committedRecords
              <<YAML::Key<<"validation_records"<<YAML::Value<<p.validationCommittedRecords
              <<YAML::EndMap;
    if (config_.trajectory.enabled)
      {
        std::ostringstream name;
        name<<(config_.trajectory.basename.empty()?"trajectory":config_.trajectory.basename)
            <<"-rank-"<<std::setw(5)<<std::setfill('0')<<rank_<<".h5";
        prefix<<YAML::BeginMap<<YAML::Key<<"file"<<YAML::Value
              <<absolutePath(joinPath(config_.trajectory.directory,name.str()))
              <<YAML::Key<<"group"<<YAML::Value<<"trajectory"
              <<YAML::Key<<"records"<<YAML::Value<<trajectoryWriter_.committedRecords()
              <<YAML::EndMap;
      }
    if (config_.energyLedger.enabled && rank_==0)
      prefix<<YAML::BeginMap<<YAML::Key<<"file"<<YAML::Value
            <<absolutePath(joinPath(config_.energyLedger.directory,config_.energyLedger.filename))
            <<YAML::Key<<"group"<<YAML::Value<<"energy_ledger"
            <<YAML::Key<<"records"<<YAML::Value<<energyLedgerWriter_.committedRecords()
            <<YAML::EndMap;
    prefix<<YAML::EndSeq<<YAML::EndMap;
    writeText(file,"metadata/diagnostic_prefixes",prefix.c_str());
    file.commit();
    const unsigned long long localBytes=file.localTransferredBytes();
    unsigned long long totalBytes=0;
    MPI_Allreduce(&localBytes,&totalBytes,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,communicator_);
    const double elapsed=MPI_Wtime()-start;
    double maximum=0;
    MPI_Allreduce(&elapsed,&maximum,1,MPI_DOUBLE,MPI_MAX,communicator_);
    disarmShutdownWatchdog();
    // Keep checkpoint pauses separate from advancing-step throughput reports.
    if (loopWallStart_>0) loopWallStart_+=elapsed;
    if (performanceWallStart_>0) performanceWallStart_+=elapsed;
    std::ostringstream message;
    message<<std::setprecision(10)<<"[checkpoint] phase=write step="<<step_
           <<" time_box_s="<<timeBoxSI_<<" bytes="<<totalBytes
           <<" seconds="<<maximum<<" gib_per_s="
           <<(maximum>0?totalBytes/(maximum*1073741824.0):0)
           <<" path="<<absolutePath(path);
    logRoot(communicator_,message.str());
  }
}
