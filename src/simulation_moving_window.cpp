#include "simulation.h"
#include "runtime_util.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#ifdef APRL_HAS_MOVING_WINDOW
#include "blocked_cowan_reference.h"
#endif

namespace aprl {
SimulationMovingWindowState::~SimulationMovingWindowState() = default;
namespace {
// The same executable/ABI runs on every rank. Byte packets contain no pointers.
template<class T> void exchangeVector(const std::vector<T>& send,
    std::vector<T>& receive, int destination, int source, int tag, MPI_Comm comm) {
  static_assert(std::is_trivially_copyable<T>::value, "MPI packet must be trivial");
  unsigned long long n=send.size(), m=0;
  MPI_Sendrecv(&n,1,MPI_UNSIGNED_LONG_LONG,destination,tag,
      &m,1,MPI_UNSIGNED_LONG_LONG,source,tag,comm,MPI_STATUS_IGNORE);
  if(n>INT_MAX/sizeof(T)||m>INT_MAX/sizeof(T))
    throw std::overflow_error("Moving-window adjacent particle packet exceeds MPI count");
  receive.resize(static_cast<std::size_t>(m));
  MPI_Sendrecv(send.empty()?nullptr:send.data(),static_cast<int>(n*sizeof(T)),MPI_BYTE,destination,tag+1,
      receive.empty()?nullptr:receive.data(),static_cast<int>(m*sizeof(T)),MPI_BYTE,source,tag+1,
      comm,MPI_STATUS_IGNORE);
}
template<class T,class Owner> void migrateWindowObjects(std::vector<T>& values,
    Owner owner, int rank, int ranks, int tag, MPI_Comm comm) {
  std::vector<T> lower,upper,incoming;
  std::size_t retained=0;
  for(const T& p:values) {
    const int target=owner(p);
    if(target==rank) values[retained++]=p;
    else if(target==rank-1&&target>=0) lower.push_back(p);
    else if(target==rank+1&&target<ranks) upper.push_back(p);
    else throw std::runtime_error("Moving-window remap particle crossed a non-adjacent slab or outer boundary");
  }
  values.resize(retained);
  const int lo=rank?rank-1:MPI_PROC_NULL,hi=rank+1<ranks?rank+1:MPI_PROC_NULL;
  exchangeVector(upper,incoming,hi,lo,tag,comm);
  values.insert(values.end(),incoming.begin(),incoming.end());
  exchangeVector(lower,incoming,lo,hi,tag+2,comm);
  values.insert(values.end(),incoming.begin(),incoming.end());
}
}

void Simulation::initializeMovingWindow() {
#ifndef APRL_HAS_MOVING_WINDOW
  throw std::runtime_error("This executable was built without APRL_ENABLE_MOVING_WINDOW");
#else
  movingWindow_.reset(new SimulationMovingWindowState);
  auto& w=*movingWindow_;
  w.referenceGeometry=globalGeometry_;w.referenceOrigin=globalOriginBox_;
  for(unsigned a=0;a<3;++a)w.guards[a]=config_.boundary.cpml.cells[a];
  if(config_.mesh.movingWindow.referenceThroughCpml) {
    w.guards[0]=w.guards[1]=0;
    logRoot(communicator_,
      "Reference transverse boundary: extends through CPML to outer termination; active field geometry and CPML update unchanged.");
  }
  long double local[4]={},global[4]={};
  for(const auto& p:particles_) {
    const double beta=frame_.heldLongitudinalBetaZ(p.properVelocity);
    if(!std::isfinite(beta)||beta<0.||beta>=1.)
      throw std::runtime_error("Moving-window first integration requires nonnegative subluminal held longitudinal velocities");
    local[0]+=p.mass; local[1]+=static_cast<long double>(p.mass)*beta;
    local[2]+=beta;local[3]+=1.;
  }
  MPI_Allreduce(local,global,4,MPI_LONG_DOUBLE,MPI_SUM,communicator_);
  if(!(global[0]>0.)||!(global[3]>0.))throw std::runtime_error("Moving reference requires particles with positive mass");
  w.solverBeta=static_cast<double>(global[1]/global[0]);
  w.translationBeta=static_cast<double>(global[2]/global[3]);
  const auto& cfg=config_.mesh.movingWindow;
  const double release=frame_.boxZFromLabZAndBoxT(config_.reference.activationPlaneZ,0.);
  const double last=frame_.boxZFromLabZAndBoxT(lastBeamlineInteractionExitLab(),0.);
  const double durationSteps=std::ceil(totalTimeBoxSI_/globalGeometry_.dt)+1.;
  const double startCell=std::floor((release-w.referenceOrigin[2])/globalGeometry_.dz)-cfg.lowerGuardCells;
  const double minimumCells=std::ceil((last-release)/globalGeometry_.dz)+cfg.lowerGuardCells+cfg.upperGuardCells+2.;
  if(!std::isfinite(startCell)||!std::isfinite(minimumCells)||!std::isfinite(durationSteps)||
      startCell<1.||minimumCells<8.||durationSteps>std::numeric_limits<std::ptrdiff_t>::max()/4)
    throw std::runtime_error("Moving-window release/last diagnostic/reference geometry is invalid");
  const std::size_t cells=cfg.cells?cfg.cells:static_cast<std::size_t>(minimumCells);
  if(cells<minimumCells||cells<8*static_cast<std::size_t>(size_))
    throw std::runtime_error("Moving-window cells must cover release, all interactions/detectors, 64/128 guards and at least eight cells per rank");
  if(cfg.lowerGuardCells<4||cfg.upperGuardCells<8)
    throw std::runtime_error("Moving-window particle stencils require lower/upper guards of at least 4/8 cells");
  w.initialOffset=w.offset=static_cast<std::ptrdiff_t>(startCell);
  // Exact range needed by the tested translated finite-domain reference. Never
  // silently extrapolate, insert zeros, or enlarge the simulation field grid.
  const double lastQuery=startCell-(1.+w.translationBeta)*durationSteps-1.;
  const double maximumQuery=std::max(startCell+cells,startCell-1.-w.translationBeta*durationSteps);
  if(std::min(startCell-1.,lastQuery)<0.||maximumQuery>static_cast<double>(w.referenceGeometry.nz-1)) {
    std::ostringstream error;error<<"Prepared reference domain does not cover moving inlet through mesh.duration: required native z range ["
      <<std::min(startCell-1.,lastQuery)<<", "<<maximumQuery<<"], available [0, "<<w.referenceGeometry.nz-1
      <<"]. Set the explicit reference/source mesh from this range; the moving E/B window remains compact.";
    throw std::runtime_error(error.str());
  }
  const std::size_t nativeCollar=config_.mesh.wideLongitudinalSource()?8:4;
  if(cfg.particleCollarCells&&cfg.particleCollarCells<nativeCollar)
    throw std::runtime_error("Explicit particle collar is smaller than the source support");
  const double rightMargin=static_cast<double>(cfg.particleCollarCells?cfg.particleCollarCells:nativeCollar)+std::ceil((1.-frame_.beta())*durationSteps);
  const double windowOrigin=std::fma(startCell,globalGeometry_.dz,w.referenceOrigin[2]);
  const double safeExit=std::fma(static_cast<double>(cells)-rightMargin,globalGeometry_.dz,windowOrigin);
  w.exitLabZ=frame_.labZFromBoxZT(safeExit,0.);
  if(!(safeExit>last+2.*globalGeometry_.dz))
    throw std::runtime_error("Moving-window right causal particle collar overlaps a requested interaction/detector");
  for(const auto& p:particles_) {
    if(p.position[2]>=safeExit||(!p.frozen&&p.position[2]<windowOrigin+2.*globalGeometry_.dz))
      throw std::runtime_error("Initial active particle is outside the moving-window trusted region");
  }
  // Only source particles gather. No full-domain rho/phi/E/B allocation occurs
  // in Simulation; the reference generator streams its scalar blocks to disk.
  std::vector<unsigned long long> counts(static_cast<std::size_t>(size_));
  const unsigned long long localCount=particles_.size();
  MPI_Gather(&localCount,1,MPI_UNSIGNED_LONG_LONG,counts.data(),1,MPI_UNSIGNED_LONG_LONG,0,communicator_);
  std::vector<RelativisticParticleSI> source;
  static_assert(std::is_trivially_copyable<RelativisticParticleSI>::value,"Particle packet must be trivial");
  const std::size_t batch=static_cast<std::size_t>(64*1024*1024)/sizeof(RelativisticParticleSI);
  if(rank_==0) {
    std::size_t total=0;for(auto count:counts){if(count>source.max_size()-total)throw std::overflow_error("Reference particle gather overflow");total+=count;}
    source.resize(total);std::copy(particles_.begin(),particles_.end(),source.begin());
    std::size_t offset=counts[0];
    for(int peer=1;peer<size_;++peer) {
      for(std::size_t n=0;n<counts[peer];n+=batch) {
        const auto take=std::min(batch,static_cast<std::size_t>(counts[peer])-n);
        MPI_Recv(source.data()+offset+n,static_cast<int>(take*sizeof(RelativisticParticleSI)),MPI_BYTE,peer,810,communicator_,MPI_STATUS_IGNORE);
      }offset+=counts[peer];
    }
    // Stable identity order makes the reference independent of slab ownership.
    std::sort(source.begin(),source.end(),[](const RelativisticParticleSI& a,const RelativisticParticleSI& b){return a.id<b.id;});
    w.reference.reset(new verification::BlockedCowanReference(w.referenceGeometry,w.referenceOrigin,
        source,w.solverBeta,w.guards,cfg.blockLayers,cfg.referenceDirectory,MPI_COMM_SELF,config_.mesh.longitudinalQuintic,config_.mesh.longitudinalNotch));
    w.reference->writeReport(cfg.referenceDirectory+"/preparation.json");
  } else {
    for(std::size_t n=0;n<particles_.size();n+=batch) {
      const auto take=std::min(batch,particles_.size()-n);
      MPI_Send(particles_.data()+n,static_cast<int>(take*sizeof(RelativisticParticleSI)),MPI_BYTE,0,810,communicator_);
    }
  }
  MPI_Barrier(communicator_);
  if(rank_!=0)w.reference=verification::BlockedCowanReference::openPrepared(w.referenceGeometry,w.referenceOrigin,
      w.solverBeta,w.guards,cfg.blockLayers,cfg.referenceDirectory,MPI_COMM_SELF,config_.mesh.longitudinalQuintic,config_.mesh.longitudinalNotch);
  if(cfg.referenceSurfaceCurrent) {
    const std::string currentFile=cfg.referenceDirectory+"/transverse-current.bin";
    if(rank_==0)TransverseReferenceCurrent::prepare(currentFile,*w.reference,w.referenceGeometry,w.guards,w.solverBeta);
    MPI_Barrier(communicator_);
    std::ostringstream workName;
    if(config_.energyLedger.enabled)workName<<cfg.referenceDirectory<<"/transverse-source-work-rank-"<<std::setw(5)<<std::setfill('0')<<rank_<<".csv";
    w.referenceCurrent.reset(new TransverseReferenceCurrent(currentFile,w.referenceGeometry,w.guards,w.solverBeta,workName.str()));
    logRoot(communicator_,"Reference transverse surface current: prescribed from reference B at four interfaces; outgoing perturbations keep original CPML; work is included in boundary_source_work_J.");
  }
  config_.mesh.cells[2]=cells;config_.mesh.extent[2]=cells*globalGeometry_.dz;
  config_.mesh.center[2]=std::fma(.5,config_.mesh.extent[2],windowOrigin);
  initializeGeometry();
  // Keep world-to-lab frame origins fixed. Integer reference indices, not
  // repeated additions, determine the translated field-grid origin.
  globalOriginBox_=w.referenceOrigin;globalOriginBox_[2]=windowOrigin;
  localOriginBox_=globalOriginBox_;localOriginBox_[2]=std::fma(static_cast<double>(localZOffset_),globalGeometry_.dz,windowOrigin);
  particleCPML_.reset(new ParticleCPMLRegion(globalGeometry_,globalOriginBox_,config_.boundary.cpml));
  // The compact window assigns the complete pending upstream source to rank
  // zero. That destination can exceed MPI_Alltoallv's int byte displacements
  // even when every original slab fit. Reuse the already gathered ID-ordered
  // source and scatter bounded packets without another all-particle gather.
  std::vector<unsigned long long> targetCounts(static_cast<std::size_t>(size_),0);
  std::vector<std::vector<RelativisticParticleSI> > buckets;
  if(rank_==0) {
    buckets.resize(static_cast<std::size_t>(size_));
    for(const auto& particle:source) {
      const int target=ownerRank(particle.position[2]);
      if(target<0||target>=size_)
        throw std::out_of_range("Reference source is outside the compact particle ownership range");
      ++targetCounts[static_cast<std::size_t>(target)];
    }
    for(int peer=0;peer<size_;++peer) {
      auto& bucket=buckets[static_cast<std::size_t>(peer)];
      const auto count=targetCounts[static_cast<std::size_t>(peer)];
      if(count>bucket.max_size())
        throw std::overflow_error("Compact particle ownership count exceeds vector capacity");
      bucket.resize(static_cast<std::size_t>(count));
    }
    std::vector<std::size_t> cursor(static_cast<std::size_t>(size_),0);
    for(const auto& particle:source) {
      const auto target=static_cast<std::size_t>(ownerRank(particle.position[2]));
      buckets[target][cursor[target]++]=particle;
    }
    std::vector<RelativisticParticleSI>().swap(source);
  }
  unsigned long long ownedCount=0;
  MPI_Scatter(targetCounts.data(),1,MPI_UNSIGNED_LONG_LONG,
      &ownedCount,1,MPI_UNSIGNED_LONG_LONG,0,communicator_);
  std::vector<RelativisticParticleSI>().swap(particles_);
  if(ownedCount>particles_.max_size())
    throw std::overflow_error("Compact particle receive count exceeds vector capacity");
  if(rank_==0) {
    particles_.swap(buckets[0]);
    for(int peer=1;peer<size_;++peer) {
      auto& bucket=buckets[static_cast<std::size_t>(peer)];
      for(std::size_t n=0;n<bucket.size();n+=batch) {
        const auto take=std::min(batch,bucket.size()-n);
        MPI_Send(bucket.data()+n,static_cast<int>(take*sizeof(RelativisticParticleSI)),
            MPI_BYTE,peer,812,communicator_);
      }
      std::vector<RelativisticParticleSI>().swap(bucket);
    }
  } else {
    particles_.resize(static_cast<std::size_t>(ownedCount));
    for(std::size_t n=0;n<particles_.size();n+=batch) {
      const auto take=std::min(batch,particles_.size()-n);
      MPI_Recv(particles_.data()+n,static_cast<int>(take*sizeof(RelativisticParticleSI)),
          MPI_BYTE,0,812,communicator_,MPI_STATUS_IGNORE);
    }
  }
  // Keep diagnostics upstream of the physical inlet on the original path:
  // they can observe held trajectories before those trajectories enter PIC.
  bool schedulePending = true;
  for (const auto& plane : detectorConfig_.particlePlanes)
    if (frame_.boxZFromLabZAndBoxT(plane.z, 0.) < windowOrigin)
      schedulePending = false;
  for (const auto& plane : detectorConfig_.energyPlanes)
    if (frame_.boxZFromLabZAndBoxT(plane.z, 0.) < windowOrigin)
      schedulePending = false;
  for (const auto& plane : detectorConfig_.fieldPlanes)
    if (plane.particleBackgroundReference) schedulePending = false;
  if (rank_ == 0 && schedulePending) {
    PendingWindowGeometry geometry;
    geometry.referenceOriginZ = w.referenceOrigin[2];
    geometry.initialOffset = w.initialOffset;
    geometry.dz = globalGeometry_.dz;
    geometry.heldDistancePerBeta = SI::c * globalGeometry_.dt;
    geometry.maximumStep = static_cast<std::size_t>(durationSteps);
    w.pending.reset(new PendingParticleQueue(geometry));
    w.pending->extract(particles_);
    // Do not leave a second full-bunch capacity attached to the active list.
    std::vector<RelativisticParticleSI>(particles_).swap(particles_);
    if (config_.energyLedger.enabled)
      w.pendingStatistics.reset(new PendingParticleStatistics(
        *w.pending, frame_, geometry.heldDistancePerBeta));
    std::ostringstream message;
    message << "Pending admission queue: waiting=" << w.pending->size()
      << " resident=" << particles_.size()
      << "; exact repeated-addition held coordinates at admission; no per-step waiting-particle scan; "
         "pending axial ledger moments use an affine diagnostic trajectory.";
    logRoot(communicator_, message.str());
  } else if (rank_ == 0) {
    logRoot(communicator_, "Pending admission queue disabled: upstream particle/reference diagnostics require the original held traversal.");
  }
  hasLongitudinalLayoutReport_=false; // The old full/source-domain report is not the moving mesh.
  std::ostringstream name;name<<cfg.referenceDirectory<<"/window-rank-"<<rank_<<".csv";
  w.ledger.open(name.str());if(!w.ledger)throw std::runtime_error("Cannot create moving-volume ledger");
  w.ledger<<"step,origin_box_z_m,energy_measured,inserted_field_J,cropped_field_J,particle_exits,exit_charge_C,cropped_carriers,cropped_charge_C,shift_s,migration_s\n"<<std::setprecision(17);
  if(rank_==0) {
    std::ostringstream message;message<<std::setprecision(17)<<"Moving window: reference_nz="<<w.referenceGeometry.nz
      <<" field_nz="<<cells<<" initial_offset="<<w.initialOffset<<" ranks="<<size_<<" motion=-c; fixed_lab_particle_exit_m="<<w.exitLabZ
      <<" reference_solver_beta="<<w.solverBeta<<" translation_beta="<<w.translationBeta
      <<"; all source particles retained, source preparation on rank 0, per-rank bounded disk readers; no full-domain E/B donor.";
    logRoot(communicator_,message.str());
    if(config_.energyLedger.enabled)logRoot(communicator_,"Moving-volume energy: HDF5 accounts same-position six-face propagation, inlet overwrite and remap (including interface corrections); longitudinal_boundary_outward_J combines longitudinal exchanges. The residual remains a numerical/source convergence diagnostic, not a pure-radiation estimate.");
  }
#endif
}

std::size_t Simulation::pendingParticleCount() const {
  return movingWindow_ && movingWindow_->pending ? movingWindow_->pending->size() : 0;
}
std::size_t Simulation::storedParticleCount() const {
  return particles_.size() + pendingParticleCount();
}
void Simulation::admitMovingWindowParticles() {
  if (movingWindow_->pending)
    movingWindow_->admittedParticles += movingWindow_->pending->admit(step_, particles_);
}

Double Simulation::movingWindowFieldValue(void* context,unsigned a,std::size_t i,std::size_t j,
    std::ptrdiff_t k,std::size_t step) {
#ifdef APRL_HAS_MOVING_WINDOW
  auto& w=*static_cast<SimulationMovingWindowState*>(context);
  const double z=static_cast<double>(k)-w.translationBeta*step;
  const std::size_t count=w.referenceGeometry.nz+((a==0||a==1||a==5)?1:0);
  if(a>=6||!std::isfinite(z)||z<0.||z>count-1)throw std::out_of_range("Moving reference query outside prepared native field");
  const auto lo=static_cast<std::size_t>(std::floor(z)),hi=std::min(lo+1,count-1);
  const double f=z-lo,lower=w.reference->field(a,i,j,lo);
  return f==0.?lower:(1.-f)*lower+f*w.reference->field(a,i,j,hi);
#else
  (void)context;(void)a;(void)i;(void)j;(void)k;(void)step;
  throw std::logic_error("Moving reference not compiled");
#endif
}
void Simulation::installMovingWindowBoundary() {
  MovingWindowReference ref;ref.context=movingWindow_.get();ref.value=&Simulation::movingWindowFieldValue;
  movingWindow_->fields.reset(new MovingWindowFields(localGeometry_,globalGeometry_.nz,localZOffset_,communicator_,config_.boundary.cpml,ref));
  movingWindow_->fields->installBoundary(*fields_);
  movingWindow_->fields->enableEnergyLedger(config_.energyLedger.enabled);
  if(config_.mesh.movingWindow.inletChargeMatch) {
    movingWindow_->inletCompletion.reset(new ChargeMatchedInlet(localGeometry_,config_.boundary.cpml.cells[0],config_.boundary.cpml.cells[1],config_.mesh.longitudinalQuintic,config_.mesh.longitudinalNotch));
    movingWindow_->fields->setNewLayerCompletion(&Simulation::completeMovingInletLayer,this);
    if(rank_==0&&config_.energyLedger.enabled) {
      movingWindow_->inletCompletionLedger.open(config_.mesh.movingWindow.referenceDirectory+"/inlet-charge-match.csv");
      if(!movingWindow_->inletCompletionLedger)throw std::runtime_error("Cannot create inlet charge-match ledger");
      movingWindow_->inletCompletionLedger<<std::setprecision(21)<<"step,source_particles,max_gauss_before_C_per_m3,max_gauss_after_C_per_m3,max_delta_Ez_V_per_m,native_layer_energy_change_J\n";
    }
    logRoot(communicator_,"Inlet charge completion: approximate background retained; new normal-E layer matches actual CIC charge and adjacent evolved fields; all exchange included in remap energy.");
  }
  movingWindow_->particleHalo.reset(new EBParticleFieldHalo(localGeometry_,communicator_));
}
void Simulation::completeMovingInletLayer(void* context,EBFieldGrid& f,std::size_t step,std::ptrdiff_t offset) {
  auto& sim=*static_cast<Simulation*>(context);auto& w=*sim.movingWindow_;
  auto origin=w.referenceOrigin;origin[2]=std::fma(static_cast<double>(offset),f.geometry().dz,w.referenceOrigin[2]);
  w.inletCompletion->begin(origin);
  for(const auto& p:sim.particles_)w.inletCompletion->particle(p);
  if(w.pending){const auto visit=[&](const RelativisticParticleSI& p){w.inletCompletion->particle(p);};if(sim.config_.mesh.wideLongitudinalSource())w.pending->visitNear(step,8,visit);else w.pending->visitDue(step,visit);}
  const auto report=w.inletCompletion->complete(f,sim.config_.energyLedger.enabled);
  if(w.inletCompletionLedger.is_open()) {
    w.inletCompletionLedger<<step<<','<<report.particles<<','<<report.maxBefore<<','<<report.maxAfter<<','<<report.maxDeltaEz<<','<<report.nativeWorkJ<<'\n';
    if(!w.inletCompletionLedger)throw std::runtime_error("Cannot write inlet charge-match ledger");
  }
}
void Simulation::initializeMovingWindowFields() {
#ifdef APRL_HAS_MOVING_WINDOW
  movingWindow_->reference->fill(*fields_,static_cast<std::size_t>(movingWindow_->offset)+localZOffset_);
  movingWindow_->reference->clearCache();
#else
  throw std::logic_error("Moving reference not compiled");
#endif
}
EBElectricHaloView Simulation::beginMovingWindowStep() {
  auto& w=*movingWindow_;w.fields->beginStep(*fields_,step_,w.offset);
  return w.fields->physicalElectricHalos(*fields_,step_);
}
void Simulation::finishMovingWindowStep() {
  auto& w=*movingWindow_;
  w.offset=w.initialOffset-static_cast<std::ptrdiff_t>(step_);
  const auto report=w.fields->shift(*fields_,step_,w.offset);
  w.shiftSeconds+=report.totalSeconds;w.movedBytes+=report.fieldMovedBytes+report.historyMovedBytes;
  w.communicatedBytes+=report.sentBytes+report.receivedBytes;
  if(config_.energyLedger.enabled) {
    w.insertedFieldEnergy=w.fields->energyLedger().insertedEnergyJ;
    w.croppedFieldEnergy=w.fields->energyLedger().croppedEnergyJ;
  }
  globalOriginBox_[2]=std::fma(static_cast<double>(w.offset),globalGeometry_.dz,w.referenceOrigin[2]);
  localOriginBox_[2]=std::fma(static_cast<double>(localZOffset_),globalGeometry_.dz,globalOriginBox_[2]);
  unsigned long long counts[6];Double charges[6];particleBoundary_->checkpointCumulative(counts,charges);
  particleBoundary_.reset(new ParticleOpenBoundary(localGeometry_,localOriginBox_,globalGeometry_,globalOriginBox_,localZOffset_));
  particleBoundary_->restoreCumulative(counts,charges);
  particleCPML_.reset(new ParticleCPMLRegion(globalGeometry_,globalOriginBox_,config_.boundary.cpml));
  const double start=MPI_Wtime();
  migrateWindowObjects(particles_,[&](const RelativisticParticleSI& p){
    if(p.position[2]<globalOriginBox_[2]&&!p.frozen)throw std::runtime_error("Active particle fell behind moving inlet");
    return ownerRank(p.position[2]);},rank_,size_,820,communicator_);
  const double upper=std::fma(static_cast<double>(globalGeometry_.nz),globalGeometry_.dz,globalOriginBox_[2]);
  std::size_t retained=0;
  for(const auto& p:pmlCarriers_) {
    if(p.position[2]>=upper){++w.croppedCarriers;w.croppedCharge+=p.charge*p.currentWeight;}
    else if(p.position[2]<globalOriginBox_[2])throw std::runtime_error("Current carrier fell behind moving inlet");
    else pmlCarriers_[retained++]=p;
  }pmlCarriers_.resize(retained);
  migrateWindowObjects(pmlCarriers_,[&](const ParticleCPMLCarrier& p){return ownerRank(p.position[2]);},rank_,size_,824,communicator_);
  if(!retirementCarriers_.empty())throw std::logic_error("Retirement carriers are unsupported by moving window");
  w.migrationSeconds+=MPI_Wtime()-start;
  if(detectors_)detectors_->updateGridOrigins(globalOriginBox_,localOriginBox_);
  // Intermediate rows inherit the cadence of an explicitly enabled energy
  // ledger. With it disabled, only the final row is written, marked unmeasured.
  if(config_.energyLedger.enabled&&
      step_%std::max<std::size_t>(1,config_.energyLedger.sampleIntervalSteps)==0) {
    w.ledger<<step_<<','<<globalOriginBox_[2]<<",1,"<<w.insertedFieldEnergy<<','<<w.croppedFieldEnergy<<','
      <<w.exits<<','<<w.exitCharge<<','<<w.croppedCarriers<<','<<w.croppedCharge<<','<<w.shiftSeconds<<','<<w.migrationSeconds<<'\n';
  }
}
bool Simulation::firstMovingWindowEntry(const FieldVector<Double>& start,const FieldVector<Double>& end,
    Double startFraction,ParticleBoundaryHit& hit) const {
  bool present=particleCPML_->firstEntry(start,end,hit);
  if(!movingWindow_)return present;
  const auto& w=*movingWindow_;
  const double a=frame_.labZFromBoxZT(start[2],timeBoxSI_+startFraction*globalGeometry_.dt)-w.exitLabZ;
  const double b=frame_.labZFromBoxZT(end[2],timeBoxSI_+globalGeometry_.dt)-w.exitLabZ;
  if(a>=0.)throw std::runtime_error("Active particle started beyond trusted moving-window exit");
  if(b>=0.&&b>a) {
    const double fraction=-a/(b-a);
    if(!present||fraction<hit.fraction) {
      for(unsigned d=0;d<3;++d)hit.position[d]=start[d]+fraction*(end[d]-start[d]);
      hit.fraction=fraction;hit.face=ParticleBoundaryFace::UpperZ;present=true;
    }
  }return present;
}
void Simulation::recordMovingWindowExit(Double charge) {
  ++movingWindow_->exits;movingWindow_->exitCharge+=charge;
}
void Simulation::applyMovingTransverseReferenceCurrent() {
#ifdef APRL_HAS_MOVING_WINDOW
  auto& w=*movingWindow_;
  w.referenceCurrent->apply(*fields_,w.offset+static_cast<std::ptrdiff_t>(localZOffset_),step_,w.translationBeta,config_.energyLedger.enabled);
#endif
}
void Simulation::finishMovingTransverseReferenceCurrent() {
#ifdef APRL_HAS_MOVING_WINDOW
  auto& w=*movingWindow_;const auto work=w.referenceCurrent->finish(*fields_);
  energyLedgerBoundarySourceWork_+=work.signedJ;
  energyLedgerBoundarySourceAbsWork_+=work.absoluteJ;
#endif
}
void Simulation::finalizeMovingWindow() {
  auto& w=*movingWindow_;
  if(w.inletCompletionLedger.is_open()) {w.inletCompletionLedger.flush();if(!w.inletCompletionLedger)throw std::runtime_error("Cannot flush inlet charge-match ledger");w.inletCompletionLedger.close();}
#ifdef APRL_HAS_MOVING_WINDOW
  if(w.referenceCurrent)w.referenceCurrent->close();
#endif
  w.ledger<<step_<<','<<globalOriginBox_[2]<<','<<(config_.energyLedger.enabled?1:0)<<','
    <<w.insertedFieldEnergy<<','<<w.croppedFieldEnergy<<','
    <<w.exits<<','<<w.exitCharge<<','<<w.croppedCarriers<<','<<w.croppedCharge<<','<<w.shiftSeconds<<','<<w.migrationSeconds<<'\n';
  w.ledger.close();if(!w.ledger)throw std::runtime_error("Moving-volume ledger close failed");
  unsigned long long local[5]={w.exits,w.croppedCarriers,w.movedBytes,w.communicatedBytes,
    static_cast<unsigned long long>(w.particleHalo->memoryBytes())},global[5]={};
  MPI_Reduce(local,global,5,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,communicator_);
  double timing[3]={w.shiftSeconds,w.migrationSeconds,w.particleHalo->communicationSeconds()},maxima[3]={};
  MPI_Reduce(timing,maxima,3,MPI_DOUBLE,MPI_MAX,0,communicator_);
  if(rank_==0){std::ostringstream message;message<<std::setprecision(10)<<"Moving window completed: steps="<<step_
    <<" particle_exits="<<global[0]<<" cropped_carriers="<<global[1]<<" moved_bytes="<<global[2]<<" shift_communicated_bytes="<<global[3]
    <<" max_shift_s="<<maxima[0]<<" max_particle_remap_s="<<maxima[1]
    <<" sampling_halo_bytes_global="<<global[4]<<" max_sampling_halo_communication_s="<<maxima[2];logRoot(communicator_,message.str());}
}
}
