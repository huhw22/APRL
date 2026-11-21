#include "simulation.h"
#include "particle_balance.h"
#include "runtime_util.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace aprl { namespace {
// Ownership changes at complete-step boundaries only. Packets carry the same
// native ABI as the existing moving-window particle migration.
template<class T,class Owner>
std::size_t redistribute(std::vector<T>& values,Owner owner,int ranks,MPI_Comm comm) {
  static_assert(std::is_trivially_copyable<T>::value,"MPI balance packets must be trivial");
  std::vector<int> counts(ranks,0),received(ranks,0),starts(ranks,0),offsets(ranks,0);
  int invalid=0;
  for(const auto& value:values){const int peer=owner(value);if(peer<0||peer>=ranks)invalid=1;
    else if(counts[peer]==INT_MAX/static_cast<int>(sizeof(T)))invalid=1;else ++counts[peer];}
  int any=0;MPI_Allreduce(&invalid,&any,1,MPI_INT,MPI_MAX,comm);
  if(any)throw std::overflow_error("Invalid or oversized particle balance destination");
  MPI_Alltoall(counts.data(),1,MPI_INT,received.data(),1,MPI_INT,comm);
  std::size_t incoming=0;
  for(int p=0;p<ranks;++p){if(incoming>static_cast<std::size_t>(INT_MAX)/sizeof(T))invalid=1;
    offsets[p]=static_cast<int>(std::min(incoming,static_cast<std::size_t>(INT_MAX)/sizeof(T)));
    incoming+=received[p];if(p)starts[p]=starts[p-1]+counts[p-1];}
  if(incoming>static_cast<std::size_t>(INT_MAX)/sizeof(T)||values.size()>static_cast<std::size_t>(INT_MAX)/sizeof(T))invalid=1;
  MPI_Allreduce(&invalid,&any,1,MPI_INT,MPI_MAX,comm);
  if(any)throw std::overflow_error("Particle balance packet exceeds MPI int byte range");
  std::vector<T> send(values.size()),next(incoming);auto cursor=starts;
  for(const auto& value:values)send[cursor[owner(value)]++]=value;
  for(int p=0;p<ranks;++p){counts[p]*=sizeof(T);received[p]*=sizeof(T);starts[p]*=sizeof(T);offsets[p]*=sizeof(T);}
  MPI_Alltoallv(send.empty()?nullptr:send.data(),counts.data(),starts.data(),MPI_BYTE,
      next.empty()?nullptr:next.data(),received.data(),offsets.data(),MPI_BYTE,comm);
  values.swap(next);return send.size()*sizeof(T);
}
}

bool Simulation::maybeBalanceParticles() {
  const auto& config=config_.mesh.particleBalance;
  if(!config.enabled || size_==1 || step_%config.intervalSteps!=0)return false;
  const double started=MPI_Wtime();
  const std::size_t cells=globalGeometry_.nz;
  if(cells>static_cast<std::size_t>(INT_MAX))throw std::overflow_error("MPI balance histogram exceeds MPI count");
  std::vector<double> local(cells,0.),work(cells,0.);
  const auto cellOf=[&](const FieldVector<Double>& position){
    const double z=(position[2]-globalOriginBox_[2])/globalGeometry_.dz;
    return z<=0.?std::size_t(0):(z>=cells?cells-1:static_cast<std::size_t>(std::floor(z)));};
  for(const auto& p:particles_)local[cellOf(p.position)]+=p.frozen?config.frozenWeight:1.;
  for(const auto& p:pmlCarriers_)local[cellOf(p.position)]+=config.frozenWeight;
  MPI_Allreduce(local.data(),work.data(),static_cast<int>(cells),MPI_DOUBLE,MPI_SUM,communicator_);
  // A nonzero field cost distributes otherwise empty regions as well. This
  // dimensionless cost affects scheduling only, never any physical equation.
  const double field=config.fieldCellWeight>0.?config.fieldCellWeight:
      static_cast<double>(globalGeometry_.nx)*globalGeometry_.ny/64.;
  for(auto& value:work)value+=field;
  std::vector<std::size_t> current(size_+1);
  for(int p=0;p<size_;++p)current[p]=slabForRank(p).offset;current.back()=cells;
  std::vector<unsigned long long> wire(size_+1);
  if(rank_==0){const auto next=particle_balance::plan(work,size_,config.minimumCells,0,current,config.minimumImprovement);
    std::copy(next.begin(),next.end(),wire.begin());}
  MPI_Bcast(wire.data(),size_+1,MPI_UNSIGNED_LONG_LONG,0,communicator_);
  std::vector<std::size_t> next(wire.begin(),wire.end());
  const auto before=particle_balance::measure(work,current),after=particle_balance::measure(work,next);
  if(next==current){particleBalanceSeconds_+=MPI_Wtime()-started;return false;}
  if(!movingWindow_||!sourceCoupling_)throw std::logic_error("MPI balance outside conservative moving window");
  const auto newOffset=next[rank_],newCells=next[rank_+1]-newOffset;
  EBGridGeometry geometry=globalGeometry_;geometry.nz=newCells;
  std::unique_ptr<EBFieldGrid> replacement(new EBFieldGrid(geometry));
  auto newWindow=movingWindow_->fields->repartition(*fields_,*replacement,newOffset);
  // Release every object holding references to the old native grid before it
  // is destroyed. The coupling CSV and cumulative energy diagnostics persist.
  sourceCoupling_->reconfigurePartition();
  movingWindow_->particleHalo.reset();
  auto previousWindow=std::move(movingWindow_->fields);
  auto previousFields=std::move(fields_);
  movingWindow_->fields=std::move(newWindow);fields_=std::move(replacement);
  previousFields.reset();previousWindow.reset();
  slabCuts_=next;localZOffset_=newOffset;localGeometry_=geometry;
  localOriginBox_=globalOriginBox_;localOriginBox_[2]=std::fma(static_cast<double>(newOffset),geometry.dz,globalOriginBox_[2]);
  unsigned long long boundaryCounts[6];Double boundaryCharges[6];
  particleBoundary_->checkpointCumulative(boundaryCounts,boundaryCharges);
  particleBoundary_.reset(new ParticleOpenBoundary(localGeometry_,localOriginBox_,globalGeometry_,globalOriginBox_,localZOffset_));
  particleBoundary_->restoreCumulative(boundaryCounts,boundaryCharges);
  movingWindow_->particleHalo.reset(new EBParticleFieldHalo(localGeometry_,communicator_));
  if(config_.mesh.movingWindow.inletChargeMatch)
    movingWindow_->inletCompletion.reset(new ChargeMatchedInlet(localGeometry_,config_.boundary.cpml.cells[0],config_.boundary.cpml.cells[1],
      config_.mesh.longitudinalQuintic,config_.mesh.longitudinalNotch,config_.mesh.longitudinalBinomial3,config_.mesh.transverseXBinomial3));
  if(detectors_)detectors_->updateGridPartition(slabCuts_,globalOriginBox_,localOriginBox_);
  const unsigned long long localCounts[2]={static_cast<unsigned long long>(particles_.size()),static_cast<unsigned long long>(pmlCarriers_.size())};
  unsigned long long countsBefore[2]={},countsAfter[2]={};
  MPI_Allreduce(localCounts,countsBefore,2,MPI_UNSIGNED_LONG_LONG,MPI_SUM,communicator_);
  std::size_t bytes=redistribute(particles_,[&](const RelativisticParticleSI& p){return ownerRank(p.position[2]);},size_,communicator_);
  bytes+=redistribute(pmlCarriers_,[&](const ParticleCPMLCarrier& p){return ownerRank(p.position[2]);},size_,communicator_);
  const unsigned long long finalCounts[2]={static_cast<unsigned long long>(particles_.size()),static_cast<unsigned long long>(pmlCarriers_.size())};
  MPI_Allreduce(finalCounts,countsAfter,2,MPI_UNSIGNED_LONG_LONG,MPI_SUM,communicator_);
  if(countsBefore[0]!=countsAfter[0]||countsBefore[1]!=countsAfter[1])throw std::logic_error("MPI balancing changed particle counts");
  std::sort(particles_.begin(),particles_.end(),[](const RelativisticParticleSI& a,const RelativisticParticleSI& b){return a.id<b.id;});
  ++particleBalanceEvents_;const double elapsed=MPI_Wtime()-started;particleBalanceSeconds_+=elapsed;
  unsigned long long localBytes=bytes,totalBytes=0;double slowest=0.;
  MPI_Reduce(&localBytes,&totalBytes,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,communicator_);
  MPI_Reduce(&elapsed,&slowest,1,MPI_DOUBLE,MPI_MAX,0,communicator_);
  if(rank_==0){std::ostringstream message;message<<std::setprecision(10)
    <<"[particle-balance] step="<<step_<<" event="<<particleBalanceEvents_
    <<" modeled_max_before="<<static_cast<double>(before.maximum)<<" modeled_max_after="<<static_cast<double>(after.maximum)
    <<" modeled_imbalance_before="<<static_cast<double>(before.imbalance)<<" modeled_imbalance_after="<<static_cast<double>(after.imbalance)
    <<" repartition_s="<<slowest<<" particle_packet_bytes="<<totalBytes<<" cuts=";
    for(auto cut:next)message<<cut<<':' ;logRoot(communicator_,message.str());}
  return true;
}
}
