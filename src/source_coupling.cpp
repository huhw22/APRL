#include "source_coupling.h"
#include "source_coupling_math.h"
#include "longitudinal_shape.h"
#include "particle_momentum_probe.h"
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <sstream>
#include <cstdlib>
#ifdef APRL_HAS_OPENMP
#include <omp.h>
#endif
namespace aprl {
using namespace source_coupling_math;
struct ConservativeSourceCoupling::Impl {
 struct Slot {RelativisticParticleSI first,last;Vec guess=Vec(0.),usedE=Vec(0.),usedB=Vec(0.),nextE=Vec(0.),nextB=Vec(0.);bool seen=false,initialized=false,replayed=false,nextReady=false;unsigned pushes=0;Double time=0.,dt=0.;std::vector<RelativisticParticleSI> states;std::vector<BorisWork> works;};
 struct Constraint {Vec begin,end;Double charge;bool carrier;};
 bool quintic=false;std::unique_ptr<longitudinal_shape::Current> shapedCurrent;
 Vec electric(const Vec& a,const Vec& b,bool open=true){
  if(!quintic)return open?openOrbitElectric(*fields,&old,origin,a,b,&newHalo,&oldHalo):orbitElectric(*fields,&old,origin,a,b,&newHalo,&oldHalo);
  const auto path=open?clipOpenOrbit(fields->geometry(),origin,a,b,&newHalo):OpenOrbit();
  const Vec& first=open?path.begin:a;const Vec& last=open?path.end:b;if(open&&!path.valid)throw std::out_of_range("Quintic active orbit outside domain");
  Vec result(0.);const YeeComponent* e[3]={&fields->ex(),&fields->ey(),&fields->ez()};
  longitudinal_shape::weights(fields->geometry(),origin,first,last,[&](unsigned d,long long i,long long j,long long k,double w){const auto& v=*e[d];if(i<0||j<0||i>=static_cast<long long>(v.nx())||j>=static_cast<long long>(v.ny()))return;if((k<0&&!newHalo.hasLower())||(k>=static_cast<long long>(v.nz())&&!newHalo.hasUpper()))return;result[d]+=.5*w*(newHalo.value(v,d,i,j,k)+oldHalo.value(old[d],d,i,j,k));});
  if(open)result*=path.last-path.first;return result;
 }
 EBFieldGrid* fields=nullptr;Vec origin=Vec(0.);Electric old;NativeHalo oldHalo,newHalo,magneticHalo;MPI_Comm comm;int rank=0,size=1;bool replay=false;
 std::unordered_map<std::uint64_t,std::size_t> index;std::vector<Slot> slots;std::vector<Constraint> constraints,activeSegments,longitudinalHandoffs;unsigned cx=0,cy=0;std::size_t completedStep=0;std::ofstream exitAudit;std::string auditPath;
 std::ofstream log;std::string path;unsigned long long trialPushes=0,committedPushes=0;
 explicit Impl(MPI_Comm c,const std::string& dir):comm(c){MPI_Comm_rank(c,&rank);MPI_Comm_size(c,&size);std::ostringstream p;p<<dir<<"/source-coupling-rank-"<<std::setw(5)<<std::setfill('0')<<rank;path=p.str();log.open(path+".csv");if(!log)throw std::runtime_error("Cannot open source-coupling log");log<<std::setprecision(18)<<"step,iterations,self_force_defect,relative_fixed_point_work_defect,trial_pushes,committed_pushes\n";const char* audit=std::getenv("APRL_SOURCE_EXIT_AUDIT");if(audit&&std::string(audit)!="0"){auditPath=dir+"/source-exit-rank-"+p.str().substr(p.str().size()-5);exitAudit.open(auditPath+".csv");if(!exitAudit)throw std::runtime_error("Cannot open source exit audit");exitAudit<<std::setprecision(21)<<"step,active_deposited_chord_work_J,frozen_deposited_chord_work_J,carrier_deposited_chord_work_J,held_global_z_node_work_J,handoff_reclassified_work_J,active_segments\n";}}
 Slot& slot(const RelativisticParticleSI& p){const auto& map=index;auto i=map.find(p.id);if(i!=map.end())return slots[i->second];
#ifdef APRL_HAS_OPENMP
 if(omp_in_parallel())throw std::logic_error("Unseeded parallel source-coupling particle");
#endif
 index.emplace(p.id,slots.size());slots.emplace_back();return slots.back();}
};
ConservativeSourceCoupling::ConservativeSourceCoupling(MPI_Comm c,const std::string& directory,bool quintic):impl_(new Impl(c,directory)){impl_->quintic=quintic;}
ConservativeSourceCoupling::~ConservativeSourceCoupling()=default;
void ConservativeSourceCoupling::begin(EBFieldGrid& f,const Vec& origin,const std::vector<RelativisticParticleSI>& p,unsigned cx,unsigned cy){auto& s=*impl_;s.cx=cx;s.cy=cy;s.fields=&f;s.origin=origin;s.old[0]=f.ex();s.old[1]=f.ey();s.old[2]=f.ez();s.oldHalo.exchange(f,s.comm,false,s.quintic?8:1);if(s.quintic&&!s.shapedCurrent)s.shapedCurrent.reset(new longitudinal_shape::Current(f,s.comm));s.magneticHalo.exchange(f,s.comm,true);s.index.clear();s.slots.clear();s.slots.reserve(p.size()+32);for(const auto& q:p){s.index.emplace(q.id,s.slots.size());s.slots.emplace_back();}s.constraints.clear();s.activeSegments.clear();s.longitudinalHandoffs.clear();s.trialPushes=s.committedPushes=0;s.replay=false;}
longitudinal_shape::Current* ConservativeSourceCoupling::longitudinalCurrent(){return impl_->shapedCurrent.get();}
void ConservativeSourceCoupling::completeCurrent(){if(impl_->shapedCurrent)impl_->shapedCurrent->complete();}
void ConservativeSourceCoupling::beginTrial(){if(impl_->shapedCurrent)impl_->shapedCurrent->clear();for(auto& q:impl_->slots)q.seen=false;}
void ConservativeSourceCoupling::setReplay(bool enabled){impl_->replay=enabled;if(enabled&&impl_->shapedCurrent)impl_->shapedCurrent->clear();}
void ConservativeSourceCoupling::exchangeCandidate(){impl_->newHalo.exchange(*impl_->fields,impl_->comm,false,impl_->quintic?8:1);}
void ConservativeSourceCoupling::swapElectric(){auto& s=*impl_;using std::swap;swap(s.fields->ex(),s.old[0]);swap(s.fields->ey(),s.old[1]);swap(s.fields->ez(),s.old[2]);}
void ConservativeSourceCoupling::applyCurrentSource(Double sign){auto& s=*impl_;auto& f=*s.fields;const auto& g=f.geometry();YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};const YeeComponent* j[3]={&f.jx(),&f.jy(),&f.jz()};
 for(unsigned a=0;a<3;++a){auto* dst=e[a]->data();const auto* current=j[a]->data();const auto nx=e[a]->nx(),ny=e[a]->ny(),nz=e[a]->nz();
#pragma omp parallel for collapse(2) schedule(static)
  for(std::size_t k=0;k<nz;++k)for(std::size_t y=0;y<ny;++y)for(std::size_t x=0;x<nx;++x){const bool applied=a==0?(y>0&&y<g.ny&&!(k==0&&s.rank==0)&&!(k==g.nz&&s.rank+1==s.size)):(a==1?(x>0&&x<g.nx&&!(k==0&&s.rank==0)&&!(k==g.nz&&s.rank+1==s.size)):(x>0&&x<g.nx&&y>0&&y<g.ny));if(applied){const auto n=(k*ny+y)*nx+x;if(current[n]!=0.)dst[n]+=sign*g.dt/SI::epsilon0*current[n];}}
 }
}
void ConservativeSourceCoupling::push(RelativisticParticleSI& p,const SIFieldSourceSet& sources,const BoostFrameTransform& frame,Double time,Double dt,unsigned steps,Double* prescribed,BorisWork* work,ParticleSubstepObserver* observer){auto& s=*impl_;auto& q=s.slot(p);
 if(s.replay){if(!q.seen||q.time!=time||q.dt!=dt||q.states.size()!=steps+1)throw std::logic_error("Missing converged particle replay");for(unsigned a=0;a<3;++a)if(q.first.position[a]!=p.position[a]||q.first.properVelocity[a]!=p.properVelocity[a])throw std::logic_error("Release state changed during coupled replay");}
 else {
  if(!q.initialized){q.guess=p.position;const Double gamma=BoostFrameTransform::gammaFromProperVelocity(p.properVelocity);q.guess.pmv(dt*SI::c/gamma,p.properVelocity);q.initialized=true;}
  q.first=p;q.time=time;q.dt=dt;q.seen=true;++q.pushes;
  if(q.nextReady){q.usedE=q.nextE;q.usedB=q.nextB;}
  else{q.usedE=s.electric(p.position,q.guess);q.usedB=openOrbitMagnetic(*s.fields,s.origin,p.position,q.guess,&s.magneticHalo);}
  q.states.clear();q.works.clear();q.states.reserve(steps+1);q.works.reserve(steps);q.states.push_back(p);const Double h=dt/steps;
  for(unsigned n=0;n<steps;++n){const auto before=p;const Double g0=BoostFrameTransform::gammaFromProperVelocity(p.properVelocity);Vec point(p.position);point.pmv(.5*h*SI::c/g0,p.properVelocity);Vec ee(0.),bb(0.);sources.addPrescribedBox(point,time+(n+.5)*h,frame,ee,bb);Vec e=q.usedE;e+=ee;Vec b=q.usedB;b+=bb;momentum(p.properVelocity,e,b,p.charge,p.mass,h,g0);const Double g1=BoostFrameTransform::gammaFromProperVelocity(p.properVelocity);Vec mean(before.properVelocity);mean+=p.properVelocity;p.position.pmv(h*SI::c/(g0+g1),mean);
   BorisWork w;w.gridWork=chargeWork(p.charge,q.usedE,before.position,p.position);w.prescribedWork=chargeWork(p.charge,ee,before.position,p.position);w.totalDeltaK=(long double)p.mass*SI::c*SI::c*deltaGamma(before.properVelocity,p.properVelocity);w.roundoffResidual=w.totalDeltaK-w.gridWork-w.prescribedWork;w.absoluteWork=std::abs(w.gridWork)+std::abs(w.prescribedWork);q.works.push_back(w);q.states.push_back(p);
  }q.last=p;q.guess=p.position;
 }
 if(s.replay){q.replayed=true;p=q.last;for(unsigned n=0;n<steps;++n){if(work)work->add(q.works[n]);if(prescribed)*prescribed+=static_cast<Double>(q.works[n].prescribedWork);if(observer)observer->observe(q.states[n],q.states[n+1],time+n*(dt/steps),time+(n+1)*(dt/steps),&q.works[n]);}if(ParticleMomentumProbe::selected(time))ParticleMomentumProbe::record(q.first,p,q.usedE,q.usedB,time,dt);}
}
std::array<Double,4> ConservativeSourceCoupling::localResidual(){auto& s=*impl_;Double force=0.,defect=0.,scale=0.,calls=0.;for(auto& q:s.slots)if(q.seen){calls+=1;const Vec e=s.electric(q.first.position,q.last.position);const Vec b=openOrbitMagnetic(*s.fields,s.origin,q.first.position,q.last.position,&s.magneticHalo);q.nextE=e;q.nextB=b;q.nextReady=true;Vec de=e;de-=q.usedE;Vec db=b;db-=q.usedB;force=std::max(force,std::abs(q.first.charge/q.first.mass)*q.dt/SI::c*(std::sqrt(de.norm2())+SI::c*std::sqrt(db.norm2())));defect+=std::abs(chargeWork(q.first.charge,de,q.first.position,q.last.position));for(const auto& w:q.works)scale+=std::abs(w.gridWork);scale+=std::abs(chargeWork(q.first.charge,e,q.first.position,q.last.position));}return {{force,defect,scale,calls}};}
void ConservativeSourceCoupling::recordConstraint(const Vec& a,const Vec& b,Double charge,bool carrier){impl_->constraints.push_back({a,b,charge,carrier});}
void ConservativeSourceCoupling::recordActiveSegment(const Vec& a,const Vec& b,Double charge){auto& s=*impl_;if(s.exitAudit.is_open())s.activeSegments.push_back({a,b,charge,false});}
void ConservativeSourceCoupling::recordLongitudinalHandoff(std::uint64_t id,const Vec& crossing,bool carrier){
 if(!carrier)return;auto& s=*impl_;const auto found=s.index.find(id);if(found==s.index.end())throw std::logic_error("Longitudinal handoff has no coupled orbit");const auto& q=s.slots[found->second];s.longitudinalHandoffs.push_back({crossing,q.last.position,q.first.charge,true});
}
std::array<long double,8> ConservativeSourceCoupling::finalizeConstraints(){auto& s=*impl_;s.newHalo.exchange(*s.fields,s.comm,false,s.quintic?8:1);std::array<long double,8> out{};for(const auto& q:s.constraints){const Vec e=s.electric(q.begin,q.end,false);const auto work=-chargeWork(q.charge,e,q.begin,q.end);const unsigned i=q.carrier?3:0;out[i]+=work;out[i+1]+=std::abs(work);out[i+2]+=1;}long double held=0.;
 for(const auto& q:s.longitudinalHandoffs){const auto path=clipOpenOrbit(s.fields->geometry(),s.origin,q.begin,q.end,&s.newHalo);if(path.valid){const Vec e=s.electric(path.begin,path.end,false);out[7]+=chargeWork(q.charge,e,path.begin,path.end);}}
 long double active=0.;if(s.exitAudit.is_open())for(const auto& q:s.activeSegments){const Vec e=s.electric(q.begin,q.end,false);active+=chargeWork(q.charge,e,q.begin,q.end);}
 const auto& f=*s.fields;const auto& g=f.geometry();const long double v=(long double)g.dt*g.dx*g.dy*g.dz;
 const YeeComponent* field[2]={&f.ex(),&f.ey()};const YeeComponent* current[2]={&f.jx(),&f.jy()};
 for(unsigned face=0;face<2;++face)if((face==0&&s.rank==0)||(face==1&&s.rank+1==s.size)) {
  const auto k=face==0?0:g.nz;
  for(unsigned a=0;a<2;++a)for(std::size_t j=s.cy;j<g.ny-s.cy+(a==0);++j)for(std::size_t i=s.cx;i<g.nx-s.cx+(a==1);++i) {
   long double weight=1.;if(a==0&&(j==s.cy||j==g.ny-s.cy))weight*=.5;if(a==1&&(i==s.cx||i==g.nx-s.cx))weight*=.5;
   held+=v*weight*(*current[a])(i,j,k)*((long double)s.old[a](i,j,k)+(*field[a])(i,j,k))/2;
  }
 }
 out[6]=held;
 if(s.exitAudit.is_open())s.exitAudit<<s.completedStep<<','<<active<<','<<-out[0]<<','<<-out[3]<<','<<held<<','<<out[7]<<','<<s.activeSegments.size()<<'\n';if(s.exitAudit.is_open()&&!s.exitAudit)throw std::runtime_error("Cannot write source exit audit");s.activeSegments.clear();s.longitudinalHandoffs.clear();s.constraints.clear();return out;}
void ConservativeSourceCoupling::finish(std::size_t step,unsigned iterations,Double force,Double work){auto& s=*impl_;s.completedStep=step;unsigned long long tries=0,commits=0;for(const auto& q:s.slots){tries+=q.pushes;commits+=q.replayed?1:0;}s.log<<step<<','<<iterations<<','<<force<<','<<work<<','<<tries<<','<<commits<<'\n';if(!s.log)throw std::runtime_error("Cannot write source coupling step");}
void ConservativeSourceCoupling::close(std::size_t step){auto& s=*impl_;s.log.flush();if(!s.log)throw std::runtime_error("Cannot flush source coupling log");s.log.close();std::ofstream f(s.path+".complete");f<<"1 "<<step<<'\n';if(!f)throw std::runtime_error("Cannot complete source coupling log");if(s.exitAudit.is_open()){s.exitAudit.flush();if(!s.exitAudit)throw std::runtime_error("Cannot flush source exit audit");s.exitAudit.close();std::ofstream mark(s.auditPath+".complete");mark<<"1 "<<step<<'\n';if(!mark)throw std::runtime_error("Cannot complete source exit audit");}}
}
