#include "source_coupling.h"
#include "longitudinal_shape.h"
#include "pending_particle_queue.h"
#include "eb_deposition.h"
#include "eb_mpi.h"
#include "runtime_util.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace aprl;
using Vec = FieldVector<Double>;
namespace {
YeeComponent& current(EBFieldGrid& f,unsigned a) {
  return a==0?f.jx():(a==1?f.jy():f.jz());
}
void require(bool yes,const char* what) { if(!yes)throw std::runtime_error(what); }
void fillElectric(EBFieldGrid& f,std::size_t globalOffset) {
  YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};
  for(unsigned a=0;a<3;++a)for(std::size_t k=0;k<e[a]->nz();++k)
    for(std::size_t j=0;j<e[a]->ny();++j)for(std::size_t i=0;i<e[a]->nx();++i)
      (*e[a])(i,j,k)=1e5*(.2+.013*i+.017*j+.023*(k+globalOffset)+.07*a);
}
long double meshWork(EBFieldGrid& f) {
  const YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};
  long double out=0.;
  for(unsigned a=0;a<3;++a)for(std::size_t i=0;i<e[a]->size();++i)
    out+=(long double)e[a]->data()[i]*current(f,a).data()[i];
  const auto& g=f.geometry();return out*g.dt*g.dx*g.dy*g.dz;
}
double nodeWeight(double z,long long k) {
  const auto b=longitudinal_shape::basis(z);
  return k<b.first||k>=b.first+b.count?0.:b.v[k-b.first];
}
long double chargeNode(const EBGridGeometry& g,const Vec& origin,const Vec& p,
                       long long i,long long j,long long k,double q) {
  const double x=(p[0]-origin[0])/g.dx,y=(p[1]-origin[1])/g.dy,
    z=(p[2]-origin[2])/g.dz;
  return (long double)q/(g.dx*g.dy*g.dz)*std::max(0.,1.-std::abs(x-i))*
    std::max(0.,1.-std::abs(y-j))*nodeWeight(z,k);
}
long double retainedContinuity(EBFieldGrid& f,const Vec& origin,const Vec& begin,
                               const Vec& end,double q) {
  const auto& g=f.geometry();long double out=0.;
  // Node zero is the open boundary. Every retained physical node is tested.
  for(std::size_t k=1;k<g.nz;++k)for(std::size_t j=1;j<g.ny;++j)
    for(std::size_t i=1;i<g.nx;++i) {
      const long double div=((long double)f.jx()(i,j,k)-f.jx()(i-1,j,k))/g.dx+
        ((long double)f.jy()(i,j,k)-f.jy()(i,j-1,k))/g.dy+
        ((long double)f.jz()(i,j,k)-f.jz()(i,j,k-1))/g.dz;
      const long double r=chargeNode(g,origin,end,i,j,k,q)-
        chargeNode(g,origin,begin,i,j,k,q)+g.dt*div;
      out+=std::abs(r)*g.dx*g.dy*g.dz/std::abs(q);
    }
  return out;
}
double currentError(EBFieldGrid& candidate,EBFieldGrid& reference,std::size_t offset=0) {
  double error=0.,scale=0.;
  for(unsigned a=0;a<3;++a) {
    auto& v=current(candidate,a);auto& r=current(reference,a);
    for(std::size_t k=0;k<v.nz();++k)for(std::size_t j=0;j<v.ny();++j)
      for(std::size_t i=0;i<v.nx();++i) {
        error=std::max(error,std::abs(v(i,j,k)-r(i,j,k+offset)));
        scale=std::max(scale,std::abs(r(i,j,k+offset)));
      }
  }
  return error/std::max(scale,1e-300);
}
struct Metrics {
  unsigned phases=0,pendingSteps=0,negativeControls=0;
  double currentRelative=0.,continuity=0.,oldClipRetained=0.,workRelative=0.;
  bool linearUnchanged=true,trialReplayOnce=true;
};

void admissionCase(MPI_Comm comm,const std::string& output,double startCells,
                   double heldBeta,unsigned caseIndex,Metrics& metrics) {
  int rank=0,size=1;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
  const double dz=8.293564128124092e-6;
  EBGridGeometry g(8,9,16,4e-5,4e-5,dz,dz/SI::c,EBMaxwellSolver::CowanZ),gg=g;
  gg.nz*=size;
  EBFieldGrid f(g),reference(gg),old(g);
  fillElectric(f,rank*g.nz);fillElectric(reference,0);
  Vec origin(0.),local(0.);local[2]=rank*g.nz*dz;
  const double q=-SI::elementaryCharge;
  Vec begin(0.);begin[0]=3.27*g.dx;begin[1]=4.31*g.dy;begin[2]=startCells*dz;
  Vec end(begin);end[2]+=heldBeta*dz;
  Vec entry(begin);entry[2]=0.;
  const std::string dir=output+"/phase-"+std::to_string(caseIndex);
  createDirectories(dir,comm);
  ConservativeSourceCoupling coupling(comm,dir,true);
  coupling.begin(f,local,std::vector<RelativisticParticleSI>());
  // Repeated trial sources must not leave work records or duplicate current.
  for(unsigned trial=0;trial<3;++trial) {
    f.clearCurrent();coupling.beginTrial();
    if(rank==0) {
      coupling.depositFrozenSourceSupport(begin,entry,q,false);
      ChargeConservingCurrentDepositor(f,local,coupling.longitudinalCurrent()).depositSegment(entry,end,q);
    }
    coupling.completeCurrent();
  }
  f.clearCurrent();coupling.setReplay(true);
  if(rank==0) {
    coupling.depositFrozenSourceSupport(begin,entry,q,true);
    ChargeConservingCurrentDepositor(f,local,coupling.longitudinalCurrent()).depositSegment(entry,end,q);
    coupling.recordConstraint(entry,end,q,false);
    longitudinal_shape::Current oldDeposit(old,MPI_COMM_SELF);
    oldDeposit.clear();oldDeposit.deposit(local,entry,end,q);
    metrics.oldClipRetained=std::max(metrics.oldClipRetained,
      (double)retainedContinuity(old,local,begin,end,q));
    ++metrics.negativeControls;
  }
  coupling.completeCurrent();
  EBZSlabHaloExchange exchange(comm);exchange.sumSharedTangentialCurrent(f);
  longitudinal_shape::Current whole(reference,MPI_COMM_SELF);
  whole.clear();whole.deposit(origin,begin,end,q);
  metrics.currentRelative=std::max(metrics.currentRelative,currentError(f,reference,rank*g.nz));
  if(rank==0)metrics.continuity=std::max(metrics.continuity,
      (double)retainedContinuity(f,local,begin,end,q));
  const auto work=coupling.finalizeConstraints();
  const long double mesh=meshWork(f),error=std::abs(mesh+work[0]);
  metrics.workRelative=std::max(metrics.workRelative,
    (double)(error/std::max(std::abs(mesh),1e-100L)));
  metrics.trialReplayOnce=metrics.trialReplayOnce&&(work[2]==(rank==0?2:0));
  coupling.close(1);++metrics.phases;
}

void linearGuard(const std::string& output,Metrics& metrics) {
  const double dz=8.293564128124092e-6;
  EBGridGeometry g(8,9,16,4e-5,4e-5,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
  EBFieldGrid f(g),reference(g);Vec origin(0.);
  createDirectories(output,MPI_COMM_SELF);
  ConservativeSourceCoupling coupling(MPI_COMM_SELF,output,false);
  coupling.begin(f,origin,std::vector<RelativisticParticleSI>());
  Vec begin(0.);begin[0]=3.27*g.dx;begin[1]=4.31*g.dy;begin[2]=-.2*dz;
  Vec entry(begin);entry[2]=0.;Vec end(begin);end[2]+=.51479762393803552*dz;
  bool rejected=false;
  try {coupling.depositFrozenSourceSupport(begin,entry,-SI::elementaryCharge,false);}
  catch(const std::logic_error&) {rejected=true;}
  require(rejected,"Linear mode accepted high-order support API");
  ChargeConservingCurrentDepositor(f,origin).depositSegment(entry,end,-SI::elementaryCharge);
  ChargeConservingCurrentDepositor(reference,origin).depositSegment(entry,end,-SI::elementaryCharge);
  for(unsigned a=0;a<3;++a) {
    auto& v=current(f,a);auto& r=current(reference,a);
    metrics.linearUnchanged=metrics.linearUnchanged&&
      std::memcmp(v.data(),r.data(),v.size()*sizeof(double))==0;
  }
  coupling.close(0);
}

void pendingHandoff(const std::string& output,Metrics& metrics) {
  const double dz=8.293564128124092e-6,beta=.5151304962827177;
  EBGridGeometry g(8,9,16,4e-5,4e-5,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
  PendingWindowGeometry geometry;geometry.dz=dz;geometry.heldDistancePerBeta=dz;
  geometry.referenceOriginZ=0.;geometry.initialOffset=0;geometry.maximumStep=20;
  RelativisticParticleSI initial;initial.position=Vec(0.);
  initial.position[0]=3.27*g.dx;initial.position[1]=4.31*g.dy;
  initial.position[2]=(-8.*(1.+beta)-.215322880587)*dz;
  initial.charge=-SI::elementaryCharge;initial.mass=SI::electronMass;
  initial.frozen=true;initial.frozenBetaZ=beta;
  PendingParticleQueue queue(geometry);std::vector<RelativisticParticleSI> active(1,initial);
  queue.extract(active);require(active.empty()&&queue.nextDueStep()==8,"Pending admission setup differs");
  unsigned admitted=0;
  for(std::size_t step=0;step<=10;++step) {
    admitted+=queue.admit(step,active);unsigned pendingVisited=0;
    EBFieldGrid f(g),reference(g);fillElectric(f,0);fillElectric(reference,0);
    Vec origin(0.);origin[2]=queue.origin(step);
    const std::string dir=output+"/step-"+std::to_string(step);createDirectories(dir,MPI_COMM_SELF);
    ConservativeSourceCoupling coupling(MPI_COMM_SELF,dir,true);
    coupling.begin(f,origin,std::vector<RelativisticParticleSI>());
    queue.visitNear(step,8,[&](const RelativisticParticleSI& p){
      ++pendingVisited;Vec end(p.position);end[2]+=dz*beta;
      coupling.depositFrozenSourceSupport(p.position,end,p.charge,true);
    });
    if(step==8)require(pendingVisited==0&&active.size()==1,"Pending source duplicated at admission");
    for(auto& p:active) {
      const Vec begin(p.position);p.position[2]+=dz*beta;
      if(begin[2]<origin[2]) {
        Vec entry(begin);entry[2]=origin[2];
        coupling.depositFrozenSourceSupport(begin,entry,p.charge,true);
        ChargeConservingCurrentDepositor(f,origin,coupling.longitudinalCurrent()).depositSegment(entry,p.position,p.charge);
        coupling.recordConstraint(entry,p.position,p.charge,false);
      } else {
        ChargeConservingCurrentDepositor(f,origin,coupling.longitudinalCurrent()).depositSegment(begin,p.position,p.charge);
        coupling.recordConstraint(begin,p.position,p.charge,false);
      }
    }
    Vec baseline(initial.position);baseline[2]=repeatHeldAddition(initial.position[2],dz*beta,step);
    Vec baselineEnd(baseline);baselineEnd[2]+=dz*beta;
    longitudinal_shape::Current direct(reference,MPI_COMM_SELF);direct.clear();
    direct.deposit(origin,baseline,baselineEnd,initial.charge);
    coupling.completeCurrent();
    metrics.currentRelative=std::max(metrics.currentRelative,currentError(f,reference));
    metrics.continuity=std::max(metrics.continuity,
      (double)retainedContinuity(f,origin,baseline,baselineEnd,initial.charge));
    const auto work=coupling.finalizeConstraints();const auto mesh=meshWork(f);
    metrics.workRelative=std::max(metrics.workRelative,
      (double)(std::abs(mesh+work[0])/std::max(std::abs(mesh),1e-100L)));
    coupling.close(step+1);++metrics.pendingSteps;
  }
  require(admitted==1&&queue.empty(),"Pending event count changed");
}
}

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);int rank=0,size=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  try {
    require(argc==2,"Need an independent output directory");createDirectories(argv[1],MPI_COMM_WORLD);
    Metrics m;const double beta=.51479762393803552;
    unsigned n=0;for(double phase:{.01,.1,.2,.37,.5})admissionCase(MPI_COMM_WORLD,argv[1],-phase,beta,n++,m);
    admissionCase(MPI_COMM_WORLD,argv[1],-.215322880587,.5151304962827177,n++,m);
    if(rank==0) {linearGuard(std::string(argv[1])+"/linear",m);pendingHandoff(std::string(argv[1])+"/pending",m);}
    double local[4]={m.currentRelative,m.continuity,m.oldClipRetained,m.workRelative},global[4]={};
    MPI_Allreduce(local,global,4,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    require(global[0]<2e-10,"Admission current differs from full trajectory");
    require(global[1]<2e-10,"Admission retained-node continuity failed");
    require(global[2]>.1,"Old clipped-source negative control did not reproduce");
    require(global[3]<2e-10,"Source constraint does not pair with actual mesh work");
    int flags=(m.linearUnchanged&&m.trialReplayOnce)?1:0,all=0;
    MPI_Allreduce(&flags,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);require(all==1,"Linear or trial/replay contract failed");
    if(rank==0)std::cout<<std::setprecision(18)<<"{\"passed\":true,\"physical_acceptance\":false,\"ranks\":"<<size
      <<",\"admission_phases\":"<<m.phases<<",\"pending_steps\":"<<m.pendingSteps
      <<",\"full_current_relative_error\":"<<global[0]<<",\"retained_continuity_L1_over_abs_q\":"<<global[1]
      <<",\"old_clip_retained_L1_over_abs_q\":"<<global[2]<<",\"constraint_mesh_work_relative_error\":"<<global[3]
      <<",\"linear_bitwise_unchanged\":true,\"trial_replay_work_recorded_once\":true}\n";
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
  MPI_Finalize();return 0;
}
