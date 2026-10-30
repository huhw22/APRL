#include "charge_matched_inlet.h"
#include "pending_particle_queue.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace aprl;
void require(bool p,const char* message){if(!p)throw std::runtime_error(message);}
bool same(const FieldVector<Double>& a,const FieldVector<Double>& b){for(unsigned i=0;i<3;++i)if(a[i]!=b[i])return false;return true;}
int main(){try {
  PendingWindowGeometry pg;pg.referenceOriginZ=0.;pg.initialOffset=100;pg.dz=1.;pg.heldDistancePerBeta=1.;pg.maximumStep=8;
  RelativisticParticleSI incoming;incoming.position=FieldVector<Double>(0.);incoming.position[2]=99.;incoming.frozen=true;incoming.frozenBetaZ=.25;incoming.sourceId=3;
  std::vector<RelativisticParticleSI> active(1,incoming);PendingParticleQueue queue(pg);queue.extract(active);require(active.empty(),"Fixture must wait upstream");
  const auto count=queue.size(),due=queue.nextDueStep();RelativisticParticleSI peek;std::size_t visits=0;
  queue.visitDue(1,[&](const RelativisticParticleSI& p){peek=p;++visits;});require(visits==1&&queue.size()==count&&queue.nextDueStep()==due&&queue.admittedCount()==0,"Read-only inlet view changed queue");
  queue.admit(1,active);require(active.size()==1&&active[0].position[2]==peek.position[2]&&same(active[0].properVelocity,peek.properVelocity),"Inlet view does not match actual admission");
  EBGridGeometry g(12,12,8,1.,1.,1.,1./SI::c,EBMaxwellSolver::CowanZ);EBFieldGrid f(g);FieldVector<Double> origin(0.);
  for(unsigned a=0;a<3;++a){YeeComponent* v=a==0?&f.ex():a==1?&f.ey():&f.ez();for(std::size_t k=0;k<v->nz();++k)for(std::size_t j=0;j<v->ny();++j)for(std::size_t i=0;i<v->nx();++i)(*v)(i,j,k)=.03*i+.02*j+.04*k;}
  const auto retained=f.ez()(6,6,1),tangent=f.ex()(6,6,0);ChargeMatchedInlet completion(g,2,2);completion.begin(origin);
  RelativisticParticleSI p;p.position=FieldVector<Double>(6.25);p.position[1]=5.75;p.position[2]=1.2;p.charge=-1e-12;const auto oldPosition=p.position,oldMomentum=p.properVelocity;
  completion.particle(p);const auto r=completion.complete(f,true);
  require(r.maxBefore>1e-14&&r.maxAfter<1e-26,"Manufactured charge must satisfy Gauss after layer completion");
  require(f.ez()(6,6,1)==retained&&f.ex()(6,6,0)==tangent,"Completion changed retained/tangential fields");
  require(same(p.position,oldPosition)&&same(p.properVelocity,oldMomentum),"Completion changed source phase space");
  require(std::isfinite(r.nativeWorkJ)&&r.maxDeltaEz>0.,"Layer energy report invalid");
  std::cout<<"CHARGE_MATCHED_INLET_PASS\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}return 0;}
