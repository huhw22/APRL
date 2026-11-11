#include "particle_energy_audit.h"
#include <cassert>
#include <iostream>
using namespace aprl;
int main(int argc,char**argv){
 try{
  if(argc!=2)throw std::runtime_error("side unit directory missing");
  BoostFrameTransform frame;frame.setOriginsFromGamma(3.,SI::c,0.,0.,0.);
  LabSideParticleAudit::open(argv[1],0,1.,1.);
  RelativisticParticleSI a,b;a.mass=b.mass=9.1093837015e-31;a.sourceId=b.sourceId=7;
  a.position=FieldVector<Double>(0.);b.position=FieldVector<Double>(0.);b.position[0]=2.;
  a.properVelocity=FieldVector<Double>(0.);b.properVelocity=FieldVector<Double>(0.);a.properVelocity[2]=1.;b.properVelocity[2]=2.;
  std::vector<LabSideParticleEvent> events;
  LabSideParticleAudit::preview(a,b,0.,2e-10,frame,events);
  if(events.size()!=1||events[0].face!=1||events[0].direction!=1)throw std::runtime_error("outward side event");
  const auto ka=ParticleEnergyMoment::sample(a,frame).kineticLab;
  const auto kb=ParticleEnergyMoment::sample(b,frame).kineticLab;
  if(std::abs(events[0].kineticLab-(ka+kb)/2)>1e-14L*events[0].kineticLab)throw std::runtime_error("scalar endpoint energy");
  LabSideParticleAudit::commit(events,.5e-10); // excluded: true crossing at 1e-10
  LabSideParticleAudit::commit(events,1.5e-10);
  events.clear();LabSideParticleAudit::preview(b,a,2e-10,4e-10,frame,events);
  if(events.size()!=1||events[0].direction!=-1)throw std::runtime_error("incoming sign");
  LabSideParticleAudit::commit(events,4e-10);
  b.position[1]=3.;events.clear();LabSideParticleAudit::preview(a,b,0.,2e-10,frame,events);
  if(events.size()!=1||events[0].face!=3)throw std::runtime_error("other transverse coordinate clipping");
  LabSideParticleAudit::close();
  std::cout<<"PASS outward/inward signs, scalar endpoint interpolation, terminal clipping, rectangular face clipping\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}return 0;
}
