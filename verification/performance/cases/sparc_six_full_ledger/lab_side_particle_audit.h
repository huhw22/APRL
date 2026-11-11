#ifndef APRL_LAB_SIDE_PARTICLE_AUDIT_H
#define APRL_LAB_SIDE_PARTICLE_AUDIT_H
#include "eb_particles.h"
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace aprl {
struct LabSideParticleEvent {
  unsigned face;int direction;std::uint64_t source;
  double timeBox,timeLab,zLab;long double kineticLab;
};
// Only the serial, terminal-clipped commit writes records. Parallel replay
// callbacks merely append to the existing particle-local audit trace.
class LabSideParticleAudit {
  struct State {bool enabled=false;double half[2]={};std::ofstream out;};
  static State& state(){static State s;return s;}
  static long double kinetic(const RelativisticParticleSI& p,const BoostFrameTransform& frame){
    long double u2=0.;for(unsigned a=0;a<3;++a)u2+=(long double)p.properVelocity[a]*p.properVelocity[a];
    const long double gp=std::sqrt(1.L+u2),uz=frame.gamma()*(long double)p.properVelocity[2]+frame.gammaBeta()*gp;
    const long double v2=(long double)p.properVelocity[0]*p.properVelocity[0]+(long double)p.properVelocity[1]*p.properVelocity[1]+uz*uz;
    return (long double)p.mass*SI::c*SI::c*v2/(std::sqrt(1.L+v2)+1.L);
  }
public:
  static bool enabled(){return state().enabled;}
  static void open(const std::string& directory,int rank,double hx,double hy){
    auto& s=state();if(s.enabled)throw std::runtime_error("Side-particle observer already enabled");
    s.half[0]=hx;s.half[1]=hy;std::ostringstream name;
    name<<directory<<"/particle-sides-rank-"<<std::setw(5)<<std::setfill('0')<<rank<<".csv";
    s.out.open(name.str());if(!s.out)throw std::runtime_error("Cannot open side-particle observer");
    s.out<<std::setprecision(21)<<"face,outward_direction,source_id,time_box_s,time_lab_s,z_lab_m,kinetic_lab_J\n";s.enabled=true;
  }
  static void preview(const RelativisticParticleSI& a,const RelativisticParticleSI& b,
      double t0,double t1,const BoostFrameTransform& frame,std::vector<LabSideParticleEvent>& out){
    if(!enabled()||!(t1>t0))return;const auto& s=state();
    bool haveK=false;long double ka=0.,kb=0.;
    for(unsigned face=0;face<4;++face){const unsigned axis=face/2,other=1-axis;const int normal=face%2?1:-1;const double boundary=normal*s.half[axis];
      const double x0=a.position[axis],x1=b.position[axis];
      const int crossing=x0<boundary&&x1>=boundary?1:(x0>boundary&&x1<=boundary?-1:0);if(!crossing)continue;
      const long double fraction=((long double)boundary-x0)/((long double)x1-x0);
      const long double transverse=a.position[other]+fraction*(b.position[other]-a.position[other]);
      // Face edges use the same half-open rectangular convention as the field aperture.
      if(transverse < -s.half[other]||transverse >= s.half[other])continue;
      if(!haveK){ka=kinetic(a,frame);kb=kinetic(b,frame);haveK=true;}
      const double time=t0+fraction*(t1-t0),z=a.position[2]+fraction*(b.position[2]-a.position[2]);
      LabSideParticleEvent event={face,crossing*normal,b.sourceId,time,
        frame.labTFromBoxZT(z,time),frame.labZFromBoxZT(z,time),ka+fraction*(kb-ka)};
      out.push_back(event);
    }
  }
  static void commit(const std::vector<LabSideParticleEvent>& events,double terminal){
    if(!enabled())return;auto& out=state().out;
    for(const auto& v:events)if(v.timeBox<=terminal)out<<v.face<<','<<v.direction<<','<<v.source<<','<<v.timeBox<<','<<v.timeLab<<','<<v.zLab<<','<<v.kineticLab<<'\n';
    if(!out)throw std::runtime_error("Side-particle observer write failed");
  }
  static void close(){auto& s=state();if(s.enabled){s.out.flush();if(!s.out)throw std::runtime_error("Side-particle observer flush failed");s.out.close();s.enabled=false;}}
};
}
#endif
