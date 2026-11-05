#ifndef APRL_PARTICLE_ENERGY_AUDIT_H
#define APRL_PARTICLE_ENERGY_AUDIT_H

#include <algorithm>
#include <cmath>
#include <vector>
#include "eb_particles.h"
#include "lab_detectors.h"
#include "particle_work_geometry.h"

namespace aprl {

// Diagnostic four-momentum at a particle event. cP has energy units. These
// values never feed the particle pusher, current deposition or grid fields.
struct ParticleEnergyMoment {
  long double kineticBox, totalBox, cPzBox, kineticLab, totalLab, cPzLab;
  ParticleEnergyMoment() : kineticBox(0), totalBox(0), cPzBox(0),
    kineticLab(0), totalLab(0), cPzLab(0) {}
  static ParticleEnergyMoment sample(const RelativisticParticleSI& p,
                                    const BoostFrameTransform& frame) {
    ParticleEnergyMoment r;
    const long double rest = static_cast<long double>(p.mass) * SI::c * SI::c;
    long double u2 = 0;
    for (unsigned a=0; a<3; ++a)
      u2 += static_cast<long double>(p.properVelocity[a]) * p.properVelocity[a];
    const long double gamma = std::sqrt(1 + u2);
    r.kineticBox = rest * u2 / (gamma + 1);
    r.totalBox = rest * gamma;
    r.cPzBox = rest * p.properVelocity[2];
    FieldVector<Double> lab(0.);
    frame.properVelocityBoxToLab(p.properVelocity, lab);
    u2 = 0;
    for (unsigned a=0; a<3; ++a) u2 += static_cast<long double>(lab[a]) * lab[a];
    const long double gammaLab = std::sqrt(1 + u2);
    r.kineticLab = rest * u2 / (gammaLab + 1);
    r.totalLab = rest * gammaLab;
    r.cPzLab = rest * lab[2];
    return r;
  }
  void add(const ParticleEnergyMoment& v) {
    kineticBox+=v.kineticBox; totalBox+=v.totalBox; cPzBox+=v.cPzBox;
    kineticLab+=v.kineticLab; totalLab+=v.totalLab; cPzLab+=v.cPzLab;
  }
  static ParticleEnergyMoment interpolate(const ParticleEnergyMoment& a,
                                         const ParticleEnergyMoment& b,
                                         long double f) {
    ParticleEnergyMoment r;
    r.kineticBox=a.kineticBox+f*(b.kineticBox-a.kineticBox);
    r.totalBox=a.totalBox+f*(b.totalBox-a.totalBox);
    r.cPzBox=a.cPzBox+f*(b.cPzBox-a.cPzBox);
    r.kineticLab=a.kineticLab+f*(b.kineticLab-a.kineticLab);
    r.totalLab=a.totalLab+f*(b.totalLab-a.totalLab);
    r.cPzLab=a.cPzLab+f*(b.cPzLab-a.cPzLab);
    return r;
  }
};

// Per-preparation-batch storage, never a whole-bunch history. Work is computed
// from actual electric kicks. A terminal event inside a substep is explicitly
// approximated by a linear prefix, rather than scaling the entire field step.
class ParticleEnergyTrace : public ParticleSubstepObserver {
public:
  struct Segment {
    Double begin, end;
    ParticleEnergyMoment before, after;
    BorisWork work;
  };
  const BoostFrameTransform* frame = NULL;
  const LabDetectorManager* detectors = NULL;
  bool recordWork = false;
  bool recordSubstepParticles = false;
  unsigned observedSubsteps = 0;
  std::vector<Segment> segments;
  std::vector<EnergyPlaneCrossing> crossings;
  std::vector<SubstepParticleCrossing> substepCrossings;
  ParticleWorkGeometry geometry;
  void reset(const BoostFrameTransform& f, const LabDetectorManager* d, bool work) {
    frame=&f; detectors=d; recordWork=work;
    recordSubstepParticles=d&&d->hasSubstepParticlePlanes();
    segments.clear(); crossings.clear(); substepCrossings.clear(); observedSubsteps=0;
    geometry.fields=NULL;
  }
  void observe(const RelativisticParticleSI& before,
               const RelativisticParticleSI& after, Double begin, Double end,
               const BorisWork* work) override {
    ++observedSubsteps;
    if (geometry.fields) geometry.observe(before,after,work);
    if (detectors)
      {
        detectors->previewEnergyParticleStep(before, after, begin, end, crossings);
        if (recordSubstepParticles)
          detectors->previewSubstepParticleStep(before, after, begin, end, substepCrossings);
      }
    if (recordWork) {
      Segment s;
      s.begin=begin; s.end=end;
      s.before=ParticleEnergyMoment::sample(before,*frame);
      s.after=ParticleEnergyMoment::sample(after,*frame);
      if (work) s.work=*work;
      segments.push_back(s);
    }
  }
  BorisWork prefix(Double stop, unsigned long long* partial=NULL) const {
    BorisWork result;
    for (const auto& s:segments) {
      if (stop<=s.begin) break;
      const long double f=std::min(1.L,
        (static_cast<long double>(stop)-s.begin)/(s.end-s.begin));
      result.add(s.work,f);
      if (partial && f<1) ++*partial;
    }
    return result;
  }
  ParticleEnergyMoment at(Double time, const RelativisticParticleSI& fallback,
                         const BoostFrameTransform& f) const {
    if (segments.empty()) return ParticleEnergyMoment::sample(fallback,f);
    for (const auto& s:segments) {
      if (time<=s.begin) return s.before;
      if (time<=s.end) return ParticleEnergyMoment::interpolate(s.before,s.after,
        (static_cast<long double>(time)-s.begin)/(s.end-s.begin));
    }
    return segments.back().after;
  }
};
}
#endif
