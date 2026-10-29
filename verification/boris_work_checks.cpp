#include "eb_particles.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using namespace aprl;
typedef FieldVector<Double> Vec;
unsigned momentumComparisons = 0, trajectoryComparisons = 0;
long double maximumRelativeRounding = 0;

void require(bool pass, const char* message) {
  if (!pass) throw std::runtime_error(message);
}
bool sameBits(Double a, Double b) {
  return std::memcmp(&a, &b, sizeof a) == 0;
}
void sameVector(const Vec& a, const Vec& b, const char* message) {
  for (unsigned i=0; i<3; ++i) require(sameBits(a[i], b[i]), message);
}
void sameParticle(const RelativisticParticleSI& a, const RelativisticParticleSI& b) {
  sameVector(a.position,b.position,"Diagnostic changed particle position bits");
  sameVector(a.properVelocity,b.properVelocity,"Diagnostic changed particle momentum bits");
  require(sameBits(a.charge,b.charge) && sameBits(a.mass,b.mass) &&
    sameBits(a.weight,b.weight) && a.id==b.id && a.sourceId==b.sourceId &&
    a.frozen==b.frozen && sameBits(a.frozenBetaZ,b.frozenBetaZ),
    "Diagnostic changed particle metadata");
}
Vec vector(Double x, Double y, Double z) {
  Vec v(0.); v[0]=x; v[1]=y; v[2]=z; return v;
}
long double gamma(const Vec& u) {
  // Independent norm evaluation; no shared audit helper or difference formula.
  return std::hypot(std::hypot(1.L,static_cast<long double>(u[0])),
                    std::hypot(static_cast<long double>(u[1]),
                               static_cast<long double>(u[2])));
}
long double restEnergy(Double mass) {
  return static_cast<long double>(mass)*SI::c*SI::c;
}
void near(long double a, long double b, long double bound, const char* message) {
  if (!(std::isfinite(a) && std::isfinite(b) && std::abs(a-b)<=bound)) {
    std::cerr.precision(20);
    std::cerr << message << " actual=" << a << " expected=" << b
              << " bound=" << bound << '\n';
    throw std::runtime_error(message);
  }
}
void audit(const Vec& before, const Vec& after, Double mass, const BorisWork& w) {
  const long double energyScale=restEnergy(mass)*(gamma(before)+gamma(after));
  const long double eps=std::numeric_limits<long double>::epsilon();
  near(w.totalDeltaK,restEnergy(mass)*(gamma(after)-gamma(before)),
       16*eps*energyScale,"Stable delta K disagrees with independent energy difference");
  const long double roundingScale=energyScale+std::abs(w.gridWork)+
    std::abs(w.prescribedWork);
  const long double relative=std::abs(w.roundoffResidual)/roundingScale;
  maximumRelativeRounding=std::max(maximumRelativeRounding,relative);
  require(relative<=64*std::numeric_limits<Double>::epsilon(),
    "Electric kick attribution exceeded floating-point error scale");
  near(w.gridWork+w.prescribedWork+w.magneticRotationDeltaK+w.roundoffResidual,
       w.totalDeltaK,16*eps*roundingScale,"Boris work components do not reconstruct delta K");
}

void momentumChecks() {
  const Vec initial[] = {vector(0.,0.,0.),vector(1e-9,-2e-9,3e-9),
    vector(.3,-.7,.23),vector(12.,-25.,999.),vector(0.,0.,1e6)};
  const Vec grids[] = {Vec(0.),vector(2e7,-3e7,7e7),vector(-8e8,0.,0.)};
  const Vec prescribed[] = {Vec(0.),vector(-1e7,5e7,-2e7),vector(8e8,0.,0.)};
  const Vec magnets[] = {Vec(0.),vector(.1,2.,-.3)};
  for (const Vec& start:initial) for (const Vec& eg:grids)
    for (const Vec& ep:prescribed) for (const Vec& b:magnets)
      for (Double weight:{1.,4096.}) for (Double sign:{-1.,1.}) {
        Vec old(start),observed(start),e(eg); e+=ep;
        const Double mass=weight*SI::electronMass;
        const Double charge=sign*weight*SI::elementaryCharge;
        for (unsigned n=0;n<5;++n) {
          const Vec before(observed);
          BorisWork work;
          RelativisticBorisPusher::pushMomentum(old,e,b,charge,mass,2e-13);
          RelativisticBorisPusher::pushMomentumWithWork(
            observed,e,b,eg,ep,charge,mass,2e-13,work);
          sameVector(old,observed,"Audited Boris momentum is not bitwise identical");
          audit(before,observed,mass,work);
          ++momentumComparisons;
        }
      }

  // Analytic straight acceleration and deceleration with independently known
  // source fractions; the grid work cannot be recovered by subtracting an
  // arbitrary global conservation residual.
  for (Double uz:{-2.,0.,2.}) {
    Vec u=vector(0.,0.,uz),before(u),eg=vector(0.,0.,6e7);
    const Vec ep=vector(0.,0.,2e7),b(0.);
    Vec e(eg); e+=ep;
    BorisWork w;
    RelativisticBorisPusher::pushMomentumWithWork(u,e,b,eg,ep,
      SI::elementaryCharge,SI::electronMass,1e-12,w);
    const long double independent=restEnergy(SI::electronMass)*(gamma(u)-gamma(before));
    const long double bound=64*std::numeric_limits<Double>::epsilon()*
      restEnergy(SI::electronMass)*(gamma(before)+gamma(u));
    near(w.gridWork,.75L*independent,bound,"Known grid work fraction failed");
    near(w.prescribedWork,.25L*independent,bound,"Known prescribed work fraction failed");
    require(w.magneticRotationDeltaK==0.,"Zero B introduced magnetic work");
  }

  // Opposing E contributions do zero net work, but the source work channels
  // must remain equal and opposite rather than both being reported as zero.
  Vec cancel=vector(.2,.1,.3),eg=vector(2e8,-1e8,4e8),ep(eg); ep*=-1.;
  const Vec zero(0.);
  BorisWork cancellation;
  RelativisticBorisPusher::pushMomentumWithWork(cancel,zero,zero,eg,ep,
    -SI::elementaryCharge,SI::electronMass,1e-12,cancellation);
  require(cancellation.gridWork!=0. && cancellation.prescribedWork==-cancellation.gridWork &&
          cancellation.totalDeltaK==0.,"Canceling source work was not retained");

  // Exchange attribution labels with the same exact total field and require
  // identical momentum plus interchanged source work.
  Vec mixed=vector(.23,-.11,.67),swapped(mixed),e(eg); ep=vector(7e7,2e7,-4e7); e+=ep;
  const Vec b=vector(.3,.7,-.1);
  BorisWork first,second;
  RelativisticBorisPusher::pushMomentumWithWork(mixed,e,b,eg,ep,
    -SI::elementaryCharge,SI::electronMass,2e-12,first);
  RelativisticBorisPusher::pushMomentumWithWork(swapped,e,b,ep,eg,
    -SI::elementaryCharge,SI::electronMass,2e-12,second);
  sameVector(mixed,swapped,"Attribution labels changed state");
  require(first.gridWork==second.prescribedWork && first.prescribedWork==second.gridWork,
    "Exchanged source labels did not exchange work");

  // A long pure magnetic rotation must expose its actual numerical energy
  // defect independently of both electric work channels.
  Vec gyro=vector(.5,-.3,.7),ordinary(gyro);
  BorisWork cumulative;
  for (unsigned n=0;n<1000;++n) {
    const Vec before(gyro); BorisWork w;
    RelativisticBorisPusher::pushMomentum(ordinary,zero,b,
      -SI::elementaryCharge,SI::electronMass,1e-12);
    RelativisticBorisPusher::pushMomentumWithWork(gyro,zero,b,zero,zero,
      -SI::elementaryCharge,SI::electronMass,1e-12,w);
    sameVector(ordinary,gyro,"Pure B audited state changed");
    require(w.gridWork==0. && w.prescribedWork==0. &&
      w.totalDeltaK==w.magneticRotationDeltaK && w.roundoffResidual==0.,
      "Pure B work attribution failed");
    audit(before,gyro,SI::electronMass,w); cumulative.add(w);
  }
}

struct Recorder : ParticleSubstepObserver {
  struct Record { RelativisticParticleSI before,after; Double begin,end;
    bool hasWork; BorisWork work; };
  std::vector<Record> records;
  void observe(const RelativisticParticleSI& before,
      const RelativisticParticleSI& after,Double begin,Double end,const BorisWork* work) override {
    Record r; r.before=before; r.after=after; r.begin=begin; r.end=end;
    r.hasWork=work!=0; if (work) r.work=*work; records.push_back(r);
  }
};

void sameWork(const BorisWork& a,const BorisWork& b,const char* message) {
  require(a.gridWork==b.gridWork && a.prescribedWork==b.prescribedWork &&
    a.magneticRotationDeltaK==b.magneticRotationDeltaK &&
    a.totalDeltaK==b.totalDeltaK && a.roundoffResidual==b.roundoffResidual,message);
}

void subcycleChecks() {
  const Double dt=4e-15;
  EBGridGeometry geometry(8,8,8,5e-6,5e-6,1.25e-6,dt,EBMaxwellSolver::CowanZ);
  EBFieldGrid grid(geometry);
  grid.ex().fill(1e5); grid.ey().fill(-3e5); grid.ez().fill(2e5);
  grid.bx().fill(1e-3); grid.by().fill(2e-3); grid.bz().fill(-3e-3);
  const Vec origin=vector(-20e-6,-20e-6,-5e-6);
  BoostFrameTransform frame; frame.setOriginsFromGamma(929.1305492858543,SI::c,0.,0.,0.);
  SIFieldSourceSet sources;
  sources.addPrescribedLabMagnet(SIMagneticElement::planarUndulatorFromK(.6,.0234,0.,8,.1));
  RelativisticParticleSI initial;
  initial.position=vector(0.,0.,.5e-6); initial.properVelocity=vector(.02,-.01,.23);
  initial.weight=1024.; initial.charge=-SI::elementaryCharge*initial.weight;
  initial.mass=SI::electronMass*initial.weight; initial.id=71; initial.sourceId=731;
  for (bool includeGrid:{false,true}) for (unsigned steps:{1u,5u,17u})
    for (Double duration:{dt,.37*dt}) {
      const Double start=.31*dt;
      RelativisticParticleSI plain(initial),workOnly(initial),observerOnly(initial),both(initial);
      Double oldLegacy=0.,workLegacy=0.,observerLegacy=0.,bothLegacy=0.;
      BorisWork a,c; Recorder observer,collector;
      RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
        plain,grid,origin,sources,frame,start,duration,steps,&oldLegacy,includeGrid);
      RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
        workOnly,grid,origin,sources,frame,start,duration,steps,&workLegacy,includeGrid,&a);
      RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
        observerOnly,grid,origin,sources,frame,start,duration,steps,&observerLegacy,includeGrid,0,&observer);
      RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
        both,grid,origin,sources,frame,start,duration,steps,&bothLegacy,includeGrid,&c,&collector);
      sameParticle(plain,workOnly); sameParticle(plain,observerOnly); sameParticle(plain,both);
      require(sameBits(oldLegacy,workLegacy) && sameBits(oldLegacy,observerLegacy) &&
        sameBits(oldLegacy,bothLegacy),"Legacy prescribed work changed");
      sameWork(a,c,"Observer changed Boris work");
      require(observer.records.size()==steps && collector.records.size()==steps,
        "Wrong observer event count");
      BorisWork summed,prefix,manualPrefix;
      const unsigned stop=steps/2;
      for (unsigned n=0;n<steps;++n) {
        const Recorder::Record& o=observer.records[n];
        const Recorder::Record& r=collector.records[n];
        require(!o.hasWork && r.hasWork,"Observer work optionality was not respected");
        require(r.end>r.begin,"Substep time not increasing");
        near(r.begin,std::fma(static_cast<Double>(n),duration/steps,start),0.,
          "Observer start time differs from source sample time");
        sameParticle(o.before,r.before); sameParticle(o.after,r.after);
        if (n) sameParticle(collector.records[n-1].after,r.before);
        else sameParticle(initial,r.before);
        audit(r.before.properVelocity,r.after.properVelocity,initial.mass,r.work);
        summed.add(r.work);
        if (n<=stop) prefix.add(r.work,n==stop?.25L:1.L);
        if (n<stop) manualPrefix.add(r.work);
      }
      sameParticle(collector.records.back().after,both);
      sameWork(summed,c,"Observer sum differs from subcycle work");
      BorisWork quarter=collector.records[stop].work;
      quarter.gridWork*=.25L; quarter.prescribedWork*=.25L;
      quarter.magneticRotationDeltaK*=.25L; quarter.totalDeltaK*=.25L;
      quarter.roundoffResidual*=.25L; manualPrefix.add(quarter);
      sameWork(prefix,manualPrefix,"Fractional terminal substep accumulation failed");
      if (!includeGrid) require(c.gridWork==0.,"Disabled grid fields produced grid work");
      ++trajectoryComparisons;
    }
}
}

int main() {
  try {
    momentumChecks(); subcycleChecks();
    std::cout.precision(18);
    std::cout << "PASS boris_work momentum_comparisons=" << momentumComparisons
      << " trajectory_comparisons=" << trajectoryComparisons
      << " maximum_local_roundoff_over_energy_scale=" << maximumRelativeRounding << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FAIL boris_work: " << error.what() << '\n'; return 1;
  }
}
