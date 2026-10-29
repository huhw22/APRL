// Diagnostic-only free-drift coupling probe. No Maxwell advance, CPML, MPI,
// moving-window supply, or external magnetic field is exercised here.
// Link the existing eb_field/eb_deposition/eb_particles/eb_sources/boostframe
// implementation; do not substitute a new pusher or current depositor.
#include "eb_deposition.h"
#include "eb_particles.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
using namespace aprl;
using Vec = FieldVector<Double>;
const double baseDx = 5e-6;
const double baseDz = 1.2543262415359034e-6;
const double boostGamma = 929.1305492858543;
const double amplitude = 1e6; // Synthetic E in V/m, not a sampled PIC field.
const double gradient = .30; // Synthetic fractional change per BASE dz.
const unsigned baseSteps = 8;

Vec vector(double x, double y, double z) {
  Vec a(0.); a[0]=x; a[1]=y; a[2]=z; return a;
}
void require(bool good, const char* message) {
  if (!good) throw std::runtime_error(message);
}
void near(long double a, long double b, long double bound, const char* message) {
  if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a-b)>bound) {
    std::cerr << std::setprecision(21) << message << " actual=" << a
              << " expected=" << b << " bound=" << bound << '\n';
    throw std::runtime_error(message);
  }
}
long double gamma(const Vec& u) {
  return std::hypot(std::hypot(1.L, static_cast<long double>(u[0])),
    std::hypot(static_cast<long double>(u[1]), static_cast<long double>(u[2])));
}
void sameParticle(const RelativisticParticleSI& a,
                  const RelativisticParticleSI& b) {
  for (unsigned axis=0; axis<3; ++axis) {
    require(std::memcmp(&a.position[axis], &b.position[axis], sizeof(double))==0,
      "Work observer changed particle position");
    require(std::memcmp(&a.properVelocity[axis], &b.properVelocity[axis],
                        sizeof(double))==0,
      "Work observer changed particle momentum");
  }
}
long double edgeDot(const YeeComponent& j, const YeeComponent& e) {
  require(j.nx()==e.nx() && j.ny()==e.ny() && j.nz()==e.nz(),
          "Current/electric lattice mismatch");
  long double sum=0.;
  for (std::size_t k=0;k<j.nz();++k)
    for (std::size_t y=0;y<j.ny();++y)
      for (std::size_t x=0;x<j.nx();++x)
        sum += static_cast<long double>(j(x,y,k))*e(x,y,k);
  return sum;
}
long double currentWork(const EBFieldGrid& f) {
  const auto& g=f.geometry();
  return static_cast<long double>(g.dt)*g.dx*g.dy*g.dz*
    (edgeDot(f.jx(),f.ex())+edgeDot(f.jy(),f.ey())+edgeDot(f.jz(),f.ez()));
}

// Integral from zero to s of floor(s)+1/2; independent of the depositor's
// transverse shape integration. The selected field has only Ez(z), so those
// transverse shape weights must sum to one for each segment.
long double centrePrimitive(long double s) {
  const long double cell=std::floor(s);
  return .5L*cell*cell+(s-cell)*(cell+.5L);
}
long double analyticDepositedWork(double start, double end, double q,
                                  const EBGridGeometry& g, bool varying) {
  const long double a=start, b=end, dz=g.dz;
  if (!varying) return static_cast<long double>(q)*amplitude*(b-a);
  const long double anchor=6.23L*baseDz;
  const long double piecewiseIntegral=
    (1.L-gradient*anchor/baseDz)*(b-a)+
    gradient/baseDz*dz*dz*(centrePrimitive(b/dz)-centrePrimitive(a/dz));
  return static_cast<long double>(q)*amplitude*piecewiseIntegral;
}
void setElectric(EBFieldGrid& f, bool varying) {
  f.ex().fill(0.); f.ey().fill(0.);
  f.bx().fill(0.); f.by().fill(0.); f.bz().fill(0.);
  const double anchor=6.23*baseDz;
  auto& e=f.ez();
  for (std::size_t k=0;k<e.nz();++k)
    for (std::size_t j=0;j<e.ny();++j)
      for (std::size_t i=0;i<e.nx();++i) {
        const double z=(static_cast<double>(k)+.5)*f.geometry().dz;
        e(i,j,k)=amplitude*(varying ? 1.+gradient*(z-anchor)/baseDz : 1.);
      }
}

struct Result {
  long double boris=0., path=0., deposited=0.;
  long double analyticDeposit=0., analyticVelocityDefect=0., totalDeltaK=0.;
  long double maximumContinuity=0., maximumWorkRounding=0.;
  unsigned segments=0;
};

Result run(bool varying, bool matchedMesh, unsigned refinement,
           const Vec& initialU, const BoostFrameTransform& frame) {
  const double baseDt=baseDz/SI::c;
  const double particleDt=baseDt/refinement;
  // fixed_mesh keeps geometry.dt at baseDt because the Cowan geometry
  // contract is c*dt=dz. J is normalized by that enclosing time, and its work
  // is multiplied by the SAME time. Only the represented path is shorter.
  // matched_dz_dt refines both to retain the physical domain and c*dt=dz.
  const double dz=matchedMesh ? baseDz/refinement : baseDz;
  const unsigned nz=matchedMesh ? 16*refinement : 16;
  const EBGridGeometry geometry(16,16,nz,baseDx,baseDx,dz,dz/SI::c,
                                EBMaxwellSolver::CowanZ);
  EBFieldGrid fields(geometry); setElectric(fields,varying);
  const Vec origin(0.);
  ChargeConservingCurrentDepositor deposit(fields,origin);
  YeeComponent rho0(17,17,nz+1),rho1(17,17,nz+1);
  SIFieldSourceSet emptySources;
  RelativisticParticleSI particle;
  particle.position=vector(6.21*baseDx,6.37*baseDx,6.23*baseDz);
  particle.properVelocity=initialU;
  particle.weight=1.;
  particle.charge=-SI::elementaryCharge;
  particle.mass=SI::electronMass;
  particle.id=particle.sourceId=1;
  const long double rest=static_cast<long double>(particle.mass)*SI::c*SI::c;
  const long double eps=std::numeric_limits<double>::epsilon();
  Result result;
  for (unsigned step=0;step<baseSteps*refinement;++step) {
    const RelativisticParticleSI before(particle);
    Vec e(0.),b(0.);
    fields.sampleFieldsPosition(before.position,origin,e,b);
    require(b.norm2()==0. && e[0]==0. && e[1]==0.,"Unexpected field component");
    const long double analyticGather=amplitude*(varying ?
      1.L+gradient*(static_cast<long double>(before.position[2])-6.23L*baseDz)/baseDz
      : 1.L);
    near(e[2],analyticGather,256*eps*std::abs(analyticGather),
         "CIC does not reproduce the linear synthetic field");
    RelativisticParticleSI plain(before);
    RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
      plain,fields,origin,emptySources,frame,step*particleDt,particleDt,1);
    BorisWork work;
    RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
      particle,fields,origin,emptySources,frame,step*particleDt,particleDt,1,
      nullptr,true,&work);
    sameParticle(plain,particle);
    require(work.prescribedWork==0. && work.magneticRotationDeltaK==0.,
            "Unrequested prescribed or magnetic work");
    for (unsigned axis=0;axis<3;++axis) {
      const double spacing=axis==2 ? dz : baseDx;
      const double extent=(axis==2 ? nz : 16)*spacing;
      require(particle.position[axis]>2*spacing &&
              particle.position[axis]<extent-2*spacing,
              "Probe entered a boundary collar");
    }

    fields.clearCurrent(); rho0.fill(0.); rho1.fill(0.);
    deposit.depositCharge(before.position,before.charge,rho0);
    const auto deposited=deposit.depositSegment(before.position,particle.position,
                                                before.charge);
    deposit.depositCharge(particle.position,particle.charge,rho1);
    const long double wj=currentWork(fields);
    const long double wp=static_cast<long double>(before.charge)*e[2]*
      (static_cast<long double>(particle.position[2])-before.position[2]);
    const long double wa=analyticDepositedWork(before.position[2],
      particle.position[2],before.charge,geometry,varying);
    const long double ga=gamma(before.properVelocity);
    const long double gb=gamma(particle.properVelocity);
    const long double conjugateVelocity=SI::c*
      (static_cast<long double>(before.properVelocity[2])+particle.properVelocity[2])/
      (ga+gb);
    const long double afterVelocity=SI::c*particle.properVelocity[2]/gb;
    const long double analyticBoris=static_cast<long double>(before.charge)*
      particleDt*e[2]*conjugateVelocity;
    const long double analyticDefect=static_cast<long double>(before.charge)*
      particleDt*e[2]*(afterVelocity-conjugateVelocity);
    const long double arithmeticBound=256*eps*rest*(ga+gb);
    near(work.gridWork,analyticBoris,arithmeticBound,
         "B=0 Boris work violates independent conjugate-velocity identity");
    near(wp-work.gridWork,analyticDefect,arithmeticBound,
         "Position/Boris defect violates analytic velocity formula");
    near(wj,wa,4096*eps*std::max(std::abs(wj),std::abs(wa)),
         "Deposited J work differs from independent piecewise-field integral");
    if (!varying)
      near(wj,wp,4096*eps*std::max(std::abs(wj),std::abs(wp)),
           "Uniform E deposition does not give q E dot displacement");
    const double continuity=ChargeConservingCurrentDepositor::maxContinuityResidual(
      fields,rho0,rho1);
    double continuityScale=0.;
    for (std::size_t k=0;k<rho0.nz();++k)
      for (std::size_t j=0;j<rho0.ny();++j)
        for (std::size_t i=0;i<rho0.nx();++i)
          continuityScale=std::max(continuityScale,
            std::abs(rho1(i,j,k)-rho0(i,j,k))/geometry.dt);
    require(continuityScale>0. && continuity/continuityScale<2e-11,
            "Discrete continuity failed");
    result.maximumContinuity=std::max(result.maximumContinuity,
      static_cast<long double>(continuity/continuityScale));
    result.maximumWorkRounding=std::max(result.maximumWorkRounding,
      std::abs(work.roundoffResidual));
    result.boris+=work.gridWork; result.path+=wp; result.deposited+=wj;
    result.analyticDeposit+=wa; result.analyticVelocityDefect+=analyticDefect;
    result.totalDeltaK+=work.totalDeltaK; result.segments+=deposited.subsegments;
  }
  require(std::isfinite(result.boris) && std::isfinite(result.deposited),
          "Nonfinite cumulative result");
  return result;
}

void print(const Result& r, bool varying, bool matched, unsigned n) {
  std::cout << "{\"type\":\"result\",\"field\":\""
    << (varying ? "linear_Ez_z" : "uniform_Ez")
    << "\",\"refinement_mode\":\""
    << (matched ? "matched_dz_dt" : "fixed_mesh_shorter_particle_interval")
    << "\",\"refinement\":" << n
    << ",\"steps\":" << baseSteps*n
    << ",\"particle_dt_s\":" << baseDz/SI::c/n
    << ",\"current_normalization_dt_s\":" << baseDz/SI::c/(matched?n:1)
    << ",\"W_Boris_J\":" << r.boris
    << ",\"q_Egather_dx_J\":" << r.path
    << ",\"dt_sum_J_E_J\":" << r.deposited
    << ",\"Wpath_minus_WBoris_J\":" << r.path-r.boris
    << ",\"Wdeposit_minus_Wpath_J\":" << r.deposited-r.path
    << ",\"Wdeposit_minus_WBoris_J\":" << r.deposited-r.boris
    << ",\"analytic_deposit_work_J\":" << r.analyticDeposit
    << ",\"analytic_velocity_defect_J\":" << r.analyticVelocityDefect
    << ",\"sum_delta_K_J\":" << r.totalDeltaK
    << ",\"max_continuity_relative\":" << r.maximumContinuity
    << ",\"max_Boris_roundoff_J\":" << r.maximumWorkRounding
    << ",\"deposition_subsegments\":" << r.segments
    << ",\"bitwise_trajectory\":true,\"pass\":true}\n";
}
}

int main() {
  try {
    std::cout << std::setprecision(21);
    BoostFrameTransform frame;
    frame.setOriginsFromGamma(boostGamma,SI::c,0.,0.,0.);
    // Exact first-record momentum from the existing short-case input. The
    // initial position is relocated, and charge/mass represent one electron;
    // this is a coupling probe, not that macroparticle's physical trajectory.
    const Vec labU=vector(-0.0017307752502591376,-0.040745133968771524,
                          1171.4907206638086);
    Vec boxU(0.); frame.properVelocityLabToBox(labU,boxU);
    std::cout << "{\"type\":\"metadata\",\"scope\":\"No Maxwell evolution;"
      " original gather, Boris and deposition algebra only; synthetic static E,"
      " B=0; not a PIC energy or physical convergence acceptance\","
      "\"momentum_source\":\"campaign-r3/inputs/held.h5 particles/records[0]"
      " source_id=1\",\"input_sha256\":\""
      "58f4cc2d750d2b00b665b16ec81988824d04a3403209c3841c490d90846717a4\","
      "\"boost_gamma\":" << boostGamma
      << ",\"initial_box_u\":[" << boxU[0] << ',' << boxU[1] << ',' << boxU[2]
      << "],\"E0_V_per_m\":" << amplitude
      << ",\"fractional_gradient_per_base_dz\":" << gradient
      << ",\"base_dx_m\":" << baseDx << ",\"base_dz_m\":" << baseDz
      << ",\"base_dt_s\":" << baseDz/SI::c
      << ",\"physical_duration_s\":" << baseSteps*baseDz/SI::c
      << ",\"current_metric\":\"native edge work, interior support;"
      " not the collocated production ledger field energy\"}\n";
    for (bool matched:{false,true})
      for (bool varying:{false,true})
        for (unsigned refinement:{1u,2u,4u,8u,16u})
          print(run(varying,matched,refinement,boxU,frame),
                varying,matched,refinement);
    std::cout << "{\"type\":\"complete\",\"rows\":20,\"pass\":true}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FAIL deposition_work_probe: " << error.what() << '\n';
    return 1;
  }
}
