// Isolated self-field startup experiment. Every branch releases the identical
// generated-gaussian x/p at t=0. There is no held interval, source injection,
// moving window or external device. Zero/scaled fields are deliberately
// constraint-inconsistent negative/sensitivity controls, not physical models.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <mpi.h>
#include "compute_threads.h"
#include "eb_deposition.h"
#include "eb_initial_field.h"
#include "particle_initializer.h"
#include "cowan_reference_field.h"

namespace {
using aprl::Double;
using aprl::FieldVector;
using aprl::RelativisticParticleSI;
using aprl::YeeComponent;
const double sigma = 4.e-6, totalCharge = -1.e-12;
const std::size_t count = 4096;
const unsigned seed = 710271;

std::ofstream fresh(const std::string& path) {
  std::ifstream existing(path);
  if (existing.good()) throw std::runtime_error("Refusing to overwrite " + path);
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Cannot open " + path);
  out << std::setprecision(17);
  return out;
}

void particlesCsv(const std::string& path,
                  const std::vector<RelativisticParticleSI>& particles) {
  auto out=fresh(path);
  out << "id,x_m,y_m,z_m,ux,uy,uz,q_C,m_kg\n";
  for(const auto& p:particles) {
    out << p.sourceId;
    for(unsigned a=0;a<3;++a) out << ',' << p.position[a];
    for(unsigned a=0;a<3;++a) out << ',' << p.properVelocity[a];
    out << ',' << p.charge << ',' << p.mass << '\n';
  }
}

long double gamma(const FieldVector<double>& u) {
  long double u2=0;
  for(unsigned a=0;a<3;++a) u2+=static_cast<long double>(u[a])*u[a];
  return std::sqrt(1.L+u2);
}

struct Moments {
  long double mass=0, charge=0, kinetic=0, deltaK=0;
  long double x[3]={}, u[3]={}, varianceX[3]={}, varianceU[3]={}, xu[3]={};
};

Moments moments(const std::vector<RelativisticParticleSI>& p,
                const std::vector<RelativisticParticleSI>& initial) {
  Moments s;
  const long double c=aprl::SI::c;
  for(std::size_t i=0;i<p.size();++i) {
    const long double m=p[i].mass, g=gamma(p[i].properVelocity);
    const long double g0=gamma(initial[i].properVelocity);
    long double u2=0, du2=0;
    for(unsigned a=0;a<3;++a) {
      const long double u=p[i].properVelocity[a], u0=initial[i].properVelocity[a];
      u2+=u*u; du2+=(u-u0)*(u+u0);
      s.x[a]+=m*p[i].position[a]; s.u[a]+=m*u;
    }
    s.mass+=m; s.charge+=p[i].charge;
    s.kinetic+=m*c*c*u2/(g+1.L);
    s.deltaK+=m*c*c*du2/(g+g0);
  }
  for(unsigned a=0;a<3;++a) {s.x[a]/=s.mass;s.u[a]/=s.mass;}
  for(const auto& particle:p) for(unsigned a=0;a<3;++a) {
    const long double dx=particle.position[a]-s.x[a];
    const long double du=particle.properVelocity[a]-s.u[a];
    s.varianceX[a]+=particle.mass*dx*dx;
    s.varianceU[a]+=particle.mass*du*du;
    s.xu[a]+=particle.mass*dx*du;
  }
  for(unsigned a=0;a<3;++a) {
    s.varianceX[a]/=s.mass;s.varianceU[a]/=s.mass;s.xu[a]/=s.mass;
  }
  return s;
}

// Trapezoidal dual-volume weights: E has two nodal axes, B one nodal axis.
double nodeWeight(std::size_t index,std::size_t cells) {
  return index==0||index==cells ? .5 : 1.;
}
double edgeWeight(unsigned a,std::size_t i,std::size_t j,std::size_t k,
                  const aprl::EBGridGeometry& g) {
  return (a==0?1.:nodeWeight(i,g.nx)) *
         (a==1?1.:nodeWeight(j,g.ny)) *
         (a==2?1.:nodeWeight(k,g.nz));
}
double faceWeight(unsigned a,std::size_t i,std::size_t j,std::size_t k,
                  const aprl::EBGridGeometry& g) {
  return a==0?nodeWeight(i,g.nx):(a==1?nodeWeight(j,g.ny):nodeWeight(k,g.nz));
}

struct FieldEnergy {long double electric=0,magneticCentered=0,magneticNative=0;};
FieldEnergy energy(aprl::EBFieldGrid& f,const std::array<std::vector<double>,3>& oldB) {
  const auto& g=f.geometry();
  const long double dv=static_cast<long double>(g.dx)*g.dy*g.dz;
  const YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};
  const YeeComponent* b[3]={&f.bx(),&f.by(),&f.bz()};
  FieldEnergy out;
  for(unsigned a=0;a<3;++a) {
    for(std::size_t k=0;k<e[a]->nz();++k)
      for(std::size_t j=0;j<e[a]->ny();++j)
        for(std::size_t i=0;i<e[a]->nx();++i) {
          const long double v=(*e[a])(i,j,k);
          out.electric+=edgeWeight(a,i,j,k,g)*v*v;
        }
    for(std::size_t k=0;k<b[a]->nz();++k)
      for(std::size_t j=0;j<b[a]->ny();++j)
        for(std::size_t i=0;i<b[a]->nx();++i) {
          const std::size_t n=(k*b[a]->ny()+j)*b[a]->nx()+i;
          const long double value=(*b[a])(i,j,k);
          const long double centered=.5L*(value+oldB[a][n]);
          const double w=faceWeight(a,i,j,k,g);
          out.magneticCentered+=w*centered*centered;
          out.magneticNative+=w*value*value;
        }
  }
  out.electric*=.5L*aprl::SI::epsilon0*dv;
  out.magneticCentered*=.5L/aprl::SI::mu0*dv;
  out.magneticNative*=.5L/aprl::SI::mu0*dv;
  return out;
}

struct Gauss {double max=0,relativeMax=0,relativeL2=0;};
Gauss gauss(aprl::EBFieldGrid& f,const YeeComponent& rho) {
  const auto& g=f.geometry();Gauss out;long double err2=0,source2=0;double sourceMax=0;
  for(std::size_t k=1;k<g.nz;++k) for(std::size_t j=1;j<g.ny;++j)
    for(std::size_t i=1;i<g.nx;++i) {
      const double source=rho(i,j,k)/aprl::SI::epsilon0;
      const double residual=(f.ex()(i,j,k)-f.ex()(i-1,j,k))/g.dx+
        (f.ey()(i,j,k)-f.ey()(i,j-1,k))/g.dy+
        (f.ez()(i,j,k)-f.ez()(i,j,k-1))/g.dz-source;
      out.max=std::max(out.max,std::abs(residual));
      sourceMax=std::max(sourceMax,std::abs(source));
      err2+=static_cast<long double>(residual)*residual;
      source2+=static_cast<long double>(source)*source;
    }
  out.relativeMax=out.max/sourceMax;
  out.relativeL2=std::sqrt(err2/source2);return out;
}

std::array<std::vector<double>,3> copyComponents(YeeComponent& a,YeeComponent& b,YeeComponent& c) {
  return {{std::vector<double>(a.data(),a.data()+a.size()),
           std::vector<double>(b.data(),b.data()+b.size()),
           std::vector<double>(c.data(),c.data()+c.size())}};
}

long double currentWork(aprl::EBFieldGrid& f,const std::array<std::vector<double>,3>& oldE) {
  const auto& g=f.geometry();
  const YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};
  const YeeComponent* current[3]={&f.jx(),&f.jy(),&f.jz()};long double work=0;
  for(unsigned a=0;a<3;++a) for(std::size_t k=0;k<e[a]->nz();++k)
    for(std::size_t j=0;j<e[a]->ny();++j) for(std::size_t i=0;i<e[a]->nx();++i) {
      const std::size_t n=(k*e[a]->ny()+j)*e[a]->nx()+i;
      work+=edgeWeight(a,i,j,k,g)*.5L*(oldE[a][n]+(*e[a])(i,j,k))*(*current[a])(i,j,k);
    }
  return work*static_cast<long double>(g.dt)*g.dx*g.dy*g.dz;
}

std::array<double,3> forceMoment(const std::vector<RelativisticParticleSI>& particles,
    aprl::EBFieldGrid& f,const FieldVector<double>& origin,bool analytic,double beta) {
  const auto stats=moments(particles,particles);long double out[3]={};
  const double beamGamma=1/std::sqrt(1-beta*beta),pi=std::acos(-1.);
  for(const auto& p:particles) {
    FieldVector<double> e(0.),b(0.);
    if(analytic) {
      double r[3]={p.position[0],p.position[1],beamGamma*p.position[2]};
      const double radius=std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
      const double a=radius/sigma;
      const double enclosed=std::erf(a/std::sqrt(2.))-std::sqrt(2./pi)*a*std::exp(-.5*a*a);
      const double scale=radius>0?totalCharge*enclosed/(4*pi*aprl::SI::epsilon0*radius*radius*radius):0;
      e[0]=beamGamma*scale*r[0];e[1]=beamGamma*scale*r[1];e[2]=scale*r[2];
      b[0]=-beta*e[1]/aprl::SI::c;b[1]=beta*e[0]/aprl::SI::c;
    } else f.sampleFieldsPosition(p.position,origin,e,b);
    const double inverseGamma=1/static_cast<double>(gamma(p.properVelocity));
    double v[3];for(unsigned a=0;a<3;++a)v[a]=aprl::SI::c*p.properVelocity[a]*inverseGamma;
    const double force[3]={e[0]+v[1]*b[2]-v[2]*b[1],e[1]+v[2]*b[0]-v[0]*b[2],e[2]+v[0]*b[1]-v[1]*b[0]};
    for(unsigned a=0;a<3;++a)
      out[a]+=(p.position[a]-stats.x[a])*p.charge*force[a]/aprl::SI::c;
  }
  return {{static_cast<double>(out[0]/stats.mass),static_cast<double>(out[1]/stats.mass),static_cast<double>(out[2]/stats.mass)}};
}

void run(const std::string& root,const std::string& mode,double beta,const std::string& grid) {
  if(beta!=0.&&std::abs(beta-.23)>1e-12)throw std::runtime_error("beta must be 0 or .23");
  const bool cowanReference=mode=="cowan"||mode=="cowan-fieldonly";
  const bool staticSourceOnly=mode=="fieldonly"||mode=="cowan-fieldonly";
  if(staticSourceOnly&&beta!=0.)throw std::runtime_error("fieldonly requires beta=0 and J=0");
  const std::size_t n=grid=="base"?48:(grid=="wide"?64:(grid=="fine"?96:0));
  if(!n)throw std::runtime_error("grid must be base, wide, or fine");
  const double spacing=grid=="fine"?.5e-6:1e-6;
  const unsigned steps=grid=="fine"?30:15;
  const unsigned heldSteps=mode=="held50"?steps:(mode=="held100"?2*steps:0);
  double scale=1.;
  if(mode=="zero")scale=0.;else if(mode=="scale0.9")scale=.9;else if(mode=="scale1.1")scale=1.1;
  else if(mode!="qs"&&!staticSourceOnly&&!heldSteps&&!cowanReference)throw std::runtime_error("unknown initialization mode");
  const aprl::EBGridGeometry g(n,n,n,spacing,spacing,spacing,spacing/aprl::SI::c,aprl::EBMaxwellSolver::CowanZ);
  aprl::EBFieldGrid fields(g);fields.clearFields();fields.clearCurrent();
  FieldVector<double> origin(-.5*n*spacing);
  aprl::SimulationConfig config;
  config.beam.type=aprl::BeamInputType::GeneratedGaussian;
  config.beam.macroparticles=count;config.beam.electrons=-totalCharge/aprl::SI::elementaryCharge;
  config.beam.randomSeed=seed;config.beam.gamma=1/std::sqrt(1-beta*beta);
  config.beam.direction=FieldVector<double>(0.);config.beam.direction[2]=1.;
  config.beam.center=FieldVector<double>(0.);config.beam.positionOffset=FieldVector<double>(0.);
  config.beam.sigmaPosition=FieldVector<double>(sigma);config.beam.sigmaPosition[2]=sigma/config.beam.gamma;
  config.beam.sigmaProperVelocity=FieldVector<double>(0.);
  aprl::ParticleInitializationReport particleReport;
  auto particles=aprl::ParticleInitializer::create(config,MPI_COMM_WORLD,particleReport);
  for(auto& p:particles) {
    p.id=p.sourceId;p.frozen=false;
    for(unsigned a=0;a<3;++a)
      if(!(p.position[a]>origin[a]+spacing&&p.position[a]<-origin[a]-spacing))
        throw std::runtime_error("Particle lacks one-cell initial guard");
  }
  const auto initial=particles;particlesCsv(root+"/initial-particles.csv",initial);
  // Prescribed negative-time trajectories use the same fixed release state.
  // At the last prehistory step remaining==0, so the actual deposited endpoint
  // is exactly p0; there is no post-deposition position reset or momentum kick.
  const auto prescribedPosition = [&](std::size_t index,unsigned remaining) {
    FieldVector<double> position(initial[index].position);
    const double inverseGamma=1/static_cast<double>(gamma(initial[index].properVelocity));
    for(unsigned a=0;a<3;++a)position[a]-=
      remaining*g.dt*aprl::SI::c*initial[index].properVelocity[a]*inverseGamma;
    return position;
  };
  if(heldSteps) {
    for(std::size_t i=0;i<particles.size();++i) {
      particles[i].position=prescribedPosition(i,heldSteps);
      for(unsigned a=0;a<3;++a)
        if(!(particles[i].position[a]>origin[a]+spacing&&particles[i].position[a]<-origin[a]-spacing))
          throw std::runtime_error("Held prehistory reaches the particle guard");
    }
    particlesCsv(root+"/held-start-particles.csv",particles);
  }
  const std::size_t guards[3]={0,0,0};
  aprl::EBGaussInitializationReport initialReport;
  aprl::verification::CowanReferenceReport cowanReport={};
  const auto initialBegin=std::chrono::steady_clock::now();
  if(cowanReference)cowanReport=aprl::verification::initializeCowanReference(fields,origin,particles,beta);
  else if(scale!=0.)initialReport=aprl::EBGaussFieldInitializer::initialize(fields,g,0,origin,particles,
    aprl::InitialSelfFieldModel::RelativisticPoisson,1e-8,guards,1e-10,3000,MPI_COMM_WORLD);
  const double initialSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-initialBegin).count();
  for(YeeComponent* c:{&fields.ex(),&fields.ey(),&fields.ez(),&fields.bx(),&fields.by(),&fields.bz()})
    for(std::size_t i=0;i<c->size();++i)c->data()[i]*=scale;
  aprl::ChargeConservingCurrentDepositor depositor(fields,origin);
  YeeComponent rho(g.nx+1,g.ny+1,g.nz+1),nextRho(g.nx+1,g.ny+1,g.nz+1);
  rho.fill(0.);for(const auto& p:particles)depositor.depositCharge(p.position,p.charge,rho);
  long double heldGridWork=0;
  double heldContinuityMax=0;
  const auto heldBegin=std::chrono::steady_clock::now();
  if(heldSteps) {
    auto hold=fresh(root+"/held-history.csv");
    hold<<"step,time_s,UE_J,UB_native_minus_half_J,Gauss_relative_max,Gauss_relative_L2,continuity_relative_rho_over_dt,JdotE_cumulative_J\n";
    double continuityRelative=0;
    for(unsigned h=0;h<=heldSteps;++h) {
      const auto beforeB=copyComponents(fields.bx(),fields.by(),fields.bz());
      const auto en=energy(fields,beforeB);const auto gc=gauss(fields,rho);
      hold<<h<<','<<-static_cast<double>(heldSteps-h)*g.dt<<','<<en.electric<<','<<en.magneticNative
          <<','<<gc.relativeMax<<','<<gc.relativeL2<<','<<continuityRelative<<','<<heldGridWork<<'\n';
      if(h==heldSteps)break;
      const auto oldE=copyComponents(fields.ex(),fields.ey(),fields.ez());
      fields.advanceMagnetic();fields.clearCurrent();nextRho.fill(0.);
      for(std::size_t i=0;i<particles.size();++i) {
        const auto before=particles[i].position;
        particles[i].position=prescribedPosition(i,heldSteps-h-1);
        depositor.depositSegment(before,particles[i].position,particles[i].charge);
        depositor.depositCharge(particles[i].position,particles[i].charge,nextRho);
      }
      const double residual=aprl::ChargeConservingCurrentDepositor::maxContinuityResidual(fields,rho,nextRho);
      double sourceScale=0;for(std::size_t i=0;i<rho.size();++i)
        sourceScale=std::max(sourceScale,std::max(std::abs(rho.data()[i]),std::abs(nextRho.data()[i]))/g.dt);
      continuityRelative=residual/sourceScale;heldContinuityMax=std::max(heldContinuityMax,continuityRelative);
      fields.advanceElectric();heldGridWork+=currentWork(fields,oldE);
      std::copy(nextRho.data(),nextRho.data()+nextRho.size(),rho.data());
    }
    for(std::size_t i=0;i<particles.size();++i)for(unsigned a=0;a<3;++a)
      if(particles[i].position[a]!=initial[i].position[a]||particles[i].properVelocity[a]!=initial[i].properVelocity[a])
        throw std::runtime_error("Held trajectory did not reproduce exact release x/p");
    if(!hold)throw std::runtime_error("Held history output failed");
    fields.clearCurrent();
  }
  const double heldSeconds=heldSteps?std::chrono::duration<double>(std::chrono::steady_clock::now()-heldBegin).count():0.;
  particlesCsv(root+"/release-particles.csv",particles);
  const aprl::SIFieldSourceSet emptySources;aprl::BoostFrameTransform identity;
  identity.setFromGamma(1.,aprl::SI::c,0.);
  auto csv=fresh(root+"/history.csv");
  csv<<"step,time_s,K_J,delta_K_J,total_charge_C,mass_kg";
  for(const auto* name:{"mean_x_m","mean_y_m","mean_z_m","mean_ux","mean_uy","mean_uz","var_x_m2","var_y_m2","var_z_m2","var_ux","var_uy","var_uz","cov_x_ux_m","cov_y_uy_m","cov_z_uz_m"})csv<<','<<name;
  csv<<",UE_J,UB_centered_J,UB_native_plus_half_J,Gauss_max_V_m2,Gauss_relative_max,Gauss_relative_L2,continuity_max_C_m3_s,continuity_relative_rho_over_dt,JdotE_step_J,JdotE_cumulative_J\n";
  long double cumulativeWork=0,stepWork=0;double continuity=0,continuityRelative=0;
  std::array<double,3> analyticMoment={{0,0,0}},gridMoment={{0,0,0}};
  aprl::verification::CowanMagneticDivergenceReport initialMagneticDivergence={},finalMagneticDivergence={};
  for(unsigned step=0;step<=steps;++step) {
    auto oldB=copyComponents(fields.bx(),fields.by(),fields.bz());
    fields.advanceMagnetic(); // At the final step this is diagnostic centering only.
    const auto s=moments(particles,initial);const auto en=energy(fields,oldB);const auto gc=gauss(fields,rho);
    if(step==0){analyticMoment=forceMoment(particles,fields,origin,true,beta);gridMoment=forceMoment(particles,fields,origin,false,beta);initialMagneticDivergence=aprl::verification::cowanMagneticDivergence(fields);}
    if(step==steps)finalMagneticDivergence=aprl::verification::cowanMagneticDivergence(fields);
    csv<<step<<','<<step*g.dt<<','<<s.kinetic<<','<<s.deltaK<<','<<s.charge<<','<<s.mass;
    for(const long double* group:{s.x,s.u,s.varianceX,s.varianceU,s.xu})for(unsigned a=0;a<3;++a)csv<<','<<group[a];
    csv<<','<<en.electric<<','<<en.magneticCentered<<','<<en.magneticNative<<','<<gc.max<<','<<gc.relativeMax<<','<<gc.relativeL2<<','<<continuity<<','<<continuityRelative<<','<<stepWork<<','<<cumulativeWork<<'\n';
    if(step==steps)break;
    auto oldE=copyComponents(fields.ex(),fields.ey(),fields.ez());
    fields.clearCurrent();nextRho.fill(0.);
    for(auto& p:particles) {
      const auto before=p.position;
      if(!staticSourceOnly)
        aprl::RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(p,fields,origin,emptySources,identity,step*g.dt,g.dt,1);
      for(unsigned a=0;a<3;++a)if(!(p.position[a]>origin[a]+spacing&&p.position[a]<-origin[a]-spacing))
        throw std::runtime_error("Particle reached guard; this probe forbids boundary crossing");
      depositor.depositSegment(before,p.position,p.charge);depositor.depositCharge(p.position,p.charge,nextRho);
    }
    continuity=aprl::ChargeConservingCurrentDepositor::maxContinuityResidual(fields,rho,nextRho);
    double sourceScale=0;for(std::size_t i=0;i<rho.size();++i)sourceScale=std::max(sourceScale,std::max(std::abs(rho.data()[i]),std::abs(nextRho.data()[i]))/g.dt);
    continuityRelative=continuity/sourceScale;
    fields.advanceElectric();stepWork=currentWork(fields,oldE);cumulativeWork+=stepWork;
    std::copy(nextRho.data(),nextRho.data()+nextRho.size(),rho.data());
  }
  particlesCsv(root+"/final-particles.csv",particles);
  auto json=fresh(root+"/metadata.json");
  const double pi=std::acos(-1.);
  const double gaussianEnergy=totalCharge*totalCharge/(8*std::pow(pi,1.5)*aprl::SI::epsilon0*sigma);
  const double restMass=moments(initial,initial).mass;
  json<<"{\n\"mode\":\""<<mode<<"\",\"grid\":\""<<grid<<"\",\"beta\":"<<beta<<",\"gamma\":"<<config.beam.gamma<<",\"particles\":"<<count<<",\"seed\":"<<seed
      <<",\"charge_C\":"<<totalCharge<<",\"rest_sigma_m\":"<<sigma<<",\"nx\":"<<n<<",\"spacing_m\":"<<spacing<<",\"dt_s\":"<<g.dt<<",\"steps\":"<<steps<<",\"initial_field_scale\":"<<scale
      <<",\"poisson_iterations\":"<<initialReport.iterations<<",\"poisson_relative_residual\":"<<initialReport.relativeResidual
      <<",\"poisson_unscaled_collocated_UE_J\":"<<initialReport.initialElectricEnergy<<",\"poisson_unscaled_collocated_UB_J\":"<<initialReport.initialMagneticEnergy
      <<",\"gaussian_rest_self_energy_J\":"<<gaussianEnergy<<",\"continuum_cov_xu_slope_m_per_s\":"<<gaussianEnergy/(3*restMass*aprl::SI::c*config.beam.gamma)
      <<",\"analytic_sampled_force_covariance_m_per_s\":["<<analyticMoment[0]<<','<<analyticMoment[1]<<','<<analyticMoment[2]<<"]"
      <<",\"grid_sampled_force_covariance_m_per_s\":["<<gridMoment[0]<<','<<gridMoment[1]<<','<<gridMoment[2]<<"]"
      <<",\"static_source_only\":"<<(staticSourceOnly?"true":"false")
      <<",\"held_prehistory_steps\":"<<heldSteps<<",\"held_prehistory_duration_s\":"<<heldSteps*g.dt
      <<",\"held_grid_work_J\":"<<heldGridWork<<",\"held_max_continuity_relative\":"<<heldContinuityMax
      <<",\"initialization_wall_s\":"<<initialSeconds<<",\"held_prehistory_wall_s\":"<<heldSeconds
      <<",\"cowan_reference\":"<<(cowanReference?"true":"false")
      <<",\"cowan_transforms\":"<<cowanReport.transforms
      <<",\"cowan_relative_gauss_L2\":"<<cowanReport.relativeGaussL2
      <<",\"cowan_maximum_gauss_residual\":"<<cowanReport.maximumGaussResidual
      <<",\"cowan_minimum_denominator\":"<<cowanReport.minimumDenominator
      <<",\"cowan_minimum_smoothing_symbol\":"<<cowanReport.minimumSmoothingSymbol
      <<",\"cowan_minimum_longitudinal_symbol\":"<<cowanReport.minimumLongitudinalSymbol
      <<",\"cowan_source_charge_C\":"<<cowanReport.sourceCharge
      <<",\"magnetic_divergence_samples\":"<<initialMagneticDivergence.samples
      <<",\"release_ordinary_divB_max_T_per_m\":"<<initialMagneticDivergence.ordinaryMaximum
      <<",\"release_modified_divB_max_T_per_m\":"<<initialMagneticDivergence.modifiedMaximum
      <<",\"release_ordinary_divB_relative_L2\":"<<initialMagneticDivergence.ordinaryRelativeL2
      <<",\"release_modified_divB_relative_L2\":"<<initialMagneticDivergence.modifiedRelativeL2
      <<",\"final_ordinary_divB_max_T_per_m\":"<<finalMagneticDivergence.ordinaryMaximum
      <<",\"final_modified_divB_max_T_per_m\":"<<finalMagneticDivergence.modifiedMaximum
      <<",\"final_ordinary_divB_relative_L2\":"<<finalMagneticDivergence.ordinaryRelativeL2
      <<",\"final_modified_divB_relative_L2\":"<<finalMagneticDivergence.modifiedRelativeL2
      <<",\"initial_field_model\":\""<<(cowanReference?"Cowan-compatible semidiscrete rigid reference":
          (heldSteps?"ordinary relativistic-Poisson followed by prescribed held history":
          (scale==0.?"zero-field negative control":"ordinary relativistic-Poisson")))<<"\""
      <<",\"held_reference_scope\":\"prescribed straight trajectories and production current deposition; synchronized t0 release, not production lab-plane crossing events; held source work excluded from release deltaK\""
      <<",\"completed\":true,\"physical_acceptance\":\"not_assessed_by_probe\",\"boundary\":\"PEC, zero-potential initialization; assess domain convergence\",\"release\":\"all at t0 except fieldonly, which fixes stationary sources with J=0; no held lab-plane interval\",\"energy_note\":\"Yee dual-volume E energy and centered B=(Bminus+Bplus)/2; not an exact Cowan discrete invariant\",\"control_note\":\"zero and scaled fields deliberately violate initial Gauss, not admissible candidates\",\"continuity_normalization\":\"max(abs(rho_before),abs(rho_after))/dt\",\"final_B_update\":\"diagnostic time centering only; particles/E stop at stated time\"}\n";
  if(!csv||!json)throw std::runtime_error("Output write failed");
  std::cout<<"startup_self_field_probe completed mode="<<mode<<" beta="<<beta<<" grid="<<grid<<" steps="<<steps<<"; physical acceptance requires independent review\n";
}
}

int main(int argc,char** argv) {
  int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int status=0,ranks;MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  try {
    if(ranks!=1||argc!=5)throw std::runtime_error("Use one MPI rank: OUT MODE BETA GRID");
    aprl::configureComputeThreads(provided,1);run(argv[1],argv[2],std::stod(argv[3]),argv[4]);
  } catch(const std::exception& e){std::cerr<<e.what()<<'\n';status=1;}
  MPI_Finalize();return status;
}
