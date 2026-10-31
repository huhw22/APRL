// Diagnostic executable linked to the unchanged accepted particle/source core.
// No Maxwell advance, grid feedback, charge deposit, CPML or energy correction.
#include "eb_particles.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace aprl;
using Vec=FieldVector<Double>;
long double gamma(const Vec& u) {
  return std::sqrt(1.L+(long double)u[0]*u[0]+(long double)u[1]*u[1]+(long double)u[2]*u[2]);
}
long double difference(const Vec& a,const Vec& b) {
  long double v=0;for(unsigned k=0;k<3;++k)v+=((long double)a[k]-b[k])*((long double)a[k]+b[k]);
  return v/(gamma(a)+gamma(b));
}
int main(int argc,char** argv) {
  try {
    if(argc!=7)throw std::runtime_error("Need sample index, x, y, ux, uy, uz");
    const unsigned long sample=std::strtoul(argv[1],nullptr,10);
    Vec initialLab(0.);for(unsigned k=0;k<3;++k)initialLab[k]=std::strtod(argv[k+4],nullptr);
    const Double dz=1.2543262415359034e-6,dt=dz/SI::c;
    BoostFrameTransform frame;frame.setOriginsFromGamma(929.1305492858543,SI::c,0.,-.05,0.);
    EBGridGeometry geometry(8,8,8,5e-6,5e-6,dz,dt,EBMaxwellSolver::CowanZ);
    EBFieldGrid dummyGrid(geometry);const Vec origin(0.);
    std::cout<<std::setprecision(21);
    for(bool magnet:{false,true}) for(unsigned substeps:{1u,5u,10u,20u}) {
      SIFieldSourceSet sources;
      if(magnet) {
        auto element=SIMagneticElement::planarUndulatorFromK(1.092714657005,.0234,0.,170,std::acos(-1.)/2);
        element.gaussianFringe=true;element.fringeRelativeCutoff=1e-9;element.prepare();
        sources.addPrescribedLabMagnet(element);
      }
      RelativisticParticleSI particle;
      particle.position=Vec(0.);particle.position[0]=std::strtod(argv[2],nullptr);particle.position[1]=std::strtod(argv[3],nullptr);
      frame.properVelocityLabToBox(initialLab,particle.properVelocity);
      particle.charge=-SI::elementaryCharge;particle.mass=SI::electronMass;particle.weight=1;particle.id=particle.sourceId=sample+1;
      Double time=0.;unsigned steps=0;long double maximum=0.;Vec lab(initialLab);
      while(frame.labZFromBoxZT(particle.position[2],time)<4.5) {
        if(steps>20000)throw std::runtime_error("Probe step bound exceeded");
        RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
          particle,dummyGrid,origin,sources,frame,time,dt,substeps,nullptr,false,nullptr,nullptr);
        time+=dt;++steps;frame.properVelocityBoxToLab(particle.properVelocity,lab);
        maximum=std::max(maximum,std::abs(difference(lab,initialLab)));
      }
      const long double mc2eV=(long double)SI::electronMass*SI::c*SI::c/SI::elementaryCharge;
      std::cout<<"{\"sample_index\":"<<sample<<",\"magnet\":"<<(magnet?"true":"false")
        <<",\"particle_substeps\":"<<substeps<<",\"field_dt_s\":"<<dt
        <<",\"steps\":"<<steps<<",\"lab_z_final_m\":"<<frame.labZFromBoxZT(particle.position[2],time)
        <<",\"lab_energy_change_eV\":"<<difference(lab,initialLab)*mc2eV
        <<",\"max_sampled_absolute_lab_energy_change_eV\":"<<maximum*mc2eV
        <<",\"scope\":\"Unchanged external-only pusher; static lab magnetic field requires constant lab energy; no Maxwell/Cowan advance\"}\n";
    }
    return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
