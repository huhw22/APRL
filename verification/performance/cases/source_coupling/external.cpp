#include "source_coupling.h"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace source_trial;
int main(int argc,char** argv) {
 try {
  if(argc!=2)throw std::runtime_error("Need output CSV");
  std::ofstream out(argv[1]);out<<std::setprecision(21)<<"sample,method,substeps,field_steps,final_lab_energy_change_eV,max_lab_energy_change_eV,x_at_4p5m,y_at_4p5m,lab_ux,lab_uy,lab_uz\n";
  const Double dz=1.2543262415359034e-6,dt=dz/SI::c;
  BoostFrameTransform frame;frame.setOriginsFromGamma(929.1305492858543,SI::c,0.,-.05,0.);
  EBGridGeometry g(8,8,8,5e-6,5e-6,dz,dt,EBMaxwellSolver::CowanZ);EBFieldGrid dummy(g);Vec origin(0.);
  SIFieldSourceSet sources;auto magnet=SIMagneticElement::planarUndulatorFromK(1.092714657005,.0234,0.,170,std::acos(-1.)/2);magnet.gaussianFringe=true;magnet.fringeRelativeCutoff=1e-9;magnet.prepare();sources.addPrescribedLabMagnet(magnet);
  for(unsigned sample=0;sample<3;++sample)for(bool coupled:{false,true})for(unsigned substeps:{5u,10u,20u}) {
   Vec initialLab(0.);initialLab[0]=sample==1?.002:0.;initialLab[1]=sample==2?-.003:0.;initialLab[2]=std::sqrt(1174.171072271572*1174.171072271572-1-initialLab[0]*initialLab[0]-initialLab[1]*initialLab[1]);
   RelativisticParticleSI p;p.position=Vec(0.);p.position[0]=(sample==1?30e-6:0.);p.position[1]=(sample==2?-40e-6:0.);frame.properVelocityLabToBox(initialLab,p.properVelocity);const Vec initial(p.properVelocity);p.charge=-SI::elementaryCharge;p.mass=SI::electronMass;p.weight=1;
   Double t=0.;unsigned steps=0;long double maximum=0.;Vec eventX(0.),finalLab(0.);
   while(true) {
    const auto before=p;const Double z0=frame.labZFromBoxZT(p.position[2],t);
    if(coupled)push(p,Vec(0.),Vec(0.),sources,frame,t,dt,substeps);
    else RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(p,dummy,origin,sources,frame,t,dt,substeps,nullptr,false,nullptr,nullptr);
    t+=dt;++steps;if(steps>20000)throw std::runtime_error("External path exceeded bound");
    const auto change=(long double)SI::electronMass*SI::c*SI::c/SI::elementaryCharge*frame.gamma()*(deltaGamma(initial,p.properVelocity)+(long double)frame.beta()*((long double)p.properVelocity[2]-initial[2]));maximum=std::max(maximum,std::abs(change));
    const Double z1=frame.labZFromBoxZT(p.position[2],t);
    if(z1>=4.5){const Double fraction=(4.5-z0)/(z1-z0);for(unsigned a=0;a<3;++a)eventX[a]=before.position[a]+fraction*(p.position[a]-before.position[a]);frame.properVelocityBoxToLab(p.properVelocity,finalLab);out<<sample<<','<<(coupled?"mean_velocity":"original_boris")<<','<<substeps<<','<<steps<<','<<change<<','<<maximum<<','<<eventX[0]<<','<<eventX[1]<<','<<finalLab[0]<<','<<finalLab[1]<<','<<finalLab[2]<<'\n';break;}
   }
  }
  return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
