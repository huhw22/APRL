#ifndef APRL_LONGITUDINAL_COUPLING_DIAGNOSTIC_H
#define APRL_LONGITUDINAL_COUPLING_DIAGNOSTIC_H
#include <limits>
// Test-only observations. No diagnostic below changes E/B, p, or particle positions.
namespace longitudinal_diagnostic {
using LD=long double;
std::size_t index(long long i,long long j,long long k){return (wrapIndex(k,nz)-padding)*nx*ny+(wrapIndex(j,ny)-padding)*nx+wrapIndex(i,nx)-padding;}
struct Nodes {
 std::vector<LD> rho,gauss;
 std::vector<Double> force;
 LD chargeScale=0.,gaussL1=0.,forceGauss=0.,forceCharge=0.;
};
Nodes nodes(EBFieldGrid& topology,const Components& E,const std::vector<RelativisticParticleSI>& particles){
 const auto& g=topology.geometry();const LD dv=(LD)g.dx*g.dy*g.dz;Nodes n;n.rho.assign(nx*ny*nz,0.);n.gauss.resize(n.rho.size());n.force.resize(n.rho.size());
 YeeComponent rho(g.nx+1,g.ny+1,g.nz+1);ChargeConservingCurrentDepositor dep(topology,Vec(0.));
 for(const auto& p:particles){dep.depositCharge(p.position,p.charge,rho);n.chargeScale+=std::abs((LD)p.charge);}
 for(std::size_t k=0;k<rho.nz();++k)for(std::size_t j=0;j<rho.ny();++j)for(std::size_t i=0;i<rho.nx();++i)if(rho(i,j,k)!=0.)n.rho[index(i,j,k)]+=rho(i,j,k);
 const auto& m=matchedMetric(g);auto raw=[&](unsigned a,long long i,long long j,long long k){return E[a](wrapIndex(i,nx),wrapIndex(j,ny),wrapIndex(k,nz));};
 for(unsigned k=padding;k<padding+nz;++k)for(unsigned j=padding;j<padding+ny;++j)for(unsigned i=padding;i<padding+nx;++i){
  const auto at=index(i,j,k);n.gauss[at]=SI::epsilon0*((raw(0,i,j,k)-raw(0,i-1,j,k))/g.dx+(raw(1,i,j,k)-raw(1,i,j-1,k))/g.dy+(raw(2,i,j,k)-raw(2,i,j,k-1))/g.dz);
  n.force[at]=forceMode==0?.5*(m.apply(m.electric[2],raw,2,i,j,k-1)+m.apply(m.electric[2],raw,2,i,j,k)):.5*(raw(2,i,j,k-1)+raw(2,i,j,k));n.gaussL1+=std::abs(n.gauss[at]-n.rho[at])*dv;n.forceGauss+=n.gauss[at]*n.force[at]*dv;n.forceCharge+=n.rho[at]*n.force[at]*dv;
 }return n;
}
Double sample(const std::vector<Double>& values,const Vec& p,const EBGridGeometry& g){
 const Double ds[3]={g.dx,g.dy,g.dz};long long cell[3];Double w[3][2];for(unsigned a=0;a<3;++a){const Double u=p[a]/ds[a];cell[a]=(long long)std::floor(u);w[a][1]=u-cell[a];w[a][0]=1-w[a][1];}Double out=0.;
 for(unsigned k=0;k<2;++k)for(unsigned j=0;j<2;++j)for(unsigned i=0;i<2;++i)out+=w[0][i]*w[1][j]*w[2][k]*values[index(cell[0]+i,cell[1]+j,cell[2]+k)];return out;
}
Double orbit(const std::vector<Double>& values,const Vec& before,const Vec& after,const EBGridGeometry& g){
 const Double ds[3]={g.dx,g.dy,g.dz};Double cuts[32]={0.,1.};unsigned count=2;Vec dx(after);dx-=before;
 for(unsigned a=0;a<3;++a){const Double x=before[a]/ds[a],d=dx[a]/ds[a];if(d==0.)continue;for(long long k=(long long)std::floor(std::min(x,x+d))+1;k<std::max(x,x+d);++k){if(count==32)throw std::runtime_error("Diagnostic orbit split bound");const Double t=(k-x)/d;if(t>0.&&t<1.)cuts[count++]=t;}}
 std::sort(cuts,cuts+count);Double value=0.;for(unsigned t=1;t<count;++t){const Double mid=.5*(cuts[t-1]+cuts[t]),width=cuts[t]-cuts[t-1];Vec a(before),b(before);a.pmv(mid-width/std::sqrt(12.),dx);b.pmv(mid+width/std::sqrt(12.),dx);value+=.5*width*(sample(values,a,g)+sample(values,b,g));}return value;
}
struct Step {
 LD space=0.,time=0.,gauss=0.,adjoint=0.,continuityL1=0.,gaussL1=0.;
 LD endForce=0.,pathForce=0.,particleForce=0.;
 LD impossible=0.,trials=0.,maxTransverseKick=0.,maxRequiredTransverseField=0.,maxAbsDeltaUz=0.;
};
Step observe(EBFieldGrid& f,const Components& old,const std::vector<RelativisticParticleSI>& before,const std::vector<RelativisticParticleSI>& after,const std::vector<Vec>& usedE,LD actualFieldImpulse){
 const auto& g=f.geometry();const auto a=nodes(f,old,before),b=nodes(f,copyElectric(f),after);std::vector<Double> mean=a.force;for(std::size_t k=0;k<mean.size();++k)mean[k]=.5*(a.force[k]+b.force[k]);Step d;
 for(std::size_t t=0;t<before.size();++t){
  const auto& p=before[t];const auto& q=after[t];const LD factor=(LD)p.charge*SI::c*g.dt;const Double endpoint=.5*(sample(a.force,p.position,g)+sample(b.force,q.position,g)),path=orbit(mean,p.position,q.position,g);
  d.space+=factor*(usedE[t][2]-path);d.time+=factor*(path-endpoint);d.endForce+=factor*endpoint;d.pathForce+=factor*path;d.particleForce+=factor*usedE[t][2];
  // Necessary one-step feasibility only: fixed original end energy, prescribed
  // change in longitudinal impulse, smallest possible transverse momentum change.
  // This is NOT a proposed post-push correction and is never applied.
  const LD du=(LD)p.charge*g.dt/(p.mass*SI::c)*(endpoint-usedE[t][2]);const LD uz=q.properVelocity[2],perp=(LD)q.properVelocity[0]*q.properVelocity[0]+(LD)q.properVelocity[1]*q.properVelocity[1];const LD exchange=-du*(2*uz+du),target=perp+exchange,tol=128*std::numeric_limits<Double>::epsilon()*(perp+std::abs(exchange));d.trials+=1.;d.maxAbsDeltaUz=std::max(d.maxAbsDeltaUz,std::abs(du));
  if(target < -tol)d.impossible+=1.;else d.maxTransverseKick=std::max(d.maxTransverseKick,std::abs(std::sqrt(std::max(0.L,target))-std::sqrt(perp)));
  const LD dx=(LD)q.position[0]-p.position[0],dy=(LD)q.position[1]-p.position[1],dz=(LD)q.position[2]-p.position[2],dr=std::sqrt(dx*dx+dy*dy);if(dr>0.)d.maxRequiredTransverseField=std::max(d.maxRequiredTransverseField,std::abs(dz*(endpoint-usedE[t][2]))/dr);
 }
 d.gauss=d.endForce-actualFieldImpulse;d.adjoint=d.endForce-.5L*SI::c*g.dt*(a.forceCharge+b.forceCharge);d.gaussL1=std::max(a.gaussL1/a.chargeScale,b.gaussL1/b.chargeScale);
 auto current=[&](unsigned axis,long long i,long long j,long long k){const YeeComponent* v[3]={&f.jx(),&f.jy(),&f.jz()};return (*v[axis])(wrapIndex(i,nx),wrapIndex(j,ny),wrapIndex(k,nz));};const LD dv=(LD)g.dx*g.dy*g.dz;
 for(unsigned k=padding;k<padding+nz;++k)for(unsigned j=padding;j<padding+ny;++j)for(unsigned i=padding;i<padding+nx;++i){const auto at=index(i,j,k);const LD div=(current(0,i,j,k)-current(0,i-1,j,k))/g.dx+(current(1,i,j,k)-current(1,i,j-1,k))/g.dy+(current(2,i,j,k)-current(2,i,j,k-1))/g.dz;d.continuityL1+=std::abs(b.rho[at]-a.rho[at]+g.dt*div)*dv/a.chargeScale;}return d;
}
struct Accumulator {
 LD space=0.,time=0.,gauss=0.,maxGaussL1=0.,maxContinuityL1=0.,maxAdjointJ=0.,impossible=0.,trials=0.,maxTransverseKick=0.,maxRequiredTransverseField=0.,maxAbsDeltaUz=0.;
 void add(const Step& s){space+=s.space;time+=s.time;gauss+=s.gauss;maxGaussL1=std::max(maxGaussL1,s.gaussL1);maxContinuityL1=std::max(maxContinuityL1,s.continuityL1);maxAdjointJ=std::max(maxAdjointJ,std::abs(s.adjoint));impossible+=s.impossible;trials+=s.trials;maxTransverseKick=std::max(maxTransverseKick,s.maxTransverseKick);maxRequiredTransverseField=std::max(maxRequiredTransverseField,s.maxRequiredTransverseField);maxAbsDeltaUz=std::max(maxAbsDeltaUz,s.maxAbsDeltaUz);}
 void row(std::ostream& out,const Step& s)const{out<<','<<space<<','<<time<<','<<gauss<<','<<s.gaussL1<<','<<s.continuityL1<<','<<s.adjoint<<','<<s.impossible<<','<<s.trials<<','<<s.maxTransverseKick<<','<<s.maxRequiredTransverseField<<','<<s.maxAbsDeltaUz<<'\n';}
 void summary(std::ostream& out)const{out<<std::setprecision(21)<<"LONGITUDINAL_DIAGNOSTICS {\"space_gap_J\":"<<space<<",\"time_gap_J\":"<<time<<",\"Gauss_gap_J\":"<<gauss<<",\"max_Gauss_charge_L1_over_total_abs_charge\":"<<maxGaussL1<<",\"max_continuity_L1_over_total_abs_charge\":"<<maxContinuityL1<<",\"max_charge_force_adjoint_error_J\":"<<maxAdjointJ<<",\"frozen_energy_transverse_redirection_impossible_fraction\":"<<impossible/trials<<",\"max_minimum_required_transverse_du\":"<<maxTransverseKick<<",\"max_work_neutral_transverse_delta_E_V_per_m\":"<<maxRequiredTransverseField<<",\"max_required_longitudinal_du\":"<<maxAbsDeltaUz<<",\"particle_force_or_state_changed\":false}\n";}
};
void initializeSpaceCharge(EBFieldGrid& f,const std::vector<RelativisticParticleSI>& particles){
 for(auto* v:{&f.ex(),&f.ey(),&f.ez(),&f.bx(),&f.by(),&f.bz()})std::fill(v->data(),v->data()+v->size(),0.);const auto& g=f.geometry();const auto n=nodes(f,copyElectric(f),particles);std::vector<LD> rho(nz),E(nz);LD total=0.;
 for(unsigned z=0;z<nz;++z){for(unsigned y=0;y<ny;++y)for(unsigned x=0;x<nx;++x)rho[z]+=n.rho[(z*ny+y)*nx+x]/(nx*ny);total+=rho[z];}
 LD value=0.,mean=0.;for(unsigned z=0;z<nz;++z){value+=(rho[z]-total/nz)*g.dz/SI::epsilon0;E[z]=value;mean+=value/nz;}
 for(unsigned z=padding;z<padding+nz;++z)for(unsigned y=padding;y<padding+ny;++y)for(unsigned x=padding;x<padding+nx;++x)f.ez()(x,y,z)=E[z-padding]-mean;wrapE(f);wrapB(f);
}
std::vector<Double> oldChargeForce,newChargeForce;
Vec testOrbitElectric(const EBFieldGrid& f,const Components* old,const Vec& origin,const Vec& begin,const Vec& end){
 if(forceMode==0)return matchedOrbitElectric(f,old,origin,begin,end);
 if(forceMode==1)return orbitElectric(f,old,origin,begin,end);
 Vec value(0.);value[2]=.5*(sample(oldChargeForce,begin,f.geometry())+sample(newChargeForce,end,f.geometry()));return value;
}
}
#endif
