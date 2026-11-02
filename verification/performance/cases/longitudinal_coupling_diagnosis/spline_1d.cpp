// Independent 1-D reduction, not an APRL production depositor/pusher change.
// Original mean-velocity momentum kernel; pairwise-neutral periodic sheets.
#include "source_coupling_math.h"
#include "compute_threads.h"
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace aprl;using namespace source_coupling_math;
using LD=long double;using Array=std::vector<Double>;
const Double Area=60e-6*50e-6,Length=24*1.2543262415359034e-6;
struct Particle {Double z,u,q,m;};
Double cardinal(Double x,unsigned m){const Double t=.5*(m+1)-std::abs(x);if(t<=0.)return 0.;Double sum=0.,choose=1.,factorial=1.;for(unsigned k=2;k<=m;++k)factorial*=k;for(unsigned k=0;k<=m+1&&k<t;++k){sum+=(k%2?-1.:1.)*choose*std::pow(t-k,(int)m);choose*=Double(m+1-k)/(k+1);}return sum/factorial;}
Double shape(Double x,unsigned order){if(order==5)return cardinal(x,5);x=std::abs(x);if(order==1)return std::max(0.,1.-x);if(x<1.)return (4.-6*x*x+3*x*x*x)/6.;if(x<2.)return std::pow(2.-x,3)/6.;return 0.;}
Double cdf(Double x,unsigned order){
 if(order==0)return std::max(0.,std::min(1.,x+.5));
 if(x<=-1.5)return 0.;if(x<-.5)return std::pow(x+1.5,3)/6.;if(x<=.5)return .5+.75*x-x*x*x/3.;if(x<1.5)return 1.-std::pow(1.5-x,3)/6.;return 1.;
}
Double currentShape(Double x,unsigned order){if(order==4)return cardinal(x,4);if(order==0)return x>=-.5&&x<.5?1.:0.;x=std::abs(x);if(x<=.5)return .75-x*x;if(x<1.5)return .5*std::pow(1.5-x,2);return 0.;}
Array compensate(const Array& x,unsigned order,unsigned degree){
 Array y=x;const unsigned n=x.size();if(!order)return y;
 const Double power=order>=3?degree:degree-1,a=power/6.,b=order==1?0.:3*power/40.+power*(power-1)/72.;
 const Double c=order==4?5*power/112.+power*(power-1)/80.+power*(power-1)*(power-2)/1296.:0.;
 const Double center=1+a/2+3*b/8+5*c/16,first=-a/4-b/4-15*c/64,second=b/16+3*c/32,third=-c/64;
 for(unsigned k=0;k<n;++k)y[k]=center*x[k]+first*(x[(k+1)%n]+x[(k+n-1)%n])+second*(x[(k+2)%n]+x[(k+n-2)%n])+third*(x[(k+3)%n]+x[(k+n-3)%n]);return y;
}
Array weights(Double a,Double b,unsigned cells,unsigned degree,Double dz){
 Array out(cells,0.);const Double x0=a/dz,x1=b/dz,d=x1-x0;Double cuts[8]={0.,1.};unsigned count=2;
 if(d!=0.)for(long long k=(long long)std::floor(std::min(x0,x1))+1;k<std::max(x0,x1);++k){if(count==8)throw std::runtime_error("Sheet orbit split bound");cuts[count++]=(k-x0)/d;}
 std::sort(cuts,cuts+count);
 for(unsigned t=1;t<count;++t){const Double mid=.5*(cuts[t-1]+cuts[t]),width=cuts[t]-cuts[t-1],u=x0+(mid-width/std::sqrt(12.))*d,v=x0+(mid+width/std::sqrt(12.))*d;
  for(unsigned k=0;k<cells;++k)for(int shift:{-1,0,1}){const Double center=k+.5+shift*(Double)cells;if(degree<5)out[k]+=.5*width*(currentShape(u-center,degree-1)+currentShape(v-center,degree-1));else {const Double m=x0+mid*d,offset=.5*width*std::sqrt(3./5.)*d;out[k]+=width*((5./18.)*(currentShape(m-offset-center,4)+currentShape(m+offset-center,4))+(4./9.)*currentShape(m-center,4));}}
 }return out;
}
Array chargePathWeights(Double a,Double b,unsigned cells,unsigned degree,Double dz){
 Array out(cells,0.);const Double x0=a/dz,x1=b/dz,d=x1-x0;Double cuts[8]={0.,1.};unsigned count=2;
 if(d!=0.)for(long long k=(long long)std::floor(std::min(x0,x1))+1;k<std::max(x0,x1);++k){if(count==8)throw std::runtime_error("Sheet orbit split bound");cuts[count++]=(k-x0)/d;}
 std::sort(cuts,cuts+count);
 for(unsigned t=1;t<count;++t){const Double mid=.5*(cuts[t-1]+cuts[t]),width=cuts[t]-cuts[t-1],u=x0+(mid-width/std::sqrt(12.))*d,v=x0+(mid+width/std::sqrt(12.))*d;
  for(unsigned k=0;k<cells;++k)for(int shift:{-1,0,1}){const Double center=k+shift*(Double)cells;if(degree<3)out[k]+=.5*width*(shape(u-center,degree)+shape(v-center,degree));else {const Double m=x0+mid*d,offset=.5*width*std::sqrt(3./5.)*d;out[k]+=width*((5./18.)*(shape(m-offset-center,degree)+shape(m+offset-center,degree))+(4./9.)*shape(m-center,degree));}}
 }return out;
}
Array charge(const std::vector<Particle>& p,unsigned cells,unsigned degree,Double dz,unsigned filter){
 Array rho(cells,0.);for(const auto& a:p)for(unsigned k=0;k<cells;++k)for(int s:{-1,0,1})rho[k]+=a.q/(Area*dz)*shape(a.z/dz-k-s*(Double)cells,degree);return compensate(rho,filter,degree);
}
std::pair<LD,LD> fundamental(const Array& E){LD a=0.,b=0.;const Double pi=std::acos(-1.);for(unsigned k=0;k<E.size();++k){const Double phase=2*pi*(k+.5)/E.size();a+=2.L/E.size()*E[k]*std::cos(phase);b-=2.L/E.size()*E[k]*std::sin(phase);}return {a,b};}
LD fieldEnergy(const Array& E,Double dz){LD out=0.;for(const auto& e:E)out+=(LD)e*e;return .5L*SI::epsilon0*Area*dz*out;}
Double dot(const Array& a,const Array& b){Double out=0.;for(unsigned k=0;k<a.size();++k)out+=a[k]*b[k];return out;}
int main(int argc,char** argv){
 MPI_Init(&argc,&argv);configureComputeThreads(MPI_THREAD_SINGLE,1);
 if(argc<8)throw std::runtime_error("output degree compensation longitudinal_refinement sheet_pairs phase initial_uz");
 const unsigned degree=std::stoul(argv[2]),filter=std::stoul(argv[3]),refine=std::stoul(argv[4]),pairs=std::stoul(argv[5]),cells=24*refine,duration=argc>8?std::stoul(argv[8]):1,steps=160*refine*duration;
 if((degree!=1&&degree!=3&&degree!=5)||filter>4)throw std::runtime_error("Unimplemented sheet basis");
 const Double phase=std::stod(argv[6]),u0=std::stod(argv[7]),amplitude=argc>9?std::stod(argv[9]):.15,dz=Length/cells,dt=dz/SI::c,represented=192e5/(2*pairs),pi=std::acos(-1.);std::vector<Particle> particles;
 for(unsigned k=0;k<pairs;++k)for(int sign:{-1,1}){const Double relative=(k+phase)*Length/pairs;Particle p;p.z=relative+sign*amplitude*(Length/24)*std::sin(2*pi*relative/Length);p.u=u0;p.q=sign*SI::elementaryCharge*represented;p.m=SI::electronMass*represented;particles.push_back(p);}
 const bool diagnostics=argc>10&&std::stoul(argv[10]);LD cumulativeSpace=0.,cumulativeTime=0.,cumulativeEndpoint=0.,cumulativeIteration=0.;
 const auto initialParticles=particles;auto rho=charge(particles,cells,degree,dz,filter);Array E(cells);LD rmean=0.;for(auto r:rho)rmean+=r/cells;Double integral=0.,mean=0.;for(unsigned k=0;k<cells;++k){integral+=(rho[k]-rmean)*dz/SI::epsilon0;E[k]=integral;mean+=integral/cells;}for(auto& x:E)x-=mean;
 const auto initialMode=fundamental(E);const LD omega=std::sqrt(192e5L*SI::elementaryCharge*SI::elementaryCharge/(Area*Length*SI::electronMass*SI::epsilon0)),linearAmplitude=192e5L*SI::elementaryCharge*amplitude*(Length/24)/(Area*Length*SI::epsilon0);const LD initial=fieldEnergy(E,dz),G=929.1305492858543L,beta=std::sqrt(1-1/(G*G)),volume=(LD)Area*dz;LD dK=0.,dP=0.,worstEnergy=0.,worstGauss=0.,worstContinuity=0.,worstWork=0.,maxP=0.,maximumActivity=0.;unsigned maxIterations=0;
 std::ofstream rows(argv[1]);rows<<std::setprecision(21)<<"step,particle_delta_Kprime_J,particle_delta_cPzprime_J,field_delta_U_J,energy_residual_J,mean_delta_uz,lab_energy_from_momentum_eV_per_particle,gauss_L1_relative,continuity_L1_relative,work_gap_J,individual_momentum_change_activity_J,iterations,fundamental_E_real_V_per_m,fundamental_E_imag_V_per_m,spatial_impulse_gap_J,time_impulse_gap_J,endpoint_charge_impulse_J,iteration_impulse_gap_J\n";
 for(unsigned step=1;step<=steps;++step){
  const auto old=particles;const auto E0=E;const auto rho0=rho;std::vector<Double> guess;for(const auto& p:old)guess.push_back(p.z+dt*SI::c*p.u/std::sqrt(1+p.u*p.u));std::vector<Particle> candidate;Array used(old.size()),J(cells);unsigned iterations=0;bool converged=false;
  for(unsigned trial=0;trial<12;++trial){
   Array meanE(cells);for(unsigned k=0;k<cells;++k)meanE[k]=.5*(E0[k]+E[k]);const auto filtered=compensate(meanE,filter,degree);candidate=old;J.assign(cells,0.);
   for(unsigned i=0;i<old.size();++i){
    const auto w=weights(old[i].z,guess[i],cells,degree,dz);used[i]=dot(w,filtered);Vec u(0.),ef(0.),bf(0.);u[2]=old[i].u;ef[2]=used[i];momentum(u,ef,bf,old[i].q,old[i].m,dt,std::sqrt(1+old[i].u*old[i].u));auto& p=candidate[i];p.u=u[2];p.z=old[i].z+dt*SI::c*(old[i].u+p.u)/(std::sqrt(1+old[i].u*old[i].u)+std::sqrt(1+p.u*p.u));const auto actual=weights(old[i].z,p.z,cells,degree,dz);for(unsigned k=0;k<cells;++k)J[k]+=p.q*(p.z-old[i].z)/(Area*dz*dt)*actual[k];
   }
   J=compensate(J,filter,degree);for(unsigned k=0;k<cells;++k)E[k]=E0[k]-dt/SI::epsilon0*J[k];for(unsigned k=0;k<cells;++k)meanE[k]=.5*(E0[k]+E[k]);const auto next=compensate(meanE,filter,degree);Double force=0.;LD defect=0.,scale=0.;
   for(unsigned i=0;i<old.size();++i){const Double f=dot(weights(old[i].z,candidate[i].z,cells,degree,dz),next),delta=f-used[i];force=std::max(force,std::abs(old[i].q/old[i].m)*dt/SI::c*std::abs(delta));defect+=std::abs((LD)old[i].q*delta*(candidate[i].z-old[i].z));scale+=std::abs((LD)old[i].q*f*(candidate[i].z-old[i].z));guess[i]=candidate[i].z;}
   iterations=trial+1;if(trial>=1&&force<1e-13&&defect<1e-10L*std::max(scale,1e-40L)){converged=true;break;}
  }
  if(!converged)throw std::runtime_error("Sheet particle/field iteration failed");maxIterations=std::max(maxIterations,iterations);LD particleWork=0.,activity=0.;
  for(unsigned i=0;i<old.size();++i){Vec a(0.),b(0.);a[2]=old[i].u;b[2]=candidate[i].u;dK+=(LD)old[i].m*SI::c*SI::c*deltaGamma(a,b);dP+=(LD)old[i].m*SI::c*SI::c*((LD)candidate[i].u-old[i].u);particleWork+=(LD)old[i].q*used[i]*(candidate[i].z-old[i].z);activity+=(LD)old[i].m*SI::c*SI::c*std::abs((LD)candidate[i].u-u0);while(candidate[i].z<0.)candidate[i].z+=Length;while(candidate[i].z>=Length)candidate[i].z-=Length;}
  particles=candidate;rho=charge(particles,cells,degree,dz,filter);LD gauss=0.,continuity=0.,Q=0.,meshWork=0.;for(const auto& p:particles)Q+=std::abs((LD)p.q);
  for(unsigned k=0;k<cells;++k){const unsigned prev=(k+cells-1)%cells;gauss+=std::abs((LD)SI::epsilon0*(E[k]-E[prev])/dz-rho[k])*volume/Q;continuity+=std::abs((LD)rho[k]-rho0[k]+(LD)dt*(J[k]-J[prev])/dz)*volume/Q;meshWork+=(LD)dt*volume*J[k]*.5L*(E0[k]+E[k]);}
  if(diagnostics){Array meanE(cells),nodeE(cells);for(unsigned k=0;k<cells;++k){meanE[k]=.5*(E0[k]+E[k]);nodeE[k]=.25*(E0[k]+E[k]+E0[(k+cells-1)%cells]+E[(k+cells-1)%cells]);}const auto filteredE=compensate(meanE,filter,degree),filteredNode=compensate(nodeE,filter,degree);LD actual=0.,predicted=0.,pathCharge=0.,endpoint=0.;for(unsigned i=0;i<old.size();++i){Double end=particles[i].z;while(end-old[i].z>Length/2)end-=Length;while(end-old[i].z<-Length/2)end+=Length;actual+=(LD)old[i].m*SI::c*SI::c*((LD)particles[i].u-old[i].u);predicted+=(LD)old[i].q*dt*SI::c*dot(weights(old[i].z,end,cells,degree,dz),filteredE);pathCharge+=(LD)old[i].q*dt*SI::c*dot(chargePathWeights(old[i].z,end,cells,degree,dz),filteredNode);}for(unsigned k=0;k<cells;++k)endpoint+=.5L*((LD)rho0[k]+rho[k])*nodeE[k]*volume*dt*SI::c;cumulativeSpace+=predicted-pathCharge;cumulativeTime+=pathCharge-endpoint;cumulativeEndpoint+=endpoint;cumulativeIteration+=actual-predicted;}
  const LD R=dK+fieldEnergy(E,dz)-initial;worstEnergy=std::max(worstEnergy,std::abs(R)/initial);worstGauss=std::max(worstGauss,gauss);worstContinuity=std::max(worstContinuity,continuity);worstWork=std::max(worstWork,std::abs(particleWork-meshWork)/initial);maxP=std::max(maxP,std::abs(dP));maximumActivity=std::max(maximumActivity,activity);
  rows<<step<<','<<dK<<','<<dP<<','<<fieldEnergy(E,dz)-initial<<','<<R<<','<<dP/(192e5L*SI::electronMass*SI::c*SI::c)<<','<<G*beta*dP/(192e5L*SI::elementaryCharge)<<','<<gauss<<','<<continuity<<','<<particleWork-meshWork<<','<<activity<<','<<iterations<<','<<fundamental(E).first<<','<<fundamental(E).second<<','<<cumulativeSpace<<','<<cumulativeTime<<','<<cumulativeEndpoint<<','<<cumulativeIteration<<'\n';
 }
 rows.close();std::ofstream states(std::string(argv[1])+".particles.csv");states<<std::setprecision(21)<<"id,charge_C,mass_kg,initial_z_m,initial_uz,final_z_m,final_uz\n";for(unsigned i=0;i<particles.size();++i)states<<i<<','<<particles[i].q<<','<<particles[i].m<<','<<initialParticles[i].z<<','<<initialParticles[i].u<<','<<particles[i].z<<','<<particles[i].u<<'\n';states.close();std::cout<<std::setprecision(21)<<"{\"complete\":true,\"degree\":"<<degree<<",\"compensation_order\":"<<filter<<",\"cells\":"<<cells<<",\"sheet_pairs\":"<<pairs<<",\"initial_uz\":"<<u0<<",\"displacement_amplitude_m\":"<<amplitude*(Length/24)<<",\"linear_plasma_omega_rad_s\":"<<omega<<",\"continuum_linear_initial_E_amplitude_V_per_m\":"<<linearAmplitude<<",\"initial_E_fundamental_real\":"<<initialMode.first<<",\"initial_E_fundamental_imag\":"<<initialMode.second<<",\"initial_field_energy_J\":"<<initial<<",\"max_energy_residual_over_initial_field\":"<<worstEnergy<<",\"max_work_gap_over_initial_field\":"<<worstWork<<",\"max_gauss_L1_relative\":"<<worstGauss<<",\"max_continuity_L1_relative\":"<<worstContinuity<<",\"max_cPz_residual_J\":"<<maxP<<",\"max_individual_cPz_change_activity_J\":"<<maximumActivity<<",\"max_Pz_residual_over_individual_change_activity\":"<<maxP/maximumActivity<<",\"max_iterations\":"<<maxIterations<<",\"production_modified\":false,\"physical_acceptance\":false}\n";
 if(worstEnergy>1e-7L||worstWork>1e-8L||worstGauss>1e-10L||worstContinuity>1e-10L)throw std::runtime_error("Independent sheet conservation failed");MPI_Finalize();return 0;
}
