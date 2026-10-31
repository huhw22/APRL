#include "source_coupling.h"
#include "eb_cpml.h"
#include "eb_mpi.h"
#include "compute_threads.h"
#include "discrete_field_ledger.h"
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sys/resource.h>
#include <omp.h>
using namespace source_trial;
using Clock=std::chrono::steady_clock;
Double now(){return std::chrono::duration<Double>(Clock::now().time_since_epoch()).count();}
struct Boundary:EBBoundaryOperator {
 EBConvolutionalPML cpml;DiscreteFieldLedger* ledger;unsigned magnetic=0,electric=0;
 Boundary(const EBGridGeometry& g,const EBCPMLParameters& p,DiscreteFieldLedger* l):cpml(g,g.nz,0,p,true,true),ledger(l){}
 void afterMagneticUpdate(EBFieldGrid& f,const EBElectricHaloView& h) override {
  ++magnetic;if(ledger)DiscreteFieldLedger::observe(ledger,"magnetic_curl",f,&h);
  cpml.afterMagneticUpdate(f,h);if(ledger)DiscreteFieldLedger::observe(ledger,"magnetic_cpml",f,nullptr);
 }
 void afterElectricUpdate(EBFieldGrid& f) override {
  ++electric;if(ledger)DiscreteFieldLedger::observe(ledger,"electric_curl_local",f,nullptr);
  cpml.afterElectricUpdate(f);if(ledger){DiscreteFieldLedger::observe(ledger,"electric_cpml",f,nullptr);DiscreteFieldLedger::observe(ledger,"inlet_overwrite",f,nullptr);}
 }
 std::size_t memoryBytes()const override{return cpml.memoryBytes();}
};
long double meshWork(const EBFieldGrid& f,const Electric& old,bool absolute=false) {
 const YeeComponent* es[3]={&f.ex(),&f.ey(),&f.ez()};const YeeComponent* js[3]={&f.jx(),&f.jy(),&f.jz()};long double sum=0;
 for(unsigned a=0;a<3;++a) {
  const auto count=es[a]->size();const auto* e=es[a]->data();const auto* j=js[a]->data();const auto* o=old[a].data();
#pragma omp parallel for reduction(+:sum) schedule(static)
  for(std::ptrdiff_t i=0;i<(std::ptrdiff_t)count;++i){const long double v=(long double)j[i]*(o[i]+(long double)e[i])/2;sum+=absolute?std::abs(v):v;}
 }
 const auto& g=f.geometry();return sum*(long double)g.dt*g.dx*g.dy*g.dz;
}
// Undo the preceding trial's source, then install the newly deposited source.
// Zero-current nodes are never rewritten. Curl and CPML are not reevaluated.
long double currentSource(EBFieldGrid& f,Double sign,const Electric* old,long double* absolute) {
 YeeComponent* es[3]={&f.ex(),&f.ey(),&f.ez()};const YeeComponent* js[3]={&f.jx(),&f.jy(),&f.jz()};const auto& g=f.geometry();long double sum=0,absSum=0;
 for(unsigned a=0;a<3;++a) {
  auto* e=es[a]->data();const auto* j=js[a]->data();const Double* o=old?(*old)[a].data():nullptr;const auto count=es[a]->size();
#pragma omp parallel for reduction(+:sum,absSum) schedule(static)
  for(std::ptrdiff_t i=0;i<(std::ptrdiff_t)count;++i)if(j[i]!=0.) {
   e[i]+=sign*(g.dt/SI::epsilon0)*j[i];
   if(o){const long double work=(long double)j[i]*(o[i]+(long double)e[i])/2;sum+=work;absSum+=std::abs(work);}
  }
 }
 const long double scale=(long double)g.dt*g.dx*g.dy*g.dz;if(absolute)*absolute=absSum*scale;return sum*scale;
}
Vec continuumE(const Vec& point,Double time,Double gammaD,Double betaD,Double physicalN,bool uniform,Double radius) {
 Vec r(point);r[2]=gammaD*(point[2]-betaD*SI::c*time);const Double length=std::sqrt(r.norm2());if(length==0.)return Vec(0.);
 Double enclosed=1.;
 if(uniform)enclosed=std::min(1.,std::pow(length/radius,3));
 else {const Double s=length/(std::sqrt(2.)*radius);enclosed=std::erf(s)-2*s/std::sqrt(std::acos(-1.))*std::exp(-s*s);}
 r*=(-physicalN*SI::elementaryCharge)*enclosed/(4*std::acos(-1.)*SI::epsilon0*length*length*length);
 r[0]*=gammaD;r[1]*=gammaD;return r;
}
void initializeFields(EBFieldGrid& f,const Vec& origin,Double physicalN,bool uniform,Double radius,Double gammaD,Double betaD,bool vacuum) {
 const auto& g=f.geometry();YeeComponent* fields[6]={&f.ex(),&f.ey(),&f.ez(),&f.bx(),&f.by(),&f.bz()};
 const Double offsets[6][3]={{.5,0,0},{0,.5,0},{0,0,.5},{0,.5,.5},{.5,0,.5},{.5,.5,0}};
 for(unsigned a=0;a<6;++a) {
  auto& v=*fields[a];
#pragma omp parallel for collapse(2) schedule(static)
  for(std::size_t k=0;k<v.nz();++k)for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i) {
   Vec x(0.);x[0]=origin[0]+(i+offsets[a][0])*g.dx;x[1]=origin[1]+(j+offsets[a][1])*g.dy;x[2]=origin[2]+(k+offsets[a][2])*g.dz;
   if(vacuum) {
    const Double z=x[2]-(a<3?0.:-g.dt*SI::c/2);const Double pulse=1e5*std::exp(-std::pow(z/(8*g.dz),2)-std::pow(x[0]/(20*g.dx),2)-std::pow(x[1]/(20*g.dy),2));
    v(i,j,k)=a==1?pulse:(a==3?-pulse/SI::c:0.);continue;
   }
   Vec e=continuumE(x,a<3?0.:-g.dt/2,gammaD,betaD,physicalN,uniform,radius);
   v(i,j,k)=a<3?e[a]:(a==3?-betaD*e[1]/SI::c:(a==4?betaD*e[0]/SI::c:0.));
  }
 }
}
std::vector<RelativisticParticleSI> generate(unsigned count,Double physicalN,bool uniform,Double radius,const Vec& proper) {
 if(count%8)throw std::runtime_error("Quiet-start count must be a multiple of eight");
 std::vector<RelativisticParticleSI> p;std::mt19937 random(20261031);std::normal_distribution<Double> normal(0.,1.);
 const Double gammaD=BoostFrameTransform::gammaFromProperVelocity(proper);
 for(unsigned n=0;n<count/8;++n) {
  Vec x(0.);
  if(uniform) {
   const Double r=radius*std::cbrt((n+.5)/(count/8.));const Double z=std::fmod((n+.5)*.6180339887498949,1.);const Double phi=.5*std::acos(-1.)*std::fmod((n+.5)*.4142135623730951,1.);
   x[0]=r*std::sqrt(1-z*z)*std::cos(phi);x[1]=r*std::sqrt(1-z*z)*std::sin(phi);x[2]=r*z;
  }else do{for(unsigned a=0;a<3;++a)x[a]=radius*normal(random);}while(x.norm2()>16*radius*radius);
  for(unsigned signs=0;signs<8;++signs) {
   RelativisticParticleSI q;for(unsigned a=0;a<3;++a)q.position[a]=(signs&(1<<a)?1.:-1.)*std::abs(x[a]);q.position[2]/=gammaD;
   q.properVelocity=proper;q.weight=physicalN/count;q.mass=q.weight*SI::electronMass;q.charge=-q.weight*SI::elementaryCharge;q.id=q.sourceId=p.size()+1;q.frozen=false;p.push_back(q);
  }
 }
 return p;
}
std::uint64_t stateHash(const EBFieldGrid& f) {
 std::uint64_t h=1469598103934665603ULL;const YeeComponent* v[6]={&f.ex(),&f.ey(),&f.ez(),&f.bx(),&f.by(),&f.bz()};
 for(auto* c:v){const auto* p=reinterpret_cast<const unsigned char*>(c->data());for(std::size_t i=0;i<c->bytes();++i){h^=p[i];h*=1099511628211ULL;}}return h;
}
int main(int argc,char** argv) {
 int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
 try {
  if(argc!=10)throw std::runtime_error("Need output, kind, mode, particles, steps, audit, threads, dt divisor, transverse cells");
  const std::string out=argv[1],kind=argv[2],mode=argv[3];unsigned count=std::stoul(argv[4]),steps=std::stoul(argv[5]);const bool audit=std::stoi(argv[6]);const int threads=std::stoi(argv[7]);const Double divide=std::stod(argv[8]);const unsigned nx=std::stoul(argv[9]);
  if(mode!="baseline"&&mode!="coupled")throw std::runtime_error("Unknown coupling mode");
  if(kind!="sphere"&&kind!="drift"&&kind!="undulator"&&kind!="vacuum")throw std::runtime_error("Unknown physical case");
  configureComputeThreads(provided,threads);omp_set_num_threads(threads);const bool coupled=mode=="coupled",uniform=kind=="sphere",vacuum=kind=="vacuum";
  if(vacuum)count=0;
  const Double dz=1.2543262415359034e-6,dt=dz/SI::c/divide,physicalN=vacuum?0.:(uniform?1e8:21300.),radius=uniform?60e-6:30e-6;
  EBGridGeometry g(nx,nx,256,5e-6,5e-6,dz,dt,EBMaxwellSolver::CowanZ);EBFieldGrid f(g);Vec origin(0.);origin[0]=-g.nx*g.dx/2;origin[1]=-g.ny*g.dy/2;origin[2]=-g.nz*g.dz/2;
  BoostFrameTransform frame;frame.setOriginsFromGamma(929.1305492858543,SI::c,0.,.0117,0.);
  Vec proper(0.);if(!uniform&&!vacuum){Vec lab(0.);lab[2]=std::sqrt(1174.171072271572*1174.171072271572-1);frame.properVelocityLabToBox(lab,proper);}
  const Double gammaD=BoostFrameTransform::gammaFromProperVelocity(proper),betaD=proper[2]/gammaD;
  auto particles=generate(count,physicalN,uniform,radius,proper);const auto initial=particles;auto trial=particles;
  SIFieldSourceSet sources;if(kind=="undulator"){auto magnet=SIMagneticElement::planarUndulatorFromK(1.092714657005,.0234,-.234,170,std::acos(-1.)/2);magnet.gaussianFringe=true;magnet.fringeRelativeCutoff=1e-9;magnet.prepare();sources.addPrescribedLabMagnet(magnet);}
  std::unique_ptr<DiscreteFieldLedger> ledger;if(audit)ledger.reset(new DiscreteFieldLedger(out,g,g.nz,0,8,8,MPI_COMM_WORLD));
  EBCPMLParameters cp;cp.cells[0]=cp.cells[1]=8;cp.cells[2]=0;auto* boundary=new Boundary(g,cp,ledger.get());f.setBoundary(std::unique_ptr<EBBoundaryOperator>(boundary));
  EBZSlabHaloExchange halo(MPI_COMM_WORLD);halo.installPhysicalBoundaryMask(f);initializeFields(f,origin,vacuum?0:physicalN,uniform,radius,gammaD,betaD,vacuum);
  ChargeConservingCurrentDepositor deposit(f,origin);Electric old;
  std::ofstream csv(out+"/history.csv");csv<<std::setprecision(21)<<"step,time_s,delta_particle_box_J,particle_grid_work_J,prescribed_work_J,mesh_mid_J,coupling_gap_J,pusher_residual_J,field_operator_residual_J,rms_radius_m,mean_uz,mean_lab_loss_eV,std_lab_loss_eV,outer_iterations,state_residual,work_residual,physical_step_s\n";
  long double cumulativeGrid=0,cumulativePrescribed=0,cumulativeMesh=0,cumulativeField=0,maximumContinuity=0;Double physicalSeconds=0;unsigned maxIterations=0,totalIterations=0,maxRoots=0;long double finalDelta=0,finalGap=0,finalPusher=0;Double finalMeanLoss=0,finalStdLoss=0;
  const auto depositAll=[&](){f.clearCurrent();for(unsigned i=0;i<count;++i){const auto& x=trial[i].position;const auto& p=particles[i].position;for(unsigned a=0;a<3;++a){const Double s=a==2?g.dz:g.dx,n=a==2?g.nz:g.nx;const Double lo=a==2?2.:10.,hi=n-lo;if((x[a]-origin[a])/s<=lo||(x[a]-origin[a])/s>=hi)throw std::runtime_error("Particle entered excluded boundary collar");}deposit.depositSegment(p,x,trial[i].charge);}};
  for(unsigned step=0;step<steps;++step) {
   const Double started=now();if(ledger)ledger->begin(f,step);if(coupled||audit)old=copyElectric(f);
   f.clearCurrent();halo.advanceMagnetic(f);long double grid=0,prescribed=0,absGrid=0,mid=0,absMid=0;unsigned iterations=1;Double stateResidual=0,workResidual=0;
   if(coupled&&count) {
    halo.advanceElectric(f); // The only electric curl/CPML update this step.
    for(unsigned i=0;i<count;++i){trial[i]=particles[i];trial[i].position.pmv(g.dt*SI::c/BoostFrameTransform::gammaFromProperVelocity(particles[i].properVelocity),particles[i].properVelocity);}
    bool converged=false;
    for(iterations=1;iterations<=8;++iterations) {
     grid=prescribed=absGrid=0;stateResidual=0;unsigned roots=0;
#pragma omp parallel for reduction(+:grid,prescribed,absGrid) reduction(max:stateResidual,roots) schedule(static)
     for(std::ptrdiff_t index=0;index<(std::ptrdiff_t)count;++index) {
      const unsigned i=index;const auto previous=trial[i];auto q=particles[i];Vec middle(particles[i].position);middle+=previous.position;middle*=.5;
      Vec unused(0.),b(0.);f.sampleFieldsPosition(middle,origin,unused,b);
      const Vec e=orbitElectric(f,&old,origin,particles[i].position,previous.position);
      const auto work=push(q,e,b,sources,frame,step*g.dt,g.dt,5);grid+=work.grid;prescribed+=work.prescribed;absGrid+=work.absoluteGrid;roots=std::max(roots,work.rootIterations);
      for(unsigned a=0;a<3;++a){stateResidual=std::max(stateResidual,std::abs(q.properVelocity[a]-previous.properVelocity[a])/(1+std::abs(q.properVelocity[a])));stateResidual=std::max(stateResidual,std::abs(q.position[a]-previous.position[a])/g.dz);}
      trial[i]=q;
     }
     maxRoots=std::max(maxRoots,roots);if(iterations>1)currentSource(f,+1,nullptr,nullptr);
     depositAll();mid=currentSource(f,-1,&old,&absMid);
     workResidual=(Double)(std::abs(grid-mid)/std::max(absGrid+absMid,1e-100L));
     if(stateResidual<1e-9&&workResidual<1e-8){converged=true;break;}
    }
    if(!converged)throw std::runtime_error("Source/particle fixed point failed declared tolerances");
   }else {
    if(!vacuum) {
#pragma omp parallel for reduction(+:grid,prescribed) schedule(static)
     for(std::ptrdiff_t i=0;i<(std::ptrdiff_t)count;++i) {
      trial[i]=particles[i];BorisWork work;
      RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(trial[i],f,origin,sources,frame,step*g.dt,g.dt,5,nullptr,true,audit?&work:nullptr,nullptr);
      if(audit){grid+=work.gridWork;prescribed+=work.prescribedWork;}
     }
     depositAll();
    }
    halo.advanceElectric(f);if(audit)mid=meshWork(f,old);
   }
   if(ledger){ledger->finishElectric(f);ledger->afterRemap(f);const auto& r=ledger->lastRecord();cumulativeField+=r.residual;if(std::abs(mid-r.jMid)>1e-11L*std::max(std::abs(mid),1e-30L))throw std::runtime_error("Independent mesh-work disagrees with original ledger");}
   physicalSeconds+=now()-started;maxIterations=std::max(maxIterations,iterations);totalIterations+=iterations;
   if(audit&&count&&(step==0||step+1==steps||step==steps/2)) {
    YeeComponent rho0(g.nx+1,g.ny+1,g.nz+1),rho1(g.nx+1,g.ny+1,g.nz+1);rho0.fill(0.);rho1.fill(0.);
    for(unsigned i=0;i<count;++i){deposit.depositCharge(particles[i].position,particles[i].charge,rho0);deposit.depositCharge(trial[i].position,trial[i].charge,rho1);}
    Double scale=0;for(std::size_t i=0;i<rho0.size();++i)scale=std::max(scale,std::max(std::abs(rho0.data()[i]),std::abs(rho1.data()[i]))/g.dt);
    maximumContinuity=std::max(maximumContinuity,(long double)ChargeConservingCurrentDepositor::maxContinuityResidual(f,rho0,rho1)/std::max(scale,1e-100));
   }
   particles.swap(trial);cumulativeGrid+=grid;cumulativePrescribed+=prescribed;cumulativeMesh+=mid;
   long double delta=0,r2=0,uz=0,labLoss=0,labLoss2=0;
   for(unsigned i=0;i<count;++i) {
    const auto& q=particles[i];const auto dg=deltaGamma(initial[i].properVelocity,q.properVelocity);delta+=(long double)q.mass*SI::c*SI::c*dg;r2+=q.position.norm2();uz+=q.properVelocity[2];
    const long double loss=-(long double)SI::electronMass*SI::c*SI::c/SI::elementaryCharge*frame.gamma()*(dg+(long double)frame.beta()*((long double)q.properVelocity[2]-initial[i].properVelocity[2]));labLoss+=loss;labLoss2+=loss*loss;
   }
   finalMeanLoss=count?(Double)(labLoss/count):0.;finalStdLoss=count?(Double)std::sqrt(std::max(0.L,labLoss2/count-labLoss*labLoss/count/count)):0.;
   finalDelta=delta;finalGap=cumulativeGrid-cumulativeMesh;finalPusher=delta-cumulativeGrid-cumulativePrescribed;
   csv<<step+1<<','<<(step+1)*g.dt<<','<<delta<<','<<cumulativeGrid<<','<<cumulativePrescribed<<','<<cumulativeMesh<<','<<finalGap<<','<<finalPusher<<','<<cumulativeField<<','<<(count?std::sqrt(r2/count):0)<<','<<(count?uz/count:0)<<','<<finalMeanLoss<<','<<finalStdLoss<<','<<iterations<<','<<stateResidual<<','<<workResidual<<','<<now()-started<<'\n';
  }
  csv.close();if(ledger)ledger->close(true);if(boundary->magnetic!=steps||boundary->electric!=steps)throw std::runtime_error("Repeated wave/CPML update");
  if(audit&&count){std::ofstream file(out+"/particles-final.csv");file<<std::setprecision(17)<<"id,x,y,z,ux,uy,uz\n";for(const auto& q:particles)file<<q.id<<','<<q.position[0]<<','<<q.position[1]<<','<<q.position[2]<<','<<q.properVelocity[0]<<','<<q.properVelocity[1]<<','<<q.properVelocity[2]<<'\n';}
  struct rusage usage;getrusage(RUSAGE_SELF,&usage);std::ofstream result(out+"/RESULT.json");result<<std::setprecision(21)<<"{\"complete\":true,\"kind\":\""<<kind<<"\",\"mode\":\""<<mode<<"\",\"particles\":"<<count<<",\"represented_electrons\":"<<physicalN<<",\"radius_or_sigma_m\":"<<radius<<",\"dt_s\":"<<dt<<",\"steps\":"<<steps<<",\"threads\":"<<threads<<",\"audit\":"<<(audit?"true":"false")<<",\"physical_seconds\":"<<physicalSeconds<<",\"seconds_per_step\":"<<physicalSeconds/steps<<",\"peak_RSS_KiB\":"<<usage.ru_maxrss<<",\"maximum_outer_iterations\":"<<maxIterations<<",\"mean_outer_iterations\":"<<(Double)totalIterations/steps<<",\"maximum_local_scalar_iterations\":"<<maxRoots<<",\"relative_continuity_residual_max\":"<<maximumContinuity<<",\"final_grid_coupling_gap_J\":"<<finalGap<<",\"final_pusher_residual_J\":"<<finalPusher<<",\"final_particle_delta_box_J\":"<<finalDelta<<",\"field_operator_residual_J\":"<<cumulativeField<<",\"mean_lab_loss_eV\":"<<finalMeanLoss<<",\"std_lab_loss_eV\":"<<finalStdLoss<<",\"field_hash\":\""<<stateHash(f)<<"\",\"magnetic_CPML_calls\":"<<boundary->magnetic<<",\"electric_CPML_calls\":"<<boundary->electric<<",\"production_method_changed\":false}\n";result.close();
  std::cout<<"COMPLETE "<<kind<<' '<<mode<<" seconds="<<physicalSeconds<<" iterations="<<maxIterations<<" gap="<<(double)finalGap<<'\n';MPI_Finalize();return 0;
 }catch(const std::exception& e){std::cerr<<"PROTOTYPE FAILED: "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,2);return 2;}
}
