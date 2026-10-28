/* Verification-only held inlet / causal right cut with production grid,
 * boost, transverse CPML and original particle/Maxwell kernels.
 * Independent held donor and frozen prewarmed reference are explicitly timed.
 */
#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <mpi.h>
#include "boostframe.h"
#include "config.h"
#include "compute_threads.h"
#include "eb_deposition.h"
#include "eb_initial_field.h"
#include "eb_particles.h"
#include "particle_initializer.h"
#include "held_cut_boundary.h"
#include "moving_window_remap.h"
#include "eb_cpml.h"
namespace {
using aprl::Double;using aprl::EBFieldGrid;using aprl::EBGridGeometry;
using aprl::FieldVector;using aprl::RelativisticParticleSI;using aprl::YeeComponent;
using aprl::verification::MovingWindowRemap;
using aprl::verification::HeldCutBoundary;
using aprl::verification::HeldCutReference;
using aprl::verification::copyHeldCutFields;
using aprl::verification::copyHeldCutHistories;
using Clock=std::chrono::steady_clock;
const std::size_t NX=160,NY=160,SIDE=10,LOWER_GUARD=64,UPPER_GUARD=128;
const double DX=5e-6,DY=5e-6,DZ=1.2543262415359034e-6;
const double BOOST_GAMMA=929.1305492858543,LAB_GAMMA=1174.185661478854;
const double LAB_GAMMA_SIGMA=2.568638506144439;
const double LAB_U_MEAN[3]={3.526198929811919e-6,2.294658317104036e-5,1174.1852350615736};
const double LAB_U_SIGMA[3]={.02709500523622301,.02543628834980338,2.568639437850893};
const double TRANSVERSE_CORRELATION[2]={-.7334153049192772,-.6546423205475509};
const double SIGMA_Z_CELLS=8.,SOURCE_SIGMA_Z_CELLS=10540.26945;
const double TOTAL_CHARGE=-6.686327698e-12*SIGMA_Z_CELLS/SOURCE_SIGMA_Z_CELLS;
const double RELEASE_LAB=-.05,TRUE_DETECTOR_LAB=4.5;
const unsigned SEED=710271;
const std::size_t LOCAL_PLANE_CELLS=16,WIDTH=LOWER_GUARD+LOCAL_PLANE_CELLS+UPPER_GUARD;
void require(bool b,const std::string& s){if(!b)throw std::runtime_error(s);}
double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
std::ofstream fresh(const std::string& p){require(!std::ifstream(p).good(),"Refusing overwrite "+p);std::ofstream f(p);require(f.good(),"Cannot create "+p);f<<std::setprecision(17);return f;}
void directory(const std::string& p){if(mkdir(p.c_str(),0755)!=0)require(errno==EEXIST,"Cannot create directory "+p);}
YeeComponent& component(EBFieldGrid& f,unsigned a){if(a==0)return f.ex();if(a==1)return f.ey();if(a==2)return f.ez();if(a==3)return f.bx();if(a==4)return f.by();return f.bz();}
const YeeComponent& component(const EBFieldGrid& f,unsigned a){return component(const_cast<EBFieldGrid&>(f),a);}
EBGridGeometry geometry(std::size_t nz){return EBGridGeometry(NX,NY,nz,DX,DY,DZ,DZ/aprl::SI::c,aprl::EBMaxwellSolver::CowanZ);}
FieldVector<double> origin(std::size_t offset){FieldVector<double> o(0.);o[0]=-.5*NX*DX;o[1]=-.5*NY*DY;o[2]=offset*DZ;return o;}
long double gamma(const FieldVector<double>& u){return std::sqrt(1.L+static_cast<long double>(u[0])*u[0]+static_cast<long double>(u[1])*u[1]+static_cast<long double>(u[2])*u[2]);}
struct Segment {FieldVector<double> start,end;double first,last;FieldVector<double> startU,endU;Segment()=default;Segment(const FieldVector<double>& a,const FieldVector<double>& b,double f,double l,const FieldVector<double>& ua,const FieldVector<double>& ub):start(a),end(b),first(f),last(l),startU(ua),endU(ub){}};
bool inside(const FieldVector<double>& p,const FieldVector<double>& lo,
            const EBGridGeometry& g) {
  return p[0]>=lo[0]&&p[0]<=lo[0]+g.nx*g.dx&&
    p[1]>=lo[1]&&p[1]<=lo[1]+g.ny*g.dy&&
    p[2]>=lo[2]&&p[2]<=lo[2]+g.nz*g.dz;
}
bool clip(const Segment& segment,const FieldVector<double>& lo,
          const EBGridGeometry& g,FieldVector<double>& start,FieldVector<double>& end) {
  const double extent[3]={g.nx*g.dx,g.ny*g.dy,g.nz*g.dz};
  double first=0,last=1;
  for(unsigned a=0;a<3;++a) {
    const double d=segment.end[a]-segment.start[a];
    if(d==0.) {if(segment.start[a]<lo[a]||segment.start[a]>lo[a]+extent[a])return false;}
    else {
      double low=(lo[a]-segment.start[a])/d,high=(lo[a]+extent[a]-segment.start[a])/d;
      if(low>high)std::swap(low,high);
      first=std::max(first,low);last=std::min(last,high);
      if(last<=first)return false;
    }
  }
  for(unsigned a=0;a<3;++a) {
    start[a]=std::max(lo[a],std::min(lo[a]+extent[a],segment.start[a]+first*(segment.end[a]-segment.start[a])));
    end[a]=std::max(lo[a],std::min(lo[a]+extent[a],segment.start[a]+last*(segment.end[a]-segment.start[a])));
  }
  return true;
}

void cropCurrent(const EBFieldGrid& ext,EBFieldGrid& out,std::size_t start) {
  const YeeComponent* from[]={&ext.jx(),&ext.jy(),&ext.jz()};
  YeeComponent* to[]={&out.jx(),&out.jy(),&out.jz()};
  for(unsigned a=0;a<3;++a) {
    const std::size_t plane=to[a]->nx()*to[a]->ny();
    require(start+to[a]->nz()<=from[a]->nz(),"Current crop left deposition collar");
    std::copy(from[a]->data()+start*plane,
      from[a]->data()+(start+to[a]->nz())*plane,to[a]->data());
  }
}

long double interiorCharge(const YeeComponent& rho) {
  long double result=0;
  for(std::size_t k=1;k+1<rho.nz();++k)
    for(std::size_t j=1;j+1<rho.ny();++j)
      for(std::size_t i=1;i+1<rho.nx();++i)result+=rho(i,j,k);
  return result*DX*DY*DZ;
}
long double outwardCurrent(const EBFieldGrid& f) {
  const auto& g=f.geometry();long double result=0;
  for(std::size_t k=1;k<g.nz;++k)for(std::size_t j=1;j<g.ny;++j)
    result+=(f.jx()(g.nx-1,j,k)-f.jx()(0,j,k))*DY*DZ;
  for(std::size_t k=1;k<g.nz;++k)for(std::size_t i=1;i<g.nx;++i)
    result+=(f.jy()(i,g.ny-1,k)-f.jy()(i,0,k))*DX*DZ;
  for(std::size_t j=1;j<g.ny;++j)for(std::size_t i=1;i<g.nx;++i)
    result+=(f.jz()(i,j,g.nz-1)-f.jz()(i,j,0))*DX*DY;
  return result;
}
double continuity(const EBFieldGrid& f,const YeeComponent& before,
                  const YeeComponent& after) {
  double maximum=0,scale=0;const auto& g=f.geometry();
  for(std::size_t k=1;k<g.nz;++k)for(std::size_t j=1;j<g.ny;++j)
    for(std::size_t i=1;i<g.nx;++i) {
      const double value=(after(i,j,k)-before(i,j,k))/g.dt+
        (f.jx()(i,j,k)-f.jx()(i-1,j,k))/DX+
        (f.jy()(i,j,k)-f.jy()(i,j-1,k))/DY+
        (f.jz()(i,j,k)-f.jz()(i,j,k-1))/DZ;
      maximum=std::max(maximum,std::abs(value));
      scale=std::max(scale,std::max(std::abs(before(i,j,k)),std::abs(after(i,j,k)))/g.dt);
    }
  return scale?maximum/scale:maximum;
}
double sampleComponent(const EBFieldGrid& f,std::size_t offset,unsigned a,
                       double gx,double gy,double gz) {
  static const double stagger[6][3]={{.5,0,0},{0,.5,0},{0,0,.5},
    {0,.5,.5},{.5,0,.5},{.5,.5,0}};
  const double q[3]={gx-stagger[a][0],gy-stagger[a][1],gz-stagger[a][2]};
  std::ptrdiff_t index[3];double frac[3];
  for(unsigned d=0;d<3;++d){index[d]=static_cast<std::ptrdiff_t>(std::floor(q[d]));frac[d]=q[d]-index[d];}
  const auto& c=component(f,a);double value=0;
  for(unsigned k=0;k<2;++k)for(unsigned j=0;j<2;++j)for(unsigned i=0;i<2;++i) {
    const double weight=(i?frac[0]:1-frac[0])*(j?frac[1]:1-frac[1])*(k?frac[2]:1-frac[2]);
    if(weight==0)continue;
    const std::ptrdiff_t x=index[0]+i,y=index[1]+j,z=index[2]+k-offset;
    require(x>=0&&y>=0&&z>=0&&x<static_cast<std::ptrdiff_t>(c.nx())&&
      y<static_cast<std::ptrdiff_t>(c.ny())&&z<static_cast<std::ptrdiff_t>(c.nz()),"Detector sample outside active lattice");
    value+=weight*c(x,y,z);
  }
  require(std::isfinite(value),"Non-finite detector field");return value;
}

void particleHeader(std::ofstream& out) {
  out<<"id,step,fraction,time_box_s,time_lab_s,x_m,y_m,z_box_m,ux,uy,uz,charge_C,delta_energy_box_eV,delta_energy_lab_eV,frozen\n";
}
void particleRow(std::ofstream& out,const RelativisticParticleSI& p,
                 const RelativisticParticleSI& initial,unsigned step,double fraction,
                 const aprl::BoostFrameTransform& frame) {
  const double time=(step+fraction)*DZ/aprl::SI::c;
  long double numerator=0;
  for(unsigned a=0;a<3;++a)numerator+=(static_cast<long double>(p.properVelocity[a])-initial.properVelocity[a])*
    (static_cast<long double>(p.properVelocity[a])+initial.properVelocity[a]);
  const long double dg=numerator/(gamma(p.properVelocity)+gamma(initial.properVelocity));
  const long double rest=static_cast<long double>(aprl::SI::electronMass)*aprl::SI::c*aprl::SI::c/aprl::SI::elementaryCharge;
  const long double lab=frame.gamma()*(dg+static_cast<long double>(frame.beta())*
    (static_cast<long double>(p.properVelocity[2])-initial.properVelocity[2]));
  out<<p.id<<','<<step<<','<<fraction<<','<<time<<','<<frame.labTFromBoxZT(p.position[2],time)
    <<','<<p.position[0]<<','<<p.position[1]<<','<<p.position[2]
    <<','<<p.properVelocity[0]<<','<<p.properVelocity[1]<<','<<p.properVelocity[2]
    <<','<<p.charge<<','<<dg*rest<<','<<lab*rest<<','<<(p.frozen?1:0)<<'\n';
}


struct Timing {double initialize=0,prewarm=0,donor=0,full=0,cut=0,feed=0,diagnostics=0;};
struct Fixture {
  unsigned steps,particles,warm,warmHorizon;std::size_t offset0,nz;
  double release,detector,betaMean=0,betaMin=1,betaMax=-1;
  unsigned substeps=1;
  aprl::BoostFrameTransform frame;
  aprl::EBCPMLParameters cpml;
  aprl::SIFieldSourceSet sources;
  std::vector<RelativisticParticleSI> seed;
  Fixture(unsigned s,unsigned p,unsigned w,unsigned wh=64)
    :steps(s),particles(p),warm(w),warmHorizon(wh) {
    require(s>=16&&s<=128&&p>=16&&p<=1024&&wh>=64&&wh<=256&&w<=wh,
      "Outside bounded fixture");
    // The explicit preparation horizon fixes one common reference geometry
    // for all prewarm durations in a convergence series. The default retains
    // the previous 16-versus-64 fixture exactly. Native Cowan plus transverse-
    // CPML dependency bounds are checked
    // independently before this fixture is used: lower-end contamination
    // reaches at most n+1 native layers; upper-end contamination n layers.
    // The existing production-shaped cut remains 208 cells. For horizon 90,
    // offset=248,nz=624 also reproduce the earlier 64-step/16-warm geometry.
    const unsigned horizon=std::max(90u,steps);
    offset0=2*horizon+warmHorizon+4;
    nz=((offset0+WIDTH+horizon+warmHorizon+4+15)/16)*16;
    release=(offset0+LOWER_GUARD)*DZ;detector=release+LOCAL_PLANE_CELLS*DZ;
    frame.setOriginsFromGamma(BOOST_GAMMA,aprl::SI::c,0.,RELEASE_LAB,release);
    cpml.cells[0]=cpml.cells[1]=SIDE;cpml.cells[2]=0;
    cpml.polynomialOrder=3.;cpml.targetReflection=1e-8;cpml.kappaMax=8.;cpml.alphaFraction=0.;
    auto magnet=aprl::SIMagneticElement::planarUndulatorFromK(1.092714657005,.0234,0.,170,1.5707963267948966);
    magnet.gaussianFringe=true;magnet.fringeRelativeCutoff=1e-9;magnet.prepare();
    sources.addPrescribedLabMagnet(magnet);
    aprl::SimulationConfig config;config.beam.type=aprl::BeamInputType::GeneratedGaussian;
    config.beam.macroparticles=particles;config.beam.electrons=-TOTAL_CHARGE/aprl::SI::elementaryCharge;
    config.beam.randomSeed=SEED;config.beam.gamma=LAB_GAMMA;
    config.beam.direction=FieldVector<double>(0.);config.beam.direction[2]=1.;
    config.beam.center=FieldVector<double>(0.);config.beam.positionOffset=FieldVector<double>(0.);
    config.beam.sigmaPosition=FieldVector<double>(0.);
    config.beam.sigmaPosition[0]=12.78982675*DX;config.beam.sigmaPosition[1]=12.2478603*DY;
    config.beam.sigmaPosition[2]=SIGMA_Z_CELLS*DZ;
    // The generator supplies independent unit Gaussians in transverse u.
    // Correlate them with x/y using measured original-source moments before
    // applying the production Lorentz transform. No sample recentering occurs.
    config.beam.sigmaProperVelocity=FieldVector<double>(1.);config.beam.sigmaProperVelocity[2]=LAB_U_SIGMA[2];
    aprl::ParticleInitializationReport report;seed=aprl::ParticleInitializer::create(config,MPI_COMM_WORLD,report);
    for(auto& q:seed) {
      auto lab=q.properVelocity;
      for(unsigned a=0;a<2;++a) {
        const double correlation=TRANSVERSE_CORRELATION[a];
        const double gaussianPosition=q.position[a]/config.beam.sigmaPosition[a];
        lab[a]=LAB_U_MEAN[a]+LAB_U_SIGMA[a]*(correlation*gaussianPosition+
          std::sqrt(1.-correlation*correlation)*lab[a]);
      }
      lab[2]+=LAB_U_MEAN[2]-aprl::BoostFrameTransform::gammaBetaFromGamma(LAB_GAMMA);
      frame.properVelocityLabToBox(lab,q.properVelocity);
      q.id=q.sourceId;q.frozen=true;q.frozenBetaZ=frame.heldLongitudinalBetaZ(q.properVelocity);
      betaMean+=q.frozenBetaZ;betaMin=std::min(betaMin,q.frozenBetaZ);betaMax=std::max(betaMax,q.frozenBetaZ);
    }
    betaMean/=seed.size();
    const double center=release-(betaMean+frame.beta())*60.*DZ;
    for(auto& q:seed) {
      q.position[2]+=center;
      require(q.position[2]<release,"Initial source already beyond release plane");
      require(q.position[0]>origin(0)[0]+(SIDE+2)*DX&&q.position[0]<origin(0)[0]+(NX-SIDE-2)*DX&&
              q.position[1]>origin(0)[1]+(SIDE+2)*DY&&q.position[1]<origin(0)[1]+(NY-SIDE-2)*DY,
              "Generated particle touches transverse CPML");
    }
    const double maxLabAdvance=frame.gamma()*DZ*(frame.beta()+betaMax);
    const double ratio=80.*maxLabAdvance/.0234;
    substeps=std::max(1.,std::ceil(ratio-128*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(ratio))));
    require(substeps<=4096,"Production substep bound exceeded");
    require(offset0>2*steps+warm+2,"Fixed lower endpoint enters inlet dependency envelope");
    require(nz-(offset0+WIDTH)>steps+warm+1,"Fixed upper endpoint enters retained dependency envelope");
  }
  double moving(double z,double n)const{return z-frame.beta()*DZ*n;}
  std::size_t offset(unsigned n)const {
    const auto offset=static_cast<std::size_t>(std::floor(static_cast<double>(offset0)-frame.beta()*n));
    require(offset==offset0-n,"Production-beta schedule ceased shifting one cell each step");return offset;
  }
  RelativisticParticleSI held(std::size_t i,double n)const {
    auto q=seed[i];q.position[2]+=q.frozenBetaZ*DZ*n;return q;
  }
};

struct Reference {
  EBFieldGrid* fields=nullptr;const MovingWindowRemap* layout=nullptr;
  bool frozen=false;double beta=0.;unsigned state=0;
  static double query(void* raw,unsigned a,std::size_t i,std::size_t j,std::ptrdiff_t k,unsigned n) {
    auto& r=*static_cast<Reference*>(raw);require(a<6,"Invalid reference field");
    if(!r.frozen)require(n==r.state,"Held donor time mismatch");
    const auto& c=component(*r.fields,a);const double z=static_cast<double>(k)-(r.frozen?r.beta*n:0.);
    require(z>=0.&&z<=c.nz()-1,"Reference query outside prepared domain");
    const auto lo=static_cast<std::size_t>(std::floor(z)),hi=std::min(lo+1,c.nz()-1);const double f=z-lo;
    return f==0.?c(i,j,lo):(1-f)*c(i,j,lo)+f*c(i,j,hi);
  }
  void histories(EBFieldGrid& target,const MovingWindowRemap& targetLayout,std::size_t offset,unsigned n) {
    if(!frozen) {copyHeldCutHistories(target,targetLayout,*fields,*layout,offset,0,1);return;}
    auto& to=dynamic_cast<aprl::EBConvolutionalPML&>(target.boundary());
    const auto& from=dynamic_cast<const aprl::EBConvolutionalPML&>(fields->boundary());
    for(unsigned h=0;h<12;++h) {
      const auto& tl=targetLayout.historyLayout(h);const auto& sl=layout->historyLayout(h);
      if(tl.derivativeAxis==2){require(tl.values==0,"Unexpected longitudinal history");continue;}
      const auto q=tl.dimensions[tl.derivativeAxis==0?1:0];
      const double z=offset-beta*n;require(z>=0.&&z<=sl.dimensions[2]-1,"Frozen side history outside cache");
      const auto lo=static_cast<std::size_t>(std::floor(z)),hi=std::min(lo+1,sl.dimensions[2]-1);const double f=z-lo;
      auto& tv=to.checkpointHistory(h);const auto& sv=from.checkpointHistory(h);
      for(std::size_t slot=0;slot<tl.activeAxisIndices.size();++slot)for(std::size_t i=0;i<q;++i) {
        const auto a=(slot*sl.dimensions[2]+lo)*q+i,b=(slot*sl.dimensions[2]+hi)*q+i;
        tv[(slot*tl.dimensions[2])*q+i]=f==0.?sv[a]:(1-f)*sv[a]+f*sv[b];
      }
    }
  }
};

void advanceHeld(EBFieldGrid& grid,const Fixture& f,int n) {
  grid.advanceMagnetic();grid.clearCurrent();
  aprl::ChargeConservingCurrentDepositor d(grid,origin(0));
  for(std::size_t i=0;i<f.seed.size();++i) {
    const auto before=f.held(i,n),after=f.held(i,n+1);
    require(inside(before.position,origin(0),grid.geometry())&&inside(after.position,origin(0),grid.geometry()),"Donor source escaped full domain");
    d.depositSegment(before.position,after.position,before.charge);
  }
  grid.advanceElectric();
}

struct World {
  const Fixture& fixture;bool cut;std::size_t offset;
  std::unique_ptr<EBFieldGrid> fields,deposit;
  std::unique_ptr<MovingWindowRemap> layout;
  HeldCutBoundary* boundary=nullptr;
  std::vector<RelativisticParticleSI> particles;
  std::vector<std::size_t> active,queue;
  std::vector<bool> admitted,released,detected;
  std::size_t cursor=0,initialActive=0,peakActive=0,admissionCount=0,discardCount=0;
  std::uint64_t depositedSegments=0,shiftBytes=0,historyBytes=0;
  double minPicClearance=std::numeric_limits<double>::infinity(),maxContinuity=0,maxChargeBalance=0;
  double maximumLabZ=-1e100,referenceSeconds=0,preparationSeconds=0,advanceSeconds=0;
  std::ofstream initial,releases,crossings,final,plane,history;
  YeeComponent rhoBefore,rhoAfter;
  World(const std::string& root,const Fixture& f,bool isCut,EBFieldGrid& donor,
        const MovingWindowRemap& donorLayout,Reference& ref)
    :fixture(f),cut(isCut),offset(cut?f.offset0:0),particles(f.seed),admitted(f.seed.size(),!cut),
     released(f.seed.size(),false),detected(f.seed.size(),false),
     rhoBefore(NX+1,NY+1,(cut?WIDTH+4:f.nz)+1),rhoAfter(NX+1,NY+1,(cut?WIDTH+4:f.nz)+1) {
    directory(root);initial=fresh(root+"/initial-particles.csv");releases=fresh(root+"/release-particles.csv");
    crossings=fresh(root+"/detector-crossings.csv");final=fresh(root+"/final-particles.csv");
    plane=fresh(root+"/field-plane.csv");history=fresh(root+"/history.csv");
    particleHeader(initial);particleHeader(releases);particleHeader(crossings);particleHeader(final);
    for(const auto& p:particles)particleRow(initial,p,p,0,0.,f.frame);
    plane<<"step,pixel_i,pixel_j,time_E_box_s,time_B_box_s,time_E_lab_s,x_m,y_m,z_box_m,Ex_V_m,Ey_V_m,Ez_V_m,Bx_T,By_T,Bz_T\n";
    history<<"step,origin_z_m,active,held,admitted,released,detected,mean_delta_energy_lab_eV,sigma_energy_lab_eV,rms_x_m,rms_y_m,mean_x_m,mean_y_m,continuity_max,charge_balance_max_C\n";
    fields.reset(new EBFieldGrid(geometry(cut?WIDTH:f.nz)));
    if(cut) {
      HeldCutReference provider;provider.context=&ref;provider.field=&Reference::query;
      boundary=new HeldCutBoundary(fields->geometry(),f.cpml,offset,provider,0.);
      fields->setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(boundary));
      deposit.reset(new EBFieldGrid(geometry(WIDTH+4)));
    } else fields->setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(new aprl::EBConvolutionalPML(fields->geometry(),f.nz,0,f.cpml,true,true)));
    layout.reset(new MovingWindowRemap(fields->geometry(),f.cpml,{-1}));
    copyHeldCutFields(*fields,donor,offset);copyHeldCutHistories(*fields,*layout,donor,donorLayout,offset);
    queue.resize(particles.size());std::iota(queue.begin(),queue.end(),0);
    std::sort(queue.begin(),queue.end(),[&](std::size_t a,std::size_t b){return entryTime(a)<entryTime(b);});
    if(!cut){active=queue;cursor=queue.size();admissionCount=active.size();}
    else admit(0);
    initialActive=active.size();peakActive=active.size();
  }
  double entryTime(std::size_t i)const {
    return ((fixture.offset0-2)*DZ-fixture.seed[i].position[2])/
      ((fixture.seed[i].frozenBetaZ+1.)*DZ);
  }
  void admit(unsigned n) {
    if(!cut)return;
    while(cursor<queue.size()&&entryTime(queue[cursor])<=n+1.) {
      const auto i=queue[cursor++];require(!admitted[i],"Duplicate source admission");
      particles[i]=fixture.held(i,n);admitted[i]=true;active.push_back(i);++admissionCount;
      require(particles[i].position[2]<fixture.moving(fixture.release,n),"Source released before admission");
    }
    peakActive=std::max(peakActive,active.size());
  }
  EBFieldGrid& depositGrid(){return cut?*deposit:*fields;}
  FieldVector<double> depositOrigin()const{return origin(cut?offset-2:0);}
  void rho(YeeComponent& target) {
    target.fill(0.);auto& dgrid=depositGrid();const auto o=depositOrigin();aprl::ChargeConservingCurrentDepositor d(dgrid,o);
    for(auto i:active)if(inside(particles[i].position,o,dgrid.geometry()))d.depositCharge(particles[i].position,particles[i].charge,target);
  }
  void prepare(unsigned n) {
    const auto begin=Clock::now();
    admit(n);if(boundary)boundary->beginStep(*fields,n);
    preparationSeconds+=elapsed(begin);
  }
  void advance(unsigned n) {
    const auto begin=Clock::now();const auto& f=fixture;const double dt=DZ/aprl::SI::c;
    const auto fieldOrigin=origin(offset),depOrigin=depositOrigin();auto& dep=depositGrid();
    const bool audit=n%8==0||n+1==f.steps;
    if(audit)rho(rhoBefore);
    dep.clearCurrent();
    const aprl::EBElectricHaloView halo=boundary?boundary->electricHalos(*fields,n):aprl::EBElectricHaloView();
    fields->advanceMagnetic(halo);
    aprl::ChargeConservingCurrentDepositor depositor(dep,depOrigin);
    for(const auto i:active) {
      auto& p=particles[i];const auto& seed=f.seed[i];const auto before=p.position;
      Segment segments[2];unsigned used=0;double fraction=0.;
      if(p.frozen) {
        auto end=f.held(i,n+1).position;
        const double a=before[2]-f.moving(f.release,n),b=end[2]-f.moving(f.release,n+1);
        if(b>=0.) {
          fraction=a>=0.?0.:std::max(0.,std::min(1.,-a/(b-a)));
          p.position=before;p.position[2]+=fraction*(end[2]-before[2]);
          if(fraction>0.)segments[used++]=Segment(before,p.position,0.,fraction,p.properVelocity,p.properVelocity);
          p.frozen=false;require(!released[i],"Duplicate release");released[i]=true;
          for(unsigned axis=0;axis<3;++axis)require(p.properVelocity[axis]==seed.properVelocity[axis],"Held release changed momentum");
          require(p.position[0]==seed.position[0]&&p.position[1]==seed.position[1],"Held release changed x/y");
          particleRow(releases,p,seed,n,fraction,f.frame);
        }else {p.position=end;segments[used++]=Segment(before,end,0.,1.,p.properVelocity,p.properVelocity);}
      }
      if(!p.frozen&&fraction<1.) {
        const auto start=p.position,startU=p.properVelocity;
        // After B and before remap, the right trusted native range has lost
        // one layer. Keep the full trilinear sample strictly inside that range.
        const double clearance=std::min(start[2]-fieldOrigin[2],fieldOrigin[2]+(fields->geometry().nz-2)*DZ-start[2]);
        minPicClearance=std::min(minPicClearance,clearance);
        if(cut)require(clearance>=2*DZ,"PIC sample touched transient right sacrifice collar");
        aprl::RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(p,*fields,fieldOrigin,f.sources,f.frame,
          (n+fraction)*dt,(1-fraction)*dt,f.substeps);
        segments[used++]=Segment(start,p.position,fraction,1.,startU,p.properVelocity);
      }
      require(used>0&&used<=2,"Invalid source segment split");
      for(unsigned a=0;a<3;++a)require(std::isfinite(p.position[a])&&std::isfinite(p.properVelocity[a]),"Nonfinite particle");
      require(p.position[0]>fieldOrigin[0]+(SIDE+2)*DX&&p.position[0]<fieldOrigin[0]+(NX-SIDE-2)*DX&&
              p.position[1]>fieldOrigin[1]+(SIDE+2)*DY&&p.position[1]<fieldOrigin[1]+(NY-SIDE-2)*DY,"Particle reached transverse CPML");
      maximumLabZ=std::max(maximumLabZ,f.frame.labZFromBoxZT(p.position[2],(n+1)*dt));
      for(unsigned s=0;s<used;++s) {
        const auto& leg=segments[s];
        if(!detected[i]) {
          const double a=leg.start[2]-f.moving(f.detector,n+leg.first),b=leg.end[2]-f.moving(f.detector,n+leg.last);
          if(a<=0.&&b>=0.&&b>a) {
            const double t=-a/(b-a);auto event=p;
            for(unsigned axis=0;axis<3;++axis){event.position[axis]=leg.start[axis]+t*(leg.end[axis]-leg.start[axis]);event.properVelocity[axis]=leg.startU[axis]+t*(leg.endU[axis]-leg.startU[axis]);}
            require(released[i],"Auxiliary plane crossed before release");
            particleRow(crossings,event,seed,n,leg.first+t*(leg.last-leg.first),f.frame);detected[i]=true;
          }
        }
        FieldVector<double> a(0.),b(0.);if(clip(leg,depOrigin,dep.geometry(),a,b))depositedSegments+=depositor.depositSegment(a,b,p.charge).subsegments;
      }
    }
    if(audit) {
      rho(rhoAfter);maxContinuity=std::max(maxContinuity,continuity(dep,rhoBefore,rhoAfter));
      const long double error=interiorCharge(rhoAfter)-interiorCharge(rhoBefore)+dt*outwardCurrent(dep);
      maxChargeBalance=std::max(maxChargeBalance,static_cast<double>(std::abs(error)));
    }
    if(cut)cropCurrent(dep,*fields,2);
    fields->advanceElectric();advanceSeconds+=elapsed(begin);
  }
  void shift(unsigned n,Reference& reference) {
    if(!cut)return;
    const auto begin=Clock::now();const auto next=fixture.offset(n);require(offset==next+1,"Cut must discard exactly one cell");
    struct Incoming {Reference* r;std::size_t offset;unsigned n;} context{&reference,next,n};
    aprl::verification::MovingWindowIncoming incoming;incoming.context=&context;
    incoming.value=[](void* raw,unsigned a,std::size_t i,std::size_t j,std::size_t k,double)->double {
      auto& v=*static_cast<Incoming*>(raw);return a<6?Reference::query(v.r,a,i,j,v.offset+k,v.n):0.;
    };
    const auto report=layout->shift(*fields,-1,origin(next)[2],incoming,aprl::verification::MovingWindowCurrentPolicy::Reset);
    reference.histories(*fields,*layout,next,n);offset=next;
    shiftBytes+=report.fieldMovedBytes;historyBytes+=report.historyMovedBytes;
    // Removed source trajectories must first be fully deposited through the
    // open collar. No particle reaches it in the initial bounded fixture.
    active.erase(std::remove_if(active.begin(),active.end(),[&](std::size_t i){
      if(particles[i].position[2]>origin(offset)[2]+(WIDTH+2)*DZ){require(detected[i],"Discarded unobserved particle");++discardCount;return true;}return false;
    }),active.end());
    referenceSeconds+=elapsed(begin);
  }
  void diagnose(unsigned n) {
    const auto& f=fixture;const double dt=DZ/aprl::SI::c;
    long double dx=0,dy=0,x2=0,y2=0,energy=0,e2=0,delta=0;
    std::size_t held=0;const long double rest=static_cast<long double>(aprl::SI::electronMass)*aprl::SI::c*aprl::SI::c/aprl::SI::elementaryCharge;
    for(std::size_t i=0;i<particles.size();++i) {
      const auto q=admitted[i]?particles[i]:f.held(i,n);held+=q.frozen;
      const long double x=q.position[0],y=q.position[1];dx+=x;dy+=y;x2+=x*x;y2+=y*y;
      const auto& u=q.properVelocity;const auto& u0=f.seed[i].properVelocity;
      long double numerator=0;for(unsigned a=0;a<3;++a)numerator+=(static_cast<long double>(u[a])-u0[a])*(static_cast<long double>(u[a])+u0[a]);
      const long double dg=numerator/(gamma(u)+gamma(u0));
      const long double de=f.frame.gamma()*(dg+static_cast<long double>(f.frame.beta())*(static_cast<long double>(u[2])-u0[2]))*rest;delta+=de;
      const long double e=f.frame.gamma()*(gamma(u)+static_cast<long double>(f.frame.beta())*u[2])*rest;
      // Center energy on the initial mean scale to avoid subtracting squares
      // of ~600 MeV numbers when computing a narrow spread.
      const long double centered=e-LAB_GAMMA*rest;energy+=centered;e2+=centered*centered;
    }
    const long double N=particles.size(),mx=dx/N,my=dy/N,me=energy/N;
    history<<n<<','<<origin(offset)[2]<<','<<active.size()<<','<<held<<','<<admissionCount<<','<<std::count(released.begin(),released.end(),true)<<','<<std::count(detected.begin(),detected.end(),true)
      <<','<<delta/N<<','<<std::sqrt(std::max(0.L,e2/N-me*me))<<','<<std::sqrt(std::max(0.L,x2/N-mx*mx))<<','<<std::sqrt(std::max(0.L,y2/N-my*my))
      <<','<<mx<<','<<my<<','<<maxContinuity<<','<<maxChargeBalance<<'\n';
    const double z=f.moving(f.detector,n),k=z/DZ;
    for(std::size_t j=SIDE+2;j+SIDE+2<NY;j+=4)for(std::size_t i=SIDE+2;i+SIDE+2<NX;i+=4) {
      plane<<n<<','<<i<<','<<j<<','<<n*dt<<','<<(n-.5)*dt<<','<<f.frame.labTFromBoxZT(z,n*dt)
        <<','<<origin(0)[0]+(i+.5)*DX<<','<<origin(0)[1]+(j+.5)*DY<<','<<z;
      for(unsigned a=0;a<6;++a)plane<<','<<sampleComponent(*fields,offset,a,i+.5,j+.5,k);
      plane<<'\n';
    }
  }
  void finish(unsigned n) {
    for(std::size_t i=0;i<particles.size();++i)particleRow(final,admitted[i]?particles[i]:fixture.held(i,n),fixture.seed[i],n,0.,fixture.frame);
    require(initial.good()&&releases.good()&&crossings.good()&&final.good()&&plane.good()&&history.good(),"Scientific write failed");
  }
};

void run(const std::string& out,unsigned steps,unsigned particles,unsigned warm,
         const std::string& mode,unsigned warmHorizon=64,
         const std::string& initialModel="relativistic-poisson") {
  require(mode=="live"||mode=="frozen","Reference mode must be live or frozen");
  require(initialModel=="relativistic-poisson"||initialModel=="cowan-relativistic-poisson",
    "Initial model must be relativistic-poisson or cowan-relativistic-poisson");
  const auto model=initialModel=="cowan-relativistic-poisson"?
    aprl::InitialSelfFieldModel::CowanRelativisticPoisson:
    aprl::InitialSelfFieldModel::RelativisticPoisson;
  directory(out);Fixture f(steps,particles,warm,warmHorizon);Timing timing;const auto overall=Clock::now();
  auto donor=std::unique_ptr<EBFieldGrid>(new EBFieldGrid(geometry(f.nz)));
  donor->setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(new aprl::EBConvolutionalPML(donor->geometry(),f.nz,0,f.cpml,true,true)));
  MovingWindowRemap donorLayout(donor->geometry(),f.cpml,{-1});
  std::vector<RelativisticParticleSI> before;for(std::size_t i=0;i<f.seed.size();++i)before.push_back(f.held(i,-static_cast<int>(warm)));
  const std::size_t guards[3]={SIDE,SIDE,0};auto begin=Clock::now();
  const auto initialization=aprl::EBGaussFieldInitializer::initialize(*donor,donor->geometry(),0,origin(0),before,
    model,3e-5,guards,1e-8,10000,MPI_COMM_WORLD,false,true,&f.frame);
  timing.initialize=elapsed(begin);std::cout<<"Initialized production-matched fixture: full_nz="<<f.nz<<" cut_nz="<<WIDTH<<" iterations="<<initialization.iterations<<std::endl;
  begin=Clock::now();for(int n=-static_cast<int>(warm);n<0;++n)advanceHeld(*donor,f,n);timing.prewarm=elapsed(begin);
  // This records the actually prepared native state without changing it.
  // The norms are raw stored-value norms, not a field-energy invariant.
  begin=Clock::now();double preparedNorm[6]={},preparedMaximum[6]={};
  double preparedHistoryNorm[12]={};
  for(unsigned a=0;a<6;++a) {
    const auto& values=component(*donor,a);long double sum=0.;
    for(std::size_t i=0;i<values.size();++i) {
      const double value=values.data()[i];require(std::isfinite(value),"Nonfinite prepared reference field");
      sum+=static_cast<long double>(value)*value;
      preparedMaximum[a]=std::max(preparedMaximum[a],std::abs(value));
    }
    preparedNorm[a]=std::sqrt(sum);
  }
  const auto& preparedBoundary=dynamic_cast<const aprl::EBConvolutionalPML&>(donor->boundary());
  for(unsigned h=0;h<12;++h) {
    long double sum=0.;for(double value:preparedBoundary.checkpointHistory(h)) {
      require(std::isfinite(value),"Nonfinite prepared reference history");sum+=static_cast<long double>(value)*value;
    }
    preparedHistoryNorm[h]=std::sqrt(sum);
  }
  // Compare different preparation durations at the same t=0 particle state,
  // world lattice and native time levels. No interpolation or field mutation
  // enters this record. Component staggering and physical coordinate mapping
  // are stated in metrics; these raw samples are not a field-energy measure.
  const std::size_t preparedFirst=f.offset0>LOWER_GUARD?f.offset0-LOWER_GUARD:0;
  const std::size_t preparedLast=std::min(f.nz-1,f.offset0+WIDTH+32);
  const std::size_t preparedXYFirst=SIDE+2,preparedXYStop=NX-SIDE-2,preparedStride=4;
  std::size_t preparedRows=0;
  {
    auto prepared=fresh(out+"/prepared-reference.csv");
    prepared<<"component,i,j,global_k,native_time_box_s,value\n";
    for(unsigned a=0;a<6;++a) {
      const auto& values=component(*donor,a);
      const double nativeTime=a>=3?-.5*DZ/aprl::SI::c:0.;
      for(std::size_t k=preparedFirst;k<=preparedLast;k+=preparedStride)
        for(std::size_t j=preparedXYFirst;j<NY-SIDE-2;j+=preparedStride)
          for(std::size_t i=preparedXYFirst;i<preparedXYStop;i+=preparedStride) {
            require(i<values.nx()&&j<values.ny()&&k<values.nz(),
              "Prepared diagnostic outside native lattice");
            prepared<<a<<','<<i<<','<<j<<','<<k<<','<<nativeTime<<','<<values(i,j,k)<<'\n';
            ++preparedRows;
          }
    }
    prepared.flush();require(prepared.good(),"Prepared reference write failed");
  }
  timing.diagnostics+=elapsed(begin);
  Reference reference;reference.fields=donor.get();reference.layout=&donorLayout;reference.frozen=mode=="frozen";reference.beta=f.betaMean;
  World full(out+"/full",f,false,*donor,donorLayout,reference);
  World cut(out+"/cut",f,true,*donor,donorLayout,reference);
  auto overlap=fresh(out+"/overlap.csv");overlap<<"step,E_difference_L2,E_reference_L2,B_difference_L2,B_reference_L2,E_difference_max,B_difference_max,history_difference_max\n";
  auto inlet=fresh(out+"/inlet.csv");
  inlet<<"step,relative_k,global_k,component,scope,native_time_box_s,samples,reference_L2,full_L2,reference_full_difference_L2,reference_full_difference_max,candidate_present,candidate_full_difference_L2,candidate_reference_difference_L2\n";
  auto compare=[&](unsigned n){
    const auto begin=Clock::now();long double d[2]={0,0},r[2]={0,0};double maximum[2]={0,0};
    // Compare identical world positions only. Ex/Ey/Bz end faces and the
    // left supplied endpoint are excluded from the trusted PIC comparison.
    for(unsigned a=0;a<6;++a) {
      const auto& x=component(*cut.fields,a);const auto& y=component(*full.fields,a);const unsigned t=a>=3;
      for(std::size_t k=2;k+2<x.nz();++k)for(std::size_t j=SIDE+1;j+SIDE+1<x.ny();++j)for(std::size_t i=SIDE+1;i+SIDE+1<x.nx();++i) {
        const double base=y(i,j,cut.offset+k),v=x(i,j,k)-base;require(std::isfinite(base)&&std::isfinite(v),"Nonfinite overlap field");
        d[t]+=static_cast<long double>(v)*v;r[t]+=static_cast<long double>(base)*base;maximum[t]=std::max(maximum[t],std::abs(v));
      }
    }
    // Side history values differ once physical PIC fields reach them; compare
    // candidate versus actual full, not versus the prescribed held reference.
    double historyMax=0.;const auto& cb=dynamic_cast<const aprl::EBConvolutionalPML&>(cut.fields->boundary());
    const auto& fb=dynamic_cast<const aprl::EBConvolutionalPML&>(full.fields->boundary());
    for(unsigned h=0;h<12;++h){const auto& cl=cut.layout->historyLayout(h);const auto& fl=full.layout->historyLayout(h);if(cl.derivativeAxis==2)continue;const auto q=cl.dimensions[cl.derivativeAxis==0?1:0];
      for(std::size_t slot=0;slot<cl.activeAxisIndices.size();++slot)for(std::size_t k=2;k+2<cl.dimensions[2];++k)for(std::size_t j=0;j<q;++j){
        const auto a=(slot*cl.dimensions[2]+k)*q+j,b=(slot*fl.dimensions[2]+cut.offset+k)*q+j;
        const double v=cb.checkpointHistory(h)[a]-fb.checkpointHistory(h)[b];require(std::isfinite(v),"Nonfinite CPML history");historyMax=std::max(historyMax,std::abs(v));}}
    overlap<<n<<','<<std::sqrt(d[0])<<','<<std::sqrt(r[0])<<','<<std::sqrt(d[1])<<','<<std::sqrt(r[1])<<','<<maximum[0]<<','<<maximum[1]<<','<<historyMax<<'\n';
    // Read-only comparison of the inlet information against the actual full
    // simulation at the same native lattice/time. The full state is never
    // connected to the inlet provider. k=-1 is the next electric halo;
    // k=0 is the just-created slice; k=1 has already advanced inside the cut.
    // Each raw magnetic component retains its native half-step timestamp.
    for(int localK=-1;localK<=1;++localK)for(unsigned a=0;a<6;++a) {
      const auto globalK=static_cast<std::ptrdiff_t>(cut.offset)+localK;
      const auto& fullField=component(*full.fields,a);const auto& cutField=component(*cut.fields,a);
      require(globalK>=0&&globalK<static_cast<std::ptrdiff_t>(fullField.nz()),"Inlet diagnostic outside full native lattice");
      for(unsigned scope=0;scope<2;++scope) {
        long double reference2=0,full2=0,difference2=0,candidateDifference2=0,candidateReference2=0;
        double differenceMaximum=0.;std::size_t samples=0;
        const std::size_t lower=scope?SIDE+1:0;
        const std::size_t upperI=scope?fullField.nx()-SIDE-1:fullField.nx();
        const std::size_t upperJ=scope?fullField.ny()-SIDE-1:fullField.ny();
        for(std::size_t j=lower;j<upperJ;++j)for(std::size_t i=lower;i<upperI;++i) {
          const double referenceValue=Reference::query(&reference,a,i,j,globalK,n);
          const double fullValue=fullField(i,j,static_cast<std::size_t>(globalK));
          const double difference=referenceValue-fullValue;
          require(std::isfinite(referenceValue)&&std::isfinite(fullValue),"Nonfinite inlet diagnostic field");
          reference2+=static_cast<long double>(referenceValue)*referenceValue;
          full2+=static_cast<long double>(fullValue)*fullValue;
          difference2+=static_cast<long double>(difference)*difference;
          differenceMaximum=std::max(differenceMaximum,std::abs(difference));++samples;
          if(localK>=0) {
            const double candidate=cutField(i,j,static_cast<std::size_t>(localK));
            require(std::isfinite(candidate),"Nonfinite inlet candidate field");
            const double candidateFull=candidate-fullValue,candidateReference=candidate-referenceValue;
            candidateDifference2+=static_cast<long double>(candidateFull)*candidateFull;
            candidateReference2+=static_cast<long double>(candidateReference)*candidateReference;
          }
        }
        inlet<<n<<','<<localK<<','<<globalK<<','<<a<<','<<(scope?"physical_xy":"all_xy")<<','
          <<(static_cast<double>(n)-(a>=3?.5:0.))*DZ/aprl::SI::c<<','<<samples<<','
          <<std::sqrt(reference2)<<','<<std::sqrt(full2)<<','<<std::sqrt(difference2)<<','<<differenceMaximum<<','
          <<(localK>=0?1:0)<<','<<std::sqrt(candidateDifference2)<<','<<std::sqrt(candidateReference2)<<'\n';
      }
    }
    full.diagnose(n);cut.diagnose(n);timing.diagnostics+=elapsed(begin);
  };
  compare(0);
  for(unsigned n=0;n<steps;++n) {
    cut.prepare(n);full.prepare(n);
    if(!reference.frozen){begin=Clock::now();advanceHeld(*donor,f,n);reference.state=n+1;timing.donor+=elapsed(begin);}
    full.advance(n);cut.advance(n);cut.shift(n+1,reference);
    if((n+1)%4==0||n+1==steps)compare(n+1);
    if((n+1)%16==0)std::cout<<"step="<<n+1<<" full_s="<<full.advanceSeconds<<" cut_s="<<cut.advanceSeconds<<" donor_s="<<timing.donor<<std::endl;
  }
  full.finish(steps);cut.finish(steps);
  require(full.released==cut.released&&full.detected==cut.detected,"Fixed cohort event membership differs");
  require(cut.admissionCount>cut.initialActive,"No held source actually entered through moving left boundary");
  require(std::count(full.released.begin(),full.released.end(),true)>0,"No real PIC release in bounded case");
  require(cut.maxContinuity<1e-10&&full.maxContinuity<1e-10,"Charge continuity failed");
  require(full.maximumLabZ<TRUE_DETECTOR_LAB&&cut.maximumLabZ<TRUE_DETECTOR_LAB,"Unexpected production detector reach");
  double earliest=1e100,latest=-1e100;for(const auto& p:f.seed){const double t=(f.release-p.position[2])/((p.frozenBetaZ+f.frame.beta())*DZ);earliest=std::min(earliest,t);latest=std::max(latest,t);}
  auto metrics=fresh(out+"/metrics.json");
  metrics<<"{\n\"format\":\"aprl-held-cut-bunch-v1\",\"completed\":true,\"mode\":\""<<mode<<"\",\"steps\":"<<steps<<",\"particles\":"<<particles<<",\"held_prewarm_steps\":"<<warm
    <<",\"warm_horizon_steps\":"<<f.warmHorizon
    <<",\"compute_threads\":"<<aprl::computeThreadCount()<<",\"nx\":"<<NX<<",\"ny\":"<<NY<<",\"full_nz\":"<<f.nz<<",\"cut_nz\":"<<WIDTH<<",\"dx_m\":"<<DX<<",\"dy_m\":"<<DY<<",\"dz_m\":"<<DZ<<",\"dt_s\":"<<DZ/aprl::SI::c
    <<",\"reference_offset0_cells\":"<<f.offset0<<",\"reference_causal_lower_margin_cells\":"<<f.offset0-(2*steps+warm+2)<<",\"reference_causal_upper_margin_cells\":"<<f.nz-(f.offset0+WIDTH)-(steps+warm+1)
    <<",\"boost_gamma\":"<<BOOST_GAMMA<<",\"window_beta\":"<<f.frame.beta()<<",\"particle_beta_mean\":"<<f.betaMean<<",\"particle_beta_min\":"<<f.betaMin<<",\"particle_beta_max\":"<<f.betaMax
    <<",\"lab_gamma_mean_parameter\":"<<LAB_GAMMA<<",\"lab_gamma_sigma_parameter\":"<<LAB_GAMMA_SIGMA<<",\"sigma_position_cells\":[12.78982675,12.2478603,8],\"source_sigma_z_cells\":"<<SOURCE_SIGMA_Z_CELLS<<",\"charge_C\":"<<TOTAL_CHARGE
    <<",\"lab_u_mean_parameter\":["<<LAB_U_MEAN[0]<<','<<LAB_U_MEAN[1]<<','<<LAB_U_MEAN[2]<<"],\"lab_u_sigma_parameter\":["<<LAB_U_SIGMA[0]<<','<<LAB_U_SIGMA[1]<<','<<LAB_U_SIGMA[2]<<"],\"transverse_position_u_correlation_parameter\":["<<TRANSVERSE_CORRELATION[0]<<','<<TRANSVERSE_CORRELATION[1]<<']'
    <<",\"source_moment_sample_records\":65536,\"source_moment_sample_sha256\":\"4126d7f252fe71934892e063bd9a263070d3227e718d4b70e7990129182f05d3\",\"source_moment_sample_path\":\"results/held_cut_20261028/source-moment-sample.json\""
    <<",\"cpml_cells\":[10,10,0],\"cpml_order\":3,\"cpml_kappa\":8,\"cpml_alpha\":0,\"cpml_reflection\":1e-8,\"lower_guard_cells\":64,\"upper_guard_cells\":128,\"particle_substeps\":"<<f.substeps
    <<",\"initial_self_field_model\":\""<<initialModel<<"\",\"initial_tolerance\":1e-8,\"initial_maximum_iterations\":10000,\"initial_maximum_mean_transverse_beta\":3e-5"
    <<",\"undulator_K\":1.092714657005,\"undulator_period_m\":0.0234,\"undulator_periods\":170,\"undulator_entrance_lab_m\":0,\"undulator_polarization_rad\":1.5707963267948966,\"undulator_gaussian_fringe\":true,\"undulator_fringe_cutoff\":1e-9"
    <<",\"activation_plane_lab_m\":"<<RELEASE_LAB<<",\"auxiliary_plane_lab_m\":"<<f.frame.labZFromBoxZT(f.detector,0.)<<",\"production_detector_lab_m\":4.5,\"production_detector_reached\":false"
    <<",\"source_release_step_min\":"<<earliest<<",\"source_release_step_max\":"<<latest<<",\"full_maximum_lab_z_m\":"<<full.maximumLabZ<<",\"cut_maximum_lab_z_m\":"<<cut.maximumLabZ
    <<",\"initial_cut_active\":"<<cut.initialActive<<",\"cut_admissions\":"<<cut.admissionCount<<",\"cut_peak_active\":"<<cut.peakActive<<",\"cut_discards\":"<<cut.discardCount
    <<",\"released\":"<<std::count(full.released.begin(),full.released.end(),true)<<",\"auxiliary_detected\":"<<std::count(full.detected.begin(),full.detected.end(),true)
    <<",\"full_continuity_max\":"<<full.maxContinuity<<",\"cut_continuity_max\":"<<cut.maxContinuity<<",\"full_charge_balance_max_C\":"<<full.maxChargeBalance<<",\"cut_charge_balance_max_C\":"<<cut.maxChargeBalance
    <<",\"minimum_cut_pic_clearance_m\":"<<cut.minPicClearance<<",\"shift_count\":"<<steps<<",\"field_shift_bytes\":"<<cut.shiftBytes<<",\"history_shift_bytes\":"<<cut.historyBytes
    <<",\"initial_iterations\":"<<initialization.iterations<<",\"initial_relative_residual\":"<<initialization.relativeResidual<<",\"initial_temporary_bytes\":"<<initialization.temporaryBytes
    <<",\"prepared_reference_native_component_L2\":["<<preparedNorm[0]<<','<<preparedNorm[1]<<','<<preparedNorm[2]<<','<<preparedNorm[3]<<','<<preparedNorm[4]<<','<<preparedNorm[5]<<']'
    <<",\"prepared_reference_native_component_maximum\":["<<preparedMaximum[0]<<','<<preparedMaximum[1]<<','<<preparedMaximum[2]<<','<<preparedMaximum[3]<<','<<preparedMaximum[4]<<','<<preparedMaximum[5]<<']'
    <<",\"prepared_reference_roi\":{\"file\":\"prepared-reference.csv\",\"global_k_first\":"<<preparedFirst<<",\"global_k_last_inclusive\":"<<preparedLast
    <<",\"xy_first_index\":"<<preparedXYFirst<<",\"x_stop_exclusive\":"<<preparedXYStop<<",\"y_stop_exclusive\":"<<NY-SIDE-2
    <<",\"stride_all_axes\":"<<preparedStride<<",\"rows\":"<<preparedRows
    <<",\"origin_m\":["<<origin(0)[0]<<','<<origin(0)[1]<<",0],\"native_stagger_cells\":[[0.5,0,0],[0,0.5,0],[0,0,0.5],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]],\"contract\":\"Direct native E(0), B(-dt/2) samples in the same world ROI; component numbers Ex,Ey,Ez,Bx,By,Bz; electric values V/m, magnetic values T; no resampling, physical-energy or field-collocation claim\"}"
    <<",\"prepared_reference_history_L2\":[";
  for(unsigned h=0;h<12;++h)metrics<<(h?",":"")<<preparedHistoryNorm[h];
  metrics<<']'
    <<",\"initialization_wall_s\":"<<timing.initialize<<",\"held_prewarm_wall_s\":"<<timing.prewarm<<",\"live_reference_wall_s\":"<<timing.donor<<",\"full_advance_wall_s\":"<<full.advanceSeconds<<",\"cut_advance_wall_s\":"<<cut.advanceSeconds
    <<",\"full_preparation_wall_s\":"<<full.preparationSeconds<<",\"cut_preparation_wall_s\":"<<cut.preparationSeconds
    <<",\"cut_remap_and_feed_wall_s\":"<<cut.referenceSeconds<<",\"diagnostics_wall_s\":"<<timing.diagnostics<<",\"overall_wall_s\":"<<elapsed(overall)
    <<",\"full_field_bytes\":"<<full.fields->memoryFootprint().totalBytes()<<",\"cut_field_bytes\":"<<cut.fields->memoryFootprint().totalBytes()<<",\"reference_field_bytes\":"<<donor->memoryFootprint().totalBytes()
    <<",\"reference_kind\":\""<<(reference.frozen?"short-held-prepared frozen native field and side-history translated at mean velocity; approximation with source velocity spread":"independent prescribed held-source Maxwell evolution; no actual PIC future")<<"\""
    <<",\"inlet_diagnostic_contract\":\"Read-only native comparisons at left relative k=-1,0,1; full current physical state is diagnostic only and never supplies the candidate. Per-component times E(n),B(n-1/2). Prepared field/history norms are raw stored-value norms, not physical energy.\""
    <<",\"reference_dependency_contract\":\"Fixed-domain endpoint dynamic influence excluded by independently checked native Cowan/transverse-CPML lower n+1 and upper n bounds, including preparation; this does not establish independence of the initial elliptic solve from its finite domain.\""
    <<",\"provenance\":\"Production grid/boost/CPML/guards/undulator/subcycling and transverse RMS from 18932; laboratory momentum means/spreads and x-ux/y-uy correlations from deterministic 65536-record original-input sample; Gaussian z shortened and charge scaled by sigma ratio to preserve peak line density\""
    <<",\"limitations\":[\"Verification-only single-rank full/cut comparison; not full apparatus or MPI acceptance\",\"Fixed reference domain uses transverse CPML and remote PEC z endpoints with conservative causal padding\",\"Auxiliary local plane is not replacement for production 4.5m detector\",\"B is recorded at native time n-1/2 but at the E-plane spatial position z(n); paired boost-field comparison only, not simultaneous lab E/B\",\"Gaussian surrogate matches selected measured moments and transverse correlations, not the complete original six-dimensional particle distribution\",\"All reference initialization/preparation/evolution costs must be counted; live reference retains full-domain cost\",\"Full input tables and prepared reference are retained; no linear complexity or production speed claim\"]}\n";
  require(metrics.good()&&overlap.good()&&inlet.good(),"Final output write failed");std::cout<<"held-cut probe completed mode="<<mode<<" steps="<<steps<<" admitted="<<cut.admissionCount<<" released="<<std::count(full.released.begin(),full.released.end(),true)<<std::endl;
}
}
int main(int argc,char** argv) {
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int size=0;MPI_Comm_size(MPI_COMM_WORLD,&size);
  try {aprl::configureComputeThreads(provided);require(size==1,"Held cut research requires one MPI rank");require(argc>=6&&argc<=8,"Usage: aprl_held_cut_probe OUT STEPS PARTICLES PREWARM live|frozen [WARM_HORIZON [INITIAL_MODEL]]");
    run(argv[1],std::stoul(argv[2]),std::stoul(argv[3]),std::stoul(argv[4]),argv[5],
      argc>=7?std::stoul(argv[6]):64,argc>=8?argv[7]:"relativistic-poisson");MPI_Finalize();return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,1);return 1;}
}
