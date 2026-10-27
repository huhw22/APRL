/* CPML integration probe: unchanged kernels, WorldPosition remap and existing
 * frozen-reference input. The moving CPML acts on total fields.
 * Bounded research probe: independent frozen reference feed around a moving
 * particle/PIC window. This is deliberately not a production window switch.
 * Ordinary short-window modes only read their initial full-domain Cowan solve.
 * Explicit diagnostic oracle modes replay a completed full-domain history
 * solely to isolate errors; they are not usable independent window algorithms.
 * Incoming held trajectories are selected by a sorted queue.
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
#include "eb_deposition.h"
#include "eb_initial_field.h"
#include "eb_particles.h"
#include "particle_initializer.h"
#include "streaming_reference_boundary.h"
#include "streaming_cpml_boundary.h"
#include "streaming_field_oracle.h"
#include "moving_window_remap.h"
#include "eb_cpml.h"

namespace {
using aprl::Double;
using aprl::EBFieldGrid;
using aprl::EBGridGeometry;
using aprl::FieldVector;
using aprl::RelativisticParticleSI;
using aprl::YeeComponent;
using aprl::verification::StreamingReferenceCache;
using aprl::verification::ReferenceFedCPML;
using aprl::verification::MovingWindowRemap;
using aprl::verification::MovingWindowHistoryPolicy;
using aprl::verification::NativeOracleConfig;
using aprl::verification::NativeFieldOracleReader;
using aprl::verification::NativeFieldOracleWriter;
using Clock=std::chrono::steady_clock;
const std::size_t NX=24,NY=24,NZ=424,PARTICLES=1024;
const unsigned STEPS=176;
const std::size_t CPML_Z=20,PHYSICAL_NZ=384;
const double DZ=1e-6,DX=2e-6,DY=2e-6,ZORIGIN=-84e-6;
const double RELEASE=218e-6,DETECTOR=234e-6,WINDOW_BETA=.8;
const double TOTAL_CHARGE=-.1e-12;
const unsigned SEED=710271;
double gaussSourceNorm=0.,gaussSourceMaximum=0.;

void require(bool value,const std::string& message) {
  if(!value)throw std::runtime_error(message);
}
double elapsed(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now()-start).count();
}
std::ofstream fresh(const std::string& path) {
  require(!std::ifstream(path).good(),"Refusing to overwrite "+path);
  std::ofstream result(path);
  require(result.good(),"Cannot write "+path);
  result<<std::setprecision(17);return result;
}
YeeComponent& component(EBFieldGrid& f,unsigned a) {
  if(a==0)return f.ex();
  if(a==1)return f.ey();
  if(a==2)return f.ez();
  if(a==3)return f.bx();
  if(a==4)return f.by();
  return f.bz();
}
const YeeComponent& component(const EBFieldGrid& f,unsigned a) {
  return component(const_cast<EBFieldGrid&>(f),a);
}

// Read-only, opt-in native-field snapshots for locating interface errors.
// The binary stream contains one x-z section per component; no interpolation,
// reference lookup, or field mutation occurs here. A record's state step is
// n+1 for all three stages of loop n, including the mixed-time magnetic stage.
class SectionDiagnostics {
public:
  SectionDiagnostics(const std::string& directory,unsigned maximumStep)
    :directory_(directory),maximumStep_(maximumStep) {
    static_assert(sizeof(Double)==8 && std::numeric_limits<Double>::is_iec559,
                  "Section diagnostics require IEEE binary64 fields");
    const std::uint16_t one=1;
    require(*reinterpret_cast<const unsigned char*>(&one)==1,
            "Section diagnostics currently require little-endian storage");
    const std::string path=directory_+"/field-sections.bin";
    require(!std::ifstream(path).good(),"Refusing to overwrite "+path);
    output_.open(path,std::ios::out|std::ios::binary);
    require(output_.good(),"Cannot write "+path);
    output_.write("APRLXZ01",8);
    require(output_.good(),"Cannot write section diagnostic header");
  }

  void record(const EBFieldGrid& fields,std::size_t offset,
              unsigned stateStep,unsigned stage) {
    if(stateStep>maximumStep_)return;
    require(stage<=3 && ((stage==0)==(stateStep==0)),
            "Invalid section diagnostic stage");
    const std::uint64_t header[4]={stateStep,stage,offset,6};
    output_.write(reinterpret_cast<const char*>(header),sizeof(header));
    for(unsigned a=0;a<6;++a) {
      const auto& c=component(fields,a);
      const std::size_t sectionJ=NY/2+2;
      require(sectionJ<c.ny(),"Diagnostic section outside component");
      const std::uint64_t description[5]={a,c.nx(),c.ny(),c.nz(),sectionJ};
      output_.write(reinterpret_cast<const char*>(description),sizeof(description));
      for(std::size_t k=0;k<c.nz();++k)
        output_.write(reinterpret_cast<const char*>(
          c.data()+(k*c.ny()+sectionJ)*c.nx()),c.nx()*sizeof(Double));
    }
    require(output_.good(),"Field section diagnostic write failed");
    ++records_;
  }

  void finish() {
    require(records_==1+3*maximumStep_,"Incomplete section diagnostics");
    output_.flush();require(output_.good(),"Field section flush failed");
    output_.close();require(!output_.fail(),"Field section close failed");
    auto manifest=fresh(directory_+"/field-sections.json");
    manifest<<"{\n\"format\":\"aprl-native-field-xz-sections-v1\",\"completed\":true"
      <<",\"binary_file\":\"field-sections.bin\",\"byte_order\":\"little\""
      <<",\"magic_ascii\":\"APRLXZ01\",\"maximum_state_step\":"<<maximumStep_
      <<",\"records\":"<<records_<<",\"section_j\":"<<NY/2+2
      <<",\"record_header\":\"four uint64: state_step, stage, integer_z_offset, component_count\""
      <<",\"component_header\":\"five uint64: component, nx, ny, nz, section_j\""
      <<",\"component_payload\":\"nx*nz IEEE binary64, k outer and i inner; raw native values at fixed component j\""
      <<",\"components\":[\"Ex\",\"Ey\",\"Ez\",\"Bx\",\"By\",\"Bz\"]"
      <<",\"stages\":{\"0\":\"initial E(0), B(-1/2)\","
        "\"1\":\"after magnetic including boundary hook: E(state_step-1), B(state_step-1/2)\","
        "\"2\":\"after electric including boundary hook, before remap: E(state_step), B(state_step-1/2)\","
        "\"3\":\"after remap, or unchanged if no shift: E(state_step), B(state_step-1/2)\"}"
      <<",\"world_coordinate\":\"Use integer_z_offset+k plus each component's original Yee staggering; section is not collocated\""
      <<",\"contract\":\"Read-only raw field records; no reference queries, interpolation, or particle changes; diagnostic timing is not performance evidence\"}\n";
    require(manifest.good(),"Section diagnostic manifest write failed");
  }

private:
  std::string directory_;
  unsigned maximumStep_;
  std::uint64_t records_=0;
  std::ofstream output_;
};

long double gamma(const FieldVector<double>& u) {
  return std::sqrt(1.L+static_cast<long double>(u[0])*u[0]+
    static_cast<long double>(u[1])*u[1]+static_cast<long double>(u[2])*u[2]);
}
EBGridGeometry geometry(std::size_t nz) {
  return EBGridGeometry(NX,NY,nz,DX,DY,DZ,DZ/aprl::SI::c,
    aprl::EBMaxwellSolver::CowanZ);
}
FieldVector<double> origin(std::size_t offset) {
  FieldVector<double> p(0.);p[0]=-.5*NX*DX;p[1]=-.5*NY*DY;
  p[2]=ZORIGIN+offset*DZ;return p;
}
std::size_t windowOffset(unsigned step,unsigned guard) {
  // Keep the physical release/detector planes fixed in lab. The integer
  // origin never drifts through repeated floating point subtraction.
  const long long cell=static_cast<long long>(std::floor(218.-WINDOW_BETA*step));
  require(cell+84>=static_cast<long long>(guard+CPML_Z),"Negative window offset");
  return static_cast<std::size_t>(cell+84-guard-CPML_Z);
}
double movingPlane(double initial,unsigned step,double fraction=0.) {
  return initial-WINDOW_BETA*DZ*(step+fraction);
}

struct Counters {
  std::uint64_t particleInspections=0,referenceInspections=0,diagnosticParticleInspections=0,depositedSegments=0;
  std::uint64_t shiftBytes=0,newLayerValues=0,gaussCorrectionValues=0;
  std::uint64_t fieldSamples=0;
  std::size_t activePeak=0,queueAdmissions=0,evicted=0,ballisticOutgoingSteps=0;
  std::uint64_t cpmlHistoryMovedBytes=0,cpmlHistoryZeroedBytes=0;
  double cpmlHistoryPeak=0.;bool cpmlHistorySeen=false;
  std::uint64_t heldInCpmlObservations=0,ballisticInCpmlObservations=0,picPushesInCpml=0;
  double minFrameReleaseInner=1.,minFrameDetectorInner=1.;
  double minEventReleaseInner=1.,minEventDetectorInner=1.,minPicInnerClearance=1.;
  double initialization=0,cache=0,fieldAdvance=0,particles=0,deposition=0;
  double boundary=0,shift=0,diagnostics=0,gaussCorrection=0;
  double maxContinuity=0,maxChargeBalance=0,maxRemapBalance=0;
  double maxGaussBeforeShift=0,maxGaussAfterShift=0,maxReferenceGauss=0;
  double maxEzCorrection=0;long double sumEzCorrection2=0;
};
struct ParticleState {RelativisticParticleSI p;std::size_t index;};
struct Segment {
  FieldVector<double> start,end;double first,last;
  FieldVector<double> startU,endU;
  Segment()=default;
  Segment(const FieldVector<double>& a,const FieldVector<double>& b,double f,double l,
          const FieldVector<double>& ua,const FieldVector<double>& ub)
    :start(a),end(b),first(f),last(l),startU(ua),endU(ub){}
};
struct Gauss {double relative=0,maximum=0,newLayerRelative=0;};

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

void cropRho(const YeeComponent& ext,YeeComponent& out,std::size_t start) {
  require(start+out.nz()<=ext.nz(),"Charge crop left deposition collar");
  const std::size_t plane=out.nx()*out.ny();
  std::copy(ext.data()+start*plane,ext.data()+(start+out.nz())*plane,out.data());
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
long double chargePlane(const YeeComponent& rho,std::size_t k) {
  long double result=0;
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
Gauss gauss(const EBFieldGrid& f,const YeeComponent& rho) {
  const auto& g=f.geometry();Gauss out;long double residual2=0,source2=0;
  double layer=0,scale=0;
  for(std::size_t k=CPML_Z+1;k<g.nz-CPML_Z;++k)for(std::size_t j=1;j<g.ny;++j)
    for(std::size_t i=1;i<g.nx;++i) {
      const double source=rho(i,j,k)/aprl::SI::epsilon0;
      const double value=(f.ex()(i,j,k)-f.ex()(i-1,j,k))/DX+
        (f.ey()(i,j,k)-f.ey()(i,j-1,k))/DY+
        (f.ez()(i,j,k)-f.ez()(i,j,k-1))/DZ-source;
      require(std::isfinite(value),"Non-finite total field Gauss residual");
      residual2+=static_cast<long double>(value)*value;source2+=static_cast<long double>(source)*source;
      out.maximum=std::max(out.maximum,std::abs(value));scale=std::max(scale,std::abs(source));
      if(k==CPML_Z+1)layer=std::max(layer,std::abs(value));
    }
  (void)source2;(void)scale;
  out.relative=std::sqrt(residual2)/gaussSourceNorm;
  out.newLayerRelative=layer/gaussSourceMaximum;return out;
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

std::vector<RelativisticParticleSI> initialParticles(double beta) {
  aprl::SimulationConfig c;c.beam.type=aprl::BeamInputType::GeneratedGaussian;
  c.beam.macroparticles=PARTICLES;c.beam.electrons=-TOTAL_CHARGE/aprl::SI::elementaryCharge;
  c.beam.randomSeed=SEED;c.beam.gamma=1/std::sqrt(1-beta*beta);
  c.beam.direction=FieldVector<double>(0.);c.beam.direction[2]=1.;
  c.beam.center=FieldVector<double>(0.);c.beam.center[2]=160e-6;
  c.beam.positionOffset=FieldVector<double>(0.);
  c.beam.sigmaPosition=FieldVector<double>(3e-6);c.beam.sigmaPosition[2]=12e-6;
  c.beam.sigmaProperVelocity=FieldVector<double>(0.);
  aprl::ParticleInitializationReport report;
  auto p=aprl::ParticleInitializer::create(c,MPI_COMM_WORLD,report);
  for(auto& v:p){v.id=v.sourceId;v.frozen=true;v.frozenBetaZ=beta;}
  return p;
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

void run(const std::string& directory,double beta,unsigned guard,const std::string& mode,
         unsigned diagnosticSteps=0,const std::string& suppliedOracleDirectory="") {
  const bool prescribedReference=mode=="held-record";
  const bool prescribedActualControl=mode=="residual-held-zero";
  const bool prescribedMotion=prescribedReference||prescribedActualControl;
  const bool heldReference=mode=="residual-held-reference"||prescribedActualControl;
  const bool recording=mode=="full-record"||prescribedReference,full=mode=="full"||recording;
  const bool oracleFeed=mode=="residual-oracle-feed"||mode=="residual-oracle-all"||heldReference;
  const bool oracleCpml=mode=="residual-oracle-cpml"||mode=="residual-oracle-all"||heldReference;
  const bool oracleReplay=oracleFeed||oracleCpml;
  const bool zero=mode=="zero",refresh=mode=="residual-refresh";
  const bool nativeB=mode=="residual-native-b";
  const bool residual=mode=="residual"||refresh||oracleReplay||nativeB;
  require(full||mode=="cpml"||zero||residual,"Unknown mode");
  require(beta==0.||std::abs(beta-.23)<1e-15,"Use beta 0 or .23");
  require(guard==16||guard==24||guard==40,"Use the fixed nonabsorbing guards 16/24/40");
  require(diagnosticSteps<=STEPS,"Diagnostic step limit exceeds this bounded fixture");
  require(!oracleReplay||!suppliedOracleDirectory.empty(),"Oracle replay requires its completed recording directory");
  require(suppliedOracleDirectory.empty()||recording||oracleReplay,"Oracle directory supplied to a non-oracle mode");
  const std::string oracleDirectory=recording&&suppliedOracleDirectory.empty()?
    directory+"/native-oracle":suppliedOracleDirectory;
  if(mkdir(directory.c_str(),0755)!=0)require(errno==EEXIST,"Cannot create result directory");
  auto history=fresh(directory+"/history.csv"),crossings=fresh(directory+"/detector-crossings.csv");
  auto releases=fresh(directory+"/release-particles.csv"),initialCsv=fresh(directory+"/initial-particles.csv");
  auto plane=fresh(directory+"/field-plane.csv");
  particleHeader(crossings);particleHeader(releases);particleHeader(initialCsv);
  history<<"step,time_box_s,origin_z_m,active,held,admitted,evicted,charge_interior_C,charge_balance_C,remap_balance_C,continuity_relative,total_gauss_before_shift,total_gauss_after_shift,new_layer_gauss_relative,reference_gauss_relative,reference_gauss_diagnosed,reference_first_global_k,reference_last_global_k,ez_correction_max_V_m,field_Ex_center_V_m,field_Ey_center_V_m,field_Ez_center_V_m\n";
  plane<<"step,pixel_i,pixel_j,time_E_box_s,time_B_box_s,time_E_lab_s,x_m,y_m,z_box_m,Ex_V_m,Ey_V_m,Ez_V_m,Bx_T,By_T,Bz_T\n";
  Counters count;const auto overallStart=Clock::now();
  auto initial=initialParticles(beta);
  require(initial.size()==PARTICLES,"Generated particle count differs");
  aprl::BoostFrameTransform frame;frame.setFromGamma(5./3.,aprl::SI::c,0.);
  const auto global=geometry(NZ);const double dt=global.dt;
  const NativeOracleConfig oracleConfig(global,STEPS,guard,SEED,beta,TOTAL_CHARGE,
                                        ZORIGIN,RELEASE,DETECTOR,WINDOW_BETA,
                                        (prescribedReference||heldReference)?1:0);
  std::unique_ptr<NativeFieldOracleReader> oracleReader;
  std::unique_ptr<NativeFieldOracleWriter> oracleWriter;
  if(oracleReplay)oracleReader.reset(new NativeFieldOracleReader(oracleDirectory,oracleConfig));
  if(recording)oracleWriter.reset(new NativeFieldOracleWriter(oracleDirectory,oracleConfig));
  std::vector<std::size_t> sorted(initial.size());std::iota(sorted.begin(),sorted.end(),0);
  std::sort(sorted.begin(),sorted.end(),[&](std::size_t a,std::size_t b){return initial[a].position[2]<initial[b].position[2];});
  for(const auto& p:initial) {
    require(inside(p.position,origin(0),global),"Initial source outside full reference domain");
    require(p.position[2]<RELEASE,"Fixture has an already-released initial particle");
    require(p.position[2]+beta*DZ*STEPS>movingPlane(DETECTOR,STEPS)+2*DZ,
      "176 steps do not cover every prescribed detector crossing");
    particleRow(initialCsv,p,p,0,0,frame);
  }
  const auto initializationStart=Clock::now();
  std::unique_ptr<EBFieldGrid> whole(new EBFieldGrid(global));
  aprl::EBCPMLParameters cpmlParameters;
  cpmlParameters.cells[0]=cpmlParameters.cells[1]=0;cpmlParameters.cells[2]=CPML_Z;
  whole->setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(
    new aprl::EBConvolutionalPML(global,NZ,0,cpmlParameters,true,true)));
  const std::size_t guards[3]={0,0,CPML_Z};
  const auto init=aprl::EBGaussFieldInitializer::initialize(*whole,global,0,origin(0),initial,
    aprl::InitialSelfFieldModel::CowanRelativisticPoisson,1e-8,guards,1e-10,5000,MPI_COMM_WORLD);
  count.initialization=elapsed(initializationStart);
  {
    YeeComponent initialRho(NX+1,NY+1,NZ+1);initialRho.fill(0.);
    aprl::ChargeConservingCurrentDepositor initialCharge(*whole,origin(0));
    for(const auto& p:initial)initialCharge.depositCharge(p.position,p.charge,initialRho);
    long double norm=0;
    for(std::size_t i=0;i<initialRho.size();++i) {
      const double value=initialRho.data()[i]/aprl::SI::epsilon0;
      norm+=static_cast<long double>(value)*value;
      gaussSourceMaximum=std::max(gaussSourceMaximum,std::abs(value));
    }
    gaussSourceNorm=std::sqrt(norm);require(gaussSourceNorm>0.,"Empty Gauss normalization source");
  }
  const auto cacheStart=Clock::now();StreamingReferenceCache reference(*whole,beta);
  struct NativeBInput {
    const StreamingReferenceCache* original;
    double beta;
    static double query(void* raw,unsigned a,std::size_t i,std::size_t j,double k,double state) {
      const auto& input=*static_cast<NativeBInput*>(raw);
      if((a!=3&&a!=4)||input.beta==0.||state==0.)
        return input.original->value(a,i,j,k,state);
      const auto& g=input.original->geometry();
      const double alpha=input.beta*aprl::SI::c*g.dt/g.dz;
      const double coordinate=k+.5-alpha*(state-.5);
      const double scale=input.beta*aprl::SI::c*aprl::SI::epsilon0*aprl::SI::mu0;
      return (a==3?-scale:scale)*input.original->value(a==3?1:0,i,j,coordinate,0.);
    }
  } nativeBInput{&reference,beta};
  std::unique_ptr<StreamingReferenceCache> feedReference,cpmlReference;
  if(oracleFeed) {
    feedReference.reset(new StreamingReferenceCache(*whole,beta));
    feedReference->setDiagnosticQuery(&NativeFieldOracleReader::query,oracleReader.get());
  }
  if(oracleCpml) {
    cpmlReference.reset(new StreamingReferenceCache(*whole,beta));
    cpmlReference->setDiagnosticQuery(&NativeFieldOracleReader::query,oracleReader.get());
  }
  if(nativeB) {
    feedReference.reset(new StreamingReferenceCache(*whole,beta));
    feedReference->setDiagnosticQuery(&NativeBInput::query,&nativeBInput);
  }
  count.cache=elapsed(cacheStart);
  std::size_t offset=full?0:windowOffset(0,guard);
  const std::size_t width=full?NZ:16+2*guard+2*CPML_Z;
  const std::size_t fullBytes=whole->memoryFootprint().totalBytes();
  std::unique_ptr<EBFieldGrid> grid;
  ReferenceFedCPML* boundary=nullptr;
  std::unique_ptr<MovingWindowRemap> mover;
  if(full)grid=std::move(whole);
  else {
    grid.reset(new EBFieldGrid(geometry(width)));
    require(offset>2&&offset+width+2<NZ,"Initial moving-window extent outside cache");
    require(windowOffset(STEPS,guard)>2+static_cast<std::size_t>(std::ceil(beta*(STEPS+1))),
      "Future moving-window cache halo would require extrapolation");
    for(unsigned a=0;a<6;++a) {
      auto& f=component(*grid,a);
      for(std::size_t k=0;k<f.nz();++k)for(std::size_t j=0;j<f.ny();++j)
        for(std::size_t i=0;i<f.nx();++i)f(i,j,k)=reference.value(a,i,j,offset+k,0.);
    }
    // Oracle channels are independent: replacing exterior feed must not
    // silently replace the CPML subtraction reference as well.
    boundary=new ReferenceFedCPML(grid->geometry(),cpmlParameters,
      feedReference?*feedReference:reference,offset,zero,residual,
      oracleReplay?(cpmlReference?cpmlReference.get():&reference):nullptr);
    grid->setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(boundary));
    mover.reset(new MovingWindowRemap(grid->geometry(),cpmlParameters,{-1},
                                      MovingWindowHistoryPolicy::WorldPosition));
    whole.reset();
  }
  EBFieldGrid& fields=*grid;
  if(oracleWriter)oracleWriter->record(fields,0,windowOffset(0,guard)-2);
  std::unique_ptr<SectionDiagnostics> sectionDiagnostics;
  if(diagnosticSteps) {
    sectionDiagnostics.reset(new SectionDiagnostics(directory,diagnosticSteps));
    sectionDiagnostics->record(fields,offset,0,0);
  }
  // Only the deposition J/rho use this fixed two-cell longitudinal collar.
  // Its E/B are never advanced or sampled as physical fields.
  const std::size_t collar=full?0:2;
  EBFieldGrid deposit(geometry(width+2*collar));
  YeeComponent rhoExt(NX+1,NY+1,width+2*collar+1),nextExt(NX+1,NY+1,width+2*collar+1);
  YeeComponent rho(NX+1,NY+1,width+1),nextRho(NX+1,NY+1,width+1);
  YeeComponent referenceRho(NX+1,NY+1,width+1);
  std::vector<ParticleState> active;active.reserve(PARTICLES);
  std::vector<bool> released(PARTICLES,false),detected(PARTICLES,false);
  std::size_t queue=sorted.size();
  if(full)for(std::size_t i=0;i<initial.size();++i)active.push_back({initial[i],i});
  if(full){queue=0;count.queueAdmissions=active.size();}
  const aprl::SIFieldSourceSet sources;
  double lastContinuity=0,lastPreGauss=0,lastPostGauss=0,lastCorrection=0;
  long double lastBalance=0,lastRemap=0;
  std::size_t referenceFirstK=0,referenceLastK=0;

  auto admit=[&](unsigned step,std::size_t nextOffset) {
    if(full)return;
    const double lower=origin(nextOffset-collar)[2];
    while(queue) {
      const std::size_t index=sorted[queue-1];++count.particleInspections;
      if(initial[index].position[2]+beta*DZ*(step+1)<lower)break;
      --queue;RelativisticParticleSI p=initial[index];p.position[2]+=beta*DZ*step;
      require(p.position[2]<movingPlane(RELEASE,step),"Incoming held particle already passed activation plane");
      active.push_back({p,index});++count.queueAdmissions;
    }
  };
  auto depositCharges=[&](YeeComponent& out,const FieldVector<double>& o) {
    out.fill(0.);aprl::ChargeConservingCurrentDepositor depositor(deposit,o);
    for(const auto& state:active) {
      ++count.particleInspections;
      if(inside(state.p.position,o,deposit.geometry()))depositor.depositCharge(state.p.position,state.p.charge,out);
    }
  };
  auto referenceGauss=[&](unsigned step)->double {
    referenceRho.fill(0.);const auto o=origin(offset);
    const double lo=o[2]-beta*DZ*step,hi=o[2]+width*DZ-beta*DZ*step;
    auto first=std::lower_bound(sorted.begin(),sorted.end(),lo,[&](std::size_t i,double z){return initial[i].position[2]<z;});
    auto last=std::upper_bound(first,sorted.end(),hi,[&](double z,std::size_t i){return z<initial[i].position[2];});
    aprl::ChargeConservingCurrentDepositor charge(fields,o);
    for(auto it=first;it!=last;++it) {
      ++count.referenceInspections;auto p=initial[*it];p.position[2]+=beta*DZ*step;
      charge.depositCharge(p.position,p.charge,referenceRho);
    }
    long double r2=0;
    // A translated finite initial cache has no prescribed extrapolation.
    // In the full baseline this diagnostic only uses its available overlap;
    // the actual full Maxwell/Gauss check still covers the full interior.
    const std::size_t firstK=static_cast<std::size_t>(std::max(static_cast<long double>(CPML_Z+1),
      CPML_Z+std::ceil(static_cast<long double>(beta)*step)+1.L-static_cast<long double>(offset)));
    referenceFirstK=offset+firstK;referenceLastK=offset+width-CPML_Z-1;
    for(std::size_t k=firstK;k<width-CPML_Z;++k)for(std::size_t j=1;j<NY;++j)for(std::size_t i=1;i<NX;++i) {
      const double source=referenceRho(i,j,k)/aprl::SI::epsilon0;
      const double v=(reference.value(0,i,j,offset+k,step)-reference.value(0,i-1,j,offset+k,step))/DX+
        (reference.value(1,i,j,offset+k,step)-reference.value(1,i,j-1,offset+k,step))/DY+
        (reference.value(2,i,j,offset+k,step)-reference.value(2,i,j,offset+k-1,step))/DZ-source;
      r2+=static_cast<long double>(v)*v;
    }
    return std::sqrt(r2)/gaussSourceNorm;
  };
  auto diagnostics=[&](unsigned step) {
    const auto begin=Clock::now();std::size_t held=0;
    const double innerLo=origin(offset)[2]+CPML_Z*DZ;
    const double innerHi=origin(offset)[2]+(width-CPML_Z)*DZ;
    const double releaseMargin=movingPlane(RELEASE,step)-innerLo;
    const double detectorMargin=innerHi-movingPlane(DETECTOR,step);
    count.minFrameReleaseInner=std::min(count.minFrameReleaseInner,releaseMargin);
    count.minFrameDetectorInner=std::min(count.minFrameDetectorInner,detectorMargin);
    require(releaseMargin>=2*DZ&&detectorMargin>=2*DZ,"Plane entered absorbing layer");
    auto& memory=dynamic_cast<aprl::EBConvolutionalPML&>(fields.boundary());
    for(unsigned h=0;h<12;++h)for(double value:memory.checkpointHistory(h)) {
      require(std::isfinite(value),"Nonfinite CPML history");
      count.cpmlHistoryPeak=std::max(count.cpmlHistoryPeak,std::abs(value));
      count.cpmlHistorySeen=count.cpmlHistorySeen||value!=0.;
    }
    for(const auto& state:active){held+=state.p.frozen?1:0;++count.diagnosticParticleInspections;}
    const auto g=gauss(fields,rho);
    const bool referenceCheck=step%8==0||step==STEPS;
    const double rg=referenceCheck?referenceGauss(step):0.;
    if(referenceCheck)count.maxReferenceGauss=std::max(count.maxReferenceGauss,rg);
    const double q=(movingPlane(DETECTOR,step)-ZORIGIN)/DZ;
    double central[3];for(unsigned a=0;a<3;++a)central[a]=sampleComponent(fields,offset,a,NX*.5+.5,NY*.5+.5,q);
    history<<step<<','<<step*dt<<','<<origin(offset)[2]<<','<<active.size()<<','<<held
      <<','<<count.queueAdmissions<<','<<count.evicted<<','<<interiorCharge(rho)
      <<','<<lastBalance<<','<<lastRemap<<','<<lastContinuity<<','<<lastPreGauss<<','<<g.relative
      <<','<<g.newLayerRelative<<','<<rg<<','<<(referenceCheck?1:0)
      <<','<<(referenceCheck?referenceFirstK:0)<<','<<(referenceCheck?referenceLastK:0)<<','<<lastCorrection
      <<','<<central[0]<<','<<central[1]<<','<<central[2]<<'\n';
    if(step%4==0||step==STEPS) {
      const double worldZ=movingPlane(DETECTOR,step);
      for(std::size_t j=2;j+2<NY;++j)for(std::size_t i=2;i+2<NX;++i) {
        plane<<step<<','<<i<<','<<j<<','<<step*dt<<','<<(step-.5)*dt
          <<','<<frame.labTFromBoxZT(worldZ,step*dt)<<','<<origin(offset)[0]+(i+.5)*DX
          <<','<<origin(offset)[1]+(j+.5)*DY<<','<<worldZ;
        for(unsigned a=0;a<6;++a)plane<<','<<sampleComponent(fields,offset,a,i+.5,j+.5,q);
        plane<<'\n';++count.fieldSamples;
      }
    }
    count.diagnostics+=elapsed(begin);
  };

  admit(0,full?0:windowOffset(1,guard));
  depositCharges(rhoExt,origin(offset-collar));cropRho(rhoExt,rho,collar);
  lastPreGauss=lastPostGauss=gauss(fields,rho).relative;
  count.activePeak=active.size();diagnostics(0);
  const auto evolutionStart=Clock::now();
  for(unsigned step=0;step<STEPS;++step) {
    const std::size_t nextOffset=full?0:windowOffset(step+1,guard);
    admit(step,nextOffset);count.activePeak=std::max(count.activePeak,active.size());
    const auto extOrigin=origin(offset-collar);const auto fieldOrigin=origin(offset);
    auto begin=Clock::now();depositCharges(rhoExt,extOrigin);cropRho(rhoExt,rho,collar);
    deposit.clearCurrent();nextExt.fill(0.);count.deposition+=elapsed(begin);
    aprl::EBElectricHaloView halo;
    if(boundary) {
      begin=Clock::now();boundary->beginStep(fields,step);halo=boundary->electricHalos(fields,step);
      count.boundary+=elapsed(begin);
    }
    begin=Clock::now();fields.advanceMagnetic(halo);count.fieldAdvance+=elapsed(begin);
    if(sectionDiagnostics)sectionDiagnostics->record(fields,offset,step+1,1);
    aprl::ChargeConservingCurrentDepositor depositor(deposit,extOrigin);
    for(auto& state:active) {
      ++count.particleInspections;auto& p=state.p;const auto& seed=initial[state.index];
      const auto particleStart=Clock::now();const auto before=p.position;
      Segment segments[2];unsigned used=0;double activeFraction=0.;
      if(p.frozen) {
        count.heldInCpmlObservations+=(before[2]>=fieldOrigin[2]&&before[2]<=fieldOrigin[2]+width*DZ)&&
          (before[2]<fieldOrigin[2]+CPML_Z*DZ||before[2]>fieldOrigin[2]+(width-CPML_Z)*DZ);
        FieldVector<double> end(before);end[2]=seed.position[2]+beta*DZ*(step+1);
        const double d0=before[2]-movingPlane(RELEASE,step),d1=end[2]-movingPlane(RELEASE,step+1);
        if(d1>=0.) {
          activeFraction=d0>=0.?0.:std::max(0.,std::min(1.,-d0/(d1-d0)));
          p.position=before;p.position[2]=before[2]+(end[2]-before[2])*activeFraction;
          if(activeFraction>0.)segments[used++]=Segment(before,p.position,0.,activeFraction,p.properVelocity,p.properVelocity);
          p.frozen=false;require(!released[state.index],"Particle activated twice");released[state.index]=true;
          require(p.properVelocity[0]==seed.properVelocity[0]&&p.properVelocity[1]==seed.properVelocity[1]&&
            p.properVelocity[2]==seed.properVelocity[2]&&p.position[0]==seed.position[0]&&p.position[1]==seed.position[1],
            "Prescribed release changed x/y or momentum");
          const double releaseClearance=p.position[2]-(fieldOrigin[2]+CPML_Z*DZ);
          count.minEventReleaseInner=std::min(count.minEventReleaseInner,releaseClearance);
          require(releaseClearance>=2*DZ,"Release event outside nonabsorbing PIC guard");
          particleRow(releases,p,seed,step,activeFraction,frame);
        } else {p.position=end;segments[used++]=Segment(before,end,0.,1.,p.properVelocity,p.properVelocity);}
      }
      if(!p.frozen&&activeFraction<1.) {
        const auto start=p.position,startU=p.properVelocity;const double remaining=(1-activeFraction)*dt;
        const bool sampleAllowed=full||(start[2]>=fieldOrigin[2]+(CPML_Z+2)*DZ&&
          start[2]<=fieldOrigin[2]+(width-CPML_Z-2)*DZ);
        if(prescribedMotion) {
          // The donor and explicit zero-response control use prescribed
          // sources with the original Maxwell/deposition kernels. Neither
          // samples a PIC force; ordinary candidate modes retain their push.
          p.position=seed.position;
          p.position[2]=seed.position[2]+beta*DZ*(step+1);
          p.properVelocity=seed.properVelocity;
        }
        else if(sampleAllowed) {
          const double clearance=std::min(start[2]-(fieldOrigin[2]+CPML_Z*DZ),
                                         fieldOrigin[2]+(width-CPML_Z)*DZ-start[2]);
          count.minPicInnerClearance=std::min(count.minPicInnerClearance,clearance);
          count.picPushesInCpml+=clearance<0.;
          require(clearance>=0.,"PIC push entered CPML");
          aprl::RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
            p,fields,fieldOrigin,sources,frame,(step+activeFraction)*dt,remaining,1);
        }
        else {
          require(detected[state.index],"Released particle left PIC guard before detector");
          aprl::RelativisticBorisPusher::pushPosition(p.position,p.properVelocity,remaining);
          ++count.ballisticOutgoingSteps;
          count.ballisticInCpmlObservations+=(start[2]>=fieldOrigin[2]&&start[2]<=fieldOrigin[2]+width*DZ)&&
            (start[2]<fieldOrigin[2]+CPML_Z*DZ||start[2]>fieldOrigin[2]+(width-CPML_Z)*DZ);
        }
        segments[used++]=Segment(start,p.position,activeFraction,1.,startU,p.properVelocity);
      }
      require(used>0&&used<=2,"Invalid held/active path split");
      if(prescribedMotion) {
        // Anchor every complete prescribed state to its initial source, including
        // a release exactly at the end of a step. Avoid accumulated drift.
        p.position=seed.position;
        p.position[2]=seed.position[2]+beta*DZ*(step+1);
        p.properVelocity=seed.properVelocity;
        segments[used-1].end=p.position;
        segments[used-1].endU=p.properVelocity;
      }
      for(unsigned a=0;a<3;++a)require(std::isfinite(p.position[a])&&std::isfinite(p.properVelocity[a]),"Non-finite particle state");
      require(p.position[0]>fieldOrigin[0]+DX&&p.position[0]<fieldOrigin[0]+(NX-1)*DX&&
        p.position[1]>fieldOrigin[1]+DY&&p.position[1]<fieldOrigin[1]+(NY-1)*DY,"Particle reached transverse guard");
      if(full)require(inside(p.position,origin(0),global),"Full reference particle left domain");
      for(unsigned s=0;s<used;++s)if(!detected[state.index]) {
        const auto& leg=segments[s];
        const double a=leg.start[2]-movingPlane(DETECTOR,step,leg.first);
        const double b=leg.end[2]-movingPlane(DETECTOR,step,leg.last);
        if(a<=0&&b>=0&&b>a) {
          const double f=-a/(b-a);RelativisticParticleSI crossing=p;
          for(unsigned axis=0;axis<3;++axis) {
            crossing.position[axis]=leg.start[axis]+f*(leg.end[axis]-leg.start[axis]);
            crossing.properVelocity[axis]=leg.startU[axis]+f*(leg.endU[axis]-leg.startU[axis]);
          }
          const double timeFraction=leg.first+f*(leg.last-leg.first);
          require(released[state.index]&&!p.frozen,"Detector crossed before release");
          const double detectorClearance=fieldOrigin[2]+(width-CPML_Z)*DZ-crossing.position[2];
          count.minEventDetectorInner=std::min(count.minEventDetectorInner,detectorClearance);
          require(detectorClearance>=2*DZ,"Detector event outside nonabsorbing PIC guard");
          particleRow(crossings,crossing,seed,step,timeFraction,frame);detected[state.index]=true;
        }
      }
      count.particles+=elapsed(particleStart);const auto depositStart=Clock::now();
      for(unsigned s=0;s<used;++s) {
        FieldVector<double> start(0.),end(0.);
        if(clip(segments[s],extOrigin,deposit.geometry(),start,end))
          count.depositedSegments+=depositor.depositSegment(start,end,p.charge).subsegments;
      }
      if(inside(p.position,extOrigin,deposit.geometry()))depositor.depositCharge(p.position,p.charge,nextExt);
      count.deposition+=elapsed(depositStart);
    }
    begin=Clock::now();cropCurrent(deposit,fields,collar);cropRho(nextExt,nextRho,collar);
    lastContinuity=continuity(fields,rho,nextRho);
    lastBalance=interiorCharge(nextRho)-interiorCharge(rho)+dt*outwardCurrent(fields);
    count.maxContinuity=std::max(count.maxContinuity,lastContinuity);
    count.maxChargeBalance=std::max(count.maxChargeBalance,static_cast<double>(std::abs(lastBalance)));
    count.deposition+=elapsed(begin);
    begin=Clock::now();fields.advanceElectric();count.fieldAdvance+=elapsed(begin);
    if(sectionDiagnostics)sectionDiagnostics->record(fields,offset,step+1,2);
    begin=Clock::now();lastPreGauss=gauss(fields,nextRho).relative;
    count.maxGaussBeforeShift=std::max(count.maxGaussBeforeShift,lastPreGauss);
    count.diagnostics+=elapsed(begin);
    const long double oldCharge=interiorCharge(nextRho);lastRemap=0;lastCorrection=0;
    if(nextOffset!=offset) {
      begin=Clock::now();require(offset==nextOffset+1,"Window moved more than one cell");
      struct IncomingContext {
        const StreamingReferenceCache* reference;std::size_t offset;unsigned step;bool zero;
        std::uint64_t* values;
      } context{feedReference?feedReference.get():&reference,nextOffset,step+1,zero,&count.newLayerValues};
      aprl::verification::MovingWindowIncoming incoming;
      incoming.context=&context;
      incoming.value=[](void* raw,unsigned a,std::size_t i,std::size_t j,std::size_t k,double)->double {
        auto& c=*static_cast<IncomingContext*>(raw);
        if(a>=6)return 0.; // Previous J is consumed; next deposition overwrites it.
        ++*c.values;return c.zero?0.:c.reference->value(a,i,j,c.offset+k,c.step);
      };
      const auto moved=mover->shift(fields,-1,origin(nextOffset)[2],incoming);
      count.shiftBytes+=moved.fieldMovedBytes;
      count.cpmlHistoryMovedBytes+=moved.historyMovedBytes;
      count.cpmlHistoryZeroedBytes+=moved.historyZeroedBytes;
      require(moved.historyRetainedBytes==0,"World-position policy unexpectedly retained relative histories");
      cropRho(nextExt,rho,collar-1);
      const long double added=chargePlane(nextRho,0),removed=chargePlane(nextRho,width-1);
      lastRemap=interiorCharge(rho)-oldCharge-added+removed;
      offset=nextOffset;
      if(refresh)boundary->refreshAfterRemap(fields,step+1);
      count.shift+=elapsed(begin);

    } else cropRho(nextExt,rho,collar);
    if(sectionDiagnostics)sectionDiagnostics->record(fields,offset,step+1,3);
    count.maxRemapBalance=std::max(count.maxRemapBalance,static_cast<double>(std::abs(lastRemap)));
    begin=Clock::now();lastPostGauss=gauss(fields,rho).relative;
    count.maxGaussAfterShift=std::max(count.maxGaussAfterShift,lastPostGauss);
    count.diagnostics+=elapsed(begin);
    if(!full) {
      const double high=origin(offset)[2]+(width+collar)*DZ;
      active.erase(std::remove_if(active.begin(),active.end(),[&](const ParticleState& state) {
        ++count.particleInspections;
        if(state.p.position[2]>high) {
          require(detected[state.index],"Undiagnosed particle left outgoing collar");++count.evicted;return true;
        }
        return false;
      }),active.end());
    }
    if(!full)for(const auto& state:active)if(!state.p.frozen&&!detected[state.index])
      require(state.p.position[2]>=origin(offset)[2]+(CPML_Z+2)*DZ&&
              state.p.position[2]<=origin(offset)[2]+(width-CPML_Z-2)*DZ,
              "Window shift swept an undetected PIC particle outside its interior guard");
    diagnostics(step+1);
    if(oracleWriter)oracleWriter->record(fields,step+1,windowOffset(step+1,guard)-2);
  }
  const double evolutionSeconds=elapsed(evolutionStart);
  if(sectionDiagnostics)sectionDiagnostics->finish();
  if(oracleWriter)oracleWriter->finish();
  require(std::count(released.begin(),released.end(),true)==static_cast<long>(PARTICLES),"Not all fixed-cohort particles released");
  require(std::count(detected.begin(),detected.end(),true)==static_cast<long>(PARTICLES),"Not all fixed-cohort particles crossed detector");
  require(queue==0,"Input queue not exhausted");
  require(history.good()&&crossings.good()&&releases.good()&&initialCsv.good()&&plane.good(),"Scientific output write failed");
  const std::size_t rhoBytes=(rhoExt.size()+nextExt.size()+rho.size()+nextRho.size()+referenceRho.size())*sizeof(double);
  const std::size_t particleTableBytes=initial.capacity()*sizeof(RelativisticParticleSI);
  const std::size_t activeBytes=active.capacity()*sizeof(ParticleState);
  const std::size_t queueBytes=sorted.capacity()*sizeof(std::size_t);
  const std::size_t eventFlagBytes=(released.capacity()+7)/8+(detected.capacity()+7)/8;
  const std::size_t liveBytes=fields.memoryFootprint().totalBytes()+deposit.memoryFootprint().totalBytes()+
    rhoBytes+reference.memoryBytes()+particleTableBytes+activeBytes+queueBytes+eventFlagBytes+
    (mover?mover->memoryBytes():0)+(feedReference?feedReference->memoryBytes():0)+
    (cpmlReference?cpmlReference->memoryBytes():0)+(oracleReader?oracleReader->memoryBytes():0);
  auto metrics=fresh(directory+"/metrics.json");
  metrics<<"{\n\"format\":\"aprl-streaming-cpml-bunch-probe-v1\",\"completed\":true,\"mode\":\""<<mode<<"\",\"beta_box\":"<<beta
    <<",\"guard_cells\":"<<guard<<",\"steps\":"<<STEPS<<",\"dt_s\":"<<dt<<",\"nx\":"<<NX<<",\"ny\":"<<NY
    <<",\"history_policy\":\"world-position\",\"cpml_cells\":[0,0,20],\"cpml_z_cells\":20"
    <<",\"cpml_input\":\""<<(residual?"residual":"total")<<"\""
    <<",\"residual_reference_work_bytes\":"<<(boundary?boundary->residualWorkBytes():0)
    <<",\"residual_reference_input_values\":"<<(boundary?boundary->residualInputValues():0)
    <<",\"residual_reference_input_wall_s\":"<<(boundary?boundary->residualInputSeconds():0)
    <<",\"physical_full_nz\":"<<PHYSICAL_NZ<<",\"full_origin_z_m\":"<<ZORIGIN
    <<",\"physical_window_nz\":"<<width-2*CPML_Z<<",\"inner_guard_cells\":"<<guard
    <<",\"cpml_history_peak_abs\":"<<count.cpmlHistoryPeak
    <<",\"cpml_history_nonzero_seen\":"<<(count.cpmlHistorySeen?"true":"false")
    <<",\"cpml_history_moved_bytes\":"<<count.cpmlHistoryMovedBytes
    <<",\"cpml_history_zeroed_bytes\":"<<count.cpmlHistoryZeroedBytes
    <<",\"remap_helper_bytes\":"<<(mover?mover->memoryBytes():0)
    <<",\"pic_pushes_in_cpml\":"<<count.picPushesInCpml
    <<",\"held_in_cpml_observations\":"<<count.heldInCpmlObservations
    <<",\"ballistic_in_cpml_observations\":"<<count.ballisticInCpmlObservations
    <<",\"min_frame_release_inner_m\":"<<count.minFrameReleaseInner
    <<",\"min_frame_detector_inner_m\":"<<count.minFrameDetectorInner
    <<",\"min_event_release_inner_m\":"<<count.minEventReleaseInner
    <<",\"min_event_detector_inner_m\":"<<count.minEventDetectorInner
    <<",\"min_pic_inner_clearance_m\":"<<count.minPicInnerClearance
    <<",\"full_nz\":"<<NZ<<",\"window_nz\":"<<width<<",\"particles\":"<<PARTICLES<<",\"seed\":"<<SEED
    <<",\"charge_C\":"<<TOTAL_CHARGE<<",\"release_plane_lab_m\":"<<frame.gamma()*RELEASE
    <<",\"detector_plane_lab_m\":"<<frame.gamma()*DETECTOR<<",\"window_beta\":"<<WINDOW_BETA
    <<",\"initial_gauss_relative\":"<<init.relativeResidual<<",\"initial_iterations\":"<<init.iterations
    <<",\"initial_temporary_bytes_per_rank\":"<<init.temporaryBytes
    <<",\"full_grid_bytes\":"<<fullBytes<<",\"evolving_grid_bytes\":"<<fields.memoryFootprint().totalBytes()
    <<",\"deposition_collar_grid_bytes\":"<<deposit.memoryFootprint().totalBytes()
    <<",\"rho_arrays_bytes\":"<<rhoBytes
    <<",\"reference_cache_bytes\":"<<reference.memoryBytes()<<",\"initial_particle_table_bytes\":"<<particleTableBytes
    <<",\"active_reserved_particle_bytes\":"<<activeBytes<<",\"sorted_queue_bytes\":"<<queueBytes
    <<",\"event_flag_bytes_approx\":"<<eventFlagBytes<<",\"live_storage_bytes_approx\":"<<liveBytes
    <<",\"active_peak\":"<<count.activePeak<<",\"admitted\":"<<count.queueAdmissions<<",\"evicted\":"<<count.evicted
    <<",\"particle_inspections\":"<<count.particleInspections<<",\"reference_particle_inspections\":"<<count.referenceInspections
    <<",\"diagnostic_particle_inspections\":"<<count.diagnosticParticleInspections
    <<",\"deposited_subsegments\":"<<count.depositedSegments<<",\"ballistic_outgoing_steps\":"<<count.ballisticOutgoingSteps
    <<",\"reference_sample_calls\":"<<reference.sampleCount()<<",\"boundary_values\":"<<(boundary?boundary->boundaryValues():0)
    <<",\"oracle_recording\":"<<(recording?"true":"false")
    <<",\"prescribed_reference_donor\":"<<(prescribedReference?"true":"false")
    <<",\"prescribed_actual_control\":"<<(prescribedActualControl?"true":"false")
    <<",\"prescribed_reference_input\":"<<(heldReference?"true":"false")
    <<",\"reference_provenance_kind\":\""
    <<((recording||oracleReplay)?aprl::verification::native_oracle_detail::sourceKindName(oracleConfig):
       "frozen-initial-field")<<"\""
    <<",\"native_B_reference_reconstruction\":"<<(nativeB?"true":"false")
    <<",\"oracle_replay\":"<<(oracleReplay?"true":"false")
    <<",\"oracle_feed_channel\":"<<(oracleFeed?"true":"false")
    <<",\"oracle_cpml_channel\":"<<(oracleCpml?"true":"false")
    <<",\"oracle_directory\":"<<aprl::verification::native_oracle_detail::jsonString(oracleDirectory)
    <<",\"oracle_recorded_frames\":"<<(oracleWriter?oracleWriter->frames():0)
    <<",\"oracle_recorded_frame_bytes\":"<<(oracleWriter?oracleWriter->bytes():0)
    <<",\"oracle_reader_queries\":"<<(oracleReader?oracleReader->queries():0)
    <<",\"oracle_reader_frame_loads\":"<<(oracleReader?oracleReader->loads():0)
    <<",\"oracle_reader_bytes\":"<<(oracleReader?oracleReader->bytesRead():0)
    <<",\"oracle_reader_live_bytes\":"<<(oracleReader?oracleReader->memoryBytes():0)
    <<",\"oracle_feed_sample_calls\":"<<(feedReference?feedReference->sampleCount():0)
    <<",\"oracle_cpml_sample_calls\":"<<(cpmlReference?cpmlReference->sampleCount():0)
    <<",\"oracle_extra_cache_bytes\":"<<((feedReference?feedReference->memoryBytes():0)+(cpmlReference?cpmlReference->memoryBytes():0))
    <<",\"oracle_scope\":\""
    <<((prescribedReference||heldReference)?
       "independent prescribed source from initial conditions only; full-grid donor evolution establishes feasibility, not production speed":
       "diagnostic replay requires full-domain physical history; no production applicability or speed claim")<<"\""
    <<",\"halo_values\":"<<(boundary?boundary->haloValues():0)<<",\"shift_bytes\":"<<count.shiftBytes
    <<",\"new_layer_values\":"<<count.newLayerValues<<",\"gauss_correction_values\":"<<count.gaussCorrectionValues
    <<",\"field_plane_samples\":"<<count.fieldSamples<<",\"field_cell_updates\":"<<static_cast<std::uint64_t>(NX)*NY*width*STEPS
    <<",\"max_continuity_relative\":"<<count.maxContinuity<<",\"max_open_charge_balance_C\":"<<count.maxChargeBalance
    <<",\"max_remap_charge_balance_C\":"<<count.maxRemapBalance
    <<",\"max_total_gauss_before_shift\":"<<count.maxGaussBeforeShift<<",\"max_total_gauss_after_shift\":"<<count.maxGaussAfterShift
    <<",\"max_reference_gauss_relative\":"<<count.maxReferenceGauss<<",\"max_new_Ez_correction_V_m\":"<<count.maxEzCorrection
    <<",\"reference_gauss_diagnostic_source\":\"original frozen cache retained for historical comparison; for evolving reference use its donor total_gauss diagnostics\""
    <<",\"gauss_initial_source_l2_V_m2\":"<<gaussSourceNorm<<",\"gauss_initial_source_max_V_m2\":"<<gaussSourceMaximum
    <<",\"new_Ez_correction_l2\":"<<std::sqrt(count.sumEzCorrection2)
    <<",\"initialization_wall_s\":"<<count.initialization<<",\"cache_preparation_wall_s\":"<<count.cache
    <<",\"evolution_wall_s\":"<<evolutionSeconds<<",\"field_advance_wall_s\":"<<count.fieldAdvance
    <<",\"particle_wall_s\":"<<count.particles<<",\"deposition_and_continuity_wall_s\":"<<count.deposition
    <<",\"boundary_preparation_wall_s\":"<<count.boundary<<",\"shift_wall_s\":"<<count.shift
    <<",\"diagnostics_wall_s\":"<<count.diagnostics<<",\"new_layer_gauss_correction_wall_s\":"<<count.gaussCorrection
    <<",\"overall_wall_s\":"<<elapsed(overallStart)
    <<",\"physical_acceptance\":\"measured externally; no approximation threshold imposed by this probe\""
    <<",\"field_time_contract\":\"field-plane.csv contains native box E(t_n) and B(t_n-dt/2), sampled at the fixed-lab plane position at t_n; no mixed-time lab field transform\""
    <<",\"charge_contract\":\"interior vertex control volume, open edge-current flux plus explicit moved-in/moved-out vertex plane charge; outer collar endpoints excluded\""
    <<",\"gauss_normalization\":\"Ordinary Gauss measured only on nonabsorbing interior vertices; L2 divided by fixed initial full-cohort rho/epsilon0 L2; layer maximum is first nonabsorbing interior layer, not newly inserted CPML plane; translated-reference diagnostic also uses cache-valid overlap\""
    <<",\"storage_accounting\":\"known live numerical arrays and vector capacities; vector<bool> flags estimated from bit capacity; excludes allocator/object/stream/MPI overhead; initialization CG temporary bytes and full-grid construction peak reported separately\""
    <<",\"source_contract\":\""
    <<(prescribedActualControl?
       "zero-response diagnostic control: actual short-window source follows prescribed seed positions/momenta exactly, with no PIC force; not a physical particle-evolution result":prescribedReference?
       "independent donor: initial x/y/u fixed and z=seed_z+beta*c*n*dt at every complete step; release/detector events logged but source never receives PIC force; original Cowan and charge-conserving deposition":
       "generated Gaussian input; held x/y/u fixed until common fixed-lab release; full and window use identical prescribed release paths; short-window outgoing particles become ballistic only after detector")<<"\""
    <<",\"particle_detector_contract\":\"position and proper velocity are linearly interpolated on the crossed held/active segment for diagnostics; not an exact intermediate Boris trajectory solution\""
    <<",\"limitations\":[\"One MPI rank; transverse PEC retained to isolate longitudinal CPML; no production checkpoint integration\",\""
    <<((prescribedReference||heldReference)?
       "DIAGNOSTIC ONLY: reference donor uses prescribed initial-source trajectories without physical PIC future data, but still performs a full-grid field evolution and stores it":
       oracleReplay?"DIAGNOSTIC ONLY: selected reference channels replay physical full-domain future native states; this is not an independent window algorithm":
       "Frozen reference is one initial full-domain Cowan solve plus linear native-field translation, not evolving full-domain history")
    <<"\",\"Full initial solve/cache and input table remain O(total bunch length/particle count), explicitly included in cost\",\"Total-field Cowan and particle update; CPML correction input selected by cpml_input; reference outer face/halo/new-layer feed; outgoing exterior evolution after removal is omitted\",\"No Gauss correction or reference-source mismatch correction is applied; post-detector ballistic full-charge paths cross CPML before removal\",\"Ordinary field energy is not the exact Cowan invariant; no pure radiation or full-scale speed claim\"]}\n";
  require(metrics.good(),"Metrics output write failed");
  std::cout<<"Streaming CPML bunch probe complete: "<<mode<<" beta="<<beta<<" guard="<<guard
    <<" crossings="<<PARTICLES<<" steps="<<STEPS<<" max_continuity="<<count.maxContinuity<<'\n';
}
}

int main(int argc,char** argv) {
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int rank=0,size=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  try {
    require(size==1,"Streaming bunch research probe requires one MPI rank");
    require(argc>=5&&argc<=7,"Usage: aprl_streaming_cpml_probe OUT beta inner_guard MODE [diagnostic_max_step [oracle_directory]]");
    const unsigned long diagnosticSteps=argc>=6?std::stoul(argv[5]):0;
    require(diagnosticSteps<=STEPS,"Diagnostic step limit exceeds this bounded fixture");
    run(argv[1],std::stod(argv[2]),static_cast<unsigned>(std::stoul(argv[3])),argv[4],
        static_cast<unsigned>(diagnosticSteps),argc==7?argv[6]:"");
    MPI_Finalize();return 0;
  } catch(const std::exception& error) {
    std::cerr<<"rank "<<rank<<": "<<error.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;
  }
}
