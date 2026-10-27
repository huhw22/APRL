// Integer window migration only. All physical updates call unchanged core code.
#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <mpi.h>
#include "moving_window_remap.h"
#include "moving_window_wave.h"
#include "eb_initial_field.h"
#include "eb_deposition.h"
#include "eb_particles.h"
#include "particle_initializer.h"
#include "config.h"

namespace {
using namespace aprl;
using namespace aprl::verification;
using Clock=std::chrono::steady_clock;
std::uint64_t checks=0;
void require(bool ok,const std::string& message) {
  ++checks;if(!ok)throw std::runtime_error(message);
}
double elapsed(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
bool sameDouble(const double& a,const double& b){return std::memcmp(&a,&b,sizeof(double))==0;}
YeeComponent& field(EBFieldGrid& f,unsigned a) {
  switch(a){case 0:return f.ex();case 1:return f.ey();case 2:return f.ez();
    case 3:return f.bx();case 4:return f.by();case 5:return f.bz();
    case 6:return f.jx();case 7:return f.jy();case 8:return f.jz();}
  throw std::runtime_error("Bad component");
}
EBConvolutionalPML& cpml(EBFieldGrid& f){return dynamic_cast<EBConvolutionalPML&>(f.boundary());}
std::ofstream fresh(const std::string& path) {
  struct stat st;require(stat(path.c_str(),&st)!=0&&errno==ENOENT,"Refusing existing output "+path);
  std::ofstream f(path);require(bool(f),"Cannot open "+path);f<<std::setprecision(17);return f;
}
EBGridGeometry geometry(std::size_t nx,std::size_t ny,std::size_t nz) {
  return EBGridGeometry(nx,ny,nz,2e-6,2e-6,1e-6,1e-6/SI::c,EBMaxwellSolver::CowanZ);
}
std::unique_ptr<EBFieldGrid> makeGrid(const EBGridGeometry& g,const EBCPMLParameters& p) {
  std::unique_ptr<EBFieldGrid> f(new EBFieldGrid(g));
  f->setBoundary(std::unique_ptr<EBBoundaryOperator>(new EBConvolutionalPML(g,g.nz,0,p,true,true)));
  f->clearFields();f->clearCurrent();return f;
}
void equal(EBFieldGrid& a,EBFieldGrid& b,const std::string& stage) {
  for(unsigned c=0;c<9;++c) {
    auto& x=field(a,c);auto& y=field(b,c);
    require(x.size()==y.size(),stage+" field size");
    require(std::memcmp(x.data(),y.data(),x.size()*sizeof(Double))==0,stage+" field "+std::to_string(c));
  }
  for(unsigned h=0;h<12;++h) {
    auto& x=cpml(a).checkpointHistory(h);auto& y=cpml(b).checkpointHistory(h);
    require(x.size()==y.size(),stage+" history size");
    require(x.empty()||std::memcmp(x.data(),y.data(),x.size()*sizeof(Double))==0,stage+" history "+std::to_string(h));
  }
}
void clone(EBFieldGrid& a,EBFieldGrid& b) {
  for(unsigned c=0;c<9;++c)std::copy(field(a,c).data(),field(a,c).data()+field(a,c).size(),field(b,c).data());
  for(unsigned h=0;h<12;++h) {
    const auto& v=cpml(a).checkpointHistory(h);std::copy(v.begin(),v.end(),cpml(b).checkpointHistory(h).begin());
  }
}
void finite(EBFieldGrid& f) {
  for(unsigned c=0;c<9;++c)for(std::size_t i=0;i<field(f,c).size();++i)
    require(std::isfinite(field(f,c).data()[i]),"Nonfinite field");
  for(unsigned h=0;h<12;++h)for(double x:cpml(f).checkpointHistory(h))require(std::isfinite(x),"Nonfinite history");
}

// Independent, deliberately straightforward coordinate oracle. Its buffers and
// coordinate metadata are prepared once; performance does not include per-step malloc.
struct Oracle {
  struct History {
    unsigned axis;std::size_t n[3];std::vector<long> ordinal;
    std::vector<std::array<std::size_t,3>> coordinate;std::vector<double> old;
    long index(std::size_t i,std::size_t j,long k) const {
      if(k<0||k>=static_cast<long>(n[2]))return -1;
      const std::size_t xyz[3]={i,j,static_cast<std::size_t>(k)};
      const long q=ordinal[xyz[axis]];if(q<0)return -1;
      if(axis==0)return q*n[1]*n[2]+k*n[1]+j;
      if(axis==1)return q*n[0]*n[2]+k*n[0]+i;
      return q*n[0]*n[1]+j*n[0]+i;
    }
  };
  std::array<std::vector<double>,9> oldField;
  std::array<History,12> histories;
  std::size_t memoryBytes() const {
    std::size_t n=sizeof(*this);
    for(const auto& v:oldField)n+=v.capacity()*sizeof(double);
    for(const auto& h:histories)n+=h.old.capacity()*sizeof(double)+h.ordinal.capacity()*sizeof(long)+h.coordinate.capacity()*sizeof(std::array<std::size_t,3>);
    return n;
  }
  Oracle(EBFieldGrid& f,const EBCPMLParameters& p) {
    for(unsigned a=0;a<9;++a)oldField[a].resize(field(f,a).size());
    const unsigned axes[6]={1,2,2,0,0,1},comps[6]={0,0,1,1,2,2};
    const std::size_t global[3]={f.geometry().nx,f.geometry().ny,f.geometry().nz};
    for(unsigned h=0;h<12;++h) {
      auto& out=histories[h];out.axis=axes[h%6];const bool electric=h>=6;
      auto& v=field(f,comps[h%6]+(electric?0:3));
      out.n[0]=v.nx();out.n[1]=v.ny();out.n[2]=v.nz();out.ordinal.assign(out.n[out.axis],-1);
      long q=0;
      for(std::size_t index=0;index<out.n[out.axis];++index) {
        const double position=index+(electric?0.:.5);
        const bool active=p.cells[out.axis]&&(position<p.cells[out.axis]||position>global[out.axis]-p.cells[out.axis]);
        if(!active)continue;
        out.ordinal[index]=q++;
      }
      out.old.resize(cpml(f).checkpointHistory(h).size());out.coordinate.resize(out.old.size());
      std::vector<bool> seen(out.old.size(),false);
      for(std::size_t k=0;k<out.n[2];++k)for(std::size_t j=0;j<out.n[1];++j)for(std::size_t i=0;i<out.n[0];++i) {
        const long address=out.index(i,j,k);if(address<0)continue;
        require(static_cast<std::size_t>(address)<out.old.size()&&!seen[address],"Oracle history layout collision");
        out.coordinate[address]={{i,j,k}};seen[address]=true;
      }
      require(std::find(seen.begin(),seen.end(),false)==seen.end(),"Oracle history layout hole");
    }
  }
  void shift(EBFieldGrid& f,int s,double origin,MovingWindowIncoming incoming={},
             MovingWindowCurrentPolicy policy=MovingWindowCurrentPolicy::Preserve) {
    if(!s)return;
    for(unsigned a=0;a<9;++a) {
      auto& v=field(f,a);auto& old=oldField[a];std::copy(v.data(),v.data()+v.size(),old.begin());
      for(std::size_t k=0;k<v.nz();++k)for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i) {
        const long from=static_cast<long>(k)+s;
        if(a>=6&&policy==MovingWindowCurrentPolicy::Reset)v(i,j,k)=0.;
        else if(from>=0&&from<static_cast<long>(v.nz()))v(i,j,k)=old[(from*v.ny()+j)*v.nx()+i];
        else v(i,j,k)=incoming.value?incoming.value(incoming.context,a,i,j,k,origin):0.;
      }
    }
    for(unsigned h=0;h<12;++h) {
      auto& meta=histories[h];auto& v=cpml(f).checkpointHistory(h);std::copy(v.begin(),v.end(),meta.old.begin());
      for(std::size_t n=0;n<v.size();++n) {
        const auto xyz=meta.coordinate[n];const long from=meta.index(xyz[0],xyz[1],static_cast<long>(xyz[2])+s);
        v[n]=from<0?0.:meta.old[from];
      }
    }
  }
};

double incoming(void*,unsigned a,std::size_t i,std::size_t j,std::size_t k,double origin) {
  const double value=1.+a+.03125*i+.015625*j+.0001*(origin/1e-6+k);
  return a<3?value:(a<6?value/SI::c:value*SI::epsilon0*SI::c/1e-6);
}
MovingWindowIncoming provider(){MovingWindowIncoming p;p.value=incoming;return p;}
void seedTags(EBFieldGrid& f) {
  for(unsigned a=0;a<9;++a)for(std::size_t n=0;n<field(f,a).size();++n)field(f,a).data()[n]=(a+1)*1048576.+n*.125;
  for(unsigned h=0;h<12;++h) {
    auto& v=cpml(f).checkpointHistory(h);for(std::size_t n=0;n<v.size();++n)v[n]=(h+11)*1048576.+n*.125;
  }
}
void seedWork(EBFieldGrid& f) {
  for(unsigned a=0;a<6;++a) {
    auto& v=field(f,a);
    for(std::size_t k=0;k<v.nz();++k)for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i)
      v(i,j,k)=std::sin(.09*i)*std::sin(.11*j)*std::cos(.073*k+.4*a)/(a<3?1.:SI::c);
  }
}
std::vector<const double*> addresses(EBFieldGrid& f) {
  std::vector<const double*> p;for(unsigned a=0;a<9;++a)p.push_back(field(f,a).data());
  for(unsigned h=0;h<12;++h)p.push_back(cpml(f).checkpointHistory(h).data());
  return p;
}
double maxDifference(EBFieldGrid& a,EBFieldGrid& b) {
  double out=0;
  for(unsigned c=0;c<6;++c)for(std::size_t n=0;n<field(a,c).size();++n)
    out=std::max(out,std::abs(field(a,c).data()[n]-field(b,c).data()[n])*(c<3?1.:SI::c));
  return out;
}

void verify(const std::string& root) {
  auto log=fresh(root+"/checks.csv");log<<"test,checks_after\n";
  unsigned cases=0;
  for(unsigned axes=1;axes<8;++axes) {
    auto g=geometry(11,10,19);EBCPMLParameters p;
    p.cells[0]=(axes&1)?2:0;p.cells[1]=(axes&2)?3:0;p.cells[2]=(axes&4)?4:0;
    const std::vector<int> shifts={-19,-18,-3,-1,0,1,3,18,19};
    MovingWindowRemap mover(g,p,shifts);
    for(int s:shifts)for(auto policy:{MovingWindowCurrentPolicy::Preserve,MovingWindowCurrentPolicy::Reset}) {
      auto a=makeGrid(g,p),b=makeGrid(g,p);seedTags(*a);clone(*a,*b);Oracle oracle(*b,p);
      const auto before=addresses(*a);const double origin=s*g.dz;
      const auto r=mover.shift(*a,s,origin,provider(),policy);oracle.shift(*b,s,origin,provider(),policy);
      equal(*a,*b,"tag shift");require(addresses(*a)==before,"Migration reallocated storage");
      require(r.signedCells==s,"Incorrect shift report");++cases;
    }
    auto a=makeGrid(g,p),b=makeGrid(g,p);seedTags(*a);clone(*a,*b);
    bool rejected=false;try{mover.shift(*a,2,0.);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"Unprepared shift accepted");equal(*a,*b,"Rejected shift unchanged");
    rejected=false;try{MovingWindowRemap invalid(g,p,{20});}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"Too-large planned shift accepted");
  }
  log<<"tagged_native_fields_and_all_history_supports,"<<checks<<'\n';
  const auto g=geometry(24,24,96);EBCPMLParameters p;p.cells[0]=p.cells[1]=6;p.cells[2]=8;
  // Independent equivalence over actual core updates; a missing-history negative
  // control must differ immediately and contaminate the very next core step.
  auto a=makeGrid(g,p),b=makeGrid(g,p),bad=makeGrid(g,p);seedWork(*a);a->advance();a->advance();
  clone(*a,*b);clone(*a,*bad);MovingWindowRemap mover(g,p);Oracle oracle(*b,p);
  mover.shift(*a,-1,-g.dz,provider());oracle.shift(*b,-1,-g.dz,provider());equal(*a,*b,"Warm nonzero CPML");
  mover.shift(*bad,-1,-g.dz,provider());
  double historyMax=0;for(unsigned h=0;h<12;++h)for(double x:cpml(*bad).checkpointHistory(h))historyMax=std::max(historyMax,std::abs(x));
  require(historyMax>0,"CPML negative control never excited history");
  for(unsigned h=0;h<12;++h)std::fill(cpml(*bad).checkpointHistory(h).begin(),cpml(*bad).checkpointHistory(h).end(),0.);
  a->advance();bad->advance();require(maxDifference(*a,*bad)>0,"Clearing CPML history was not detected");
  log<<"history_omission_negative_control,"<<checks<<'\n';

  SimulationConfig config;config.beam.type=BeamInputType::GeneratedGaussian;
  config.beam.macroparticles=256;config.beam.electrons=.1e-12/SI::elementaryCharge;
  config.beam.randomSeed=710271;config.beam.gamma=1./std::sqrt(1.-.23*.23);
  config.beam.direction=FieldVector<double>(0.);config.beam.direction[2]=1.;
  config.beam.center=FieldVector<double>(0.);config.beam.center[2]=48e-6;
  config.beam.positionOffset=FieldVector<double>(0.);config.beam.sigmaPosition=FieldVector<double>(2e-6);
  config.beam.sigmaProperVelocity=FieldVector<double>(0.);
  ParticleInitializationReport particleReport;auto initial=ParticleInitializer::create(config,MPI_COMM_WORLD,particleReport);
  for(auto& particle:initial){particle.id=particle.sourceId;particle.frozen=false;}
  FieldVector<double> baseOrigin(0.);baseOrigin[0]=baseOrigin[1]=-24e-6;
  auto seed=makeGrid(g,p);const std::size_t guard[3]={6,6,8};
  EBGaussFieldInitializer::initialize(*seed,g,0,baseOrigin,initial,InitialSelfFieldModel::CowanRelativisticPoisson,
                                    1e-8,guard,1e-10,3000,MPI_COMM_WORLD);
  for(unsigned schedule=0;schedule<4;++schedule) {
    a=makeGrid(g,p);b=makeGrid(g,p);clone(*seed,*a);clone(*seed,*b);Oracle reference(*b,p);
    auto pa=initial,pb=initial;long offset=0;FieldVector<double> origin=baseOrigin;
    const unsigned steps=schedule==3?64:24;
    for(unsigned step=0;step<steps;++step) {
      a->clearCurrent();b->clearCurrent();a->advanceMagnetic();b->advanceMagnetic();equal(*a,*b,"Magnetic half step");
      ChargeConservingCurrentDepositor da(*a,origin),db(*b,origin);
      for(std::size_t i=0;i<pa.size();++i) {
        const auto oldA=pa[i].position,oldB=pb[i].position;
        RelativisticBorisPusher::pushFromGrid(pa[i],*a,origin,g.dt);
        RelativisticBorisPusher::pushFromGrid(pb[i],*b,origin,g.dt);
        da.depositSegment(oldA,pa[i].position,pa[i].charge);db.depositSegment(oldB,pb[i].position,pb[i].charge);
        for(unsigned d=0;d<3;++d)require(sameDouble(pa[i].position[d],pb[i].position[d])&&sameDouble(pa[i].properVelocity[d],pb[i].properVelocity[d]),"Particle mismatch");
      }
      a->advanceElectric();b->advanceElectric();equal(*a,*b,"Electric complete step");
      int s=0;
      if(schedule==1)s=static_cast<int>(std::floor(-.8*(step+1)))-static_cast<int>(std::floor(-.8*step));
      if(schedule==2)s=static_cast<int>(std::floor(.8*(step+1)))-static_cast<int>(std::floor(.8*step));
      if(schedule==3){const int pattern[8]={-3,0,1,3,0,-1,-1,1};s=pattern[step%8];}
      offset+=s;origin[2]=baseOrigin[2]+offset*g.dz;
      mover.shift(*a,s,origin[2]);reference.shift(*b,s,origin[2]);equal(*a,*b,"Migrated complete step");
      for(std::size_t i=0;i<pa.size();++i)for(unsigned d=0;d<3;++d)
        require(sameDouble(pa[i].position[d],pb[i].position[d]),"World particle moved with array");
    }
    finite(*a);log<<"PIC_schedule_"<<schedule<<','<<checks<<'\n';
  }
  auto m=fresh(root+"/metrics.json");m<<"{\"completed\":true,\"checks\":"<<checks<<",\"bitwise_failures\":0,\"tagged_shift_cases\":"<<cases
    <<",\"particle_schedules\":4,\"particles\":256,\"status\":\"state migration and unchanged-core equivalence only\"}\n";
}

void benchmark(const std::string& root,unsigned width) {
  require(width==64||width==128||width==256,"Unsupported bounded benchmark width");
  const auto g=geometry(40,40,width);EBCPMLParameters p;p.cells[0]=p.cells[1]=6;p.cells[2]=8;
  const unsigned steps=64,repeats=3;
  {
    auto a=makeGrid(g,p),b=makeGrid(g,p);seedWork(*a);clone(*a,*b);MovingWindowRemap mover(g,p,{-1});Oracle oracle(*b,p);
    long offset=0;
    for(unsigned step=0;step<8;++step) {
      a->clearCurrent();b->clearCurrent();a->advance();b->advance();
      const long next=static_cast<long>(std::floor(-.8*(step+1)));const int s=next-offset;offset=next;
      mover.shift(*a,s,offset*g.dz,provider());oracle.shift(*b,s,offset*g.dz,provider());equal(*a,*b,"Benchmark preflight");
    }
  }
  auto csv=fresh(root+"/timing.csv");
  csv<<"repeat,mode,width,steps,kernel_s,move_s,field_move_s,history_move_s,incoming_s,total_s,setup_s,aux_bytes,field_bytes,cpml_bytes,moved_payload_bytes,incoming_bytes,history_zero_bytes,cell_updates\n";
  double largest=0;
  for(unsigned repeat=0;repeat<repeats;++repeat)for(unsigned sequence=0;sequence<3;++sequence) {
    const unsigned mode=(sequence+repeat)%3; // Rotate ordering: fixed, compact, snapshot.
    auto f=makeGrid(g,p);seedWork(*f);
    std::unique_ptr<MovingWindowRemap> fast;std::unique_ptr<Oracle> slow;
    const auto setupStart=Clock::now();
    if(mode==1)fast.reset(new MovingWindowRemap(g,p,{-1}));
    if(mode==2)slow.reset(new Oracle(*f,p));
    const double setupSeconds=mode?elapsed(setupStart):0.;
    const std::size_t auxBytes=fast?fast->memoryBytes():(slow?slow->memoryBytes():0);
    for(unsigned warm=0;warm<4;++warm){f->clearCurrent();f->advance();}
    double kernel=0,move=0,fieldTime=0,historyTime=0,fillTime=0;
    std::uint64_t bytes=0,incomingBytes=0,zeroBytes=0;long cell=0;
    const auto start=Clock::now();
    for(unsigned step=0;step<steps;++step) {
      auto t=Clock::now();f->clearCurrent();f->advance();kernel+=elapsed(t);
      t=Clock::now();
      if(mode) {
        const long next=static_cast<long>(std::floor(-.8*(step+1)));const int s=next-cell;cell=next;
        if(mode==1) {
          const auto r=fast->shift(*f,s,cell*g.dz,provider());
          fieldTime+=r.fieldWallSeconds;historyTime+=r.historyWallSeconds;fillTime+=r.incomingWallSeconds;
          bytes+=r.fieldMovedBytes+r.historyMovedBytes;incomingBytes+=r.fieldIncomingBytes;zeroBytes+=r.historyZeroedBytes;
        } else slow->shift(*f,s,cell*g.dz,provider());
      }
      move+=elapsed(t);
    }
    const double total=elapsed(start);finite(*f);largest=std::max(largest,total);
    csv<<repeat<<','<<(mode==0?"fixed":mode==1?"compact":"snapshot")<<','<<width<<','<<steps<<','<<kernel<<','<<move
       <<','<<fieldTime<<','<<historyTime<<','<<fillTime<<','<<total<<','<<setupSeconds<<','<<auxBytes
       <<','<<f->memoryFootprint().electricBytes+f->memoryFootprint().magneticBytes+f->memoryFootprint().currentBytes
       <<','<<cpml(*f).memoryBytes()<<','<<bytes<<','<<incomingBytes<<','<<zeroBytes<<','<<40ULL*40*width*steps<<'\n';
  }
  auto m=fresh(root+"/metrics.json");m<<"{\"completed\":true,\"numerical_status\":\"PASS\",\"width\":"<<width
    <<",\"steps\":"<<steps<<",\"repeats\":"<<repeats<<",\"bitwise_preflight_steps\":8,\"largest_timed_run_s\":"<<largest
    <<",\"scope\":\"synthetic kernel and migration cost; no production speedup or physical-domain acceptance\"}\n";
}
}
int main(int argc,char** argv) {
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int size=0;MPI_Comm_size(MPI_COMM_WORLD,&size);
  try {
    require(size==1,"One MPI rank required");require(argc>=3,"Usage: probe OUT verify|bench WIDTH|wave");
    require(mkdir(argv[1],0755)==0||errno==EEXIST,"Cannot create result directory");
    const std::string mode=argv[2];
    if(mode=="verify")verify(argv[1]);
    else if(mode=="bench"){require(argc==4,"Benchmark width missing");benchmark(argv[1],std::stoul(argv[3]));}
    else if(mode=="wave")aprl::verification::runMovingWindowWave(argv[1]);
    else if(mode=="wave-static")aprl::verification::runMovingWindowWave(argv[1],true);
    else if(mode=="wave20")aprl::verification::runMovingWindowWave(argv[1],false,20);
    else if(mode=="wave-static20")aprl::verification::runMovingWindowWave(argv[1],true,20);
    else throw std::runtime_error("Unknown mode");
    std::cout<<"Completed "<<mode<<"; checks="<<checks<<'\n';MPI_Finalize();return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
}
