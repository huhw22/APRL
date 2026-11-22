#include "blocked_cowan_reference.h"
#include "compute_threads.h"
#include "open_transverse_trace.h"
#include "reference_collective.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <sys/stat.h>
#include <unistd.h>

using namespace aprl;
using namespace aprl::verification;
namespace {
void require(bool ok,const char* what) {if(!ok)throw std::runtime_error(what);}
std::vector<char> contents(const std::string& path) {
  std::ifstream f(path,std::ios::binary);require(bool(f),"Missing reference comparison file");
  return std::vector<char>(std::istreambuf_iterator<char>(f),std::istreambuf_iterator<char>());
}
void clean(const std::string& path) {
  for(const auto n:{"metadata.txt","rho.bin","phi.bin","checkpoints.bin","side-potential.bin"})
    unlink((path+"/"+n).c_str());
  rmdir(path.c_str());
}
template<class Action> void rejected(MPI_Comm comm,Action action) {
  bool failed=false;
  try {action();}catch(const std::exception&) {failed=true;}
  reference_detail::collectively(comm,[&]{require(failed,"Collective failure did not reach every rank");});
}
void reference(const std::string& root,MPI_Comm comm,bool open,bool quintic,bool f3,unsigned coarse=0) {
  int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
  const double coarseDz[]={4.e-6,4.e-5,3.379676386553192e-5,2.7584815139792933e-5};
  const double coarseBeta[]={.73,.9720980709877257,.98596122023755,.99062570584060};
  const double dz=coarseDz[coarse],beta=coarseBeta[coarse];
  const EBGridGeometry g(30,26,49,4.e-5,coarse?4.e-5:3.e-5,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
  const std::size_t guards[]={3,3,0};
  FieldVector<double> origin(0.);origin[0]=-.0003;origin[2]=-.00005;
  const std::string serial=root+"/serial",parallel=root+"/parallel";
  std::vector<RelativisticParticleSI> particles;
  if(rank==0)for(unsigned n=0;n<37;++n) {
    RelativisticParticleSI p;p.id=n+1;p.mass=SI::electronMass;p.charge=-(1.+n*.01)*1.e-18;
    p.position=origin;
    p.position[0]+=(11.2+double((n*7)%41)/10.)*g.dx;
    p.position[1]+=(9.1+double((n*11)%37)/10.)*g.dy;
    p.position[2]+=(10.7+double((n*13)%251)/10.)*g.dz;
    particles.push_back(p);
  }
  std::unique_ptr<BlockedCowanReference> baseline;
  reference_detail::collectively(comm,[&] {
    if(rank==0) {
      require(mkdir(root.c_str(),0700)==0,"Cannot create test output");
      baseline.reset(new BlockedCowanReference(g,origin,particles,beta,guards,8,serial,MPI_COMM_SELF,
                                              quintic,false,f3,f3,open));
    }
  });
  auto actual=BlockedCowanReference::prepareCollective(g,origin,particles,beta,guards,8,parallel,
                                                      comm,quintic,false,f3,f3,open);
  reference_detail::collectively(comm,[&] {
    if(rank==0) {
      for(const auto name:{"rho.bin","phi.bin","checkpoints.bin","metadata.txt"})
        require(contents(serial+"/"+name)==contents(parallel+"/"+name),"Parallel reference differs bitwise");
      if(open) {
        require(contents(serial+"/side-potential.bin")==contents(parallel+"/side-potential.bin"),
                "Parallel open-wall potential differs bitwise");
        const auto& r=actual->report();
        require(r.openTraceRanks==ranks && r.openTraceSourceNodes>0 && r.openTraceMinRankEvaluations>0,
                "Open-wall work was not shared across communicator ranks");
        require(r.openTraceMaxRankEvaluations-r.openTraceMinRankEvaluations<=24,
                "Unexpected potential target imbalance");
      }
    } else require(actual->report().preparedReader && actual->report().bytesWritten==0,
                   "Worker must return a read-only reference handle");
    // Query all six native field components on every rank, including ghost
    // faces/corners and block transitions; compare against the singleton reader.
    auto oracle=BlockedCowanReference::openPrepared(g,origin,beta,guards,8,serial,MPI_COMM_SELF,
                                                   quintic,false,f3,f3,open);
    EBFieldGrid fields(g);actual->fill(fields);
    const YeeComponent* components[]={&fields.ex(),&fields.ey(),&fields.ez(),&fields.bx(),&fields.by(),&fields.bz()};
    for(unsigned c=0;c<6;++c) {
      const auto& u=*components[c];
      for(unsigned k=0;k<u.nz();++k)for(unsigned j=0;j<u.ny();++j)for(unsigned i=0;i<u.nx();++i) {
        const double expected=oracle->field(c,i,j,k),value=u(i,j,k);
        require(std::memcmp(&expected,&value,sizeof(double))==0,"Native parallel E/B differs bitwise");
      }
    }
    if(coarse) {
      // Independent full 3-D sine spectrum, expressed as native difference
      // wave numbers times transverse smoothing symbols. This does not use
      // the writer's assembled 27-point stencil or tridiagonal eigen formula.
      const auto c=EBCowanCoefficients::forZDispersion(g.dx,g.dy,g.dz);
      const double delta[]={c.deltaYZ,c.deltaZX,c.deltaXY},h[]={g.dx,g.dy,g.dz};
      const double speed=beta*SI::c,eta=speed*speed*SI::epsilon0*SI::mu0;
      double oldBound=1.,minimum=std::numeric_limits<double>::infinity();
      for(int a:{-1,1})for(int b:{-1,1})
        oldBound=std::min(oldBound,c.alpha[2]+2*c.beta[0]*a+2*c.beta[1]*b+4*delta[2]*a*b-eta);
      require(oldBound<0.,"Coarse regression must exercise the old false rejection");
      const double pi=std::acos(-1.);
      for(unsigned k=1;k<g.nz;++k)for(unsigned j=1;j<g.ny;++j)for(unsigned i=1;i<g.nx;++i) {
        const double angle[]={pi*i/g.nx,pi*j/g.ny,pi*k/g.nz};
        const double cosine[]={std::cos(angle[0]),std::cos(angle[1]),std::cos(angle[2])};
        double eigen=0.;
        for(unsigned d=0;d<3;++d) {
          const unsigned a=(d+1)%3,b=(d+2)%3;
          const double smooth=c.alpha[d]+2*c.beta[a]*cosine[a]+2*c.beta[b]*cosine[b]+
                              4*delta[d]*cosine[a]*cosine[b]-(d==2?eta:0.);
          const double wave=2*std::sin(angle[d]/2)/h[d];
          eigen+=smooth*wave*wave;
        }
        minimum=std::min(minimum,eigen);
      }
      require(minimum>0. && std::abs(actual->report().minimumModeEigenvalue/minimum-1)<2.e-11,
              "Finite-grid SPD certificate differs from full 3-D spectrum");
      long double rhs2=0.,error2=0.;
      for(unsigned k=1;k<g.nz;++k)for(unsigned j=1;j<g.ny;++j)for(unsigned i=1;i<g.nx;++i) {
        const double rhs=actual->chargeDensity(i,j,k)/SI::epsilon0;
        const double div=(fields.ex()(i,j,k)-fields.ex()(i-1,j,k))/g.dx+
                         (fields.ey()(i,j,k)-fields.ey()(i,j-1,k))/g.dy+
                         (fields.ez()(i,j,k)-fields.ez()(i,j,k-1))/g.dz;
        const double error=div-rhs;rhs2+=static_cast<long double>(rhs)*rhs;
        error2+=static_cast<long double>(error)*error;
      }
      require(rhs2>0. && std::sqrt(error2/rhs2)<5.e-11,"Coarse-grid native Gauss law failed");
      if(rank==0)std::cout<<"coarse="<<coarse<<" old_directional_bound="<<oldBound
                         <<" full_min_eigenvalue="<<minimum<<" gauss_relative_l2="
                         <<std::sqrt(error2/rhs2)<<'\n';
    }
  });
  // Failed root setup must propagate without leaving workers in a receive.
  rejected(comm,[&] {BlockedCowanReference::prepareCollective(g,origin,particles,beta,guards,8,parallel,
                                                            comm,quintic,false,f3,f3,open);});
  rejected(comm,[&] {BlockedCowanReference::prepareCollective(g,origin,particles,1.1,guards,8,root+"/invalid-beta",
                                                            comm,quintic,false,f3,f3,open);});
  actual.reset();baseline.reset();MPI_Barrier(comm);
  if(rank==0){clean(serial);clean(parallel);rmdir(root.c_str());}
}
void traces(MPI_Comm comm) {
  int rank=0;MPI_Comm_rank(comm,&rank);
  const EBGridGeometry g(128,128,1000,1.,1.,1.,1./SI::c,EBMaxwellSolver::CowanZ);
  OpenTransverseTrace::Sources source;
  if(rank==0) {
    // A signed pair with a transverse nodal zero stresses the interpolant's
    // relative-error guard and exercises native-point fallback collectively.
    source.x={63.,65.};source.y={64.,64.};source.z={430.3,430.3};source.q={1.,-1.};
  }
  OpenTransverseTrace baseline,parallel;
  reference_detail::collectively(comm,[&] {if(rank==0){baseline.configure(g);baseline.build(source,g,1.);parallel.configure(g);}});
  parallel.build(source,g,1.,comm);
  reference_detail::collectively(comm,[&] {
    if(rank==0) {
      require(baseline.fallbackTiles>0,"Signed-source negative control did not exercise fallback");
      require(baseline.fallbackTiles==parallel.fallbackTiles && baseline.values==parallel.values,
              "Distributed fallback differs from singleton trace");
    }
  });
  if(rank==0)source.x[0]=-1.;
  rejected(comm,[&]{parallel.build(source,g,1.,comm);});
  source=OpenTransverseTrace::Sources();
  if(rank==0)parallel.configure(g);
  parallel.build(source,g,1.,comm);
  reference_detail::collectively(comm,[&] {
    require(parallel.sourceNodes==0 && parallel.localPotentialEvaluations==0,"Empty-source collective did work");
  });
}
}
int main(int argc,char** argv) {
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int rank=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  try {
    configureComputeThreads(provided);
    require(argc==2,"Usage: aprl_reference_parallel_checks NEW_OUTPUT_DIRECTORY");
    const std::string root=argv[1];
    reference_detail::collectively(MPI_COMM_WORLD,[&]{if(rank==0)require(mkdir(root.c_str(),0700)==0,"Output exists");});
    reference(root+"/f3",MPI_COMM_WORLD,true,false,true);
    reference(root+"/quintic",MPI_COMM_WORLD,true,true,false);
    reference(root+"/closed",MPI_COMM_WORLD,false,false,false);
    for(unsigned c=1;c<=3;++c)reference(root+"/coarse-"+std::to_string(c),MPI_COMM_WORLD,true,false,true,c);
    traces(MPI_COMM_WORLD);
    // The helper must honor the supplied communicator, not MPI_COMM_WORLD.
    MPI_Comm subset=MPI_COMM_NULL;
    MPI_Comm_split(MPI_COMM_WORLD,rank%2,rank,&subset);
    reference(root+"/subset-"+std::to_string(rank%2),subset,true,false,true);
    MPI_Comm_free(&subset);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){rmdir(root.c_str());std::cout<<"PASS: MPI/OpenMP wall targets, native E/B, signed fallback, empty source, root/trace failures, subcommunicators\n";}
    MPI_Finalize();return 0;
  }catch(const std::exception& e){std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
  return 1;
}
