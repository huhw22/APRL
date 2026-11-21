#include "beam_cpml_start.h"
#include "blocked_cowan_reference.h"
#include <cstdlib>
#include <iostream>
#include <unistd.h>
using namespace aprl;
using namespace beam_cpml_start;

namespace {
double marker(std::size_t k,std::size_t p,std::size_t ch) {return 1.+p+1000*k+100000*ch;}
using Fields=std::array<YeeComponent*,6>;
Fields fields(EBFieldGrid& g) {return {{&g.ex(),&g.ey(),&g.ez(),&g.bx(),&g.by(),&g.bz()}};}
EBGridGeometry geometry(std::size_t nz) {return EBGridGeometry(14,12,nz,1.e-5,1.4e-5,8.e-6,8.e-6/SI::c,EBMaxwellSolver::CowanZ);}
void writePackage(const std::string& dir,bool left,SimulationConfig& cfg,std::size_t depth) {
  auto ref=geometry(depth);EBFieldGrid full(ref);
  const auto& cp=cfg.boundary.cpml;
  full.setBoundary(std::unique_ptr<EBBoundaryOperator>(new EBConvolutionalPML(ref,depth,0,cp,true,true)));
  auto& pml=dynamic_cast<EBConvolutionalPML&>(full.boundary());
  auto f=fields(full); auto keys=names(left);auto h=historyIds(left);
  auto method=cfg.mesh.movingWindow.cpmlStartup;
  YAML::Node p;
  p["version"]=1;p["method"]=method;
  p["reference_cells"]=p["window_cells"]=std::vector<std::size_t>{ref.nx,ref.ny,ref.nz};p["reference_offset"]=0;
  p["reference_origin_m"]=p["beam_position_offset_m"]=std::vector<double>{0,0,0};p["spacing_m"]=std::vector<double>{ref.dx,ref.dy,ref.dz};
  p["dt_s"]=ref.dt;p["reference_beta"]=.2;p["frame_gamma"]=cfg.mesh.boostGamma;p["beam_electrons"]=cfg.beam.electrons;
  p["beam_reference"]["initial_center_z"]=0.;
  p["cpml_cells"]=std::vector<std::size_t>{2,3,0};p["target_reflection"]=cp.targetReflection;
  p["kappa_max"]=cp.kappaMax;p["alpha_fraction"]=cp.alphaFraction;p["polynomial_order"]=cp.polynomialOrder;
  p["transverse_x_binomial3"]=true;p["longitudinal_shape"]="linear-binomial3";
  p["beam_fnv1a64"]=fingerprint(cfg.beam.file);
  for(unsigned c=0;c<keys.size();++c) {
    unsigned target=c<6?c:(h[c-6]<6?3:0)+(h[c-6]%6)/2;
    auto& u=*f[target];
    unsigned n=left?depth+1:u.nz();
    const unsigned width=left?(ref.ny+1)*2:c<6?exteriorIndices(u,2,3).size():pml.nativeZHistoryLayout(h[c-6]).valuesPerPlane;
    std::string file=dir+"/"+keys[c]+".bin";
    std::ofstream out(file,std::ios::binary|std::ios::trunc);
    for(unsigned k=0;k<n;++k)for(unsigned i=0;i<width;++i){const auto v=marker(k,i,c);out.write(reinterpret_cast<const char*>(&v),sizeof(v));}
    out.close();
    p["arrays_fnv1a64"][keys[c]]=fingerprint(file);
  }
  std::ofstream out(dir+"/startup.yaml");out<<p;
}
void one(const std::string& dir,bool left,SimulationConfig& cfg,unsigned rank,unsigned size) {
  const unsigned total=9,first=rank*total/size,stop=(rank+1)*total/size;
  auto global=geometry(total),local=geometry(stop-first);
  EBFieldGrid g(local);auto& cp=cfg.boundary.cpml;
  g.setBoundary(std::unique_ptr<EBBoundaryOperator>(new EBConvolutionalPML(local,total,first,cp,rank==0,rank+1==size)));
  auto& pml=dynamic_cast<EBConvolutionalPML&>(g.boundary());auto f=fields(g);for(auto* u:f)u->fill(7.);
  cfg.mesh.movingWindow.cpmlStartup=left?"left-20056":"osiris";
  if(rank==0)writePackage(dir,left,cfg,total);
  MPI_Barrier(MPI_COMM_WORLD);
  FieldVector<Double> origin;
  initialize(cfg,global,origin,.2,g,first,total,0,rank,MPI_COMM_WORLD);
  for(unsigned c=0;c<6;++c) {
    auto& u=*f[c];
    for(unsigned k=0;k<u.nz();++k) {
      unsigned next=0;
      for(unsigned j=0;j<u.ny();++j)for(unsigned i=0;i<u.nx();++i) {
        double expect=7.;
        if(left && i<2) expect+=marker(first+k,j*2+i,c);
        if(!left && (i<2 || i+2>=u.nx() || j<3 || j+3>=u.ny())) expect=marker(first+k,next++,c);
        require(u(i,j,k)==expect,"loaded field including corners and shared planes");
      }
    }
  }
  auto ids=historyIds(left);
  for(unsigned c=0;c<ids.size();++c) {
    const unsigned id=ids[c];const auto& u=*f[(id<6?3:0)+(id%6)/2];
    auto packed=pml.packNativeZHistory(id);
    auto layout=pml.nativeZHistoryLayout(id);
    for(unsigned k=0;k<layout.layers;++k)for(unsigned j=0;j<layout.valuesPerPlane;++j) {
      const unsigned slot=j/u.ny();
      double expected=left?(slot<2?marker(first+k,j%u.ny()*2+slot,c+6):0.):marker(first+k,j,c+6);
      require(packed[k*layout.valuesPerPlane+j]==expected,"loaded history including shared planes");
    }
  }
  cfg.beam.positionOffset[0]=1.e-5;
  bool rejected=false;
  try {initialize(cfg,global,origin,.2,g,first,total,0,rank,MPI_COMM_WORLD);}catch(const std::exception&){rejected=true;}
  require(rejected,"different physical beam rejected on every rank");
  cfg.beam.positionOffset[0]=0;
  if(rank==0) {std::fstream damaged(dir+"/Ex.bin",std::ios::binary|std::ios::out|std::ios::in);double zero=0.;damaged.write(reinterpret_cast<char*>(&zero),8);}
  MPI_Barrier(MPI_COMM_WORLD);
  rejected=false;
  try{initialize(cfg,global,origin,.2,g,first,total,0,rank,MPI_COMM_WORLD);}catch(const std::exception&){rejected=true;}
  require(rejected,"damaged array rejected on every rank");
  if(rank==0){for(auto& n:names(left))unlink((dir+"/"+n+".bin").c_str());unlink((dir+"/startup.yaml").c_str());}
}
void checkOpen(const std::string& dir,int rank) {
  const auto g=geometry(14);FieldVector<Double> origin;const std::size_t guards[]={2,2,0};const auto path=dir+"/open";
  if(!rank) {
    RelativisticParticleSI p;p.position[0]=7*g.dx;p.position[1]=6*g.dy;p.position[2]=7*g.dz;p.charge=-1.e-18;p.mass=1.e-25;
    verification::BlockedCowanReference writer(g,origin,{p},.2,guards,8,path,MPI_COMM_SELF,false,false,true,true,true);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  auto reader=verification::BlockedCowanReference::openPrepared(g,origin,.2,guards,8,path,MPI_COMM_SELF,false,false,true,true,true);
  const double beta=.2,gamma=1./std::sqrt(1.-beta*beta*SI::c*SI::c*SI::epsilon0*SI::mu0);
  double wall=0.;
  for(unsigned k=1;k<g.nz;++k)for(unsigned j=1;j<g.ny;++j)for(unsigned i=1;i<g.nx;++i) {
    const double x=i*g.dx,y=(static_cast<double>(j)-6)*g.dy,z=(static_cast<double>(k)-7)*g.dz*gamma;
    wall+=reader->chargeDensity(i,j,k)*g.dx*g.dy*g.dz*gamma/(4.*std::acos(-1.)*SI::epsilon0*std::sqrt(x*x+y*y+z*z));
  }
  require(std::abs(reader->potential(0,6,7)-wall)<2.e-11*std::abs(wall),"open reference boundary vs direct Coulomb sum");
  reader.reset();MPI_Barrier(MPI_COMM_WORLD);
  if(!rank){for(const auto n:{"metadata.txt","rho.bin","checkpoints.bin","phi.bin","side-potential.bin"})unlink((path+"/"+n).c_str());rmdir(path.c_str());}
}
}
int main(int argc,char**argv) {
  MPI_Init(&argc,&argv);int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  char dir[]="/tmp/aprl-cpml-start-XXXXXX";
  try {
    if(!rank){if(!mkdtemp(dir))throw std::runtime_error("Cannot make fixture directory");std::ofstream file(std::string(dir)+"/beam");file<<"hello";}
    MPI_Bcast(dir,sizeof(dir),MPI_CHAR,0,MPI_COMM_WORLD);
    SimulationConfig cfg;cfg.mesh.movingWindow.cpmlStartupDirectory=dir;cfg.beam.file=std::string(dir)+"/beam";
    require(fingerprint(cfg.beam.file)=="a430d84680aabd0b","published FNV1a64 check value");
    cfg.mesh.movingWindow.openTransverseReference=cfg.mesh.movingWindow.referenceThroughCpml=true;
    cfg.mesh.movingWindow.referenceSurfaceCurrent=false;cfg.mesh.transverseXBinomial3=cfg.mesh.longitudinalBinomial3=true;
    cfg.reference.inputPlaneZSet=false;cfg.reference.activationPlaneZSet=false;cfg.reference.initialCenterZ=0.;
    cfg.reference.holdTransverseUntilActivation=false;cfg.reference.frozenLabPlaneInjection=false;cfg.mesh.boostGamma=2.;
    cfg.beam.type=BeamInputType::Hdf5;cfg.beam.electrons=1000.;cfg.boundary.cpml.cells[0]=2;cfg.boundary.cpml.cells[1]=3;cfg.boundary.cpml.cells[2]=0;
    one(dir,false,cfg,rank,size);one(dir,true,cfg,rank,size);checkOpen(dir,rank);
    if(!rank){unlink(cfg.beam.file.c_str());rmdir(dir);std::cout<<"both CPML starts, MPI slabs, rejected inputs and direct open potential passed\n";}
  } catch(const std::exception& e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
  MPI_Finalize();
}
