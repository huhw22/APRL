// Manufactured timelike endpoint worldlines: detector accuracy/regression,
// not a particle-dynamics or Maxwell physical benchmark.
#include "lab_detectors.h"
#include "particle_energy_audit.h"
#include "runtime_util.h"
#include <hdf5.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace aprl;
using Vec=FieldVector<Double>;
namespace {
void require(bool good,const char* what){if(!good)throw std::runtime_error(what);}
void near(long double a,long double b,long double tolerance,const char* what){require(std::abs(a-b)<=tolerance,what);}
template<class T>std::vector<T> column(hid_t d,const char* name,hid_t native,std::size_t count){
  hid_t type=H5Tcreate(H5T_COMPOUND,sizeof(T));require(type>=0,"Cannot create column type");
  require(H5Tinsert(type,name,0,native)>=0,"Cannot bind column");std::vector<T> out(count);
  require(H5Dread(d,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,out.data())>=0,"Cannot read column");H5Tclose(type);return out;
}
std::size_t countRecords(hid_t file,const char* group){
  const std::string root(group);std::uint64_t count=0;unsigned char complete=0;
  hid_t d=H5Dopen2(file,(root+"/complete").c_str(),H5P_DEFAULT);require(d>=0,"No complete marker");
  require(H5Dread(d,H5T_NATIVE_UCHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,&complete)>=0,"Read complete");H5Dclose(d);
  d=H5Dopen2(file,(root+"/committed_records").c_str(),H5P_DEFAULT);require(d>=0,"No committed marker");
  require(H5Dread(d,H5T_NATIVE_UINT64,H5S_ALL,H5S_ALL,H5P_DEFAULT,&count)>=0,"Read committed");H5Dclose(d);
  require(complete==1,"Unclosed detector");return count;
}
bool sameLegacy(const std::string& a,const std::string& b){
  hid_t fa=H5Fopen(a.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),fb=H5Fopen(b.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  require(fa>=0&&fb>=0,"Missing legacy detector");const auto n=countRecords(fa,"/particle_plane");
  require(n==countRecords(fb,"/particle_plane"),"Legacy count changed");
  hid_t da=H5Dopen2(fa,"/particle_plane/records",H5P_DEFAULT),db=H5Dopen2(fb,"/particle_plane/records",H5P_DEFAULT);
  hid_t ta=H5Dget_type(da),tb=H5Dget_type(db);require(H5Tequal(ta,tb)>0,"Legacy dtype changed");
  std::vector<unsigned char> aa(n*H5Tget_size(ta)),bb(aa.size());
  require(H5Dread(da,ta,H5S_ALL,H5S_ALL,H5P_DEFAULT,aa.data())>=0,"Read legacy a");
  require(H5Dread(db,tb,H5S_ALL,H5S_ALL,H5P_DEFAULT,bb.data())>=0,"Read legacy b");
  const bool result=aa==bb;H5Tclose(ta);H5Tclose(tb);H5Dclose(da);H5Dclose(db);H5Fclose(fa);H5Fclose(fb);return result;
}
RelativisticParticleSI state(const BoostFrameTransform& frame,double t,double z,int rank,unsigned index,double& boxTime){
  RelativisticParticleSI p;p.id=rank+1;p.sourceId=1001+rank;p.weight=rank+1;
  p.mass=SI::electronMass*p.weight;p.charge=-SI::elementaryCharge*p.weight;
  p.position=Vec(0.);p.position[0]=.0001;p.position[1]=-.0001;
  frame.labToBox(t,z,boxTime,p.position[2]);Vec u(0.);
  const double ux[3]={.2,.4,.9},uy[3]={-.1,.1,-.3},uz[3]={.5,.7,.2};
  u[0]=ux[index];u[1]=uy[index];u[2]=uz[index];frame.properVelocityLabToBox(u,p.properVelocity);return p;
}
void check(const std::string& root,int rank,int size){
  const EBGridGeometry grid(8,8,128*size,.001,.001,.0001,.0001/SI::c);
  Vec origin(0.),local(0.);origin[0]=origin[1]=-.004;origin[2]=-.004;local=origin;local[2]+=rank*128*.0001;
  BoostFrameTransform frame;frame.setFromGamma(2.,SI::c,0.);
  DetectorConfig off;off.directory=root+"/off";
  const char* names[3]={"entry","turn","probe"};const double planes[3]={-.00025,-.00015,0.};
  for(unsigned i=0;i<3;++i){ParticleDetectorPlaneConfig p;p.name=names[i];p.z=planes[i];p.bufferRecords=2;p.compression=i%2;off.particlePlanes.push_back(p);}
  DetectorConfig on=off;on.directory=root+"/on";on.substepParticlePlanes=true;
  LabDetectorManager old(off,grid,origin,local,frame,MPI_COMM_WORLD),diagnostic(on,grid,origin,local,frame,MPI_COMM_WORLD);
  require(!old.hasSubstepParticlePlanes()&&diagnostic.hasSubstepParticlePlanes(),"Default/enable flag wrong");
  double t[3],tf;auto a=state(frame,0.,-.0002,rank,0,t[0]);auto b=state(frame,1e-12,-.00015,rank,1,t[1]);
  auto c=state(frame,3e-12,.0001,rank,2,t[2]);auto frozen=state(frame,-1e-12,-.0003,rank,0,tf);
  old.captureParticleStep(frozen,a,tf,t[0]);diagnostic.captureParticleStep(frozen,a,tf,t[0]);
  diagnostic.captureSubstepParticleStep(frozen,a,tf,t[0]);
  old.captureParticleStep(a,c,t[0],t[2]);diagnostic.captureParticleStep(a,c,t[0],t[2]);
  for(unsigned discarded=0;discarded<3;++discarded){std::vector<SubstepParticleCrossing> trial;
    diagnostic.previewSubstepParticleStep(a,b,t[0],t[1],trial);diagnostic.previewSubstepParticleStep(b,c,t[1],t[2],trial);}
  ParticleEnergyTrace trace;trace.reset(frame,&diagnostic,false);
  trace.observe(a,b,t[0],t[1],NULL);trace.observe(b,c,t[1],t[2],NULL);
  require(trace.observedSubsteps==2,"Substep observer missed segment");
  diagnostic.commitSubstepParticleCrossings(trace.substepCrossings,t[2]);
  old.collectParticleCrossings();diagnostic.collectParticleCrossings();old.close(true);diagnostic.close(true);
  // A committed terminal prefix must not leak a later substep-plane event.
  DetectorConfig terminal=on;terminal.directory=root+"/terminal";terminal.particlePlanes.erase(terminal.particlePlanes.begin(),terminal.particlePlanes.begin()+2);
  LabDetectorManager cut(terminal,grid,origin,local,frame,MPI_COMM_WORLD);
  std::vector<SubstepParticleCrossing> events;cut.previewSubstepParticleStep(b,c,t[1],t[2],events);
  require(events.size()==1,"Terminal crossing setup differs");cut.commitSubstepParticleCrossings(events,events[0].timeBox-1e-20);
  cut.collectParticleCrossings();cut.close(true);MPI_Barrier(MPI_COMM_WORLD);
  if(rank!=0)return;
  for(unsigned plane=0;plane<3;++plane){
    const std::string name=names[plane];require(sameLegacy(root+"/off/"+name+".h5",root+"/on/"+name+".h5"),"Substep enable changed legacy science");
    hid_t f=H5Fopen((root+"/on/"+name+"-substeps.h5").c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
    require(f>=0,"Missing substep file");const auto n=countRecords(f,"/particle_substep_plane");require(n==static_cast<unsigned>(size),"Duplicate or missing substep events");
    hid_t d=H5Dopen2(f,"/particle_substep_plane/records",H5P_DEFAULT);
    const auto ids=column<std::uint64_t>(d,"source_id",H5T_NATIVE_UINT64,n);
    const auto time=column<double>(d,"time_s",H5T_NATIVE_DOUBLE,n);
    for(std::size_t i=0;i<n;++i){require(ids[i]>=1001&&ids[i]<1001+size,"Source identity changed");
      near(time[i],plane==0?-5e-13L:(plane==1?1e-12L:2.2e-12L),2e-26L,"Substep crossing time wrong");}
    if(plane==2){const auto gamma=column<double>(d,"gamma_lab",H5T_NATIVE_DOUBLE,n),mass=column<double>(d,"mass_kg",H5T_NATIVE_DOUBLE,n),kinetic=column<double>(d,"kinetic_lab_J",H5T_NATIVE_DOUBLE,n);
      const long double expected=.4L*std::sqrt(1.66L)+.6L*std::sqrt(1.94L);
      for(std::size_t i=0;i<n;++i){near(gamma[i],expected,3e-15L,"Affine endpoint gamma was recomputed on shell");
        const long double energy=(long double)mass[i]*SI::c*SI::c*(expected-1.);
        near(kinetic[i],energy,3e-15L*std::abs(energy),"Stable endpoint kinetic energy mismatch");}
    }
    H5Dclose(d);H5Fclose(f);
  }
  hid_t f=H5Fopen((root+"/terminal/probe-substeps.h5").c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  require(f>=0&&countRecords(f,"/particle_substep_plane")==0,"Terminal prefix leaked crossing");H5Fclose(f);
}
}
int main(int argc,char** argv){MPI_Init(&argc,&argv);int rank=0,size=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  try{require(argc==2,"Need fresh independent output directory");createDirectories(argv[1],MPI_COMM_WORLD);check(argv[1],rank,size);
    if(rank==0)std::cout<<"{\"passed\":true,\"ranks\":"<<size<<",\"legacy_bitwise_unchanged\":true,\"substep_time_s\":2.2e-12,\"legacy_time_s\":2e-12,\"affine_four_momentum_pass\":true,\"terminal_prefix_pass\":true,\"physical_acceptance\":false}\n";
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}MPI_Finalize();return 0;}
