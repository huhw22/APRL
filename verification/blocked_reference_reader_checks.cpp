#include "blocked_cowan_reference.h"
#include "compute_threads.h"

#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <mpi.h>

namespace {
using aprl::EBGridGeometry;
using aprl::FieldVector;
using aprl::RelativisticParticleSI;
using aprl::verification::BlockedCowanReference;
const char* const fileNames[4]={"rho.bin","checkpoints.bin","phi.bin","metadata.txt"};
using Files=std::array<std::string,4>;

void require(bool condition,const std::string& message) {
  if(!condition) throw std::runtime_error(message);
}
void freshDirectory(const std::string& path) {
  require(::mkdir(path.c_str(),0700)==0,"Cannot create fresh check directory: "+path);
}
struct Descriptor {
  int value;
  explicit Descriptor(int fd):value(fd) {}
  ~Descriptor() { if(value>=0) ::close(value); }
  Descriptor(const Descriptor&)=delete;
  Descriptor& operator=(const Descriptor&)=delete;
};
void synchronize(int fd) {
  int result;
  do { result=::fsync(fd); } while(result<0 && errno==EINTR);
  require(result==0,"Cannot synchronize reader-check fixture");
}
void writeFresh(const std::string& path,const std::string& data) {
  Descriptor output(::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600));
  require(output.value>=0,"Refusing to overwrite/cannot create check file: "+path);
  std::size_t done=0;
  while(done<data.size()) {
    const ssize_t count=::write(output.value,data.data()+done,data.size()-done);
    if(count<0 && errno==EINTR) continue;
    require(count>0,"Cannot write reader-check fixture: "+path);
    done+=static_cast<std::size_t>(count);
  }
  synchronize(output.value);
}
std::string readSmall(const std::string& path) {
  std::ifstream input(path.c_str(),std::ios::binary);
  require(input.good(),"Cannot read reader-check fixture: "+path);
  input.seekg(0,std::ios::end);
  const std::streamoff length=input.tellg();
  require(length>=0 && length<=4*1024*1024,"Check fixture exceeds its fixed small size limit");
  input.seekg(0,std::ios::beg);
  std::string result(static_cast<std::size_t>(length),'\0');
  if(length) input.read(&result[0],length);
  require(!length || input.gcount()==length,"Short reader-check fixture read");
  return result;
}
Files snapshot(const std::string& directory) {
  Files result;
  for(unsigned i=0;i<4;++i) result[i]=readSmall(directory+"/"+fileNames[i]);
  return result;
}
void createFixture(const std::string& directory,const Files& data,bool includeMetadata=true) {
  freshDirectory(directory);
  for(unsigned i=0;i<(includeMetadata?4U:3U);++i)
    writeFresh(directory+"/"+fileNames[i],data[i]);
}
std::uint64_t checksum(const std::string& data) {
  std::uint64_t result=UINT64_C(14695981039346656037);
  for(unsigned char c:data) { result^=c; result*=UINT64_C(1099511628211); }
  return result;
}
std::string hex(std::uint64_t value) {
  std::ostringstream result;
  result.imbue(std::locale::classic());
  result<<std::hex<<std::setfill('0')<<std::setw(16)<<value;
  return result.str();
}
void damageBeta(std::string& metadata) {
  const std::string key="\nbeta_bits ";
  const std::size_t found=metadata.find(key);
  require(found!=std::string::npos,"Missing beta field in check metadata");
  const std::size_t digit=found+key.size()+15;
  require(digit<metadata.size(),"Incomplete beta field in check metadata");
  metadata[digit]=metadata[digit]=='0'?'1':'0';
}
void reseal(std::string& metadata) {
  const std::size_t seal=metadata.rfind("checksum ");
  require(seal!=std::string::npos,"Missing check metadata checksum");
  const std::string body=metadata.substr(0,seal);
  metadata=body+"checksum "+hex(checksum(body))+"\n";
}
bool sameBits(double a,double b) {
  std::uint64_t aa=0,bb=0;
  static_assert(sizeof(a)==sizeof(aa),"Reader checks require binary64");
  std::memcpy(&aa,&a,sizeof(a));
  std::memcpy(&bb,&b,sizeof(b));
  return aa==bb;
}

void rootAction(int rank,const char* label,const std::function<void()>& action) {
  int success=1;
  if(rank==0) {
    try { action(); }
    catch(const std::exception& error) {
      success=0;
      std::cerr<<"root "<<label<<": "<<error.what()<<'\n';
    }
  }
  MPI_Bcast(&success,1,MPI_INT,0,MPI_COMM_WORLD);
  require(success!=0,std::string("Root reader-check phase failed: ")+label);
}
void allAction(int rank,const char* label,const std::function<void()>& action) {
  int failed=0,anyFailed=0;
  try { action(); }
  catch(const std::exception& error) {
    failed=1;
    std::cerr<<"rank "<<rank<<' '<<label<<": "<<error.what()<<'\n';
  }
  MPI_Allreduce(&failed,&anyFailed,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  require(anyFailed==0,std::string("Reader-check phase failed: ")+label);
}
void expectRejected(const std::function<void()>& action,const std::string& label,unsigned& count) {
  bool rejected=false;
  try { action(); }
  catch(const std::exception&) { rejected=true; }
  require(rejected,"Reader accepted negative control: "+label);
  ++count;
}

std::vector<RelativisticParticleSI> source(const EBGridGeometry& g,
    const FieldVector<double>& origin,const std::size_t guards[3],double beta) {
  const double fraction[8]={.013,.071,.183,.207,.397,.414,.797,.983};
  const double spacing[3]={g.dx,g.dy,g.dz};
  const std::size_t cells[3]={g.nx,g.ny,g.nz};
  std::vector<RelativisticParticleSI> result;
  for(unsigned z=0;z<8;++z) for(unsigned y=0;y<3;++y) for(unsigned x=0;x<4;++x) {
    RelativisticParticleSI p;
    const unsigned index[3]={x+1,y+2,z};
    for(unsigned axis=0;axis<3;++axis) {
      const double lo=guards[axis]+1.1,hi=cells[axis]-guards[axis]-1.1;
      p.position[axis]=origin[axis]+(lo+(hi-lo)*fraction[index[axis]])*spacing[axis];
    }
    p.properVelocity=FieldVector<double>(0.);
    p.properVelocity[2]=beta/std::sqrt(1.-beta*beta);
    p.weight=500.+x+3*y+7*z;
    p.mass=aprl::SI::electronMass*p.weight;
    p.charge=-aprl::SI::elementaryCharge*p.weight;
    p.id=p.sourceId=result.size();
    result.push_back(p);
  }
  return result;
}

struct NativeShape { std::size_t nx,ny,nz; };
std::array<NativeShape,6> shapes(const EBGridGeometry& g) {
  return {{{g.nx,g.ny+1,g.nz+1},{g.nx+1,g.ny,g.nz+1},
    {g.nx+1,g.ny+1,g.nz},{g.nx+1,g.ny,g.nz},{g.nx,g.ny+1,g.nz},
    {g.nx,g.ny,g.nz+1}}};
}
std::vector<double> nativeValues(BlockedCowanReference& reference,const EBGridGeometry& g) {
  std::vector<double> result;
  const auto native=shapes(g);
  std::size_t count=0;
  for(const auto& s:native) count+=s.nx*s.ny*s.nz;
  result.reserve(count);
  for(unsigned c=0;c<6;++c) {
    const auto& s=native[c];
    for(std::size_t k=0;k<s.nz;++k) for(std::size_t j=0;j<s.ny;++j)
      for(std::size_t i=0;i<s.nx;++i) result.push_back(reference.field(c,i,j,k));
  }
  return result;
}
std::vector<double> scalarValues(BlockedCowanReference& reference,const EBGridGeometry& g) {
  std::vector<double> result;
  result.reserve((g.nx+1)*(g.ny+1)*(g.nz+1));
  for(std::size_t k=0;k<=g.nz;++k) for(std::size_t j=0;j<=g.ny;++j)
    for(std::size_t i=0;i<=g.nx;++i) result.push_back(reference.potential(i,j,k));
  return result;
}
void broadcastValues(std::vector<double>& values,int rank) {
  std::uint64_t count=values.size();
  MPI_Bcast(&count,1,MPI_UINT64_T,0,MPI_COMM_WORLD);
  require(count<=2000000,"Reader check broadcast exceeds fixed small fixture");
  if(rank!=0) values.resize(static_cast<std::size_t>(count));
  MPI_Bcast(values.data(),static_cast<int>(count),MPI_DOUBLE,0,MPI_COMM_WORLD);
}
void checkIdleReader(const BlockedCowanReference& reader,std::uint64_t particles,
                     const double history[4]) {
  const auto& r=reader.report();
  require(r.complete && r.preparedReader && r.openSeconds>=0.,"Reader mode was not reported");
  require(r.preparationSeconds==0. && r.sourceSeconds==0. && r.forwardSeconds==0. &&
    r.reconstructionSeconds==0. && r.querySeconds==0.,"Reader claimed preparation work");
  require(r.sourceIndexBytes==0 && r.sourceSortingAllowanceBytes==0 &&
    r.solverWorkspaceBytes==0 && r.residentArrayBytes==0 &&
    r.potentialCacheCapacityBytes==0 && r.sourceCacheCapacityBytes==0 &&
    r.fieldCacheCapacityBytes==0 && r.queryWorkspaceBytes==0,
    "Opening reader allocated preparation/full-field arrays");
  require(r.bytesRead==0 && r.bytesWritten==0 && r.fieldQueries==0 && r.fieldCacheBuilds==0,
    "Opening reader read or regenerated scalar data");
  require(r.particles==particles && sameBits(r.sourceCharge,history[0]) &&
    sameBits(r.depositedCharge,history[1]) && sameBits(r.minimumPivot,history[2]) &&
    sameBits(r.minimumModeEigenvalue,history[3]),"Reader lost producer history");
}
}

int main(int argc,char** argv) {
  int provided=0;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int rank=0,ranks=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  try {
    require(argc==2 && (ranks==1 || ranks==2 || ranks==4),
      "Usage: aprl_blocked_reference_reader_checks NEW_OUTPUT_DIRECTORY (1, 2 or 4 ranks)");
    aprl::configureComputeThreads(provided,1);
    const std::string root=argv[1],normal=root+"/normal";
    const double dz=1.2543262415359034e-6,beta=.23;
    const EBGridGeometry g(34,36,81,5e-6,5e-6,dz,dz/aprl::SI::c,aprl::EBMaxwellSolver::CowanZ);
    FieldVector<double> origin(0.);
    origin[0]=-17.25*g.dx; origin[1]=7.5*g.dy; origin[2]=-43.5*g.dz;
    const std::size_t guards[3]={10,10,2},block=7;
    std::unique_ptr<BlockedCowanReference> writer,reader,lateReader;
    std::vector<double> expectedField,expectedPhi;
    std::uint64_t particles=0;
    double history[4]={};
    Files normalSnapshot;

    rootAction(rank,"SELF preparation",[&] {
      freshDirectory(root);
      const auto input=source(g,origin,guards,beta);
      writer.reset(new BlockedCowanReference(g,origin,input,beta,guards,block,normal,MPI_COMM_SELF));
      const auto& r=writer->report();
      require(r.complete && !r.preparedReader && r.preparationSeconds>0.,"Writer did not prepare its reference");
      particles=r.particles;
      history[0]=r.sourceCharge; history[1]=r.depositedCharge;
      history[2]=r.minimumPivot; history[3]=r.minimumModeEigenvalue;
      expectedField=nativeValues(*writer,g);
      expectedPhi=scalarValues(*writer,g);
      require(!expectedField.empty() && !expectedPhi.empty(),"Empty reader comparison oracle");
      normalSnapshot=snapshot(normal);
      writer->writeReport(root+"/writer.json");
    });
    MPI_Bcast(&particles,1,MPI_UINT64_T,0,MPI_COMM_WORLD);
    MPI_Bcast(history,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
    broadcastValues(expectedField,rank);
    broadcastValues(expectedPhi,rank);

    allAction(rank,"open SELF reader and compare native values",[&] {
      reader=BlockedCowanReference::openPrepared(g,origin,beta,guards,block,normal);
      checkIdleReader(*reader,particles,history);
      std::size_t index=0;
      const auto native=shapes(g);
      for(unsigned c=0;c<6;++c) {
        const auto& s=native[c];
        for(std::size_t k=0;k<s.nz;++k) for(std::size_t j=0;j<s.ny;++j)
          for(std::size_t i=0;i<s.nx;++i)
            require(sameBits(reader->field(c,i,j,k),expectedField.at(index++)),
              "Reader native E/B differs bitwise from original writer");
      }
      require(index==expectedField.size(),"Native field oracle coverage mismatch");
      index=0;
      for(std::size_t k=0;k<=g.nz;++k) for(std::size_t j=0;j<=g.ny;++j)
        for(std::size_t i=0;i<=g.nx;++i)
          require(sameBits(reader->potential(i,j,k),expectedPhi.at(index++)),
            "Reader phi differs bitwise from original writer");
      require(index==expectedPhi.size(),"Potential oracle coverage mismatch");
      require(reader->report().bytesWritten==0 && reader->report().preparationSeconds==0.,
        "Reader comparison regenerated or wrote scientific data");
      reader->writeReport(root+"/reader-rank-"+std::to_string(rank)+".json");
    });

    rootAction(rank,"create separate failure fixtures",[&] {
      createFixture(root+"/fault-missing-metadata",normalSnapshot,false);
      Files data=normalSnapshot;
      data[3].resize(data[3].size()-1);
      createFixture(root+"/fault-truncated-metadata",data);
      data=normalSnapshot;
      damageBeta(data[3]);
      createFixture(root+"/fault-checksum",data);
      reseal(data[3]);
      createFixture(root+"/fault-contract",data);
      for(unsigned f=0;f<3;++f) {
        data=normalSnapshot;
        require(data[f].size()>sizeof(double),"Too-small scalar failure fixture");
        data[f].resize(data[f].size()-sizeof(double));
        createFixture(root+"/fault-length-"+std::to_string(f),data);
      }
      const double nan=std::numeric_limits<double>::quiet_NaN();
      for(unsigned f:{0U,2U}) {
        data=normalSnapshot;
        std::memcpy(&data[f][0],&nan,sizeof(nan));
        createFixture(root+"/fault-nonfinite-"+std::to_string(f),data);
      }
      createFixture(root+"/fault-after-open",normalSnapshot);
      const auto input=source(g,origin,guards,beta);
      BlockedCowanReference clamped(g,origin,input,beta,guards,128,root+"/clamped-block",MPI_COMM_SELF);
      require(clamped.report().blockLayers==g.nz-2*guards[2]-1,
        "Large requested block fixture did not exercise clamping");
    });

    rootAction(rank,"matched F3 charge with transverse smoothing",[&] {
      auto input=source(g,origin,guards,beta);
      input.resize(1);
      auto& particle=input.front();
      particle.position[0]=origin[0]+16.25*g.dx;
      particle.position[1]=origin[1]+17.75*g.dy;
      particle.position[2]=origin[2]+35.125*g.dz;
      const std::size_t sourceBlock=12;
      BlockedCowanReference linear(g,origin,input,beta,guards,sourceBlock,root+"/linear-x",
                                   MPI_COMM_SELF,false,false,true);
      BlockedCowanReference smoothed(g,origin,input,beta,guards,sourceBlock,root+"/smoothed-x",
                                     MPI_COMM_SELF,false,false,true,true);
      auto saved=BlockedCowanReference::openPrepared(g,origin,beta,guards,sourceBlock,
                            root+"/smoothed-x",MPI_COMM_SELF,false,false,true,true);
      const double scale=std::abs(particle.charge/(g.dx*g.dy*g.dz));
      for(unsigned k=30;k<=41;++k)for(unsigned j=13;j<=22;++j)for(unsigned i=12;i<=20;++i) {
        const double expected=.25*linear.chargeDensity(i-1,j,k)+.5*linear.chargeDensity(i,j,k)
                             +.25*linear.chargeDensity(i+1,j,k);
        require(std::abs(smoothed.chargeDensity(i,j,k)-expected)<1.e-14*scale,
                "Initial reference and particle source use different transverse weights");
        require(saved->chargeDensity(i,j,k)==smoothed.chargeDensity(i,j,k),
                "Prepared reader lost the x-smoothed charge density");
      }
      unsigned rejected=0;
      expectRejected([&] {
        BlockedCowanReference::openPrepared(g,origin,beta,guards,sourceBlock,root+"/smoothed-x",
                                            MPI_COMM_SELF,false,false,true,false);
      },"read smoothed charge as linear",rejected);
    });

    unsigned negatives=0;
    allAction(rank,"parameter/metadata/data negative controls",[&] {
      EBGridGeometry wrongGeometry=g;
      ++wrongGeometry.nx;
      expectRejected([&] { BlockedCowanReference::openPrepared(wrongGeometry,origin,beta,guards,block,normal); },
        "geometry",negatives);
      auto wrongOrigin=origin;
      wrongOrigin[0]=std::nextafter(wrongOrigin[0],std::numeric_limits<double>::infinity());
      expectRejected([&] { BlockedCowanReference::openPrepared(g,wrongOrigin,beta,guards,block,normal); },
        "origin by one ULP",negatives);
      expectRejected([&] { BlockedCowanReference::openPrepared(g,origin,std::nextafter(beta,1.),guards,block,normal); },
        "beta by one ULP",negatives);
      const std::size_t wrongGuards[3]={9,10,2};
      expectRejected([&] { BlockedCowanReference::openPrepared(g,origin,beta,wrongGuards,block,normal); },
        "guards",negatives);
      expectRejected([&] { BlockedCowanReference::openPrepared(g,origin,beta,guards,block+1,normal); },
        "block layers",negatives);
      expectRejected([&] { BlockedCowanReference::openPrepared(g,origin,beta,guards,block,normal,MPI_COMM_NULL); },
        "null communicator",negatives);
      if(ranks>1) {
        expectRejected([&] { BlockedCowanReference::openPrepared(g,origin,beta,guards,block,normal,MPI_COMM_WORLD); },
          "multi-rank reader communicator",negatives);
        expectRejected([&] {
          const std::vector<RelativisticParticleSI> none;
          BlockedCowanReference wrong(g,origin,none,beta,guards,block,
            root+"/must-not-prepare-rank-"+std::to_string(rank),MPI_COMM_WORLD);
        },"multi-rank writer communicator",negatives);
      }
      for(const char* name:{"fault-missing-metadata","fault-truncated-metadata","fault-checksum","fault-contract"})
        expectRejected([&] { BlockedCowanReference::openPrepared(g,origin,beta,guards,block,root+"/"+name); },
          name,negatives);
      for(unsigned f=0;f<3;++f)
        expectRejected([&] { BlockedCowanReference::openPrepared(g,origin,beta,guards,block,
          root+"/fault-length-"+std::to_string(f)); },"scalar file length",negatives);
      for(unsigned f:{0U,2U}) {
        auto bad=BlockedCowanReference::openPrepared(g,origin,beta,guards,block,
          root+"/fault-nonfinite-"+std::to_string(f));
        checkIdleReader(*bad,particles,history);
        expectRejected([&] {
          if(f==0) bad->chargeDensity(guards[0]+1,guards[1]+1,guards[2]+1);
          else bad->potential(guards[0]+1,guards[1]+1,guards[2]+1);
        },"nonfinite scalar on demand",negatives);
      }
      auto clamped=BlockedCowanReference::openPrepared(g,origin,beta,guards,128,root+"/clamped-block");
      require(clamped->report().preparedReader,"Clamped block reader failed");
      expectRejected([&] { BlockedCowanReference::openPrepared(g,origin,beta,guards,256,root+"/clamped-block"); },
        "different requested block with equal clamped size",negatives);
      lateReader=BlockedCowanReference::openPrepared(g,origin,beta,guards,block,root+"/fault-after-open");
    });

    rootAction(rank,"truncate only the dedicated already-open fixture",[&] {
      Descriptor file(::open((root+"/fault-after-open/phi.bin").c_str(),O_WRONLY));
      require(file.value>=0 && ::ftruncate(file.value,0)==0,"Cannot truncate dedicated reader failure fixture");
      synchronize(file.value);
    });
    allAction(rank,"on-demand short read after successful open",[&] {
      lateReader->clearCache();
      expectRejected([&] { lateReader->potential(guards[0]+1,guards[1]+1,guards[2]+1); },
        "short read after open",negatives);
      require(reader->report().bytesWritten==0,"Normal reader wrote scientific data");
    });

    unsigned least=0,most=0;
    MPI_Allreduce(&negatives,&least,1,MPI_UNSIGNED,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&negatives,&most,1,MPI_UNSIGNED,MPI_MAX,MPI_COMM_WORLD);
    require(least==most && least>=17,"Ranks did not execute the same negative controls");
    rootAction(rank,"confirm normal data unchanged and write result",[&] {
      require(snapshot(normal)==normalSnapshot,"Reader checks modified the normal prepared files");
      std::ostringstream result;
      result << "{\"status\":\"PASS\",\"world_ranks\":" << ranks
        << ",\"writer_comm\":\"MPI_COMM_SELF\",\"reader_comm\":\"MPI_COMM_SELF\""
        << ",\"native_values_compared_per_rank\":" << expectedField.size()
        << ",\"phi_values_compared_per_rank\":" << expectedPhi.size()
        << ",\"native_and_phi_bitwise_equal\":true,\"reader_has_no_preparation_work\":true"
        << ",\"negative_controls_per_rank\":" << negatives
        << ",\"normal_files_byte_identical\":true"
        << ",\"scope\":\"Shared-file prepared readers in a 1/2/4-rank WORLD; no moving PIC, distributed solve, peak RSS or production speed claim\"}\n";
      writeFresh(root+"/checks.json",result.str());
      std::cout<<"PASS blocked prepared readers ranks="<<ranks<<" negatives_per_rank="<<negatives<<'\n';
    });
    reader.reset(); lateReader.reset(); writer.reset();
    MPI_Finalize();
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"rank "<<rank<<": "<<error.what()<<'\n';
    MPI_Abort(MPI_COMM_WORLD,1);
    return 1;
  }
}
