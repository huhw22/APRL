#include "blocked_cowan_reference.h"
#include "transverse_source_shape.h"
#include "longitudinal_shape.h"
#include "compute_threads.h"
#include "open_transverse_trace.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <fftw3.h>
#include <iomanip>
#include <limits>
#include <locale>
#include <mpi.h>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>

namespace aprl { namespace verification {
namespace {
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now()-start).count();
}
std::size_t product(std::size_t a, std::size_t b) {
  if (b && a > std::numeric_limits<std::size_t>::max()/b)
    throw std::overflow_error("Blocked reference array size overflow");
  return a*b;
}
std::size_t bytes(std::size_t values) { return product(values,sizeof(Double)); }
std::size_t plus(std::size_t a, std::size_t b) {
  if (a > std::numeric_limits<std::size_t>::max()-b)
    throw std::overflow_error("Blocked reference dimension overflow");
  return a+b;
}
std::runtime_error ioError(const std::string& operation, const std::string& path) {
  return std::runtime_error(operation+" "+path+": "+std::strerror(errno));
}
class DiskFile {
public:
  DiskFile() : fd_(-1),writable_(false) {}
  ~DiskFile() { if (fd_>=0) ::close(fd_); }
  DiskFile(const DiskFile&)=delete;
  DiskFile& operator=(const DiskFile&)=delete;
  void create(const std::string& path) {
    if(fd_>=0) throw std::logic_error("Reference file is already open");
    path_=path;
    fd_=::open(path.c_str(),O_RDWR|O_CREAT|O_EXCL,0600);
    if (fd_<0) throw ioError("Cannot create new reference file",path_);
    writable_=true;
  }
  void openReadOnly(const std::string& path) {
    if(fd_>=0) throw std::logic_error("Reference file is already open");
    path_=path;
    fd_=::open(path.c_str(),O_RDONLY);
    if(fd_<0) throw ioError("Cannot open prepared reference file",path_);
  }
  void synchronize() const {
    if(!writable_) throw std::logic_error("Cannot synchronize a read-only reference writer");
    int result;
    do { result=::fsync(fd_); } while(result<0 && errno==EINTR);
    if(result<0) throw ioError("Cannot synchronize completed reference file",path_);
  }
  void transfer(void* data, std::size_t count, std::size_t offset,
                bool writing, std::uint64_t* counter) const {
    if(writing && !writable_)
      throw std::logic_error("Prepared reference readers cannot modify data files");
    if (offset > static_cast<std::uint64_t>(std::numeric_limits<off_t>::max()) ||
        count > static_cast<std::uint64_t>(std::numeric_limits<off_t>::max())-offset)
      throw std::overflow_error("Blocked reference file offset overflow");
    char* cursor=static_cast<char*>(data);
    while (count) {
      const std::size_t chunk=std::min(count,
        static_cast<std::size_t>(std::numeric_limits<ssize_t>::max()));
      const ssize_t done=writing ? ::pwrite(fd_,cursor,chunk,static_cast<off_t>(offset)) :
        ::pread(fd_,cursor,chunk,static_cast<off_t>(offset));
      if (done<0 && errno==EINTR) continue;
      if (done<0) throw ioError(writing?"Reference write failed":"Reference read failed",path_);
      if (done==0) throw std::runtime_error("Short reference "+
        std::string(writing?"write: ":"read: ")+path_);
      cursor+=done; offset+=static_cast<std::size_t>(done);
      count-=static_cast<std::size_t>(done);
      if (counter) *counter+=static_cast<std::uint64_t>(done);
    }
  }
  void read(Double* data,std::size_t count,std::size_t index,
            std::uint64_t& counter) const {
    transfer(data,bytes(count),bytes(index),false,&counter);
  }
  void write(const Double* data,std::size_t count,std::size_t index,
             std::uint64_t& counter) const {
    transfer(const_cast<Double*>(data),bytes(count),bytes(index),true,&counter);
  }
  std::uint64_t size() const {
    struct stat info;
    if (::fstat(fd_,&info)!=0) throw ioError("Cannot stat reference file",path_);
    if (!S_ISREG(info.st_mode)) throw std::runtime_error("Reference data must be a regular file: "+path_);
    if (info.st_size<0) throw std::runtime_error("Negative reference file length");
    return static_cast<std::uint64_t>(info.st_size);
  }
private:
  int fd_;
  bool writable_;
  std::string path_;
};
class PlaneDST {
public:
  PlaneDST(std::size_t nx,std::size_t ny) : data_(product(nx,ny)),plan_(0) {
    if (nx>static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        ny>static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::overflow_error("Blocked reference DST dimensions exceed FFTW int");
    plan_=fftw_plan_r2r_2d(static_cast<int>(ny),static_cast<int>(nx),
      data_.data(),data_.data(),FFTW_RODFT00,FFTW_RODFT00,FFTW_ESTIMATE);
    if (!plan_) throw std::runtime_error("Cannot create blocked transverse DST plan");
  }
  ~PlaneDST() { if(plan_) fftw_destroy_plan(plan_); }
  void transform(Double* plane,Double factor=1.) {
    std::copy(plane,plane+data_.size(),data_.begin());
    fftw_execute(plan_);
    for(std::size_t n=0;n<data_.size();++n) plane[n]=data_[n]*factor;
  }
  std::size_t capacityBytes() const { return bytes(data_.capacity()); }
private:
  std::vector<Double> data_;
  fftw_plan plan_;
};
Double deltaFor(const EBCowanCoefficients& c,unsigned axis) {
  return axis==0?c.deltaYZ:(axis==1?c.deltaZX:c.deltaXY);
}
std::string jsonString(const std::string& value) {
  std::ostringstream out;
  out << '"';
  for (unsigned char c : value) {
    if(c=='"' || c=='\\') out << '\\' << c;
    else if(c<32) out << "\\u" << std::hex << std::setw(4) <<
      std::setfill('0') << static_cast<unsigned>(c) << std::dec;
    else out << c;
  }
  out << '"';
  return out.str();
}
std::string hexInteger(std::uint64_t value) {
  std::ostringstream out;
  out.imbue(std::locale::classic());
  out << std::hex << std::setw(16) << std::setfill('0') << value;
  return out.str();
}
std::string doubleBits(Double value) {
  static_assert(sizeof(Double)==sizeof(std::uint64_t),"Reference metadata requires binary64 Double");
  std::uint64_t bits;
  std::memcpy(&bits,&value,sizeof(bits));
  return hexInteger(bits);
}
std::uint64_t parseHex(const std::string& value) {
  if(value.size()!=16) throw std::runtime_error("Invalid reference metadata hexadecimal field");
  std::uint64_t result=0;
  for(char c:value) {
    unsigned digit;
    if(c>='0' && c<='9') digit=static_cast<unsigned>(c-'0');
    else if(c>='a' && c<='f') digit=static_cast<unsigned>(c-'a'+10);
    else throw std::runtime_error("Invalid reference metadata hexadecimal digit");
    result=(result<<4)|digit;
  }
  return result;
}
Double parseDoubleBits(const std::string& value) {
  const std::uint64_t bits=parseHex(value);
  Double result;
  std::memcpy(&result,&bits,sizeof(result));
  if(!std::isfinite(result)) throw std::runtime_error("Nonfinite reference metadata number");
  return result;
}
std::size_t parseSize(const std::string& value) {
  if(value.empty()) throw std::runtime_error("Empty reference metadata integer");
  std::size_t result=0;
  for(char c:value) {
    if(c<'0' || c>'9') throw std::runtime_error("Invalid reference metadata integer");
    const unsigned digit=static_cast<unsigned>(c-'0');
    if(result>(std::numeric_limits<std::size_t>::max()-digit)/10)
      throw std::overflow_error("Reference metadata integer overflow");
    result=result*10+digit;
  }
  return result;
}
std::uint64_t metadataChecksum(const std::string& value) {
  // Small-manifest accidental-corruption check only; scalar files are still
  // protected by exact lengths, short-read and finite-value checks, not hashes.
  std::uint64_t hash=UINT64_C(14695981039346656037);
  for(unsigned char c:value) { hash^=c; hash*=UINT64_C(1099511628211); }
  return hash;
}
std::string byteOrder() {
  const std::uint32_t value=UINT32_C(0x01020304);
  unsigned char data[sizeof(value)];
  std::memcpy(data,&value,sizeof(value));
  if(data[0]==4 && data[3]==1) return "little";
  if(data[0]==1 && data[3]==4) return "big";
  throw std::runtime_error("Unsupported reference scalar byte order");
}
}

struct BlockedCowanReference::Impl {
  struct Cache {
    std::size_t block=std::numeric_limits<std::size_t>::max();
    std::uint64_t used=0;
    std::vector<Double> values;
  };
  bool quintic=false,linearNotch=false,linearBinomial3=false;
  const bool xBinomial;
  const bool openSides;
  OpenTransverseTrace side;
  Double ellipticStencil[3][3][3]={};
  std::size_t physicalGuards[3]={};
  EBGridGeometry g;
  FieldVector<Double> origin;
  Double beta,translation,upperWeight,lowerWeight,magneticScale;
  EBCowanCoefficients coefficient;
  std::size_t lower[3],upper[3],nodes[3],cells[3];
  std::size_t requestedBlockLayers,blockLayers,area,fullArea;
  std::size_t nx[6],ny[6],nz[6],planeSize[6];
  std::string directory;
  DiskFile rhoFile,checkpointFile,phiFile;
  BlockedCowanReferenceReport stats;
  std::array<Cache,3> phiCache;
  Cache sourceCache;
  std::uint64_t cacheClock=0;
  std::array<std::vector<Double>,6> fields;
  std::size_t fieldStart=0,fieldStop=0;

  Impl(const EBGridGeometry& geometry,const FieldVector<Double>& where,
       const std::vector<RelativisticParticleSI>& particles,Double velocity,
       const std::size_t guards[3],std::size_t block,const std::string& dir,
       MPI_Comm communicator,bool shaped,bool notch,bool binomial3,bool transverseXBinomial3,bool openTransverse,
       bool deferred=false)
    : quintic(shaped),linearNotch(notch),linearBinomial3(binomial3),xBinomial(transverseXBinomial3),openSides(openTransverse),g(geometry),origin(where),beta(velocity),translation(0.),upperWeight(0.),
      lowerWeight(0.),magneticScale(0.),requestedBlockLayers(block),blockLayers(block),
      area(0),fullArea(0),directory(dir) {
    if(static_cast<unsigned>(quintic)+static_cast<unsigned>(linearNotch)+static_cast<unsigned>(linearBinomial3)>1)throw std::invalid_argument("Ambiguous reference source shape");
    configure(guards,communicator);
    if(deferred) return;
    const auto start=Clock::now();
    std::vector<Double> diagonal(area),offDiagonal(area);
    beginPreparation(particles,diagonal,offDiagonal);
    if(openSides) {
      const auto started=Clock::now();
      const auto source=openSources();
      side.build(source,g,1./std::sqrt(1.-translation));
      finishOpenSides(started);
    }
    finishPreparation(diagonal,offDiagonal,start);
  }

  void beginPreparation(const std::vector<RelativisticParticleSI>& particles,
                        std::vector<Double>& diagonal,std::vector<Double>& offDiagonal) {
    // Validate the discrete operator before any filesystem mutation.
    diagonal.resize(area);offDiagonal.resize(area);
    makeModes(diagonal,offDiagonal);
    if (::mkdir(directory.c_str(),0700)!=0)
      throw ioError("Cannot create new blocked reference directory",directory);
    rhoFile.create(directory+"/rho.bin");
    checkpointFile.create(directory+"/checkpoints.bin");
    phiFile.create(directory+"/phi.bin");
    stats.particles=particles.size();
    if(quintic||linearNotch||linearBinomial3)depositQuintic(particles);else deposit(particles);
  }

  void finishPreparation(const std::vector<Double>& diagonal,const std::vector<Double>& offDiagonal,
                         Clock::time_point start) {
    solve(diagonal,offDiagonal);
    validateFileLengths();
    // Publish completion only after data have reached the shared filesystem.
    // Readers are released by the collective caller only after publication.
    rhoFile.synchronize();
    checkpointFile.synchronize();
    phiFile.synchronize();
    publishMetadata();
    stats.preparationSeconds=seconds(start);
    stats.complete=true;
  }

  Impl(const EBGridGeometry& geometry,const FieldVector<Double>& where,
       Double velocity,const std::size_t guards[3],std::size_t block,
       const std::string& dir,MPI_Comm communicator,bool shaped,bool notch,bool binomial3,bool transverseXBinomial3,bool openTransverse)
    : quintic(shaped),linearNotch(notch),linearBinomial3(binomial3),xBinomial(transverseXBinomial3),openSides(openTransverse),g(geometry),origin(where),beta(velocity),translation(0.),upperWeight(0.),
      lowerWeight(0.),magneticScale(0.),requestedBlockLayers(block),blockLayers(block),
      area(0),fullArea(0),directory(dir) {
    const auto start=Clock::now();
    if(static_cast<unsigned>(quintic)+static_cast<unsigned>(linearNotch)+static_cast<unsigned>(linearBinomial3)>1)throw std::invalid_argument("Ambiguous reference source shape");
    configure(guards,communicator);
    readMetadata();
    rhoFile.openReadOnly(directory+"/rho.bin");
    checkpointFile.openReadOnly(directory+"/checkpoints.bin");
    phiFile.openReadOnly(directory+"/phi.bin");
    if(openSides) {
      DiskFile file;
      file.openReadOnly(directory+"/side-potential.bin");
      if(file.size()!=bytes(side.values.size())) throw std::runtime_error("Wrong open reference boundary file length");
      file.read(side.values.data(),side.values.size(),0,stats.bytesRead);
      for(double value:side.values) if(!std::isfinite(value)) throw std::runtime_error("Nonfinite prepared wall potential");
    }
    validateFileLengths();
    stats.preparedReader=true;
    stats.openSeconds=seconds(start);
    stats.complete=true;
    // All preparation timers, source indices, solver workspace and binary I/O
    // counters stay zero: metadata history is not work performed by this rank.
  }

  void configure(const std::size_t guards[3],MPI_Comm communicator) {
    if(xBinomial && !linearBinomial3)
      throw std::invalid_argument("Transverse three-point reference requires longitudinal linear-binomial3");
    if(!guards) throw std::invalid_argument("Missing reference guards");
    for(unsigned d=0;d<3;++d) physicalGuards[d]=guards[d];
    if(openSides) {
      if(guards[2]!=0 || guards[0]==0 || guards[1]==0)
        throw std::invalid_argument("Open transverse reference requires only transverse CPML");
      for(auto n:{g.nx,g.ny,g.nz}) if(n>=std::size_t(std::numeric_limits<int>::max()/4))
        throw std::overflow_error("Open transverse reference dimensions exceed supported indices");
    }
    const std::size_t solveGuards[3]={openSides?0:guards[0],openSides?0:guards[1],guards[2]};
    validate(solveGuards,communicator);
    if(openSides) {side.configure(g);stats.openTraceBytes=bytes(side.values.size());}
    coefficient=EBCowanCoefficients::forZDispersion(g.dx,g.dy,g.dz);
    const Double speed=beta*SI::c;
    translation=speed*speed*SI::epsilon0*SI::mu0;
    upperWeight=0.5+0.5*beta*SI::c*g.dt/g.dz;
    lowerWeight=1.-upperWeight;
    magneticScale=beta*SI::c*SI::epsilon0*SI::mu0;
    if (!(upperWeight>=0. && upperWeight<=1.))
      throw std::invalid_argument("Blocked initial B half-step shift exceeds one cell");
  }

  void validateFileLengths() {
    stats.sourceFileBytes=rhoFile.size();
    stats.checkpointFileBytes=checkpointFile.size();
    stats.potentialFileBytes=phiFile.size();
    if(stats.sourceFileBytes!=bytes(product(area,nodes[2])) ||
       stats.potentialFileBytes!=stats.sourceFileBytes ||
       stats.checkpointFileBytes!=bytes(product(product(stats.blocks,area),2)))
      throw std::runtime_error("Blocked reference prepared file has incorrect length");
  }

  void validate(const std::size_t guards[3],MPI_Comm communicator) {
    int initialized=0,finished=0,size=0;
    MPI_Initialized(&initialized);
    if(initialized) MPI_Finalized(&finished);
    if(!initialized || finished)
      throw std::invalid_argument("Blocked reference requires live single-rank MPI");
    if(communicator==MPI_COMM_NULL)
      throw std::invalid_argument("Blocked reference requires a non-null singleton communicator");
    MPI_Comm_size(communicator,&size);
    if(size!=1) throw std::invalid_argument("Blocked reference communicator must contain one MPI rank");
    if(g.solver!=EBMaxwellSolver::CowanZ || !std::isfinite(beta) || std::abs(beta)>=1. ||
       blockLayers==0 || !guards || directory.empty())
      throw std::invalid_argument("Invalid blocked Cowan reference configuration");
    cells[0]=g.nx; cells[1]=g.ny; cells[2]=g.nz;
    const Double h[3]={g.dx,g.dy,g.dz};
    for(unsigned a=0;a<3;++a) {
      if(!std::isfinite(h[a]) || h[a]<=0. || !std::isfinite(origin[a]) ||
         cells[a]<3 || guards[a]>(cells[a]-3)/2 ||
         cells[a]>=static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max())-8)
        throw std::invalid_argument("Invalid blocked reference geometry or guard/padding");
      lower[a]=guards[a]; upper[a]=cells[a]-guards[a];
      nodes[a]=upper[a]-lower[a]-1;
      stats.interiorNodes[a]=nodes[a];
    }
    if(!std::isfinite(g.dt) || g.dt<=0. ||
       !std::isfinite(g.dx*g.dy*g.dz) || g.dx*g.dy*g.dz<=0.)
      throw std::invalid_argument("Invalid blocked reference time step or cell volume");
    blockLayers=std::min(blockLayers,nodes[2]);
    stats.blockLayers=blockLayers;
    stats.blocks=(nodes[2]-1)/blockLayers+1;
    stats.beta=beta;
    area=product(nodes[0],nodes[1]);
    fullArea=product(plus(g.nx,1),plus(g.ny,1));
    (void)bytes(product(area,nodes[2]));
    nx[0]=g.nx; nx[1]=g.nx+1; nx[2]=g.nx+1;
    nx[3]=g.nx+1; nx[4]=g.nx; nx[5]=g.nx;
    ny[0]=g.ny+1; ny[1]=g.ny; ny[2]=g.ny+1;
    ny[3]=g.ny; ny[4]=g.ny+1; ny[5]=g.ny;
    for(unsigned c=0;c<6;++c) {
      nz[c]=g.nz+((c==0 || c==1 || c==5)?1:0);
      planeSize[c]=product(nx[c],ny[c]);
    }
  }

  std::string metadataContract() const {
    if(!std::numeric_limits<Double>::is_iec559)
      throw std::runtime_error("Prepared reference requires IEEE-754 binary64");
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << "APRL_BLOCKED_COWAN_REFERENCE\nversion 1\n"
      << "algorithm cowan-z-dst2-blocked-thomas-v1\n"
      << "scalar_layout native-binary64-interior-x-fast-y-z\n"
      << "scalar_bytes " << sizeof(Double) << "\nsize_t_bytes " << sizeof(std::size_t)
      << "\nbyte_order " << byteOrder() << "\nsolver CowanZ\n"
      << "cells " << g.nx << ' ' << g.ny << ' ' << g.nz << '\n'
      << "spacing_bits " << doubleBits(g.dx) << ' ' << doubleBits(g.dy) << ' ' << doubleBits(g.dz) << '\n'
      << "dt_bits " << doubleBits(g.dt) << '\n'
      << "origin_bits " << doubleBits(origin[0]) << ' ' << doubleBits(origin[1]) << ' ' << doubleBits(origin[2]) << '\n'
      << "beta_bits " << doubleBits(beta) << '\n'
      << "guards " << lower[0] << ' ' << lower[1] << ' ' << lower[2] << '\n'
      << "requested_block_layers " << requestedBlockLayers << '\n'
      << "actual_block_layers " << blockLayers << '\n'
      << "blocks " << stats.blocks << '\n'
      << "interior_nodes " << nodes[0] << ' ' << nodes[1] << ' ' << nodes[2] << '\n';
    out << "cowan_alpha_bits";
    for(unsigned a=0;a<3;++a) out << ' ' << doubleBits(coefficient.alpha[a]);
    out << "\ncowan_beta_bits";
    for(unsigned a=0;a<3;++a) out << ' ' << doubleBits(coefficient.beta[a]);
    out << "\ncowan_ratio_bits";
    for(unsigned a=0;a<3;++a) out << ' ' << doubleBits(coefficient.ratio[a]);
    out << "\ncowan_delta_bits " << doubleBits(coefficient.deltaXY) << ' '
      << doubleBits(coefficient.deltaYZ) << ' ' << doubleBits(coefficient.deltaZX)
      << "\nsi_bits " << doubleBits(SI::c) << ' ' << doubleBits(SI::epsilon0) << ' ' << doubleBits(SI::mu0)
      << "\ntranslation_bits " << doubleBits(translation)
      << "\nmagnetic_interpolation_bits " << doubleBits(lowerWeight) << ' '
      << doubleBits(upperWeight) << ' ' << doubleBits(magneticScale)
      << "\nsource_bytes " << bytes(product(area,nodes[2]))
      << "\ncheckpoint_bytes " << bytes(product(product(stats.blocks,area),2))
      << "\npotential_bytes " << bytes(product(area,nodes[2]));
    if(openSides) out<<"\nopen_transverse free-rho-coulomb-outer-ghost-d7-v1\nphysical_cpml "
                     <<physicalGuards[0]<<' '<<physicalGuards[1]<<' '<<physicalGuards[2];
    if(xBinomial)out<<"\ntransverse_x_source linear-binomial3-cos2-three-point-v1";
    out<<"\nEND_CONTRACT\n";
    if(quintic)out<<"longitudinal_source quintic-continuum-seven-point-v1\n";
    if(linearNotch)out<<"longitudinal_source linear-notch-sin6-seven-point-v1\n";
    if(linearBinomial3)out<<"longitudinal_source linear-binomial3-cos2-three-point-v1\n";
    return out.str();
  }

  void publishMetadata() const {
    if(!std::isfinite(stats.sourceCharge) || !std::isfinite(stats.depositedCharge) ||
       !std::isfinite(stats.minimumPivot) || !std::isfinite(stats.minimumModeEigenvalue) ||
       !(stats.minimumPivot>0.) || !(stats.minimumModeEigenvalue>0.))
      throw std::runtime_error("Cannot publish nonfinite or incomplete prepared reference history");
    std::ostringstream payload;
    payload.imbue(std::locale::classic());
    payload << metadataContract()
      << "particles " << stats.particles
      << "\nsource_charge_bits " << doubleBits(stats.sourceCharge)
      << "\ndeposited_charge_bits " << doubleBits(stats.depositedCharge)
      << "\nminimum_pivot_bits " << doubleBits(stats.minimumPivot)
      << "\nminimum_mode_eigenvalue_bits " << doubleBits(stats.minimumModeEigenvalue)
      << "\nstatus complete\nEND\n";
    const std::string body=payload.str();
    const std::string serialized=body+"checksum "+hexInteger(metadataChecksum(body))+"\n";
    const std::string pending=directory+"/.metadata.pending",complete=directory+"/metadata.txt";
    DiskFile out;
    out.create(pending);
    out.transfer(const_cast<char*>(serialized.data()),serialized.size(),0,true,0);
    out.synchronize();
    // link publishes the complete same-directory inode atomically and refuses
    // to replace an existing marker. A leftover pending name is never accepted
    // as completion, and cleanup failure cannot invalidate a published marker.
    if(::link(pending.c_str(),complete.c_str())!=0)
      throw ioError("Cannot publish completed reference metadata",complete);
    (void)::unlink(pending.c_str());
  }

  void readMetadata() {
    DiskFile metadata;
    metadata.openReadOnly(directory+"/metadata.txt");
    const std::uint64_t length=metadata.size();
    if(length==0 || length>16384)
      throw std::runtime_error("Incomplete or oversized prepared reference metadata");
    std::string serialized(static_cast<std::size_t>(length),'\0');
    metadata.transfer(&serialized[0],serialized.size(),0,false,0);
    const std::size_t seal=serialized.rfind("checksum ");
    if(seal==std::string::npos || seal==0 || serialized[seal-1]!='\n' ||
       serialized.size()!=seal+9+16+1 || serialized.back()!='\n')
      throw std::runtime_error("Incomplete prepared reference metadata seal");
    const std::string body=serialized.substr(0,seal);
    if(parseHex(serialized.substr(seal+9,16))!=metadataChecksum(body))
      throw std::runtime_error("Prepared reference metadata checksum mismatch");
    const std::string contract=metadataContract();
    if(body.compare(0,contract.size(),contract)!=0)
      throw std::invalid_argument("Prepared reference geometry/beta/guards/block/format contract mismatch");
    std::istringstream payload(body.substr(contract.size()));
    payload.imbue(std::locale::classic());
    const auto field=[&payload](const char* expected)->std::string {
      std::string key,value;
      if(!(payload>>key>>value) || key!=expected)
        throw std::runtime_error("Invalid prepared reference metadata field order");
      return value;
    };
    stats.particles=parseSize(field("particles"));
    stats.sourceCharge=parseDoubleBits(field("source_charge_bits"));
    stats.depositedCharge=parseDoubleBits(field("deposited_charge_bits"));
    stats.minimumPivot=parseDoubleBits(field("minimum_pivot_bits"));
    stats.minimumModeEigenvalue=parseDoubleBits(field("minimum_mode_eigenvalue_bits"));
    if(!(stats.minimumPivot>0.) || !(stats.minimumModeEigenvalue>0.) || field("status")!="complete")
      throw std::runtime_error("Prepared reference metadata does not describe a complete positive operator");
    std::string end,extra;
    if(!(payload>>end) || end!="END" || (payload>>extra))
      throw std::runtime_error("Invalid prepared reference metadata terminator");
  }

  void makeModes(std::vector<Double>& diagonal,std::vector<Double>& off) {
    Double stencil[3][3][3]={};
    const Double spacing[3]={g.dx,g.dy,g.dz};
    for(unsigned axis=0;axis<3;++axis) {
      const unsigned ta=(axis+1)%3,tb=(axis+2)%3;
      const Double delta=deltaFor(coefficient,axis);
      // Individual directional terms need not be positive for their sum to
      // be SPD. The finite Dirichlet operator is diagonalized in x/y below;
      // every resulting symmetric z tridiagonal has constant diagonal/off-
      // diagonal entries. Its exact smallest eigenvalue is checked before
      // any file is created. Keep that full-operator test for coarse grids.
      const Double inverse=1./(spacing[axis]*spacing[axis]);
      for(int a=-1;a<=1;++a) for(int b=-1;b<=1;++b) {
        const Double weight=a==0 ? (b==0?coefficient.alpha[axis]:coefficient.beta[tb]) :
          (b==0?coefficient.beta[ta]:delta);
        int p[3]={1,1,1}; p[ta]+=a; p[tb]+=b;
        stencil[p[2]][p[1]][p[0]]+=2.*weight*inverse;
        --p[axis]; stencil[p[2]][p[1]][p[0]]-=weight*inverse;
        p[axis]+=2; stencil[p[2]][p[1]][p[0]]-=weight*inverse;
      }
    }
    const Double translationZ=translation/(g.dz*g.dz);
    stencil[1][1][1]-=2.*translationZ;
    stencil[0][1][1]+=translationZ; stencil[2][1][1]+=translationZ;
    for(int z=0;z<3;++z)for(int y=0;y<3;++y)for(int x=0;x<3;++x)
      ellipticStencil[z][y][x]=stencil[z][y][x];
    stats.minimumModeEigenvalue=std::numeric_limits<Double>::infinity();
    const Double pi=std::acos(-1.);
    for(std::size_t j=0;j<nodes[1];++j) for(std::size_t i=0;i<nodes[0];++i) {
      const Double tx=pi*(i+1)/(nodes[0]+1),ty=pi*(j+1)/(nodes[1]+1);
      Double symbol[3]={};
      for(unsigned z=0;z<3;++z) for(unsigned y=0;y<3;++y) for(unsigned x=0;x<3;++x)
        symbol[z]+=stencil[z][y][x]*std::cos((static_cast<int>(x)-1)*tx)*
          std::cos((static_cast<int>(y)-1)*ty);
      const Double scale=std::max(std::abs(symbol[1]),std::max(std::abs(symbol[0]),std::abs(symbol[2])));
      if(!std::isfinite(scale) || scale==0. ||
         std::abs(symbol[0]-symbol[2])>128.*std::numeric_limits<Double>::epsilon()*scale)
        throw std::runtime_error("Blocked Cowan modal matrix is not finite symmetric");
      const std::size_t m=j*nodes[0]+i;
      diagonal[m]=symbol[1]; off[m]=0.5*(symbol[0]+symbol[2]);
      const Double eigenvalue=diagonal[m]-2.*std::abs(off[m])*std::cos(pi/(nodes[2]+1));
      if(!(eigenvalue>0.) || !std::isfinite(eigenvalue))
        throw std::runtime_error("Blocked Cowan modal matrix is not positive definite");
      stats.minimumModeEigenvalue=std::min(stats.minimumModeEigenvalue,eigenvalue);
    }
  }

  bool interior(std::size_t i,std::size_t j,std::size_t k) const {
    return i>lower[0] && i<upper[0] && j>lower[1] && j<upper[1] &&
      k>lower[2] && k<upper[2];
  }
  std::size_t layers(std::size_t block) const {
    return std::min(blockLayers,nodes[2]-block*blockLayers);
  }
  void depositQuintic(const std::vector<RelativisticParticleSI>& particles) {
    const auto start=Clock::now();if(blockLayers<8)throw std::invalid_argument("Quintic reference requires block_layers >= 8");
    std::vector<std::size_t> offsets(stats.blocks+1,0),order(particles.size());long double source=0.;
    const auto bucket=[&](const RelativisticParticleSI& p){const double z=((p.position[2]-origin[2])/g.dz-lower[2]-1)/blockLayers;if(!std::isfinite(z))throw std::invalid_argument("Nonfinite quintic reference source");return static_cast<std::size_t>(std::max(0.,std::min(static_cast<double>(stats.blocks-1),std::floor(z))));};
    for(const auto& p:particles){if(!std::isfinite(p.charge))throw std::invalid_argument("Nonfinite quintic reference charge");source+=p.charge;++offsets[bucket(p)+1];}
    for(std::size_t b=1;b<offsets.size();++b)offsets[b]+=offsets[b-1];auto next=offsets;
    for(std::size_t n=0;n<particles.size();++n)order[next[bucket(particles[n])]++]=n;
    std::vector<double> block(product(blockLayers,area));long double deposited=0.;
    for(std::size_t b=0;b<stats.blocks;++b){const auto count=product(layers(b),area);std::fill(block.begin(),block.begin()+count,0.);
      const auto first=b?b-1:0,last=std::min(stats.blocks,b+2);
      for(std::size_t n=offsets[first];n<offsets[last];++n){const auto& p=particles[order[n]];const double x=(p.position[0]-origin[0])/g.dx,y=(p.position[1]-origin[1])/g.dy,z=(p.position[2]-origin[2])/g.dz;const long long ix=static_cast<long long>(std::floor(x)),iy=static_cast<long long>(std::floor(y));const double fx=x-ix,fy=y-iy,scale=p.charge/(g.dx*g.dy*g.dz);const auto weights=longitudinal_shape::basis(z,false,linearBinomial3?longitudinal_shape::Shape::LinearBinomial3:(linearNotch?longitudinal_shape::Shape::LinearNotch:longitudinal_shape::Shape::QuinticCompensated));
       const auto wx=xBinomial ? transverse_source_shape::nodeWeights(fx) :
                                  std::array<double,4>{{1.-fx,fx,0.,0.}};
       const long long firstX=ix-(xBinomial?1:0);
       const unsigned nX=xBinomial?4:2;
       for(unsigned c=0;c<weights.count;++c){
         const auto k=weights.first+c;
         if(weights.v[c]==0.)continue;
         for(unsigned j=0;j<2;++j)for(unsigned i=0;i<nX;++i){
           const double value=scale*wx[i]*(j?fy:1.-fy)*weights.v[c];
           if(value==0.)continue;
           const auto nodeX=firstX+i,nodeY=iy+j;
           if(k<0 || nodeX<0 || nodeY<0 || !interior(nodeX,nodeY,k))
             throw std::invalid_argument("Reference charge support touches Dirichlet surface");
           const std::size_t kk=k-lower[2]-1;
           if(kk/blockLayers==b)
             block[(kk%blockLayers)*area+(nodeY-lower[1]-1)*nodes[0]+nodeX-lower[0]-1]+=value;
         }
       }
      }
      for(std::size_t j=0;j<count;++j){if(!std::isfinite(block[j]))throw std::runtime_error("Nonfinite quintic reference density");deposited+=(long double)block[j]*g.dx*g.dy*g.dz;}
      rhoFile.write(block.data(),count,product(b*blockLayers,area),stats.bytesWritten);
    }
    stats.sourceCharge=source;stats.depositedCharge=deposited;stats.sourceIndexBytes=bytes(order.capacity()+offsets.capacity()+next.capacity());stats.sourceSortingAllowanceBytes=0;stats.solverWorkspaceBytes=bytes(block.capacity());stats.sourceSeconds=seconds(start);
  }

  void deposit(const std::vector<RelativisticParticleSI>& particles) {
    const auto started=Clock::now();
    struct Contribution { std::size_t block,node; Double value; };
    std::vector<Contribution> contributions;
    contributions.reserve(product(particles.size(),8));
    const Double h[3]={g.dx,g.dy,g.dz};
    const Double volume=g.dx*g.dy*g.dz,inverseVolume=1./volume;
    long double sourceCharge=0.;
    for(const auto& p:particles) {
      if(!std::isfinite(p.charge)) throw std::invalid_argument("Nonfinite blocked source charge");
      sourceCharge+=p.charge;
      std::size_t cell[3]; Double fraction[3];
      for(unsigned a=0;a<3;++a) {
        const Double relative=p.position[a]-origin[a],extent=cells[a]*h[a];
        const Double tolerance=64.*std::numeric_limits<Double>::epsilon()*
          std::max(1.,std::max(std::abs(origin[a]),extent));
        if(!std::isfinite(p.position[a]) || relative < -tolerance || relative>extent+tolerance)
          throw std::out_of_range("Blocked source particle outside global E/B domain");
        const Double normalized=std::max(0.,std::min(static_cast<Double>(cells[a]),relative/h[a]));
        cell[a]=normalized<=0.?0:(normalized>=static_cast<Double>(cells[a]) ?
          cells[a]-1:static_cast<std::size_t>(std::floor(normalized)));
        fraction[a]=std::max(0.,std::min(1.,normalized-static_cast<Double>(cell[a])));
      }
      // Same global normalization, operation order and eight-corner order as
      // ChargeConservingCurrentDepositor::depositCharge. Stable block sorting
      // below preserves the original particle accumulation order at every node.
      for(unsigned c=0;c<2;++c) for(unsigned b=0;b<2;++b) for(unsigned a=0;a<2;++a) {
        const Double value=p.charge*inverseVolume*(a?fraction[0]:1.-fraction[0])*
          (b?fraction[1]:1.-fraction[1])*(c?fraction[2]:1.-fraction[2]);
        if(!std::isfinite(value)) throw std::invalid_argument("Nonfinite blocked CIC contribution");
        if(value==0.) continue;
        const std::size_t i=cell[0]+a,j=cell[1]+b,k=cell[2]+c;
        if(!interior(i,j,k))
          throw std::invalid_argument("Blocked source CIC support touches Dirichlet surface; padding unsupported");
        const std::size_t z=k-lower[2]-1;
        contributions.push_back({z/blockLayers,(z%blockLayers)*area+
          (j-lower[1]-1)*nodes[0]+i-lower[0]-1,value});
      }
    }
    stats.sourceIndexBytes=product(contributions.capacity(),sizeof(Contribution));
    // Conservative allowance for implementation-dependent stable_sort scratch;
    // this is array accounting, not measured process RSS.
    stats.sourceSortingAllowanceBytes=stats.sourceIndexBytes;
    std::stable_sort(contributions.begin(),contributions.end(),
      [](const Contribution& a,const Contribution& b){return a.block<b.block;});
    std::vector<Double> block(product(blockLayers,area));
    std::size_t current=0;
    long double deposited=0.;
    for(std::size_t b=0;b<stats.blocks;++b) {
      const std::size_t count=product(layers(b),area);
      std::fill(block.begin(),block.begin()+count,0.);
      while(current<contributions.size() && contributions[current].block==b) {
        const auto& c=contributions[current++]; block[c.node]+=c.value;
      }
      for(std::size_t q=0;q<count;++q) {
        if(!std::isfinite(block[q])) throw std::runtime_error("Nonfinite blocked deposited charge");
        deposited+=static_cast<long double>(block[q])*volume;
      }
      rhoFile.write(block.data(),count,product(b*blockLayers,area),stats.bytesWritten);
    }
    if(current!=contributions.size()) throw std::logic_error("Unassigned blocked CIC contribution");
    stats.sourceCharge=static_cast<Double>(sourceCharge);
    stats.depositedCharge=static_cast<Double>(deposited);
    stats.solverWorkspaceBytes=bytes(block.capacity());
    stats.sourceSeconds=seconds(started);
  }

  OpenTransverseTrace::Sources openSources() {
    OpenTransverseTrace::Sources source;
    const double gamma=1./std::sqrt(1.-translation);
    const double scale=gamma*g.dx*g.dy*g.dz/(4.*std::acos(-1.)*SI::epsilon0);
    std::vector<double> rho(product(blockLayers,area));
    for(std::size_t b=0;b<stats.blocks;++b) {
      const std::size_t count=layers(b);
      rhoFile.read(rho.data(),product(count,area),b*blockLayers*area,stats.bytesRead);
      for(std::size_t p=0;p<count;++p) for(std::size_t j=0;j<nodes[1];++j) for(std::size_t i=0;i<nodes[0];++i) {
        const double density=rho[p*area+j*nodes[0]+i]; if(density==0.)continue;
        if(i+1<=physicalGuards[0] || i+1>=g.nx-physicalGuards[0] ||
           j+1<=physicalGuards[1] || j+1>=g.ny-physicalGuards[1])
          throw std::runtime_error("Open reference charge overlaps original deposition guard");
        source.x.push_back((i+1)*g.dx);source.y.push_back((j+1)*g.dy);
        source.z.push_back((b*blockLayers+p+1)*g.dz*gamma);source.q.push_back(density*scale);
      }
    }
    return source;
  }

  void finishOpenSides(Clock::time_point started) {
    stats.openTraceChecks=side.directChecks;stats.openTraceMaxRelativeError=side.maxRelativeError;
    stats.openTraceSourceNodes=side.sourceNodes;stats.openTraceCoarsePoints=side.coarsePoints;
    stats.openTraceFallbackTiles=side.fallbackTiles;stats.openTraceRanks=side.computeRanks;
    stats.openTraceThreadsPerRank=computeThreadCount();
    stats.openTraceMinRankEvaluations=stats.openTraceMaxRankEvaluations=side.localPotentialEvaluations;
    DiskFile file; file.create(directory+"/side-potential.bin");
    file.write(side.values.data(),side.values.size(),0,stats.bytesWritten);file.synchronize();
    stats.openTraceSeconds=seconds(started);
  }

  // A phi = rho/eps. Known nodal transverse-face contributions enter the
  // original sine-transform/Thomas solver as -A(interior,boundary)*phi_boundary.
  void liftOpenSides(Double* rhs,int k) const {
    const int nx=static_cast<int>(g.nx),ny=static_cast<int>(g.ny);
    const auto lift=[&](int i,int j) {
      double value=0.;
      for(int dz=-1;dz<=1;++dz)for(int dy=-1;dy<=1;++dy)for(int dx=-1;dx<=1;++dx) {
        const int x=i+dx,y=j+dy;
        if(x==0 || x==nx || y==0 || y==ny) value+=ellipticStencil[dz+1][dy+1][dx+1]*side.at(x,y,k+dz);
      }
      rhs[(j-1)*nodes[0]+i-1]-=value;
    };
    for(int i=1;i<nx;++i){lift(i,1);if(ny>2)lift(i,ny-1);}
    for(int j=2;j<ny-1;++j){lift(1,j);if(nx>2)lift(nx-1,j);}
  }

  void solve(const std::vector<Double>& diagonal,const std::vector<Double>& off) {
    PlaneDST transform(nodes[0],nodes[1]);
    std::vector<Double> rhs(product(blockLayers,area)),cprev(area,0.),dprev(area,0.);
    std::vector<Double> cprime(rhs.size()),dprime(rhs.size()),next(area,0.);
    stats.solverWorkspaceBytes=bytes(rhs.capacity()+cprev.capacity()+dprev.capacity()+
      cprime.capacity()+dprime.capacity()+next.capacity()+diagonal.capacity()+off.capacity())+
      transform.capacityBytes();
    stats.minimumPivot=std::numeric_limits<Double>::infinity();
    const auto forwardStart=Clock::now();
    for(std::size_t b=0;b<stats.blocks;++b) {
      checkpointFile.write(cprev.data(),area,product(product(b,area),2),stats.bytesWritten);
      checkpointFile.write(dprev.data(),area,product(product(b,area),2)+area,stats.bytesWritten);
      loadModes(b,rhs,transform);
      for(std::size_t k=0;k<layers(b);++k) for(std::size_t m=0;m<area;++m) {
        const Double pivot=diagonal[m]-off[m]*cprev[m];
        checkPivot(pivot);
        cprev[m]=off[m]/pivot;
        dprev[m]=(rhs[k*area+m]-off[m]*dprev[m])/pivot;
        if(!std::isfinite(dprev[m])) throw std::runtime_error("Nonfinite blocked forward summary");
      }
    }
    stats.forwardSeconds=seconds(forwardStart);
    const auto backwardStart=Clock::now();
    const Double normalization=1./(4.*(nodes[0]+1)*(nodes[1]+1));
    for(std::size_t bb=stats.blocks;bb>0;--bb) {
      const std::size_t b=bb-1,count=layers(b);
      checkpointFile.read(cprev.data(),area,product(product(b,area),2),stats.bytesRead);
      checkpointFile.read(dprev.data(),area,product(product(b,area),2)+area,stats.bytesRead);
      loadModes(b,rhs,transform);
      for(std::size_t k=0;k<count;++k) for(std::size_t m=0;m<area;++m) {
        const Double pivot=diagonal[m]-off[m]*cprev[m];
        checkPivot(pivot);
        cprev[m]=off[m]/pivot;
        dprev[m]=(rhs[k*area+m]-off[m]*dprev[m])/pivot;
        cprime[k*area+m]=cprev[m]; dprime[k*area+m]=dprev[m];
      }
      for(std::size_t kk=count;kk>0;--kk) for(std::size_t m=0;m<area;++m) {
        const std::size_t q=(kk-1)*area+m;
        dprime[q]-=cprime[q]*next[m]; next[m]=dprime[q];
      }
      for(std::size_t k=0;k<count;++k) transform.transform(dprime.data()+k*area,normalization);
      for(std::size_t q=0;q<count*area;++q)
        if(!std::isfinite(dprime[q])) throw std::runtime_error("Nonfinite reconstructed blocked potential");
      phiFile.write(dprime.data(),product(count,area),product(b*blockLayers,area),stats.bytesWritten);
    }
    stats.reconstructionSeconds=seconds(backwardStart);
  }
  void checkPivot(Double pivot) {
    if(!std::isfinite(pivot) || !(pivot>0.))
      throw std::runtime_error("Blocked Thomas matrix has nonpositive/nonfinite pivot");
    stats.minimumPivot=std::min(stats.minimumPivot,pivot);
  }
  void loadModes(std::size_t block,std::vector<Double>& rhs,PlaneDST& transform) {
    const std::size_t count=product(layers(block),area);
    rhoFile.read(rhs.data(),count,product(block*blockLayers,area),stats.bytesRead);
    for(std::size_t q=0;q<count;++q) rhs[q]/=SI::epsilon0;
    for(std::size_t k=0;k<layers(block);++k) {
      if(openSides) liftOpenSides(rhs.data()+k*area,static_cast<int>(block*blockLayers+k+1));
      transform.transform(rhs.data()+k*area);
    }
  }

  void memoryReport() {
    std::size_t potentialBytes=0,fieldBytes=0;
    for(const auto& p:phiCache) potentialBytes+=bytes(p.values.capacity());
    for(const auto& f:fields) fieldBytes+=bytes(f.capacity());
    const std::size_t sourceBytes=bytes(sourceCache.values.capacity());
    stats.potentialCacheCapacityBytes=std::max(stats.potentialCacheCapacityBytes,potentialBytes);
    stats.sourceCacheCapacityBytes=std::max(stats.sourceCacheCapacityBytes,sourceBytes);
    stats.fieldCacheCapacityBytes=std::max(stats.fieldCacheCapacityBytes,fieldBytes);
    stats.residentArrayBytes=potentialBytes+sourceBytes+fieldBytes+stats.openTraceBytes;
  }
  const std::vector<Double>& potentialBlock(std::size_t block) {
    for(auto& p:phiCache) if(p.block==block) {
      p.used=++cacheClock; ++stats.potentialCacheHits; return p.values;
    }
    auto* chosen=&phiCache[0];
    for(auto& p:phiCache)
      if(p.block==std::numeric_limits<std::size_t>::max() || p.used<chosen->used) chosen=&p;
    chosen->block=std::numeric_limits<std::size_t>::max();
    chosen->values.resize(product(layers(block),area));
    phiFile.read(chosen->values.data(),chosen->values.size(),product(block*blockLayers,area),stats.bytesRead);
    for(Double v:chosen->values) if(!std::isfinite(v))
      throw std::runtime_error("Nonfinite potential read from blocked reference disk");
    chosen->block=block; chosen->used=++cacheClock; ++stats.potentialCacheMisses;
    memoryReport();
    return chosen->values;
  }
  void vertexBounds(std::size_t i,std::size_t j,std::size_t k) const {
    if(i>g.nx || j>g.ny || k>g.nz)
      throw std::out_of_range("Blocked scalar query outside global vertex grid");
  }
  Double scalar(std::size_t i,std::size_t j,std::size_t k,bool source) {
    vertexBounds(i,j,k);
    if(!interior(i,j,k)) {
      if(!source && openSides && (i==0 || j==0 || i==g.nx || j==g.ny))
        return side.at(static_cast<int>(i),static_cast<int>(j),static_cast<int>(k));
      return 0.;
    }
    const std::size_t z=k-lower[2]-1,b=z/blockLayers;
    const std::size_t q=(z%blockLayers)*area+(j-lower[1]-1)*nodes[0]+i-lower[0]-1;
    if(!source) return potentialBlock(b)[q];
    if(sourceCache.block!=b) {
      sourceCache.block=std::numeric_limits<std::size_t>::max();
      sourceCache.values.resize(product(layers(b),area));
      rhoFile.read(sourceCache.values.data(),sourceCache.values.size(),product(b*blockLayers,area),stats.bytesRead);
      for(Double v:sourceCache.values) if(!std::isfinite(v))
        throw std::runtime_error("Nonfinite source read from blocked reference disk");
      sourceCache.block=b; memoryReport();
    }
    return sourceCache.values[q];
  }

  Double smooth(const std::vector<Double>& phi,unsigned axis,
                std::size_t i,std::size_t j,std::size_t p) const {
    const unsigned ta=(axis+1)%3,tb=(axis+2)%3;
    const std::size_t stride[3]={1,g.nx+1,fullArea};
    const std::size_t q=p*fullArea+j*(g.nx+1)+i;
    Double value=coefficient.alpha[axis]*phi[q];
    for(int sign=-1;sign<=1;sign+=2) {
      const std::size_t qa=sign<0?q-stride[ta]:q+stride[ta];
      const std::size_t qb=sign<0?q-stride[tb]:q+stride[tb];
      value+=coefficient.beta[ta]*phi[qa];
      value+=coefficient.beta[tb]*phi[qb];
    }
    const Double delta=deltaFor(coefficient,axis);
    for(int a=-1;a<=1;a+=2) for(int b=-1;b<=1;b+=2) {
      std::size_t s=a<0?q-stride[ta]:q+stride[ta];
      s=b<0?s-stride[tb]:s+stride[tb];
      value+=delta*phi[s];
    }
    return value-(axis==2?translation*phi[q]:0.);
  }
  // Evaluate the same normal-difference transverse smooth even on the outer
  // transverse faces. All physical indices, including the single external
  // ghost, use one coordinate system. Longitudinal boundaries retain zero.
  Double smoothOnWall(const std::vector<Double>& phi,unsigned axis,int i,int j,
                      std::size_t p,int k) const {
    const unsigned ta=(axis+1)%3,tb=(axis+2)%3;
    double value=0.;
    for(int a=-1;a<=1;++a) for(int b=-1;b<=1;++b) {
      int d[3]={0,0,0}; d[ta]=a; d[tb]=b;
      const int x=i+d[0],y=j+d[1],z=k+d[2];
      const auto row=static_cast<std::size_t>(static_cast<std::ptrdiff_t>(p)+d[2]);
      const double v=(x<=0 || y<=0 || x>=int(g.nx) || y>=int(g.ny)) ?
          side.at(x,y,z) : phi[row*fullArea+std::size_t(y)*(g.nx+1)+x];
      const double w=a==0?(b==0?coefficient.alpha[axis]:coefficient.beta[tb]):
                           (b==0?coefficient.beta[ta]:deltaFor(coefficient,axis));
      value+=w*v;
    }
    return value-(axis==2?translation*side.at(i,j,k):0.);
  }

  void ensureField(std::size_t requestedK) {
    if(requestedK>=fieldStart && requestedK<fieldStop) return;
    const auto started=Clock::now();
    // One extra complete native E/B plane prevents per-pixel cache thrashing
    // when a linear interpolation straddles a nominal block boundary.
    const std::size_t start=(requestedK/blockLayers)*blockLayers;
    const std::size_t stop=start+std::min(g.nz+1-start,blockLayers+1);
    const std::size_t count=stop-start;
    fieldStop=fieldStart=0; // Invalidate before reads/allocations that may fail.
    std::vector<Double> phi(product(count+3,fullArea),0.);
    std::array<std::vector<Double>,3> filtered;
    for(auto& f:filtered) f.assign(product(count+1,fullArea),0.);
    for(std::size_t p=0;p<count+3;++p) {
      const std::ptrdiff_t global=static_cast<std::ptrdiff_t>(start)+static_cast<std::ptrdiff_t>(p)-1;
      if(global<=static_cast<std::ptrdiff_t>(lower[2]) || global>=static_cast<std::ptrdiff_t>(upper[2])) continue;
      const std::size_t local=static_cast<std::size_t>(global)-lower[2]-1;
      const auto& slab=potentialBlock(local/blockLayers);
      if(openSides) {
        double* plane=phi.data()+p*fullArea;
        const int k=static_cast<int>(global),nx=static_cast<int>(g.nx),ny=static_cast<int>(g.ny);
        for(int j=0;j<=ny;++j){plane[j*(nx+1)]=side.at(0,j,k);plane[j*(nx+1)+nx]=side.at(nx,j,k);}
        for(int i=1;i<nx;++i){plane[i]=side.at(i,0,k);plane[ny*(nx+1)+i]=side.at(i,ny,k);}
      }
      for(std::size_t j=lower[1]+1;j<upper[1];++j)
        std::copy(slab.data()+(local%blockLayers)*area+(j-lower[1]-1)*nodes[0],
          slab.data()+(local%blockLayers)*area+(j-lower[1])*nodes[0],
          phi.data()+p*fullArea+j*(g.nx+1)+lower[0]+1);
    }
    const int threads=computeThreadCount();
#ifndef APRL_HAS_OPENMP
    (void)threads;
#endif
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for collapse(2) num_threads(threads) if(threads>1)
#endif
    for(std::size_t p=0;p<count+1;++p) for(std::size_t j=lower[1]+1;j<upper[1];++j) {
      if(start+p<=lower[2] || start+p>=upper[2]) continue;
      for(std::size_t i=lower[0]+1;i<upper[0];++i) {
        const std::size_t q=p*fullArea+j*(g.nx+1)+i;
        for(unsigned axis=0;axis<3;++axis) filtered[axis][q]=smooth(phi,axis,i,j,p+1);
      }
    }
    if(openSides) {
      const int nx=static_cast<int>(g.nx),ny=static_cast<int>(g.ny);
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for num_threads(threads) if(threads>1)
#endif
      for(std::size_t p=0;p<count+1;++p) {
        const int k=static_cast<int>(start+p);if(k<=0 || k>=int(g.nz)) continue;
        const auto fillWall=[&](int i,int j) {
          const std::size_t q=p*fullArea+std::size_t(j)*(g.nx+1)+i;
          for(unsigned axis=0;axis<3;++axis)filtered[axis][q]=smoothOnWall(phi,axis,i,j,p+1,k);
        };
        for(int j=0;j<=ny;++j){fillWall(0,j);fillWall(nx,j);}
        for(int i=1;i<nx;++i){fillWall(i,0);fillWall(i,ny);}
      }
    }
    const Double spacing[3]={g.dx,g.dy,g.dz};
    const std::size_t stride[3]={1,g.nx+1,fullArea};
    for(unsigned c=0;c<6;++c) {
      const std::size_t zcount=start<nz[c]?std::min(count,nz[c]-start):0;
      fields[c].resize(product(zcount,planeSize[c]));
      if(c==5) { std::fill(fields[c].begin(),fields[c].end(),0.); continue; }
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for collapse(2) num_threads(threads) if(threads>1)
#endif
      for(std::size_t p=0;p<zcount;++p) for(std::size_t j=0;j<ny[c];++j)
        for(std::size_t i=0;i<nx[c];++i) {
          const std::size_t q=p*fullArea+j*(g.nx+1)+i;
          Double value;
          if(c<3) value=-(filtered[c][q+stride[c]]-filtered[c][q])/spacing[c];
          else {
            const unsigned axis=c==3?1:0;
            const Double e0=-(filtered[axis][q+stride[axis]]-filtered[axis][q])/spacing[axis];
            const Double e1=-(filtered[axis][q+fullArea+stride[axis]]-filtered[axis][q+fullArea])/spacing[axis];
            value=(c==3?-magneticScale:magneticScale)*(lowerWeight*e0+upperWeight*e1);
          }
          fields[c][p*planeSize[c]+j*nx[c]+i]=value;
        }
    }
    for(const auto& f:fields) for(Double value:f) if(!std::isfinite(value))
      throw std::runtime_error("Nonfinite blocked native field");
    fieldStart=start; fieldStop=stop;
    ++stats.fieldCacheBuilds;
    std::size_t workspace=phi.capacity();
    for(const auto& f:filtered) workspace+=f.capacity();
    stats.queryWorkspaceBytes=std::max(stats.queryWorkspaceBytes,bytes(workspace));
    memoryReport();
    stats.querySeconds+=seconds(started);
  }
  Double field(unsigned component,std::size_t i,std::size_t j,std::size_t k) {
    if(component>=6 || i>=nx[component] || j>=ny[component] || k>=nz[component])
      throw std::out_of_range("Blocked field query outside native global Yee grid");
    ensureField(k); ++stats.fieldQueries;
    return fields[component][(k-fieldStart)*planeSize[component]+j*nx[component]+i];
  }
  void fill(EBFieldGrid& output,std::size_t offset) {
    const auto& local=output.geometry();
    if(local.nx!=g.nx || local.ny!=g.ny || local.dx!=g.dx || local.dy!=g.dy ||
       local.dz!=g.dz || local.dt!=g.dt || local.solver!=g.solver ||
       offset>g.nz || local.nz>g.nz-offset)
      throw std::invalid_argument("Blocked reference fill geometry/offset mismatch");
    YeeComponent* components[6]={&output.ex(),&output.ey(),&output.ez(),
      &output.bx(),&output.by(),&output.bz()};
    const std::size_t end=offset+local.nz+1;
    // z block outside component loop: each disk/potential/native block is
    // constructed once, then all six components are copied from it.
    for(std::size_t k=offset;k<end;) {
      ensureField(k);
      const std::size_t upperK=std::min(end,fieldStop);
      for(unsigned c=0;c<6;++c) {
        const std::size_t until=std::min(upperK,offset+components[c]->nz());
        if(until<=k) continue;
        const std::size_t count=product(until-k,planeSize[c]);
        std::copy(fields[c].data()+(k-fieldStart)*planeSize[c],
          fields[c].data()+(k-fieldStart)*planeSize[c]+count,
          components[c]->data()+(k-offset)*planeSize[c]);
      }
      k=upperK;
    }
  }
  void clear() {
    for(auto& p:phiCache) { std::vector<Double>().swap(p.values); p.block=std::numeric_limits<std::size_t>::max(); }
    std::vector<Double>().swap(sourceCache.values); sourceCache.block=std::numeric_limits<std::size_t>::max();
    for(auto& f:fields) std::vector<Double>().swap(f);
    fieldStart=fieldStop=0;
    memoryReport();
  }
};

BlockedCowanReference::BlockedCowanReference(const EBGridGeometry& g,
    const FieldVector<Double>& origin,const std::vector<RelativisticParticleSI>& particles,
    Double beta,const std::size_t guards[3],std::size_t block,const std::string& directory,
    MPI_Comm communicator,bool shaped,bool notch,bool binomial3,bool transverseXBinomial3,bool openTransverse)
  :impl_(new Impl(g,origin,particles,beta,guards,block,directory,communicator,shaped,notch,binomial3,transverseXBinomial3,openTransverse)) {}
BlockedCowanReference::BlockedCowanReference(std::unique_ptr<Impl> implementation)
  :impl_(std::move(implementation)) {}
std::unique_ptr<BlockedCowanReference> BlockedCowanReference::prepareCollective(
    const EBGridGeometry& g,const FieldVector<Double>& origin,
    const std::vector<RelativisticParticleSI>& particles,Double beta,
    const std::size_t guards[3],std::size_t block,const std::string& directory,
    MPI_Comm comm,bool shaped,bool notch,bool binomial3,bool transverseXBinomial3,bool openTransverse) {
  int initialized=0,finished=0;
  MPI_Initialized(&initialized);if(initialized)MPI_Finalized(&finished);
  if(!initialized || finished) throw std::invalid_argument("Collective reference requires live MPI");
  if(comm==MPI_COMM_NULL) throw std::invalid_argument("Null reference preparation communicator");
  int rank=0;MPI_Comm_rank(comm,&rank);
  const auto started=Clock::now();
  std::unique_ptr<Impl> writer;
  std::vector<Double> diagonal,offDiagonal;
  reference_detail::collectively(comm,[&] {
    if(rank!=0) return;
    writer.reset(new Impl(g,origin,particles,beta,guards,block,directory,MPI_COMM_SELF,
                         shaped,notch,binomial3,transverseXBinomial3,openTransverse,true));
    writer->beginPreparation(particles,diagonal,offDiagonal);
  });
  if(openTransverse) {
    const auto traceStart=Clock::now();
    OpenTransverseTrace::Sources source;
    OpenTransverseTrace worker;
    reference_detail::collectively(comm,[&] {if(rank==0)source=writer->openSources();});
    auto& trace=rank==0?writer->side:worker;
    // Use the exact native translation coefficient, including the same SI
    // constants as the singleton writer, rather than replacing it by beta^2.
    double gamma=rank==0?1./std::sqrt(1.-writer->translation):0.;
    MPI_Bcast(&gamma,1,MPI_DOUBLE,0,comm);
    trace.build(source,g,gamma,comm);
    reference_detail::collectively(comm,[&] {if(rank==0)writer->finishOpenSides(traceStart);});
    const std::uint64_t local=trace.localPotentialEvaluations;
    std::uint64_t minimum=0,maximum=0;
    MPI_Reduce(&local,&minimum,1,MPI_UINT64_T,MPI_MIN,0,comm);
    MPI_Reduce(&local,&maximum,1,MPI_UINT64_T,MPI_MAX,0,comm);
    if(rank==0) {
      writer->stats.openTraceMinRankEvaluations=minimum;
      writer->stats.openTraceMaxRankEvaluations=maximum;
    }
  }
  reference_detail::collectively(comm,[&] {
    if(rank==0)writer->finishPreparation(diagonal,offDiagonal,started);
  });
  std::unique_ptr<BlockedCowanReference> result;
  reference_detail::collectively(comm,[&] {
    if(rank==0)result.reset(new BlockedCowanReference(std::move(writer)));
    else result=openPrepared(g,origin,beta,guards,block,directory,MPI_COMM_SELF,
                             shaped,notch,binomial3,transverseXBinomial3,openTransverse);
  });
  return result;
}
std::unique_ptr<BlockedCowanReference> BlockedCowanReference::openPrepared(
    const EBGridGeometry& g,const FieldVector<Double>& origin,Double beta,
    const std::size_t guards[3],std::size_t block,const std::string& directory,
    MPI_Comm communicator,bool shaped,bool notch,bool binomial3,bool transverseXBinomial3,bool openTransverse) {
  return std::unique_ptr<BlockedCowanReference>(new BlockedCowanReference(
    std::unique_ptr<Impl>(new Impl(g,origin,beta,guards,block,directory,communicator,shaped,notch,binomial3,transverseXBinomial3,openTransverse))));
}
BlockedCowanReference::~BlockedCowanReference() = default;
Double BlockedCowanReference::field(unsigned c,std::size_t i,std::size_t j,std::size_t k) {
  return impl_->field(c,i,j,k);
}
void BlockedCowanReference::fill(EBFieldGrid& fields,std::size_t offset) { impl_->fill(fields,offset); }
Double BlockedCowanReference::chargeDensity(std::size_t i,std::size_t j,std::size_t k) {
  return impl_->scalar(i,j,k,true);
}
Double BlockedCowanReference::potential(std::size_t i,std::size_t j,std::size_t k) {
  return impl_->scalar(i,j,k,false);
}
void BlockedCowanReference::clearCache() { impl_->clear(); }
const BlockedCowanReferenceReport& BlockedCowanReference::report() const { return impl_->stats; }
void BlockedCowanReference::writeReport(const std::string& path) const {
  const auto& r=impl_->stats;
  std::ostringstream text;
  text << std::setprecision(17) << "{\n  \"scope\": \"blocked Cowan reference; open-wall potential targets use MPI/OpenMP; deposition and blocked modal solve stay on rank zero; independent singleton readers\",\n"
    << "  \"directory\": " << jsonString(impl_->directory) << ",\n"
    << "  \"files\": {\"source\": \"rho.bin\", \"checkpoints\": \"checkpoints.bin\", \"potential\": \"phi.bin\", \"completion\": \"metadata.txt\"},\n"
    << "  \"layout\": \"native-endian double, interior x-fast/y/z; checkpoint cprime,dprime at each block start\",\n"
    << "  \"fieldTimes\": \"E(0),B(-dt/2)\",\n"
    << "  \"complete\": " << (r.complete?"true":"false") << ",\n"
    << "  \"preparedReader\": " << (r.preparedReader?"true":"false") << ",\n"
    << "  \"interiorNodes\": [" << r.interiorNodes[0] << ',' << r.interiorNodes[1] << ',' << r.interiorNodes[2] << "],\n";
#define REPORT_MEMBER(name) text << "  \"" #name "\": " << r.name << ",\n"
  REPORT_MEMBER(particles); REPORT_MEMBER(blockLayers); REPORT_MEMBER(blocks); REPORT_MEMBER(beta);
  REPORT_MEMBER(openSeconds);
  if(impl_->openSides) {
    text << "  \"referenceBoundary\": \"open transverse Coulomb data / finite longitudinal zero ends\",\n";
    REPORT_MEMBER(openTraceSeconds); REPORT_MEMBER(openTraceBytes);
    REPORT_MEMBER(openTraceChecks); REPORT_MEMBER(openTraceMaxRelativeError);
    REPORT_MEMBER(openTraceSourceNodes);REPORT_MEMBER(openTraceCoarsePoints);REPORT_MEMBER(openTraceFallbackTiles);
    REPORT_MEMBER(openTraceRanks);REPORT_MEMBER(openTraceThreadsPerRank);
    REPORT_MEMBER(openTraceMinRankEvaluations);REPORT_MEMBER(openTraceMaxRankEvaluations);
  }
  REPORT_MEMBER(sourceCharge); REPORT_MEMBER(depositedCharge); REPORT_MEMBER(minimumPivot);
  REPORT_MEMBER(minimumModeEigenvalue); REPORT_MEMBER(preparationSeconds); REPORT_MEMBER(sourceSeconds);
  REPORT_MEMBER(forwardSeconds); REPORT_MEMBER(reconstructionSeconds); REPORT_MEMBER(querySeconds);
  REPORT_MEMBER(sourceFileBytes); REPORT_MEMBER(checkpointFileBytes); REPORT_MEMBER(potentialFileBytes);
  REPORT_MEMBER(bytesRead); REPORT_MEMBER(bytesWritten); REPORT_MEMBER(fieldQueries);
  REPORT_MEMBER(potentialCacheHits); REPORT_MEMBER(potentialCacheMisses); REPORT_MEMBER(fieldCacheBuilds);
  REPORT_MEMBER(sourceIndexBytes); REPORT_MEMBER(sourceSortingAllowanceBytes); REPORT_MEMBER(solverWorkspaceBytes);
  REPORT_MEMBER(potentialCacheCapacityBytes); REPORT_MEMBER(sourceCacheCapacityBytes);
  REPORT_MEMBER(fieldCacheCapacityBytes); REPORT_MEMBER(queryWorkspaceBytes); REPORT_MEMBER(residentArrayBytes);
#undef REPORT_MEMBER
  text << "  \"accounting\": \"bytesRead/Written count binary scalar data only; completion metadata and report IO excluded; capacities/workspaces are array accounting, not RSS; FFTW/MPI/runtime private allocation excluded; querySeconds measures cache builds, not cache-hit or copy time; readers inherit particle/charge history but not preparation timers or binary IO; preparation includes completion fsync/publication\",\n"
    << "  \"limitations\": \"finite Dirichlet common-beta source; open-wall potential targets use MPI, deposition and modal solve stay on rank zero; no source padding; all sources must remain available on the writer during preparation; scalar disk data and checkpoints grow with longitudinal length; global forward/backward pass remains necessary; shared filesystem required; metadata has a checksum but scalar data only have length/short-read/finite checks; Gauss/PIC/held comparison belongs to independent checks\"\n}\n";
  const std::string serialized=text.str();
  DiskFile out;
  out.create(path);
  out.transfer(const_cast<char*>(serialized.data()),serialized.size(),0,true,0);
}

}} // namespace aprl::verification
