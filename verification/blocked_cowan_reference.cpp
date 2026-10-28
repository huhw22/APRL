#include "blocked_cowan_reference.h"
#include "compute_threads.h"

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
#include <mpi.h>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

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
  DiskFile() : fd_(-1) {}
  ~DiskFile() { if (fd_>=0) ::close(fd_); }
  DiskFile(const DiskFile&)=delete;
  DiskFile& operator=(const DiskFile&)=delete;
  void create(const std::string& path) {
    path_=path;
    fd_=::open(path.c_str(),O_RDWR|O_CREAT|O_EXCL,0600);
    if (fd_<0) throw ioError("Cannot create new reference file",path_);
  }
  void transfer(void* data, std::size_t count, std::size_t offset,
                bool writing, std::uint64_t* counter) const {
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
    if (info.st_size<0) throw std::runtime_error("Negative reference file length");
    return static_cast<std::uint64_t>(info.st_size);
  }
private:
  int fd_;
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
}

struct BlockedCowanReference::Impl {
  struct Cache {
    std::size_t block=std::numeric_limits<std::size_t>::max();
    std::uint64_t used=0;
    std::vector<Double> values;
  };
  EBGridGeometry g;
  FieldVector<Double> origin;
  Double beta,translation,upperWeight,lowerWeight,magneticScale;
  EBCowanCoefficients coefficient;
  std::size_t lower[3],upper[3],nodes[3],cells[3];
  std::size_t blockLayers,area,fullArea;
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
       const std::size_t guards[3],std::size_t block,const std::string& dir)
    : g(geometry),origin(where),beta(velocity),translation(0.),upperWeight(0.),
      lowerWeight(0.),magneticScale(0.),blockLayers(block),area(0),fullArea(0),directory(dir) {
    const auto start=Clock::now();
    validate(guards);
    coefficient=EBCowanCoefficients::forZDispersion(g.dx,g.dy,g.dz);
    const Double speed=beta*SI::c;
    translation=speed*speed*SI::epsilon0*SI::mu0;
    upperWeight=0.5+0.5*beta*SI::c*g.dt/g.dz;
    lowerWeight=1.-upperWeight;
    magneticScale=beta*SI::c*SI::epsilon0*SI::mu0;
    if (!(upperWeight>=0. && upperWeight<=1.))
      throw std::invalid_argument("Blocked initial B half-step shift exceeds one cell");
    // Validate the discrete operator before any filesystem mutation.
    std::vector<Double> diagonal(area),offDiagonal(area);
    makeModes(diagonal,offDiagonal);
    if (::mkdir(directory.c_str(),0700)!=0)
      throw ioError("Cannot create new blocked reference directory",directory);
    rhoFile.create(directory+"/rho.bin");
    checkpointFile.create(directory+"/checkpoints.bin");
    phiFile.create(directory+"/phi.bin");
    stats.particles=particles.size();
    deposit(particles);
    solve(diagonal,offDiagonal);
    stats.sourceFileBytes=rhoFile.size();
    stats.checkpointFileBytes=checkpointFile.size();
    stats.potentialFileBytes=phiFile.size();
    if(stats.sourceFileBytes!=bytes(product(area,nodes[2])) ||
       stats.potentialFileBytes!=stats.sourceFileBytes ||
       stats.checkpointFileBytes!=bytes(product(product(stats.blocks,area),2)))
      throw std::runtime_error("Blocked reference prepared file has incorrect length");
    stats.preparationSeconds=seconds(start);
    stats.complete=true;
  }

  void validate(const std::size_t guards[3]) {
    int initialized=0,finished=0,size=0;
    MPI_Initialized(&initialized);
    if(initialized) MPI_Finalized(&finished);
    if(!initialized || finished)
      throw std::invalid_argument("Blocked reference requires live single-rank MPI");
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if(size!=1) throw std::invalid_argument("Blocked reference currently supports one MPI rank only");
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

  void makeModes(std::vector<Double>& diagonal,std::vector<Double>& off) {
    Double stencil[3][3][3]={};
    const Double spacing[3]={g.dx,g.dy,g.dz};
    for(unsigned axis=0;axis<3;++axis) {
      const unsigned ta=(axis+1)%3,tb=(axis+2)%3;
      const Double delta=deltaFor(coefficient,axis);
      Double minimum=std::numeric_limits<Double>::infinity();
      for(int a=-1;a<=1;a+=2) for(int b=-1;b<=1;b+=2)
        minimum=std::min(minimum,coefficient.alpha[axis]+
          2.*coefficient.beta[ta]*a+2.*coefficient.beta[tb]*b+4.*delta*a*b);
      if(axis==2) minimum-=translation;
      if(!(minimum>64.*std::numeric_limits<Double>::epsilon()))
        throw std::invalid_argument("Blocked Cowan operator fails production sufficient SPD bound");
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
    for(std::size_t k=0;k<layers(block);++k) transform.transform(rhs.data()+k*area);
  }

  void memoryReport() {
    std::size_t potentialBytes=0,fieldBytes=0;
    for(const auto& p:phiCache) potentialBytes+=bytes(p.values.capacity());
    for(const auto& f:fields) fieldBytes+=bytes(f.capacity());
    const std::size_t sourceBytes=bytes(sourceCache.values.capacity());
    stats.potentialCacheCapacityBytes=std::max(stats.potentialCacheCapacityBytes,potentialBytes);
    stats.sourceCacheCapacityBytes=std::max(stats.sourceCacheCapacityBytes,sourceBytes);
    stats.fieldCacheCapacityBytes=std::max(stats.fieldCacheCapacityBytes,fieldBytes);
    stats.residentArrayBytes=potentialBytes+sourceBytes+fieldBytes;
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
    if(!interior(i,j,k)) return 0.;
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
    Double beta,const std::size_t guards[3],std::size_t block,const std::string& directory)
  :impl_(new Impl(g,origin,particles,beta,guards,block,directory)) {}
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
  text << std::setprecision(17) << "{\n  \"scope\": \"verification-only single-rank blocked Cowan reference; not PIC acceptance or production performance\",\n"
    << "  \"directory\": " << jsonString(impl_->directory) << ",\n"
    << "  \"files\": {\"source\": \"rho.bin\", \"checkpoints\": \"checkpoints.bin\", \"potential\": \"phi.bin\"},\n"
    << "  \"layout\": \"native-endian double, interior x-fast/y/z; checkpoint cprime,dprime at each block start\",\n"
    << "  \"fieldTimes\": \"E(0),B(-dt/2)\",\n"
    << "  \"complete\": " << (r.complete?"true":"false") << ",\n"
    << "  \"interiorNodes\": [" << r.interiorNodes[0] << ',' << r.interiorNodes[1] << ',' << r.interiorNodes[2] << "],\n";
#define REPORT_MEMBER(name) text << "  \"" #name "\": " << r.name << ",\n"
  REPORT_MEMBER(particles); REPORT_MEMBER(blockLayers); REPORT_MEMBER(blocks); REPORT_MEMBER(beta);
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
  text << "  \"accounting\": \"bytesRead/Written count binary data only; capacities/workspaces are array accounting, not RSS; FFTW/MPI/runtime private allocation excluded; querySeconds measures cache builds, not cache-hit or copy time\",\n"
    << "  \"limitations\": \"finite Dirichlet common-beta source; no padding or MPI decomposition; all sources must remain available during preparation; scalar disk data and checkpoints grow with longitudinal length; global forward/backward pass remains necessary; Gauss/PIC/held comparison belongs to independent checks\"\n}\n";
  const std::string serialized=text.str();
  DiskFile out;
  out.create(path);
  out.transfer(const_cast<char*>(serialized.data()),serialized.size(),0,true,0);
}

}} // namespace aprl::verification
