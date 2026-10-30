#include "transverse_reference_current.h"
#include "blocked_cowan_reference.h"
#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
namespace aprl {namespace {
const std::uint64_t magic=0x4150524c53555246ULL,complete=0x434f4d504c455445ULL;
struct Header {std::uint64_t tag,version,nx,ny,nz,gx,gy,count;double dx,dy,dz,dt,beta;};
struct Point {std::size_t i,j;long double weight;};
std::vector<Point> points(const EBGridGeometry& g,const std::size_t guard[3]) {
  if(!guard[0]||!guard[1]||guard[2])throw std::invalid_argument("Reference surface current requires transverse-only CPML");
  const auto x0=guard[0],x1=g.nx-guard[0],y0=guard[1],y1=g.ny-guard[1];
  std::vector<Point> result;
  for(std::size_t i=x0;i<=x1;++i)for(auto j:{y0,y1})result.push_back({i,j,(i==x0||i==x1)?.25L:.5L});
  for(std::size_t j=y0+1;j<y1;++j)for(auto i:{x0,x1})result.push_back({i,j,.5L});
  return result;
}
Header header(const EBGridGeometry& g,const std::size_t guard[3],double beta) {
  return {magic,1,g.nx,g.ny,g.nz,guard[0],guard[1],points(g,guard).size(),g.dx,g.dy,g.dz,g.dt,beta};
}
void writeAll(int fd,const void* data,std::size_t size) {
  auto* p=static_cast<const char*>(data);
  while(size) {const auto n=::write(fd,p,size);if(n<0&&errno==EINTR)continue;
    if(n<=0)throw std::runtime_error("Cannot write reference surface current");
    p+=n;size-=static_cast<std::size_t>(n);}
}
void readAll(int fd,void* data,std::size_t size,std::size_t offset) {
  auto* p=static_cast<char*>(data);
  while(size) {const auto n=::pread(fd,p,size,static_cast<off_t>(offset));if(n<0&&errno==EINTR)continue;
    if(n<=0)throw std::runtime_error("Cannot read reference surface current");
    p+=n;offset+=static_cast<std::size_t>(n);size-=static_cast<std::size_t>(n);}
}
}
void TransverseReferenceCurrent::prepare(const std::string& file,
    verification::BlockedCowanReference& ref,const EBGridGeometry& g,const std::size_t guard[3],double beta) {
  const auto h=header(g,guard,beta);const auto target=points(g,guard);
  const int fd=::open(file.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);
  if(fd<0)throw std::runtime_error("Cannot create reference surface-current file");
  try {
    writeAll(fd,&h,sizeof(h));std::vector<double> values(target.size());
    for(std::size_t k=0;k<g.nz;++k) {
      for(std::size_t p=0;p<target.size();++p) {
        const auto i=target[p].i,j=target[p].j;
        values[p]=((ref.field(4,i,j,k)-ref.field(4,i-1,j,k))/g.dx
                  -(ref.field(3,i,j,k)-ref.field(3,i,j-1,k))/g.dy)/SI::mu0;
        if(!std::isfinite(values[p]))throw std::runtime_error("Nonfinite equivalent boundary current");
      }
      writeAll(fd,values.data(),values.size()*sizeof(double));
    }
    writeAll(fd,&complete,sizeof(complete));if(::fsync(fd))throw std::runtime_error("Cannot sync reference surface current");
    ::close(fd);
  }catch(...){::close(fd);throw;}
}
struct TransverseReferenceCurrent::Impl {
  struct Cache {std::size_t start=std::numeric_limits<std::size_t>::max(),stamp=0;std::vector<double> data;};
  Header h;std::vector<Point> target;int fd=-1;std::array<Cache,2> cache;std::size_t clock=0,step=0;
  std::vector<double> current,oldElectric;std::ofstream work;long double cumulative=0.;bool audit=false,applied=false;
  Impl(const std::string& file,const EBGridGeometry& g,const std::size_t guards[3],double beta,const std::string& output)
      :h(header(g,guards,beta)),target(points(g,guards)) {
    fd=::open(file.c_str(),O_RDONLY);if(fd<0)throw std::runtime_error("Cannot open reference surface current");
    try {
      Header actual;readAll(fd,&actual,sizeof(actual),0);
      if(actual.tag!=h.tag||actual.version!=h.version||actual.nx!=h.nx||actual.ny!=h.ny||actual.nz!=h.nz||actual.gx!=h.gx||actual.gy!=h.gy||actual.count!=h.count||
         actual.dx!=h.dx||actual.dy!=h.dy||actual.dz!=h.dz||actual.dt!=h.dt||actual.beta!=h.beta)
        throw std::runtime_error("Reference surface-current metadata mismatch");
      struct stat st;if(::fstat(fd,&st))throw std::runtime_error("Cannot stat reference surface current");
      const std::size_t length=sizeof(Header)+h.nz*h.count*sizeof(double)+sizeof(complete);
      if(static_cast<std::uint64_t>(st.st_size)!=length)throw std::runtime_error("Reference surface-current length mismatch");
      std::uint64_t footer;readAll(fd,&footer,sizeof(footer),length-sizeof(footer));
      if(footer!=complete)throw std::runtime_error("Reference surface-current preparation incomplete");
      if(!output.empty()) {
        std::ifstream exists(output);if(exists.good())throw std::runtime_error("Reference work file exists");
        work.open(output);if(!work)throw std::runtime_error("Cannot create reference source-work file");
        work<<std::setprecision(21)<<"step,work_J,absolute_work_J,cumulative_work_J\n";
      }
    }catch(...){::close(fd);fd=-1;throw;}
  }
  ~Impl(){if(fd>=0)::close(fd);}
  const double* row(std::size_t k) {
    const auto start=k/32*32;
    for(auto& b:cache)if(b.start==start){b.stamp=++clock;return b.data.data()+(k-start)*h.count;}
    auto& b=cache[cache[0].stamp<=cache[1].stamp?0:1];
    b.start=start;b.stamp=++clock;b.data.resize(std::min<std::size_t>(32,h.nz-start)*h.count);
    readAll(fd,b.data.data(),b.data.size()*sizeof(double),sizeof(Header)+start*h.count*sizeof(double));
    for(double x:b.data)if(!std::isfinite(x))throw std::runtime_error("Nonfinite reference current cache");
    return b.data.data()+(k-start)*h.count;
  }
};
TransverseReferenceCurrent::TransverseReferenceCurrent(const std::string& file,const EBGridGeometry& g,const std::size_t guards[3],double beta,const std::string& work)
 :impl_(new Impl(file,g,guards,beta,work)) {}
TransverseReferenceCurrent::~TransverseReferenceCurrent()=default;
void TransverseReferenceCurrent::apply(EBFieldGrid& f,std::ptrdiff_t offset,std::size_t step,double beta,bool collect) {
  auto& s=*impl_;if(s.applied)throw std::logic_error("Unfinished reference source step");s.step=step+1;s.audit=collect;
  const auto nz=f.geometry().nz;s.current.resize(collect?nz*s.target.size():0);s.oldElectric.resize(s.current.size());
  for(std::size_t k=0;k<nz;++k) {
    const double z=static_cast<double>(offset)+k-beta*static_cast<double>(step+1);
    if(!std::isfinite(z)||z<0.||z>s.h.nz-1)throw std::out_of_range("Reference surface-current query outside domain");
    const auto lo=static_cast<std::size_t>(std::floor(z)),hi=std::min<std::size_t>(lo+1,s.h.nz-1);const double fraction=z-lo;
    const auto* a=s.row(lo);const auto* b=s.row(hi);
    for(std::size_t p=0;p<s.target.size();++p) {
      const auto& q=s.target[p];const double value=fraction==0.?a[p]:(1.-fraction)*a[p]+fraction*b[p];
      if(collect){s.current[k*s.target.size()+p]=value;s.oldElectric[k*s.target.size()+p]=f.ez()(q.i,q.j,k);}
      f.jz()(q.i,q.j,k)+=value;
    }
  }
  s.applied=true;
}
TransverseReferenceCurrent::Work TransverseReferenceCurrent::finish(const EBFieldGrid& f) {
  auto& s=*impl_;if(!s.applied)throw std::logic_error("Missing reference source step");Work result;
  if(s.audit) {
    const long double scale=-.5L*f.geometry().dt*f.geometry().dx*f.geometry().dy*f.geometry().dz;
    for(std::size_t k=0;k<f.geometry().nz;++k)for(std::size_t p=0;p<s.target.size();++p) {
      const auto& q=s.target[p];const auto n=k*s.target.size()+p;
      const long double v=scale*q.weight*s.current[n]*(static_cast<long double>(s.oldElectric[n])+f.ez()(q.i,q.j,k));
      result.signedJ+=v;result.absoluteJ+=std::abs(v);
    }
    s.cumulative+=result.signedJ;
    if(s.work.is_open()){s.work<<s.step<<','<<result.signedJ<<','<<result.absoluteJ<<','<<s.cumulative<<'\n';if(!s.work)throw std::runtime_error("Cannot write reference source work");}
  }
  s.applied=false;return result;
}
void TransverseReferenceCurrent::close(){auto& s=*impl_;if(s.applied)throw std::logic_error("Reference current still applied");if(s.work.is_open()){s.work.flush();if(!s.work)throw std::runtime_error("Cannot flush reference work");s.work.close();}}
std::size_t TransverseReferenceCurrent::memoryBytes() const {
  const auto& s=*impl_;std::size_t n=s.target.capacity()*sizeof(Point)+(s.current.capacity()+s.oldElectric.capacity())*sizeof(double);
  for(const auto& c:s.cache)n+=c.data.capacity()*sizeof(double);
  return n;
}
}
