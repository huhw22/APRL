// Read-only integration of the frozen observer's bilinear spacetime currents.
// All channels are measured; no residual is used to supply a missing term.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <omp.h>
using R=long double;
struct Point{R x,y;};
using Polygon=std::vector<Point>;
Polygon clip(const Polygon& p,R a,R b,R c){
  Polygon out;if(p.empty())return out;Point prev=p.back();R vp=a+b*prev.x+c*prev.y;
  for(auto next:p){R vn=a+b*next.x+c*next.y;if((vp>=0)!=(vn>=0)){R t=vp/(vp-vn);out.push_back({prev.x+t*(next.x-prev.x),prev.y+t*(next.y-prev.y)});}if(vn>=0)out.push_back(next);prev=next;vp=vn;}return out;
}
std::array<R,4> moments(const Polygon& p){
  std::array<R,4> v{};if(p.size()<3)return v;auto a=p[0];
  for(std::size_t i=1;i+1<p.size();++i){auto b=p[i],c=p[i+1];R A=std::abs((b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x))/2,sx=a.x+b.x+c.x,sy=a.y+b.y+c.y;v[0]+=A;v[1]+=A*sx/3;v[2]+=A*sy/3;v[3]+=A*(a.x*a.y+b.x*b.y+c.x*c.y+sx*sy)/12;}return v;
}
void interval(R& l,R& h,R a,R b,R low,R high){
  if(b==0){if(!(low<=a&&a<=high)){l=1;h=0;}return;}R p=(low-a)/b,q=(high-a)/b;l=std::max(l,std::min(p,q));h=std::min(h,std::max(p,q));
}
R trace(R q00,R q10,R q01,R q11,R u,R v,R l,R h){
  R qx=q10-q00,qy=q01-q00,qxy=q11-q10-q01+q00;
  R a=q00+qy*u,b=qx+qy*v+qxy*u,c=qxy*v;
  return a*(h-l)+b*(h*h-l*l)/2+c*(h*h*h-l*l*l)/3;
}
struct Sum{R v=0,c=0;void add(R x){R y=x-c,t=v+y;c=(t-v)-y;v=t;}};
class Tube{
 public:
  std::size_t nt,nz;R g,b,eta,c,dz,zorig,boxorig,lo,hi,bz,bt;
  std::vector<R> t,o,dt,az,at,zb,tb,zs;const double* q=nullptr;std::size_t bytes=0;int fd=-1;
  Tube(const std::string& directory,const std::string& config){
    std::ifstream cfg(config);std::size_t count;
    cfg>>nt>>nz>>g>>c>>dz>>zorig>>boxorig>>lo>>hi>>count;
    if(!cfg||nt<2||nz<2||count<2)throw std::runtime_error("Invalid worldtube config");
    zs.resize(count);for(auto& z:zs)cfg>>z;if(!cfg)throw std::runtime_error("Invalid plane list");
    b=std::sqrt((g-1)*(g+1))/g;eta=1/(g*(1+b));bz=g*dz;bt=-eta*dz/c;
    std::ifstream st(directory+"/states.csv");std::string line;std::getline(st,line);
    while(std::getline(st,line)){std::istringstream s(line);std::string cell;std::getline(s,cell,',');std::getline(s,cell,',');t.push_back(std::stold(cell));std::getline(s,cell,',');o.push_back(std::stold(cell));}
    if(t.size()!=nt||o.size()!=nt)throw std::runtime_error("State length mismatch");
    dt.resize(nt-1);az.resize(nt-1);at.resize(nt-1);zb.resize(nt-1);tb.resize(nt-1);
    for(std::size_t n=0;n+1<nt;++n){dt[n]=t[n+1]-t[n];R d=o[n+1]-o[n];if(dt[n]<=0)throw std::runtime_error("Time order");az[n]=g*(d+b*c*dt[n]);at[n]=eta*(dt[n]-d/c);zb[n]=zorig+g*(o[n]+dz/2-boxorig+b*c*t[n]);tb[n]=eta*(t[n]-(o[n]+dz/2-boxorig)/c);}
    fd=open((directory+"/field_lines.bin").c_str(),O_RDONLY);if(fd<0)throw std::runtime_error("Cannot open currents");
    struct stat fs{};fstat(fd,&fs);bytes=nt*nz*6*sizeof(double);if(std::size_t(fs.st_size)!=bytes)throw std::runtime_error("Raw extent mismatch");
    void* p=mmap(nullptr,bytes,PROT_READ,MAP_SHARED,fd,0);if(p==MAP_FAILED)throw std::runtime_error("Read-only map failed");q=static_cast<const double*>(p);
    std::uint64_t bad=0;
#pragma omp parallel for reduction(+:bad) schedule(static)
    for(std::size_t i=0;i<nt*nz*6;++i)if(!std::isfinite(q[i]))++bad;
    if(bad)throw std::runtime_error("Nonfinite raw observer currents");
  }
  ~Tube(){if(q)munmap(const_cast<double*>(q),bytes);if(fd>=0)close(fd);}
  R get(std::size_t n,std::size_t k,unsigned a)const{return q[(n*nz+k)*6+a];}
  R traced(std::size_t n,std::size_t k,unsigned a,R u,R v,R l,R h)const{return trace(get(n,k,a),get(n+1,k,a),get(n,k+1,a),get(n+1,k+1,a),u,v,l,h);}
  std::array<R,2> plane(R z)const{
    Sum total,coverage;
    for(std::size_t n=0;n+1<nt;++n){R k0=(z-zb[n])/bz,k1=(z-zb[n]-az[n])/bz;
      long first=std::max(0L,long(std::floor(std::min(k0,k1)))),last=std::min(long(nz)-2,long(std::floor(std::max(k0,k1))));
      for(long k=first;k<=last;++k){R u=(z-zb[n]-bz*k)/bz,v=-az[n]/bz,l=0,h=1;interval(l,h,u,v,0,1);R a=tb[n]+bt*k+bt*u,b0=at[n]+bt*v;interval(l,h,a,b0,lo,hi);if(h<=l)continue;total.add(traced(n,k,1,u,v,l,h)*dt[n]/g);coverage.add((h-l)*dt[n]/g);}}
    return {{total.v,coverage.v}};
  }
  std::array<R,2> cap(R tau,R z1)const{
    Sum total,coverage;if(z1<=zs[0])return {{0,0}};
    for(std::size_t n=0;n+1<nt;++n){R k0=(tau-tb[n])/bt,k1=(tau-tb[n]-at[n])/bt;
      long first=std::max(0L,long(std::floor(std::min(k0,k1)))),last=std::min(long(nz)-2,long(std::floor(std::max(k0,k1))));
      for(long k=first;k<=last;++k){R u=(tau-tb[n]-bt*k)/bt,v=-at[n]/bt,l=0,h=1;interval(l,h,u,v,0,1);R a=zb[n]+bz*k+bz*u,b0=az[n]+bz*v;interval(l,h,a,b0,zs[0],z1);if(h<=l)continue;total.add(traced(n,k,0,u,v,l,h)*std::abs(b0));coverage.add((h-l)*std::abs(b0));}}
    return {{total.v,coverage.v}};
  }
  std::vector<std::array<R,5>> sides()const{
    const std::size_t bins=zs.size()-1;int threads=omp_get_max_threads();
    std::vector<std::vector<std::array<Sum,5>>> partial(threads,std::vector<std::array<Sum,5>>(bins));
#pragma omp parallel for schedule(static)
    for(std::size_t n=0;n<nt-1;++n){auto& sums=partial[omp_get_thread_num()];
      R kmin=(zs.front()-zb[n]-std::max(az[n],R(0)))/bz,kmax=(zs.back()-zb[n]-std::min(az[n],R(0)))/bz;
      long first=std::max(0L,long(std::floor(kmin))),last=std::min(long(nz)-2,long(std::floor(kmax)));
      for(long k=first;k<=last;++k){R t0=tb[n]+bt*k,z=zb[n]+bz*k;
        R tmin=t0+std::min(at[n],R(0))+bt,tmax=t0+std::max(at[n],R(0));if(tmax<lo||tmin>hi)continue;
        R zmin=z+std::min(az[n],R(0)),zmax=z+std::max(az[n],R(0))+bz;
        long b0=std::max(0L,long(std::upper_bound(zs.begin(),zs.end(),zmin)-zs.begin())-1);
        long b1=std::min(long(bins)-1,long(std::lower_bound(zs.begin(),zs.end(),zmax)-zs.begin())-1);
        for(long bin=b0;bin<=b1;++bin){std::array<R,4> m;
          if(tmin>=lo&&tmax<=hi&&zmin>=zs[bin]&&zmax<=zs[bin+1])m={{1,.5L,.5L,.25L}};
          else{Polygon p={{0,0},{1,0},{1,1},{0,1}};p=clip(p,t0-lo,at[n],bt);p=clip(p,hi-t0,-at[n],-bt);p=clip(p,z-zs[bin],az[n],bz);p=clip(p,zs[bin+1]-z,-az[n],-bz);m=moments(p);}
          if(m[0]==0)continue;R weight=dt[n]*dz;
          for(unsigned a=0;a<4;++a){R q00=get(n,k,a+2),q10=get(n+1,k,a+2),q01=get(n,k+1,a+2),q11=get(n+1,k+1,a+2);sums[bin][a].add((q00*m[0]+(q10-q00)*m[1]+(q01-q00)*m[2]+(q11-q10-q01+q00)*m[3])*weight);}
          sums[bin][4].add(m[0]*weight);
        }
      }
    }
    std::vector<std::array<R,5>> out(bins);
    for(std::size_t b0=0;b0<bins;++b0)for(unsigned a=0;a<5;++a){Sum s;for(auto& thread:partial)s.add(thread[b0][a].v);out[b0][a]=s.v;}
    return out;
  }
  void write(const std::string& path)const{
    auto side=sides();auto pin=plane(zs[0]);R duration=hi-lo;std::array<R,5> cumulative{};
    std::ofstream out(path);out<<std::setprecision(24)<<"z_lab_m,field_entrance_signed_J,field_exit_signed_J,time_cap_early_signed_J,time_cap_late_signed_J,transverse_xminus_signed_J,transverse_xplus_signed_J,transverse_yminus_signed_J,transverse_yplus_signed_J,max_coverage_relative_error\n";
    for(std::size_t j=0;j<zs.size();++j){if(j)for(unsigned a=0;a<5;++a)cumulative[a]+=side[j-1][a];auto p=plane(zs[j]),early=cap(lo,zs[j]),late=cap(hi,zs[j]);R len=zs[j]-zs[0];
      R error=std::max(std::abs(pin[1]-duration)/duration,std::abs(p[1]-duration)/duration);
      if(j)error=std::max({error,std::abs(early[1]-len)/len,std::abs(late[1]-len)/len,std::abs(cumulative[4]-len*duration)/(len*duration)});
      if(error>2e-10L)throw std::runtime_error("Uncovered worldtube surface");
      out<<zs[j]<<','<<-pin[0]<<','<<p[0]<<','<<-early[0]<<','<<late[0];for(unsigned a=0;a<4;++a)out<<','<<cumulative[a];out<<','<<error<<'\n';
      std::cout<<"FIELD_SURFACES "<<j<<' '<<zs[j]<<' '<<error<<std::endl;
    }
  }
};
int main(int argc,char** argv){try{if(argc!=4)return 2;Tube tube(argv[1],argv[2]);tube.write(argv[3]);return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
