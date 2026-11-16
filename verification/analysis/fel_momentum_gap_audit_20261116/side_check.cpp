// Independent existing-data side quadrature: affine lower/upper bounds and
// two-point Gaussian integration, without polygon clipping or polygon moments.
#include <algorithm>
#include <array>
#include <cmath>
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
struct Line { R a,b; R operator()(R x)const{return a+b*x;} };
struct Sum {R v=0,c=0;void add(R x){R y=x-c,t=v+y;c=(t-v)-y;v=t;}};
int main(int argc,char** argv){try{
  if(argc!=4)return 2;
  std::size_t nt,nz,np;R G,c,dz,zo,bo,lo,hi;
  std::ifstream cfg(argv[2]);cfg>>nt>>nz>>G>>c>>dz>>zo>>bo>>lo>>hi>>np;
  std::vector<R> zs(np);for(auto& v:zs)cfg>>v;if(!cfg)throw std::runtime_error("config");
  R z0=zs.front(),z1=zs.back(),beta=std::sqrt((G-1)*(G+1))/G,eta=1/(G*(1+beta));
  std::vector<R> t,o;std::ifstream sf(std::string(argv[1])+"/states.csv");std::string row;
  std::getline(sf,row);while(std::getline(sf,row)){std::istringstream s(row);std::string v;
    std::getline(s,v,',');std::getline(s,v,',');t.push_back(std::stold(v));
    std::getline(s,v,',');o.push_back(std::stold(v));}
  if(t.size()!=nt)throw std::runtime_error("states");
  int fd=open((std::string(argv[1])+"/field_lines.bin").c_str(),O_RDONLY);
  struct stat st{};if(fd<0||fstat(fd,&st))throw std::runtime_error("raw");
  std::size_t bytes=nt*nz*6*sizeof(double);if(bytes!=std::size_t(st.st_size))throw std::runtime_error("size");
  void* raw=mmap(nullptr,bytes,PROT_READ,MAP_SHARED,fd,0);if(raw==MAP_FAILED)throw std::runtime_error("mmap");
  const double* q=static_cast<double*>(raw);
  auto get=[&](std::size_t n,long k,unsigned a)->R{return q[(n*nz+k)*6+a];};
  const R bz=G*dz,bt=-eta*dz/c,gauss=1/std::sqrt(R(3));
  std::vector<std::array<Sum,4>> sums(omp_get_max_threads());
  unsigned long long interior=0,boundary=0;
#pragma omp parallel for reduction(+:interior,boundary) schedule(static)
  for(std::size_t n=0;n<nt-1;++n){
    auto& total=sums[omp_get_thread_num()];R dt=t[n+1]-t[n],d=o[n+1]-o[n];
    R az=G*(d+beta*c*dt),at=eta*(dt-d/c);
    R zb=zo+G*(o[n]+dz/2-bo+beta*c*t[n]),tb=eta*(t[n]-(o[n]+dz/2-bo)/c);
    long first=std::max(0L,long(std::floor((z0-zb-std::max(az,R(0)))/bz)));
    long last=std::min(long(nz)-2,long(std::floor((z1-zb-std::min(az,R(0)))/bz)));
    for(long k=first;k<=last;++k){
      R z=zb+bz*k,tau=tb+bt*k;
      R zmin=z+std::min(az,R(0)),zmax=z+std::max(az,R(0))+bz;
      R tmin=tau+std::min(at,R(0))+bt,tmax=tau+std::max(at,R(0));
      if(tmax<lo||tmin>hi)continue;
      if(zmin>=z0&&zmax<=z1&&tmin>=lo&&tmax<=hi){
        ++interior;
        for(unsigned a=0;a<4;++a)
          total[a].add((get(n,k,a+2)+get(n+1,k,a+2)+get(n,k+1,a+2)+get(n+1,k+1,a+2))*dt*dz/4);
        continue;
      }
      ++boundary;
      std::array<Line,3> lower={{{0,0},{(z0-z)/bz,-az/bz},{(hi-tau)/bt,-at/bt}}};
      std::array<Line,3> upper={{{1,0},{(z1-z)/bz,-az/bz},{(lo-tau)/bt,-at/bt}}};
      std::array<Line,6> lines;std::copy(lower.begin(),lower.end(),lines.begin());std::copy(upper.begin(),upper.end(),lines.begin()+3);
      std::vector<R> cuts={0,1};
      for(unsigned a=0;a<6;++a)for(unsigned b=a+1;b<6;++b){
        R den=lines[a].b-lines[b].b;if(den==0)continue;
        R u=(lines[b].a-lines[a].a)/den;if(u>0&&u<1)cuts.push_back(u);
      }
      std::sort(cuts.begin(),cuts.end());
      for(std::size_t j=1;j<cuts.size();++j){
        R l=cuts[j-1],h=cuts[j],m=(l+h)/2;if(h<=l)continue;
        const auto li=std::max_element(lower.begin(),lower.end(),[&](Line a,Line b){return a(m)<b(m);});
        const auto ui=std::min_element(upper.begin(),upper.end(),[&](Line a,Line b){return a(m)<b(m);});
        if((*li)(m)>=(*ui)(m))continue;
        for(unsigned a=0;a<4;++a){
          R q00=get(n,k,a+2),qx=get(n+1,k,a+2)-q00,qy=get(n,k+1,a+2)-q00;
          R qxy=get(n+1,k+1,a+2)-q00-qx-qy,value=0;
          for(int sign:{-1,1}){R u=m+sign*(h-l)*gauss/2,v0=(*li)(u),v1=(*ui)(u);
            value+=(q00+qx*u)*(v1-v0)+(qy+qxy*u)*(v1*v1-v0*v0)/2;}
          total[a].add(value*(h-l)*dt*dz/2);
        }
      }
    }
  }
  std::ofstream out(argv[3]);out<<std::setprecision(24);
  out<<"xminus_J,xplus_J,yminus_J,yplus_J,interior_cells,boundary_cells\n";
  for(unsigned a=0;a<4;++a){Sum v;for(auto& s:sums)v.add(s[a].v);if(a)out<<',';out<<v.v;}
  out<<','<<interior<<','<<boundary<<'\n';
  munmap(raw,bytes);close(fd);return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
