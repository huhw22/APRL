#ifndef APRL_OPEN_TRANSVERSE_TRACE_H
#define APRL_OPEN_TRANSVERSE_TRACE_H

#include "compute_threads.h"
#include "eb_field.h"
#include "reference_collective.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace aprl { namespace verification {

// The scalar Dirichlet data and the one external node required by the Cowan
// field stencil. Each tile fixes x or y, and spans the other coordinate and z.
// There is only one lookup convention at corners, including external corners.
class OpenTransverseTrace {
public:
  struct Sources {
    std::vector<double> x,y,z,q;
    // Coordinates relative to the grid origin, with z stretched by gamma.
    double potential(double xp,double yp,double zp) const {
      double sum=0.;
#ifdef APRL_HAS_OPENMP
#pragma omp simd reduction(+:sum)
#endif
      for(std::size_t s=0;s<q.size();++s) {
        const double dx=xp-x[s], dy=yp-y[s], dz=zp-z[s];
        sum+=q[s]/std::sqrt(dx*dx+dy*dy+dz*dz);
      }
      return sum;
    }
  };
  std::vector<double> values;
  double maxRelativeError=0.;
  std::size_t directChecks=0;
  std::size_t fallbackTiles=0,sourceNodes=0,coarsePoints=0;
  std::size_t localPotentialEvaluations=0;
  int computeRanks=1;
  void configure(const EBGridGeometry& geometry) {
    nx=static_cast<int>(geometry.nx); ny=static_cast<int>(geometry.ny);
    nz=static_cast<int>(geometry.nz);
    width=4*(std::size_t(nx)+3+ny+3);
    values.assign(width*(std::size_t(nz)+1),0.);
  }
  double at(int i,int j,int k) const {
    if(k<=0 || k>=nz) return 0.; // existing longitudinal Dirichlet end planes
    const int tx=slot(i,nx),ty=slot(j,ny);
    if(tx>=0) return values[k*width+tx*(ny+3)+j+1];
    if(ty>=0) return values[k*width+4*(ny+3)+ty*(nx+3)+i+1];
    throw std::logic_error("Open trace requested away from a transverse boundary");
  }
  // Collective when comm has several ranks. Only rank zero supplies sources
  // and retains the resulting trace. Each target is evaluated on one rank;
  // its source loop, interpolation and summation order are unchanged.
  void build(const Sources& rootSource,const EBGridGeometry& g,double gamma,
             MPI_Comm comm=MPI_COMM_SELF) {
    if(comm==MPI_COMM_NULL) throw std::invalid_argument("Null open trace communicator");
    int rank=0;
    MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&computeRanks);
    maxRelativeError=0.;directChecks=fallbackTiles=coarsePoints=localPotentialEvaluations=0;
    Sources received;
    unsigned long long count=rank==0?rootSource.q.size():0;
    reference_detail::collectively(comm,[&] {
      if(rank==0 && (rootSource.x.size()!=count || rootSource.y.size()!=count || rootSource.z.size()!=count))
        throw std::invalid_argument("Inconsistent open trace source arrays");
    });
    MPI_Bcast(&count,1,MPI_UNSIGNED_LONG_LONG,0,comm);
    sourceNodes=static_cast<std::size_t>(count);
    if(count==0) return;
    reference_detail::collectively(comm,[&] {
      if(rank!=0) {
        received.x.resize(count);received.y.resize(count);
        received.z.resize(count);received.q.resize(count);
      }
    });
    // MPI_Bcast does not modify the root buffer.
    const Sources& source=rank==0?rootSource:received;
    for(const auto* v:{&source.x,&source.y,&source.z,&source.q})
      reference_detail::broadcast(const_cast<double*>(v->data()),sourceNodes,comm);
    reference_detail::collectively(comm,[&] {
      if(rank==0) {
        if(values.size()!=4*(g.nx+3+g.ny+3)*(g.nz+1))
          throw std::invalid_argument("Open trace was not configured for this geometry");
        if(!(gamma>0.) || !std::isfinite(gamma)) throw std::invalid_argument("Invalid open trace gamma");
        for(std::size_t s=0;s<sourceNodes;++s)
          if(!std::isfinite(source.x[s]) || !std::isfinite(source.y[s]) ||
             !std::isfinite(source.z[s]) || !std::isfinite(source.q[s]))
            throw std::invalid_argument("Nonfinite open trace source");
      }
    });
    nx=static_cast<int>(g.nx);ny=static_cast<int>(g.ny);nz=static_cast<int>(g.nz);
    width=4*(std::size_t(nx)+3+ny+3);
    const auto threads=computeThreadCount();
    (void)threads;
    const double hz=g.dz*gamma;
    for(int tile=0;tile<8;++tile) {
      const bool fixedX=tile<4;
      const int t=tile%4;
      const int n=fixedX?nx:ny, ns=fixedX?ny:nx;
      const double hf=fixedX?g.dx:g.dy, hs=fixedX?g.dy:g.dx;
      const double fixed=(t<2?t-1:n+t-2)*hf;
      const auto& coords=fixedX?source.x:source.y;
      const auto extremes=std::minmax_element(coords.begin(),coords.end());
      const double dist=t<2?*extremes.first-fixed:fixed-*extremes.second;
      reference_detail::collectively(comm,[&] {
        if(!(dist>0.)) throw std::invalid_argument("Beam support must be inside each open outer face");
      });
      // Resolve the smooth Coulomb potential by at least ~13 samples per
      // closest-source distance, or use all native nodes if they are coarser.
      const int ds=std::max(1,static_cast<int>(std::min(0.075*dist/hs,(ns+2)/7.)));
      const int dz=std::max(1,static_cast<int>(std::min(0.075*dist/hz,nz/7.)));
      // Extend the table on both ends, so every native sample uses a centered stencil.
      const int cs=(ns+2+ds-1)/ds+8, cz=(nz+dz-1)/dz+8;
      std::vector<double> coarse,rows;
      coarsePoints+=std::size_t(cs)*cz;
      evaluate(std::size_t(cs)*cz,coarse,comm,[&](std::size_t index) {
        const int k=index/cs,p=index%cs;
        const double tangent=((p-3)*ds-1)*hs;
        return source.potential(fixedX?fixed:tangent,fixedX?tangent:fixed,(k-3)*dz*hz);
      });
      const std::size_t offset=fixedX?t*(ny+3):4*(ny+3)+t*(nx+3);
      reference_detail::collectively(comm,[&] {
        if(rank!=0) return;
        const auto sp=weights(ns+3,ds,cs),zp=weights(nz+1,dz,cz);
        rows.resize(std::size_t(ns+3)*cz);
        for(int k=0;k<cz;++k) for(int p=0;p<ns+3;++p)
          rows[k*(ns+3)+p]=sp[p].eval(&coarse[k*cs],1);
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for num_threads(threads) if(threads>1)
#endif
        for(int k=1;k<nz;++k) for(int p=0;p<ns+3;++p)
          values[k*width+offset+p]=zp[k].eval(&rows[p],ns+3);
      });
      double tileError=0.;
      // Sample the interpolant away from its nodes, including the near-source region.
      const auto zlimits=std::minmax_element(source.z.begin(),source.z.end());
      const auto& other=fixedX?source.y:source.x;
      const auto slimits=std::minmax_element(other.begin(),other.end());
      const auto checkIndex=[&](int a,int b) {
        const double znear=(*zlimits.first+*zlimits.second)*0.5/hz;
        const double snear=(*slimits.first+*slimits.second)*0.5/hs;
        int k=a<7?int(znear+(a-3)*dist/(3.*hz)+0.4*dz):int((a-6.4)*nz/4.2);
        int p=b<5?int(snear+(b-2)*dist/(3.*hs)+0.35*ds)+1:int((b-4.6)*(ns+2)/3.6);
        k=std::max(1,std::min(nz-1,k));p=std::max(0,std::min(ns+2,p));
        return std::array<int,2>{{k,p}};
      };
      std::vector<double> exactChecks;
      evaluate(99,exactChecks,comm,[&](std::size_t index) {
        const auto kp=checkIndex(index/9,index%9);
        const double tangent=(kp[1]-1)*hs;
        return source.potential(fixedX?fixed:tangent,fixedX?tangent:fixed,kp[0]*hz);
      });
      reference_detail::collectively(comm,[&] {
        if(rank!=0) return;
        for(int a=0;a<11;++a) for(int b=0;b<9;++b) {
          const auto kp=checkIndex(a,b);
          const double exact=exactChecks[a*9+b],actual=values[kp[0]*width+offset+kp[1]];
          if(!std::isfinite(exact)||!std::isfinite(actual))throw std::runtime_error("Nonfinite open wall trace");
          if(exact!=0.) tileError=std::max(tileError,std::abs((actual-exact)/exact));
          else if(actual!=0.) tileError=std::numeric_limits<double>::infinity();
        }
      });
      directChecks+=99;
      MPI_Bcast(&tileError,1,MPI_DOUBLE,0,comm);
      if(tileError>2.e-7) {
        ++fallbackTiles;
        std::vector<double> native;
        evaluate(std::size_t(nz-1)*(ns+3),native,comm,[&](std::size_t index) {
          const int k=1+index/(ns+3),p=index%(ns+3);
          const double tangent=(p-1)*hs;
          return source.potential(fixedX?fixed:tangent,fixedX?tangent:fixed,k*hz);
        });
        if(rank==0) for(int k=1;k<nz;++k) for(int p=0;p<ns+3;++p)
          values[k*width+offset+p]=native[std::size_t(k-1)*(ns+3)+p];
      } else maxRelativeError=std::max(maxRelativeError,tileError);
    }
  }
private:
  template<class Potential> void evaluate(std::size_t count,std::vector<double>& result,
                                         MPI_Comm comm,Potential potential) {
    int rank=0,ranks=1;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
    // Quotient/remainder partition avoids both int counts and count*rank overflow.
    const auto begin=[&](int r) {return (count/ranks)*r+std::min(count%ranks,std::size_t(r));};
    const std::size_t first=begin(rank),take=begin(rank+1)-first;
    std::vector<double> local;
    reference_detail::collectively(comm,[&] {
      local.resize(take);if(rank==0)result.resize(count);
    });
    const int threads=computeThreadCount();(void)threads;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for num_threads(threads) if(threads>1)
#endif
    for(std::size_t p=0;p<take;++p)local[p]=potential(first+p);
    localPotentialEvaluations+=take;
    const std::size_t batch=4*1024*1024;
    if(rank==0) {
      std::copy(local.begin(),local.end(),result.begin());
      for(int peer=1;peer<ranks;++peer)for(std::size_t p=begin(peer);p<begin(peer+1);p+=batch)
        MPI_Recv(result.data()+p,static_cast<int>(std::min(batch,begin(peer+1)-p)),MPI_DOUBLE,
                 peer,0,comm,MPI_STATUS_IGNORE);
    } else for(std::size_t p=0;p<take;p+=batch)
      MPI_Send(local.data()+p,static_cast<int>(std::min(batch,take-p)),MPI_DOUBLE,0,0,comm);
  }
  int nx=0,ny=0,nz=0;
  std::size_t width=0;
  static int slot(int i,int limit) {return i<=0?i+1:(i>=limit?i-limit+2:-1);}
  struct Interpolation {
    int first=0;
    std::array<double,8> w{};
    double eval(const double* v,int stride) const {
      double sum=0.;for(int j=0;j<8;++j)sum+=w[j]*v[(first+j)*stride];
      return sum;
    }
  };
  static std::vector<Interpolation> weights(int count,int stride,int coarse) {
    std::vector<Interpolation> out(count);
    for(int i=0;i<count;++i) {
      const double x=double(i)/stride+3.;
      auto& s=out[i];s.first=std::max(0,std::min(coarse-8,int(x)-3));
      for(int j=0;j<8;++j) {
        double v=1.;for(int l=0;l<8;++l)if(l!=j)v*=(x-(s.first+l))/(j-l);
        s.w[j]=v;
      }
    }
    return out;
  }
};

}} // namespace aprl::verification
#endif
