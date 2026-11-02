#ifndef APRL_VERIFICATION_YEE_SBP_OBSERVER_H
#define APRL_VERIFICATION_YEE_SBP_OBSERVER_H
#include "yee_sbp_geometry.h"
#include <mpi.h>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
namespace aprl {
class YeeReadOnlyCollar {
  std::array<const YeeComponent*,3> values_;std::array<std::vector<double>,3> below_,above_;std::size_t nz_;
public:
  YeeReadOnlyCollar(std::array<const YeeComponent*,3> values,std::size_t nz,MPI_Comm comm):values_(values),nz_(nz){
    int rank=0,size=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);const auto lower=rank?rank-1:MPI_PROC_NULL,upper=rank+1<size?rank+1:MPI_PROC_NULL;
    for(unsigned a=0;a<3;++a){const auto n=values_[a]->nx()*values_[a]->ny();below_[a].resize(n);above_[a].resize(n);
      int rc=MPI_Sendrecv(values_[a]->data()+(nz_-1)*n,n,MPI_DOUBLE,upper,23410+a,below_[a].data(),n,MPI_DOUBLE,lower,23410+a,comm,MPI_STATUS_IGNORE);if(rc!=MPI_SUCCESS)throw std::runtime_error("Yee diagnostic lower collar");
      rc=MPI_Sendrecv(values_[a]->data(),n,MPI_DOUBLE,lower,23420+a,above_[a].data(),n,MPI_DOUBLE,upper,23420+a,comm,MPI_STATUS_IGNORE);if(rc!=MPI_SUCCESS)throw std::runtime_error("Yee diagnostic upper collar");
    }
  }
  long double operator()(unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k)const{
    const auto& v=*values_[a];if(i<0||j<0||i>=(std::ptrdiff_t)v.nx()||j>=(std::ptrdiff_t)v.ny())throw std::runtime_error("Yee diagnostic transverse index");
    if(k>=0&&k<(std::ptrdiff_t)v.nz())return v(i,j,k);
    if(k==-1)return below_[a][j*v.nx()+i];
    if(k==(std::ptrdiff_t)nz_&&v.nz()==nz_)return above_[a][j*v.nx()+i];
    throw std::runtime_error("Yee diagnostic longitudinal index");
  }
};
class YeeSBPObserver {
  EBGridGeometry full_;std::ptrdiff_t start_=0;std::unique_ptr<YeeSBPGeometry> geometry_;MPI_Comm comm_;int rank_=0;std::ofstream state_,flux_;
  std::array<long double,8> reduce(const std::array<long double,8>& a){std::array<long double,8>b{};MPI_Allreduce(a.data(),b.data(),8,MPI_LONG_DOUBLE,MPI_SUM,comm_);return b;}
public:
  YeeSBPObserver(const std::string& dir,const EBGridGeometry& full,std::size_t offset,std::size_t globalNz,std::size_t gx,std::size_t gy,std::size_t gz,MPI_Comm comm):full_(full),comm_(comm){
    MPI_Comm_rank(comm,&rank_);const auto lo=std::max(offset,gz),hi=std::min(offset+full.nz,globalNz-gz);
    if(hi>lo){start_=lo-offset;EBGridGeometry sub(full.nx,full.ny,hi-lo,full.dx,full.dy,full.dz,full.dt,EBMaxwellSolver::Yee);geometry_.reset(new YeeSBPGeometry(sub,gx,gy,lo-gz,globalNz-2*gz));}
    if(rank_==0){state_.open(dir+"/sbp-state.csv");flux_.open(dir+"/sbp-flux.csv");if(!state_||!flux_)throw std::runtime_error("Cannot open Yee SBP diagnostics");state_<<std::setprecision(21)<<"step,U_sync_J,H_cross_J,B_time_variance_J,cP_sync_J,stock_identity_error_J\n";flux_<<std::setprecision(21)<<"step,Q_H_J,Q_P_J,Gauss_E_average_J,Gauss_B_J,Q_P_at_mean_J,stress_time_variance_direct_J,old_Q_S_J,old_Q_cT_J\n";}
  }
  void state(unsigned n,const EBFieldGrid& view,const std::array<YeeComponent,3>& minus,const EBFieldGrid& plus){
    std::array<long double,8> v{};
    if(geometry_){const auto& g=*geometry_;const YeeComponent* ee[3]={&view.ex(),&view.ey(),&view.ez()},*bb[3]={&view.bx(),&view.by(),&view.bz()},*bp[3]={&plus.bx(),&plus.by(),&plus.bz()};
      auto e=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k)->long double{return (*ee[a])(i,j,k+start_);};
      auto b=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k)->long double{return (*bb[a])(i,j,k+start_);};
      auto bm=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k)->long double{return minus[a](i,j,k+start_);};
      auto bn=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k)->long double{return (*bp[a])(i,j,k+start_);};
      auto d=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return bn(a,i,j,k)-bm(a,i,j,k);};
      const auto ue=.5L*SI::epsilon0*g.electricInner(e,e);v[0]=ue+.5L/SI::mu0*g.magneticInner(b,b);v[1]=ue+.5L/SI::mu0*g.magneticInner(bm,bn);v[2]=.125L/SI::mu0*g.magneticInner(d,d);v[3]=g.momentum(e,b);v[4]=v[0]-v[1]-v[2];
    }
    const auto a=reduce(v);if(!rank_){state_<<n;for(unsigned j=0;j<5;++j)state_<<','<<a[j];state_<<'\n';if(!state_)throw std::runtime_error("Yee state write");}
  }
  void transition(unsigned n,const std::array<YeeComponent,3>& old,const EBFieldGrid& f){
    std::array<const YeeComponent*,3> ea{{&old[0],&old[1],&old[2]}},eb{{&f.ex(),&f.ey(),&f.ez()}},ba{{&f.bx(),&f.by(),&f.bz()}};
    // All ranks participate even when this rank lies outside the control volume.
    YeeReadOnlyCollar oldE(ea,full_.nz,comm_),newE(eb,full_.nz,comm_),B(ba,full_.nz,comm_);std::array<long double,8> v{};
    if(geometry_){const auto& g=*geometry_;
      auto e=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return oldE(a,i,j,k+start_);};
      auto en=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return newE(a,i,j,k+start_);};
      auto b=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return B(a,i,j,k+start_);};
      auto em=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return .5L*(e(a,i,j,k)+en(a,i,j,k));};
      auto de=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return en(a,i,j,k)-e(a,i,j,k);};
      const auto qe=g.sum(g.electricStress(e)),qn=g.sum(g.electricStress(en)),qb=g.sum(g.magneticStress(b));
      v[0]=g.sum(g.energyPort(em,b));v[1]=.5L*(qe+qn)+qb;v[2]=.5L*(g.gaussElectricImpulse(e)+g.gaussElectricImpulse(en));v[3]=g.gaussMagneticImpulse(b);v[4]=g.sum(g.electricStress(em))+qb;v[5]=-.25L*g.sum(g.electricStress(de));
    }
    pending_=reduce(v);pendingStep_=n;
  }
  void oldFlux(unsigned n,long double s,long double t){
    if(n!=pendingStep_)throw std::runtime_error("Yee flux stage mismatch");if(!rank_){flux_<<n;for(unsigned j=0;j<6;++j)flux_<<','<<pending_[j];flux_<<','<<s<<','<<t<<'\n';if(!flux_)throw std::runtime_error("Yee flux write");}
  }
  void close(){if(!rank_){state_.flush();flux_.flush();if(!state_||!flux_)throw std::runtime_error("Yee SBP flush");state_.close();flux_.close();}}
private:std::array<long double,8> pending_{};unsigned pendingStep_=0;
};
}
#endif
