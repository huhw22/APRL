#ifndef APRL_DISCRETE_FIELD_LEDGER_H
#define APRL_DISCRETE_FIELD_LEDGER_H
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <mpi.h>
#include "eb_field.h"

namespace aprl {
// Read-only numerical-operator ledger. Q is the paired-curl boundary/interface
// exchange, not radiation power. Time and Cowan terms remain explicit. The
// operator residual is never used to construct/subtract a correction term.
// Physical continuous-energy estimates remain in the existing v5 HDF5 ledger.
class DiscreteFieldLedger {
public:
  struct Record {
    long double begin=0.,end=0.,exchange=0.,jOld=0.,jMid=0.,unappliedJMid=0.;
    long double timeE=0.,timeB=0.,currentPairing=0.,cowan=0.,cpmlB=0.,cpmlE=0.;
    long double inlet=0.,remap=0.,residual=0.,eError=0.,bError=0.,eErrorWork=0.,bErrorWork=0.;
    long double eScale=0.,bScale=0.;
  };
  DiscreteFieldLedger(const std::string& directory,const EBGridGeometry& g,
      std::size_t globalNz,std::size_t offset,std::size_t cx,std::size_t cy,MPI_Comm comm)
    :g_(g),globalNz_(globalNz),offset_(offset),cx_(cx),cy_(cy),comm_(comm) {
    MPI_Comm_rank(comm_,&rank_);MPI_Comm_size(comm_,&size_);
    lower_=rank_?rank_-1:MPI_PROC_NULL;upper_=rank_+1<size_?rank_+1:MPI_PROC_NULL;
    const std::size_t count=(g_.nx+1)*g_.ny+g_.nx*(g_.ny+1);
    send_.resize(count);below_.resize(count);above_.resize(count);
    std::ostringstream name;name<<directory<<"/discrete-field-rank-"<<std::setw(5)<<std::setfill('0')<<rank_;
    basename_=name.str();std::ifstream existing(basename_+".csv");
    if(existing.good())throw std::runtime_error("Discrete field ledger already exists");
    file_.open(basename_+".csv");if(!file_)throw std::runtime_error("Cannot create discrete field ledger");
    file_<<std::setprecision(21)<<"step,native_begin_J,native_end_J,discrete_curl_exchange_J,applied_current_old_J,applied_current_mid_J,unapplied_current_mid_J,time_electric_J,time_magnetic_J,current_time_pairing_J,cowan_operator_J,magnetic_cpml_J,electric_cpml_J,inlet_overwrite_J,remap_J,field_operator_residual_J,max_electric_prediction_error_V_per_m,max_magnetic_prediction_error_T,electric_prediction_error_work_J,magnetic_prediction_error_work_J,electric_prediction_scale_V_per_m,magnetic_prediction_scale_T\n";
  }
  void begin(const EBFieldGrid& f,std::size_t step) {
    if(phase_!=0)throw std::logic_error("Incomplete discrete ledger field step");
    oldE_[0]=f.ex();oldE_[1]=f.ey();oldE_[2]=f.ez();
    oldB_[0]=f.bx();oldB_[1]=f.by();oldB_[2]=f.bz();
    record_=Record();const auto u=energy(f);record_.begin=u[0]+u[1];step_=step+1;phase_=1;
  }
  static void observe(void* context,const char* stage,const EBFieldGrid& f,const EBElectricHaloView* h) {
    auto& self=*static_cast<DiscreteFieldLedger*>(context);const std::string s(stage);
    if(s=="magnetic_curl")self.magnetic(f,*h);
    else if(s=="magnetic_cpml") {self.expect(2);self.record_.cpmlB=self.energy(f)[1]-self.rawMagneticEnergy_;self.phase_=3;}
    else if(s=="electric_curl_local") {self.expect(3);self.preElectricCpml_=self.energy(f)[0];self.phase_=4;}
    else if(s=="electric_cpml") {
      self.expect(4);self.postElectricCpml_=self.energy(f)[0];
      self.record_.cpmlE=self.postElectricCpml_-self.preElectricCpml_;
      if(self.rank_==0) {const YeeComponent* e[2]={&f.ex(),&f.ey()};for(unsigned a=0;a<2;++a)
        self.inletE_[a].assign(e[a]->data(),e[a]->data()+e[a]->nx()*e[a]->ny());}
      self.phase_=5;
    } else if(s=="inlet_overwrite") {self.expect(5);self.record_.inlet=self.energy(f)[0]-self.postElectricCpml_;self.phase_=6;}
    else throw std::logic_error("Unknown discrete ledger field stage");
  }
  void finishElectric(const EBFieldGrid& f) {
    expect(6);exchangeMagnetic(f);
    const YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};
    const YeeComponent* j[3]={&f.jx(),&f.jy(),&f.jz()};
    long double ecurl=0.,jo=0.,jm=0.,ju=0.,te=0.,error=0.,errorWork=0.,scale=0.;
    const long double volume=static_cast<long double>(g_.dx)*g_.dy*g_.dz;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for reduction(+:ecurl,jo,jm,ju,te,errorWork) reduction(max:error,scale) schedule(static)
#endif
    for(std::ptrdiff_t kk=0;kk<static_cast<std::ptrdiff_t>(g_.nz);++kk) {
      const auto k=static_cast<std::size_t>(kk);
      for(std::size_t y=cy_;y<g_.ny-cy_;++y)for(std::size_t x=cx_;x<g_.nx-cx_;++x)
        for(unsigned a=0;a<3;++a)for(unsigned p=0;p<2;++p)for(unsigned q=0;q<2;++q) {
          const auto i=x+(a==0?0:p),v=y+(a==0?p:(a==1?0:q)),z=k+(a==2?0:q);
          const long double eo=oldE_[a](i,v,z);
          const long double en=(rank_==0&&a<2&&z==0)?inletE_[a][v*e[a]->nx()+i]:(*e[a])(i,v,z);
          const long double de=en-eo,jv=(*j[a])(i,v,z);
          const bool applied=!(a<2&&(offset_+z==0||offset_+z==globalNz_));
          const long double w=.25L*volume;
          const long double inc=applied?electricCurlIncrement(f,a,i,v,z):0.;
          const long double source=applied?static_cast<long double>(g_.dt)/SI::epsilon0*jv:0.;
          ecurl+=w*SI::epsilon0*eo*inc;te+=w*.5L*SI::epsilon0*de*de;
          const long double work=w*g_.dt*jv*(eo+en)/2;
          if(applied) {jo+=w*g_.dt*jv*eo;jm+=work;}else ju+=work;
          const long double defect=de-(inc-source);error=std::max(error,std::abs(defect));
          scale=std::max(scale,std::abs(eo)+std::abs(en)+std::abs(inc)+std::abs(source));
          errorWork+=w*SI::epsilon0*eo*defect;
        }
    }
    record_.jOld=jo;record_.jMid=jm;record_.unappliedJMid=ju;record_.timeE=te;
    record_.currentPairing=jm-jo;record_.exchange=-ecurl-magneticYeeWork_;
    record_.eError=error;record_.eErrorWork=errorWork;
    record_.eScale=scale;
    const auto u=energy(f);endBeforeRemap_=u[0]+u[1];phase_=7;
  }
  void afterRemap(const EBFieldGrid& f) {
    expect(7);const auto u=energy(f);record_.end=u[0]+u[1];
    record_.remap=record_.end-endBeforeRemap_;
    record_.residual=record_.end-record_.begin+record_.exchange+record_.jMid-
      record_.timeE-record_.timeB-record_.currentPairing-record_.cowan-
      record_.cpmlB-record_.cpmlE-record_.inlet-record_.remap;
    file_<<step_<<','<<record_.begin<<','<<record_.end<<','<<record_.exchange<<','<<record_.jOld<<','<<record_.jMid<<','<<record_.unappliedJMid
      <<','<<record_.timeE<<','<<record_.timeB<<','<<record_.currentPairing<<','<<record_.cowan<<','<<record_.cpmlB<<','<<record_.cpmlE
      <<','<<record_.inlet<<','<<record_.remap<<','<<record_.residual<<','<<record_.eError<<','<<record_.bError<<','<<record_.eErrorWork<<','<<record_.bErrorWork<<','<<record_.eScale<<','<<record_.bScale<<'\n';
    if(!file_)throw std::runtime_error("Cannot write discrete field ledger");
    phase_=0;
  }
  const Record& lastRecord() const {return record_;}
  void close(bool complete) {
    expect(0);file_.flush();if(!file_)throw std::runtime_error("Cannot flush discrete field ledger");file_.close();
    std::ofstream marker(basename_+".complete");marker<<(complete?1:0)<<' '<<step_<<'\n';
    if(!marker)throw std::runtime_error("Cannot complete discrete field ledger");
  }
private:
  void expect(unsigned p) const {if(phase_!=p)throw std::logic_error("Discrete ledger stage/order mismatch");}
  std::array<long double,2> energy(const EBFieldGrid& f) const {
    const YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};const YeeComponent* b[3]={&f.bx(),&f.by(),&f.bz()};
    const long double volume=static_cast<long double>(g_.dx)*g_.dy*g_.dz;long double ue=0.,ub=0.;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for reduction(+:ue,ub) schedule(static)
#endif
    for(std::ptrdiff_t kk=0;kk<static_cast<std::ptrdiff_t>(g_.nz);++kk) {
      const auto k=static_cast<std::size_t>(kk);
      for(std::size_t y=cy_;y<g_.ny-cy_;++y)for(std::size_t x=cx_;x<g_.nx-cx_;++x)for(unsigned a=0;a<3;++a) {
        for(unsigned p=0;p<2;++p)for(unsigned q=0;q<2;++q) {
          const long double v=(*e[a])(x+(a==0?0:p),y+(a==0?p:(a==1?0:q)),k+(a==2?0:q));ue+=volume*SI::epsilon0*v*v/8;
        }
        for(unsigned p=0;p<2;++p) {const long double v=(*b[a])(x+(a==0?p:0),y+(a==1?p:0),k+(a==2?p:0));ub+=volume*v*v/(4*SI::mu0);}
      }
    }
    return {{ue,ub}};
  }
  void magnetic(const EBFieldGrid& f,const EBElectricHaloView& halo) {
    expect(1);const YeeComponent* b[3]={&f.bx(),&f.by(),&f.bz()};
    const long double volume=static_cast<long double>(g_.dx)*g_.dy*g_.dz;
    long double tb=0.,wc=0.,wy=0.,error=0.,errorWork=0.,scale=0.;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for reduction(+:tb,wc,wy,errorWork) reduction(max:error,scale) schedule(static)
#endif
    for(std::ptrdiff_t kk=0;kk<static_cast<std::ptrdiff_t>(g_.nz);++kk) {
      const auto k=static_cast<std::size_t>(kk);
      for(std::size_t y=cy_;y<g_.ny-cy_;++y)for(std::size_t x=cx_;x<g_.nx-cx_;++x)
        for(unsigned a=0;a<3;++a)for(unsigned p=0;p<2;++p) {
          const auto i=x+(a==0?p:0),j=y+(a==1?p:0),z=k+(a==2?p:0);
          const long double bn=(*b[a])(i,j,z),db=bn-oldB_[a](i,j,z);
          const long double pred=f.diagnosticMagneticIncrement(a,i,j,z,halo);
          const long double yee=f.diagnosticMagneticIncrement(a,i,j,z,halo,true),w=.5L*volume/SI::mu0;
          tb-=w*db*db/2;wc+=w*bn*(pred-yee);wy+=w*bn*yee;
          error=std::max(error,std::abs(db-pred));errorWork+=w*bn*(db-pred);
          scale=std::max(scale,std::abs(bn)+std::abs(static_cast<long double>(oldB_[a](i,j,z)))+std::abs(pred));
        }
    }
    record_.timeB=tb;record_.cowan=wc;magneticYeeWork_=wy;
    record_.bError=error;record_.bErrorWork=errorWork;rawMagneticEnergy_=energy(f)[1];phase_=2;
    record_.bScale=scale;
  }
  void exchangeMagnetic(const EBFieldGrid& f) {
    if(size_==1)return;
    const std::size_t bx=(g_.nx+1)*g_.ny,by=g_.nx*(g_.ny+1);const int count=static_cast<int>(bx+by);
    auto pack=[&](std::size_t k) {std::copy(f.bx().data()+k*bx,f.bx().data()+(k+1)*bx,send_.begin());std::copy(f.by().data()+k*by,f.by().data()+(k+1)*by,send_.begin()+bx);};
    pack(g_.nz-1);int status=MPI_Sendrecv(send_.data(),count,MPI_DOUBLE,upper_,19441,below_.data(),count,MPI_DOUBLE,lower_,19441,comm_,MPI_STATUS_IGNORE);
    if(status!=MPI_SUCCESS){MPI_Abort(comm_,status);throw std::runtime_error("Discrete ledger lower B halo failed");}
    pack(0);status=MPI_Sendrecv(send_.data(),count,MPI_DOUBLE,lower_,19442,above_.data(),count,MPI_DOUBLE,upper_,19442,comm_,MPI_STATUS_IGNORE);
    if(status!=MPI_SUCCESS){MPI_Abort(comm_,status);throw std::runtime_error("Discrete ledger upper B halo failed");}
  }
  long double electricCurlIncrement(const EBFieldGrid& f,unsigned a,std::size_t i,std::size_t j,std::size_t k) const {
    const std::size_t bx=(g_.nx+1)*g_.ny;
    auto B=[&](unsigned c,std::size_t x,std::size_t y,std::ptrdiff_t z)->long double {
      if(c==2)return f.bz()(x,y,static_cast<std::size_t>(z));
      const std::size_t width=c==0?g_.nx+1:g_.nx,index=(c==0?0:bx)+y*width+x;
      if(z<0)return below_[index];
      if(z>=static_cast<std::ptrdiff_t>(g_.nz))return above_[index];
      return c==0?f.bx()(x,y,z):f.by()(x,y,z);
    };
    const long double scale=static_cast<long double>(g_.dt)/(SI::mu0*SI::epsilon0);const auto z=static_cast<std::ptrdiff_t>(k);
    if(a==0)return scale*((B(2,i,j,z)-B(2,i,j-1,z))/g_.dy-(B(1,i,j,z)-B(1,i,j,z-1))/g_.dz);
    if(a==1)return scale*((B(0,i,j,z)-B(0,i,j,z-1))/g_.dz-(B(2,i,j,z)-B(2,i-1,j,z))/g_.dx);
    return scale*((B(1,i,j,z)-B(1,i-1,j,z))/g_.dx-(B(0,i,j,z)-B(0,i,j-1,z))/g_.dy);
  }
  EBGridGeometry g_;std::size_t globalNz_,offset_,cx_,cy_,step_=0;MPI_Comm comm_;
  int rank_=0,size_=1,lower_=MPI_PROC_NULL,upper_=MPI_PROC_NULL;unsigned phase_=0;
  Record record_;YeeComponent oldE_[3],oldB_[3];std::vector<Double> inletE_[2],send_,below_,above_;
  long double rawMagneticEnergy_=0.,preElectricCpml_=0.,postElectricCpml_=0.,endBeforeRemap_=0.,magneticYeeWork_=0.;
  std::string basename_;std::ofstream file_;
};
}
#endif
