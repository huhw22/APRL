#ifndef APRL_VERIFICATION_WHOLE_COWAN_BUDGET_H
#define APRL_VERIFICATION_WHOLE_COWAN_BUDGET_H
#include "yee_sbp_observer.h"
#include <sstream>
namespace aprl {
class WholeCowanBudget {
  using Arrays=std::array<YeeComponent,3>;
  YeeSBPGeometry geometry_;MPI_Comm comm_;int rank_=0;std::size_t step_=0;std::string base_;std::ofstream out_;bool ready_=false;
  Arrays oldE_,oldB_,oldK_,newK_,rawB_,cpmlB_,plusB_,primaryE_,cpmlE_,inletE_,ghostB_;YeeComponent referenceJz_;
  long double beginH_=0.,beginP_=0.,preRemapH_=0.,preRemapP_=0.;
  static std::array<const YeeComponent*,3> electric(const EBFieldGrid& f){return {{&f.ex(),&f.ey(),&f.ez()}};}
  static std::array<const YeeComponent*,3> magnetic(const EBFieldGrid& f){return {{&f.bx(),&f.by(),&f.bz()}};}
  static void copy(std::array<const YeeComponent*,3> p,Arrays& a){for(unsigned j=0;j<3;++j)a[j]=*p[j];}
  static std::array<const YeeComponent*,3> pointers(const Arrays& a){return {{&a[0],&a[1],&a[2]}};}
  template<class E,class B,class K>std::array<long double,2> stock(E e,B b,K k){
    return {{.5L*SI::epsilon0*geometry_.electricInner(e,e)+.5L/SI::mu0*geometry_.magneticInner(b,b)+.5L/SI::mu0*geometry_.magneticInner(b,k),geometry_.momentum(e,b)+.5L*geometry_.momentum(e,k)}};
  }
  void calculateK(const EBFieldGrid& f,Arrays& k){
    YeeReadOnlyCollar collar(electric(f),geometry_.g.nz,comm_);const auto h=collar.electricHalos();const auto& g=geometry_.g;
    for(unsigned a=0;a<3;++a){
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
      for(std::size_t z=0;z<g.nz+(a==2);++z)for(std::size_t y=geometry_.cy;y<g.ny-geometry_.cy+(a==1);++y)for(std::size_t x=geometry_.cx;x<g.nx-geometry_.cx+(a==0);++x)k[a](x,y,z)=f.diagnosticMagneticIncrement(a,x,y,z,h,false);
    }
  }
  template<class B,class E,class S>void magneticSource(std::vector<long double>& v,unsigned i,B bm,E eo,S source){v[i]=geometry_.magneticInner(bm,source)/SI::mu0;v[i+1]=geometry_.momentum(eo,source);}
  template<class E,class B,class S>void electricSource(std::vector<long double>& v,unsigned i,E em,B b,S source){v[i]=SI::epsilon0*geometry_.electricInner(em,source);v[i+1]=geometry_.momentum(source,b);}
  void write(const char* stage,long double h,long double p,const std::vector<long double>& values={}){
    out_<<step_<<','<<stage<<','<<h<<','<<p;
    for(unsigned i=0;i<28;++i)out_<<','<<(i<values.size()?values[i]:0.L);
    out_<<'\n';if(!out_)throw std::runtime_error("Whole Cowan audit write failed");
  }
public:
  WholeCowanBudget(const std::string& dir,int rank,const EBGridGeometry& g,std::size_t offset,std::size_t nz,std::size_t gx,std::size_t gy,MPI_Comm comm):geometry_(g,gx,gy,offset,nz),comm_(comm),rank_(rank){
    if(g.solver!=EBMaxwellSolver::CowanZ)throw std::runtime_error("Whole Cowan audit requires original Cowan solver");
    for(unsigned a=0;a<3;++a){for(Arrays* p:{&oldK_,&newK_,&ghostB_})(*p)[a].resize(g.nx+(a==0),g.ny+(a==1),g.nz+(a==2));}
    std::ostringstream name;name<<dir<<"/whole-cowan-rank-"<<std::setw(5)<<std::setfill('0')<<rank;base_=name.str();std::ifstream existing(base_+".csv");if(existing.good())throw std::runtime_error("Whole Cowan audit already exists");out_.open(base_+".csv");if(!out_)throw std::runtime_error("Cannot open whole Cowan audit");
    out_<<std::setprecision(21)<<"step,stage,H_ref_J,cP_ref_J,Q_H_J,Q_P_J,C_H_J,C_P_J,GaussE_J,GaussB_J,WJ_particle_J,PJ_particle_J";
    for(const char* key:{"magnetic_supply","magnetic_cpml","magnetic_sync","electric_hold","electric_cpml","electric_inlet","electric_sync","reference_current"})out_<<','<<key<<"_U_J,"<<key<<"_P_J";
    out_<<",remap_H_J,remap_P_J,field_R_H_J,field_R_P_J\n";
  }
  void begin(const EBFieldGrid& f,std::size_t step){
    step_=step+1;copy(electric(f),oldE_);copy(magnetic(f),oldB_);if(!ready_)calculateK(f,oldK_);ready_=true;referenceJz_.resize(f.jz().nx(),f.jz().ny(),f.jz().nz());referenceJz_.fill(0.);
    auto e=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z)->long double{return oldE_[a](x,y,z);};auto b=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z)->long double{return oldB_[a](x,y,z);};auto k=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z)->long double{return oldK_[a](x,y,z);};const auto v=stock(e,b,k);beginH_=v[0];beginP_=v[1];write("begin",v[0],v[1]);
  }
  void magneticRaw(const EBFieldGrid& f,const EBElectricHaloView& h){
    copy(magnetic(f),rawB_);const auto& g=geometry_.g;
    for(unsigned a=0;a<3;++a){
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
      for(std::size_t z=0;z<g.nz+(a==2);++z)for(std::size_t y=geometry_.cy;y<g.ny-geometry_.cy+(a==1);++y)for(std::size_t x=geometry_.cx;x<g.nx-geometry_.cx+(a==0);++x)ghostB_[a](x,y,z)=f.diagnosticMagneticIncrement(a,x,y,z,h,false)-oldK_[a](x,y,z);
    }
  }
  void beforeReferenceCurrent(const EBFieldGrid& f){referenceJz_=f.jz();}
  void afterReferenceCurrent(const EBFieldGrid& f){for(std::size_t i=0;i<referenceJz_.size();++i)referenceJz_.data()[i]=f.jz().data()[i]-referenceJz_.data()[i];}
  void sample(const char* stage,const EBFieldGrid& f){const std::string s=stage;
    if(s=="magnetic_cpml")copy(magnetic(f),cpmlB_);
    else if(s=="magnetic_mpi")copy(magnetic(f),plusB_);
    else if(s=="electric_curl_local")copy(electric(f),primaryE_);
    else if(s=="electric_cpml")copy(electric(f),cpmlE_);
    else if(s=="inlet_overwrite")copy(electric(f),inletE_);
    else if(s=="electric_mpi")finishFixed(f);
    else if(s=="remap"){
      calculateK(f,newK_);auto p=electric(f),q=magnetic(f);
      auto e=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z)->long double{return (*p[a])(x,y,z);};auto b=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z)->long double{return (*q[a])(x,y,z);};auto k=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z)->long double{return newK_[a](x,y,z);};const auto v=stock(e,b,k);std::vector<long double> data(28,0.);data[24]=v[0]-preRemapH_;data[25]=v[1]-preRemapP_;write("remap",v[0],v[1],data);std::swap(oldK_,newK_);
    }
  }
  void finishFixed(const EBFieldGrid& f){
    calculateK(f,newK_);YeeReadOnlyCollar oldHalo(pointers(oldE_),geometry_.g.nz,comm_),newHalo(electric(f),geometry_.g.nz,comm_),bHalo(magnetic(f),geometry_.g.nz,comm_);auto eo=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return oldHalo(a,x,y,z);};auto en=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return newHalo(a,x,y,z);};auto b=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return bHalo(a,x,y,z);};const auto& m=geometry_;const auto& g=m.g;const YeeComponent* current[3]={&f.jx(),&f.jy(),&f.jz()};
    auto em=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return .5L*(eo(a,x,y,z)+en(a,x,y,z));};
    auto bm=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return .5L*((long double)oldB_[a](x,y,z)+b(a,x,y,z));};
    auto ko=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z)->long double{return oldK_[a](x,y,z);};auto kn=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z)->long double{return newK_[a](x,y,z);};
    auto fo=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return ko(a,x,y,z)-m.faraday(eo,a,x,y,z);};auto fn=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return kn(a,x,y,z)-m.faraday(en,a,x,y,z);};
    auto jr=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z)->long double{return a==2?referenceJz_(x,y,z):0.;};auto jp=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return (long double)(*current[a])(x,y,z)-jr(a,x,y,z);};auto dj=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return -(long double)g.dt/SI::epsilon0*jp(a,x,y,z);};
    std::vector<long double> v(28,0.);v[0]=m.sum(m.energyPort(em,b));v[1]=.5L*(m.sum(m.electricPortStress(eo))+m.sum(m.electricPortStress(en)))+m.sum(m.magneticStress(b));v[2]=.5L/SI::mu0*(m.magneticInner(b,fo)+m.magneticInner(b,fn));v[3]=.5L*(m.momentum(eo,fo)+m.momentum(en,fn));v[4]=.5L*(m.gaussInteriorElectricImpulse(eo)+m.gaussInteriorElectricImpulse(en));v[5]=m.gaussMagneticImpulse(b);v[6]=(long double)g.dt*m.electricInner(em,jp);v[7]=m.momentum(dj,b);
    auto supplied=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z)->long double{return ghostB_[a](x,y,z);};auto mpml=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return (long double)cpmlB_[a](x,y,z)-rawB_[a](x,y,z);};auto msync=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return b(a,x,y,z)-(long double)cpmlB_[a](x,y,z);};
    magneticSource(v,8,bm,eo,supplied);magneticSource(v,10,bm,eo,mpml);magneticSource(v,12,bm,eo,msync);
    auto held=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z)->long double{if(a<2&&(m.offset+z==0||m.offset+z==m.globalNz))return -m.ampere(b,a,x,y,z)+(long double)g.dt/SI::epsilon0*(*current[a])(x,y,z);return 0.;};
    auto epml=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return (long double)cpmlE_[a](x,y,z)-primaryE_[a](x,y,z);};auto inlet=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return (long double)inletE_[a](x,y,z)-cpmlE_[a](x,y,z);};auto esync=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return en(a,x,y,z)-(long double)inletE_[a](x,y,z);};auto eref=[&](unsigned a,std::ptrdiff_t x,std::ptrdiff_t y,std::ptrdiff_t z){return -(long double)g.dt/SI::epsilon0*jr(a,x,y,z);};
    electricSource(v,14,em,b,held);electricSource(v,16,em,b,epml);electricSource(v,18,em,b,inlet);electricSource(v,20,em,b,esync);electricSource(v,22,em,b,eref);
    const auto end=stock(en,b,kn);preRemapH_=end[0];preRemapP_=end[1];long double hu=0.,pu=0.;for(unsigned i=8;i<24;i+=2){hu+=v[i];pu+=v[i+1];}v[26]=end[0]-beginH_+v[0]+v[6]-hu-v[2];v[27]=end[1]-beginP_+v[1]-v[7]-pu-v[3]+v[4]+v[5];write("fixed_end",end[0],end[1],v);
  }
  void close(bool complete){out_.flush();if(!out_)throw std::runtime_error("Whole Cowan audit flush failed");out_.close();std::ofstream(base_+".complete")<<(complete?1:0)<<' '<<step_<<'\n';}
};
}
#endif
