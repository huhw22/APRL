/* Independent known-boundary vacuum gate for a moving Cowan-z window.
 *
 * A superposition of forward/backward discrete plane waves has an explicitly
 * known exterior field. Only boundary layers, external electric halos and
 * newly entering planes are supplied analytically. Every interior value is
 * advanced by the unmodified EBFieldGrid kernels. No full-domain simulation
 * supplies any input in this executable. This does not solve the unknown
 * exterior response of a held, self-consistent particle bunch.
 */
#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <mpi.h>

#include "eb_field.h"
#include "compute_threads.h"

namespace
{
  using aprl::Double;
  using aprl::EBFieldGrid;
  using aprl::YeeComponent;
  const std::size_t NX=8,NY=7,NZ=48;
  const double DZ=1.2543262415359034e-6,WEAK=1.e-6;
  const double PI=std::acos(-1.);
  const double VACUUM_C=1./std::sqrt(aprl::SI::mu0*aprl::SI::epsilon0);
  const double ZOFFSET[6]={0.,0.,.5,.5,.5,0.};
  const double OFFSET[6][3]={{.5,0,0},{0,.5,0},{0,0,.5},
                            {0,.5,.5},{.5,0,.5},{.5,.5,0}};

  void require(bool ok,const std::string& text)
  { if(!ok)throw std::runtime_error(text); }

  std::ofstream freshFile(const std::string& path)
  {
    struct stat status;
    require(::stat(path.c_str(),&status)!=0 && errno==ENOENT,"Refusing to overwrite "+path);
    std::ofstream stream(path.c_str());require(static_cast<bool>(stream),"Cannot create "+path);
    stream<<std::setprecision(17);return stream;
  }

  aprl::EBGridGeometry geometry()
  { return aprl::EBGridGeometry(NX,NY,NZ,4*DZ,5*DZ,DZ,DZ/aprl::SI::c,aprl::EBMaxwellSolver::CowanZ); }

  YeeComponent& component(EBFieldGrid& f,unsigned a)
  {
    switch(a){case 0:return f.ex();case 1:return f.ey();case 2:return f.ez();
      case 3:return f.bx();case 4:return f.by();default:return f.bz();}
  }
  const YeeComponent& component(const EBFieldGrid& f,unsigned a)
  { return component(const_cast<EBFieldGrid&>(f),a); }

  double wave(double z,double time,double wavelength,int direction,double phase)
  {
    const double k=2*PI/wavelength;
    // Exact one-dimensional discrete vacuum dispersion for the actual SI
    // constants used by EBFieldGrid, not an assumed rounded mu0*eps0*c^2=1.
    const double omega=2*std::asin((VACUUM_C/aprl::SI::c)*std::sin(.5*k));
    return std::cos(k*z-direction*omega*time+phase);
  }

  double analytic(unsigned a,double z,double electricStep,bool withWeak)
  {
    const double time=a<3?electricStep:electricStep-.5;
    const bool ey=(a==1 || a==3);
    if(a==2)return .17+(withWeak?.2*WEAK:0.);
    if(a==5)return (.11+(withWeak?.1*WEAK:0.))/aprl::SI::c;
    const double shift=ey?.43:0.;
    const double scale=ey?.63:1.;
    const double plus=scale*wave(z,time,61.,+1,.23+shift)+
      (withWeak?WEAK*(ey?.71:1.)*wave(z,time,13.,+1,.37+shift):0.);
    const double minus=.21*scale*wave(z,time,89.,-1,.91+shift)+
      (withWeak?.43*WEAK*(ey?.71:1.)*wave(z,time,9.,-1,.17+shift):0.);
    if(a<3)return plus+minus;
    // Bx=-Ey/v for a +z wave, By=Ex/v; both signs reverse for a -z wave.
    return (a==3?-1.:1.)*(plus-minus)/VACUUM_C;
  }

  class AnalyticBoundary:public aprl::EBBoundaryOperator
  {
  public:
    AnalyticBoundary(std::size_t& offset,bool weak,std::uint64_t& bytes)
      :offset_(offset),weak_(weak),bytes_(bytes),step_(0){}
    void setStep(unsigned before){step_=before;}
    void afterMagneticUpdate(EBFieldGrid& fields,const aprl::EBElectricHaloView&) override
    {
      // The component's outermost x/y index plane only. This corrects the
      // default parity extension at the lateral boundary, not the interior.
      for(unsigned a=3;a<6;++a)fillSurface(fields,a,false);
    }
    void afterElectricUpdate(EBFieldGrid& fields) override
    {
      // Analytic lateral walls and the two z end layers. Interior samples
      // whose indices are strictly between all component faces are untouched.
      for(unsigned a=0;a<3;++a)fillSurface(fields,a,true);
    }
    std::size_t memoryBytes() const override{return 0;}
  private:
    void fillSurface(EBFieldGrid& fields,unsigned a,bool includeZ)
    {
      YeeComponent& f=component(fields,a);
      for(std::size_t k=0;k<f.nz();++k)
        for(std::size_t j=0;j<f.ny();++j)
          for(std::size_t i=0;i<f.nx();++i)
            if(i==0 || i+1==f.nx() || j==0 || j+1==f.ny() ||
               (includeZ && (k==0 || k+1==f.nz()))) {
              f(i,j,k)=analytic(a,offset_+k+ZOFFSET[a],step_+1.,weak_);
              bytes_+=sizeof(Double);
            }
    }
    std::size_t& offset_;bool weak_;std::uint64_t& bytes_;unsigned step_;
  };

  void initialize(EBFieldGrid& fields,std::size_t offset,bool weak,std::uint64_t& bytes)
  {
    for(unsigned a=0;a<6;++a){
      YeeComponent& f=component(fields,a);
      for(std::size_t k=0;k<f.nz();++k)
        for(std::size_t j=0;j<f.ny();++j)
          for(std::size_t i=0;i<f.nx();++i){
            f(i,j,k)=analytic(a,offset+k+ZOFFSET[a],0.,weak);bytes+=sizeof(Double);
          }
    }
  }

  struct Halos
  {
    std::vector<Double> lower[3],upper[3];
    aprl::EBElectricHaloView fill(const EBFieldGrid& fields,std::size_t offset,
                                 unsigned step,bool weak,std::uint64_t& bytes)
    {
      for(unsigned a=0;a<3;++a){
        const YeeComponent& f=component(fields,a);const std::size_t count=f.nx()*f.ny();
        lower[a].assign(count,analytic(a,offset-1+ZOFFSET[a],step,weak));
        upper[a].assign(count,analytic(a,offset+f.nz()+ZOFFSET[a],step,weak));
        bytes+=2*count*sizeof(Double);
      }
      aprl::EBElectricHaloView h;
      h.lowerEx=lower[0].data();h.lowerEy=lower[1].data();h.lowerEz=lower[2].data();
      h.upperEx=upper[0].data();h.upperEy=upper[1].data();h.upperEz=upper[2].data();return h;
    }
  };

  void shift(EBFieldGrid& fields,std::size_t before,std::size_t after,
             unsigned step,bool weak,bool supply,std::uint64_t& bytes,std::uint64_t& moved)
  {
    require(after<=before && before-after<=1,"Invalid window shift");
    if(after==before)return;
    for(unsigned a=0;a<6;++a){
      YeeComponent& f=component(fields,a);const std::size_t plane=f.nx()*f.ny();
      const std::size_t n=(f.nz()-1)*plane*sizeof(Double);
      std::memmove(f.data()+plane,f.data(),n);moved+=n;
      std::fill(f.data(),f.data()+plane,supply?analytic(a,after+ZOFFSET[a],step,weak):0.);
      if(supply)bytes+=plane*sizeof(Double);
    }
  }

  double fieldSample(const EBFieldGrid& fields,unsigned a,std::size_t offset,
                      std::size_t i,std::size_t j,double globalZ)
  {
    const YeeComponent& f=component(fields,a);
    const double q[3]={i+.5-OFFSET[a][0],j+.5-OFFSET[a][1],globalZ-OFFSET[a][2]};
    std::ptrdiff_t k[3];double u[3];
    for(unsigned axis=0;axis<3;++axis){k[axis]=std::floor(q[axis]);u[axis]=q[axis]-k[axis];}
    double out=0.;
    for(unsigned z=0;z<2;++z)for(unsigned y=0;y<2;++y)for(unsigned x=0;x<2;++x){
      const double weight=(x?u[0]:1-u[0])*(y?u[1]:1-u[1])*(z?u[2]:1-u[2]);
      if(weight==0)continue;
      const std::ptrdiff_t ix=k[0]+x,iy=k[1]+y,iz=k[2]+z-offset;
      require(ix>=0 && iy>=0 && iz>=0 && ix<static_cast<std::ptrdiff_t>(f.nx()) &&
              iy<static_cast<std::ptrdiff_t>(f.ny()) && iz<static_cast<std::ptrdiff_t>(f.nz()),"Sample left window");
      out+=weight*f(ix,iy,iz);
    }
    return out;
  }

  struct Point{double value[6];};
  Point sample(const EBFieldGrid* fields,std::size_t offset,std::size_t i,std::size_t j,
                double globalZ,unsigned step,bool weak,const aprl::BoostFrameTransform& frame)
  {
    aprl::FieldVector<Double> e(0.),b(0.),el(0.),bl(0.);
    for(unsigned a=0;a<6;++a){
      double value=0.;
      if(fields)value=fieldSample(*fields,a,offset,i,j,globalZ);
      else {
        const double q=globalZ-ZOFFSET[a],floor=std::floor(q),u=q-floor;
        value=(1-u)*analytic(a,floor+ZOFFSET[a],step,weak)+
          u*analytic(a,floor+1+ZOFFSET[a],step,weak);
      }
      if(a<3)e[a]=value;else b[a-3]=value;
    }
    frame.fieldsBoxToLab(e,b,el,bl);Point out;
    for(unsigned a=0;a<3;++a){out.value[a]=el[a];out.value[a+3]=aprl::SI::c*bl[a];}
    return out;
  }

  struct Stats
  {
    unsigned steps=0;double beta=0.,weakMax=0.,weakError=0.,rawError=0.,negativeError=0.,stateError=0.,seconds=0.;
    std::uint64_t initialBytes=0,boundaryBytes=0,movedBytes=0,valuesChecked=0;
    std::size_t memory=0;std::string csv;
  };

  Stats run(const std::string& root,unsigned steps,double beta,unsigned id)
  {
    Stats out;out.steps=steps;out.beta=beta;
    const auto begin=std::chrono::steady_clock::now();
    EBFieldGrid signal(geometry()),background(geometry()),negative(geometry());
    const std::size_t initialOffset=steps+16;std::size_t offset=initialOffset;
    AnalyticBoundary* sb=new AnalyticBoundary(offset,true,out.boundaryBytes);
    AnalyticBoundary* bb=new AnalyticBoundary(offset,false,out.boundaryBytes);
    signal.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(sb));
    background.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(bb));
    initialize(signal,offset,true,out.initialBytes);initialize(background,offset,false,out.initialBytes);
    std::uint64_t ignored=0,negativeMoved=0;initialize(negative,offset,true,ignored);
    out.memory=signal.memoryFootprint().totalBytes();
    Halos sh,bh;
    aprl::BoostFrameTransform frame;frame.setFromGamma(1/std::sqrt(1-beta*beta),aprl::SI::c,0.);
    std::ostringstream name;name<<"analytic-"<<id<<"-steps-"<<steps<<".csv";out.csv=name.str();
    std::ofstream csv=freshFile(root+"/"+out.csv);
    csv<<"step,i,j,time_box_s,z_box_m";
    const char* groups[]={"reference","window","zero_boundary","reference_background","window_background"};
    const char* axes[]={"Ex","Ey","Ez","cBx","cBy","cBz"};
    for(const char* group:groups)for(const char* axis:axes)csv<<','<<group<<'_'<<axis;
    csv<<'\n';
    for(unsigned step=0;step<=steps;++step){
      if(step){
        sb->setStep(step-1);bb->setStep(step-1);
        signal.advanceMagnetic(sh.fill(signal,offset,step-1,true,out.boundaryBytes));
        background.advanceMagnetic(bh.fill(background,offset,step-1,false,out.boundaryBytes));
        negative.advanceMagnetic();
        signal.advanceElectric();background.advanceElectric();negative.advanceElectric();
        const std::size_t next=initialOffset-static_cast<std::size_t>(std::floor(beta*step));
        shift(signal,offset,next,step,true,true,out.boundaryBytes,out.movedBytes);
        shift(background,offset,next,step,false,true,out.boundaryBytes,out.movedBytes);
        shift(negative,offset,next,step,true,false,ignored,negativeMoved);
        offset=next;
      }
      for(unsigned a=0;a<6;++a){
        const YeeComponent& f=component(signal,a);
        for(std::size_t k=2;k+2<f.nz();++k)
          for(std::size_t j=1;j+1<f.ny();++j)
            for(std::size_t i=1;i+1<f.nx();++i){
              const double value=f(i,j,k),expected=analytic(a,offset+k+ZOFFSET[a],step,true);
              require(std::isfinite(value),"Nonfinite analytic window state");
              out.stateError=std::max(out.stateError,std::abs(value-expected)*(a<3?1.:aprl::SI::c));
              ++out.valuesChecked;
            }
      }
      const double q=initialOffset+.5*NZ-beta*step;
      // Exclude the one-cell lateral boundary layer from the measurement.
      // This is the entire 6x5 interior transverse plane, not a point sample.
      for(std::size_t j=1;j+1<NY;++j)for(std::size_t i=1;i+1<NX;++i){
        const Point r=sample(nullptr,0,i,j,q,step,true,frame),w=sample(&signal,offset,i,j,q,step,true,frame),
          z=sample(&negative,offset,i,j,q,step,true,frame),b=sample(nullptr,0,i,j,q,step,false,frame),
          wb=sample(&background,offset,i,j,q,step,false,frame);
        const Point values[]={r,w,z,b,wb};
        csv<<step<<','<<i<<','<<j<<','<<step*geometry().dt<<','<<q*DZ;
        for(const Point& p:values)for(double value:p.value){require(std::isfinite(value),"Nonfinite analytic plane");csv<<','<<value;}
        csv<<'\n';
        for(unsigned a=0;a<6;++a){
          out.rawError=std::max(out.rawError,std::abs(w.value[a]-r.value[a]));
          out.negativeError=std::max(out.negativeError,std::abs(z.value[a]-r.value[a]));
          out.weakMax=std::max(out.weakMax,std::abs(r.value[a]-b.value[a]));
          out.weakError=std::max(out.weakError,std::abs((w.value[a]-wb.value[a])-(r.value[a]-b.value[a])));
        }
      }
    }
    csv.close();require(static_cast<bool>(csv),"Cannot finish analytic CSV");
    out.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    require(out.weakMax>0 && out.weakError/out.weakMax<1.e-5,
            "Analytic-boundary weak-signal error exceeds the declared 1e-5 tolerance");
    require(out.rawError/out.weakMax<1.e-5,
            "Analytic-boundary total-field error exceeds 1e-5 of the weak signal");
    require(out.negativeError/out.weakMax>100.,"Analytic missing-boundary control did not fail clearly");
    return out;
  }

  void writeStats(std::ostream& json,const Stats& s)
  {
    json<<"{\"steps\":"<<s.steps<<",\"beta\":"<<s.beta<<",\"csv\":\""<<s.csv<<"\""
        <<",\"window_nz\":"<<NZ<<",\"window_bytes_per_grid\":"<<s.memory
        <<",\"state_values_checked\":"<<s.valuesChecked<<",\"max_state_error_E_units\":"<<s.stateError
        <<",\"max_lab_plane_error_E_units\":"<<s.rawError<<",\"weak_lab_plane_max_E_units\":"<<s.weakMax
        <<",\"weak_lab_plane_error_E_units\":"<<s.weakError<<",\"weak_lab_plane_relative_error\":"<<s.weakError/s.weakMax
        <<",\"raw_error_over_weak\":"<<s.rawError/s.weakMax
        <<",\"zero_boundary_error_over_weak\":"<<s.negativeError/s.weakMax
        <<",\"cell_updates_per_grid\":"<<static_cast<std::uint64_t>(NX*NY*NZ)*s.steps
        <<",\"positive_grid_count\":2,\"negative_grid_count\":1,\"full_reference_grid_count\":0"
        <<",\"analytic_initial_bytes_two_windows\":"<<s.initialBytes
        <<",\"analytic_supplied_boundary_bytes_two_windows\":"<<s.boundaryBytes
        <<",\"shift_bytes_two_windows\":"<<s.movedBytes
        <<",\"total_seconds_including_boundaries_and_csv\":"<<s.seconds<<'}';
  }
}

int main(int argc,char** argv)
{
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int size=0;MPI_Comm_size(MPI_COMM_WORLD,&size);int result=0;
  try{
    require(size==1,"Analytic window gate uses exactly one MPI rank");
    require(argc==2,"Usage: aprl_window_analytic_probe OUTPUT_DIRECTORY");
    aprl::configureComputeThreads(provided);
    require(::mkdir(argv[1],0775)==0 || errno==EEXIST,"Cannot create analytic output directory");
    const double beta[]={.9,aprl::BoostFrameTransform::betaFromGamma(929.1305492858543)};
    std::vector<Stats> results;
    for(unsigned b=0;b<2;++b)for(unsigned n:{48u,96u,192u}){
      results.push_back(run(argv[1],n,beta[b],b));
      std::cout<<"PASS analytic beta="<<std::setprecision(17)<<beta[b]<<" steps="<<n
               <<" weak_relative_error="<<results.back().weakError/results.back().weakMax<<std::endl;
    }
    std::ofstream json=freshFile(std::string(argv[1])+"/analytic-results.json");
    json<<"{\"format\":\"aprl-window-analytic-probe-v1\",\"pass\":true,\"nx\":"<<NX<<",\"ny\":"<<NY
        <<",\"diagnostic_interior_pixels\":"<<(NX-2)*(NY-2)<<",\"dz_m\":"<<DZ<<",\"dt_s\":"<<geometry().dt
        <<",\"weak_initial_scale\":"<<WEAK<<",\"declared_error_relative_to_weak_limit\":1e-5,\"cases\":[";
    for(std::size_t i=0;i<results.size();++i){if(i)json<<',';writeStats(json,results[i]);}
    json<<"],\"limitations\":[\"Independent analytic boundary supply is demonstrated only for known vacuum discrete plane waves.\","
        <<"\"One outer component-index layer in x/y is prescribed after B; E has that lateral layer plus its two z end layers prescribed. Interior is advanced by Cowan kernels.\","
        <<"\"No external full-grid history is computed; incoming field is explicitly known, unlike a general held-particle near field.\","
        <<"\"No particles, current, CPML, MPI decomposition or production checkpoint support.\","
        <<"\"Native leapfrog E/B time staggering is retained. Timing includes diagnostics and is not a production speedup claim.\"]}\n";
    json.close();require(static_cast<bool>(json),"Cannot finish analytic results");
  }catch(const std::exception& error){std::cerr<<"FAIL analytic window probe: "<<error.what()<<std::endl;result=1;}
  MPI_Finalize();return result;
}
