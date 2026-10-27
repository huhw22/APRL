/* Feasibility fixture, not a production moving-window implementation.
 *
 * The oracle cases run the unmodified Cowan-z kernels on a full reference
 * domain and a genuinely smaller, fixed-width window. Only incoming planes,
 * electric end faces and the Cowan electric halos come from the reference.
 * Interior fields are advanced, never overwritten by reference fields.
 * Reference generation and supplied boundary bytes are explicitly reported.
 * Prescribed J is a manufactured Maxwell source, not a self-consistent PIC
 * current. There is no CPML, particle injection, MPI decomposition or new
 * approximation to an unknown exterior response in this fixture.
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
  const std::size_t NX = 6, NY = 5, WIDTH = 48;
  const double DZ = 1.2543262415359034e-6;
  const double WEAK = 1.e-6;
  const double PI = std::acos(-1.0);
  const double SHIFT[6][3] = {{.5,0,0},{0,.5,0},{0,0,.5},
                            {0,.5,.5},{.5,0,.5},{.5,.5,0}};

  void require(bool ok, const std::string& message)
  { if (!ok) throw std::runtime_error(message); }

  void makeDirectory(const std::string& directory)
  {
    require(::mkdir(directory.c_str(),0775)==0 || errno==EEXIST,
            "Cannot create output directory: " + directory);
  }

  std::ofstream freshFile(const std::string& path)
  {
    struct stat status;
    require(::stat(path.c_str(),&status)!=0 && errno==ENOENT,
            "Refusing to overwrite " + path);
    std::ofstream result(path.c_str());
    require(static_cast<bool>(result),"Cannot create " + path);
    result << std::setprecision(17);
    return result;
  }

  aprl::EBGridGeometry geometry(std::size_t nz)
  {
    return aprl::EBGridGeometry(NX,NY,nz,4*DZ,5*DZ,DZ,DZ/aprl::SI::c,
                               aprl::EBMaxwellSolver::CowanZ);
  }

  YeeComponent& component(EBFieldGrid& fields,unsigned a)
  {
    switch(a) {
      case 0:return fields.ex(); case 1:return fields.ey();
      case 2:return fields.ez(); case 3:return fields.bx();
      case 4:return fields.by(); default:return fields.bz();
    }
  }

  const YeeComponent& component(const EBFieldGrid& fields,unsigned a)
  { return component(const_cast<EBFieldGrid&>(fields),a); }

  void copyPlane(YeeComponent& target,std::size_t to,
                 const YeeComponent& source,std::size_t from,
                 std::uint64_t& suppliedBytes)
  {
    require(target.nx()==source.nx() && target.ny()==source.ny() &&
            to<target.nz() && from<source.nz(),"Invalid plane copy");
    const std::size_t bytes=target.nx()*target.ny()*sizeof(Double);
    std::memcpy(target.data()+to*target.nx()*target.ny(),
                source.data()+from*source.nx()*source.ny(),bytes);
    suppliedBytes+=bytes;
  }

  void copyInitial(EBFieldGrid& window,const EBFieldGrid& reference,
                   std::size_t offset,std::uint64_t& suppliedBytes)
  {
    for(unsigned a=0;a<6;++a) {
      YeeComponent& out=component(window,a);
      for(std::size_t k=0;k<out.nz();++k)
        copyPlane(out,k,component(reference,a),offset+k,suppliedBytes);
    }
  }

  void initialize(EBFieldGrid& fields,bool weak)
  {
    for(unsigned a=0;a<6;++a) {
      YeeComponent& f=component(fields,a);
      for(std::size_t k=0;k<f.nz();++k)
        for(std::size_t j=0;j<f.ny();++j)
          for(std::size_t i=0;i<f.nx();++i) {
            const double x=i+SHIFT[a][0],y=j+SHIFT[a][1],z=k+SHIFT[a][2];
            const double background=(1.+.13*a)*(.7+.3*std::cos(.031*z+.1*a))
              *(.8+.2*std::cos(.19*x-.17*y+.11*a));
            const double signal=WEAK*(1.+.1*a)*std::cos(2*PI*z/11.+.17*x+.23*y+.3*a);
            f(i,j,k)=(background+(weak?signal:0.))/(a<3?1.:aprl::SI::c);
          }
    }
    // Apply the same transverse and external end-face PEC convention at t=0.
    fields.boundary().afterElectricUpdate(fields);
  }

  void prescribedCurrent(EBFieldGrid& fields,std::size_t offset,
                         unsigned step,bool weak)
  {
    YeeComponent* f[3]={&fields.jx(),&fields.jy(),&fields.jz()};
    const double scale=aprl::SI::epsilon0/fields.geometry().dt;
    for(unsigned a=0;a<3;++a)
      for(std::size_t k=0;k<f[a]->nz();++k)
        for(std::size_t j=0;j<f[a]->ny();++j)
          for(std::size_t i=0;i<f[a]->nx();++i) {
            const double z=offset+k+SHIFT[a][2],x=i+SHIFT[a][0],y=j+SHIFT[a][1];
            const double envelope=std::sin(PI*x/NX)*std::sin(PI*y/NY);
            const double background=.002*envelope*std::sin(.023*z-.031*step+.4*a);
            const double signal=WEAK*.03*envelope*std::cos(2*PI*z/13.-.19*step+.2*a);
            (*f[a])(i,j,k)=scale*(background+(weak?signal:0.));
          }
  }

  aprl::EBElectricHaloView referenceHalos(const EBFieldGrid& reference,
      const EBFieldGrid& window,std::size_t offset,std::uint64_t& suppliedBytes)
  {
    aprl::EBElectricHaloView halo;
    const Double* lower[3]; const Double* upper[3];
    for(unsigned a=0;a<3;++a) {
      const YeeComponent& r=component(reference,a);
      const YeeComponent& w=component(window,a);
      require(offset>0 && offset+w.nz()<r.nz(),"Window halo outside reference");
      lower[a]=r.data()+(offset-1)*r.nx()*r.ny();
      upper[a]=r.data()+(offset+w.nz())*r.nx()*r.ny();
      suppliedBytes+=2*r.nx()*r.ny()*sizeof(Double);
    }
    halo.lowerEx=lower[0];halo.lowerEy=lower[1];halo.lowerEz=lower[2];
    halo.upperEx=upper[0];halo.upperEy=upper[1];halo.upperEz=upper[2];
    return halo;
  }

  void electricFaces(EBFieldGrid& window,const EBFieldGrid& reference,
                     std::size_t offset,std::uint64_t& suppliedBytes)
  {
    // advanceElectric does not update Ex/Ey at the two local z faces.
    for(unsigned a=0;a<2;++a) {
      YeeComponent& out=component(window,a);
      copyPlane(out,0,component(reference,a),offset,suppliedBytes);
      copyPlane(out,out.nz()-1,component(reference,a),offset+out.nz()-1,suppliedBytes);
    }
  }

  void shiftWindow(EBFieldGrid& window,const EBFieldGrid& reference,
      std::size_t oldOffset,std::size_t newOffset,bool replay,
      std::uint64_t& suppliedBytes,std::uint64_t& movedBytes)
  {
    require(newOffset<=oldOffset && oldOffset-newOffset<=1,
            "Fixture only supports zero/one-cell shifts");
    if(oldOffset==newOffset) return;
    for(unsigned a=0;a<6;++a) {
      YeeComponent& f=component(window,a);
      const std::size_t plane=f.nx()*f.ny();
      const std::size_t bytes=(f.nz()-1)*plane*sizeof(Double);
      std::memmove(f.data()+plane,f.data(),bytes);
      movedBytes+=bytes;
      if(replay) copyPlane(f,0,component(reference,a),newOffset,suppliedBytes);
      else std::fill(f.data(),f.data()+plane,0.);
    }
  }

  // Global interpolation weights are shared between full and local lattices.
  // This prevents coordinate-subtraction roundoff from obscuring field-state
  // equality. It is diagnostic-only; the actual kernels remain unmodified.
  double sampleComponent(const EBFieldGrid& fields,unsigned a,std::size_t offset,
                         std::size_t i,std::size_t j,double globalZ)
  {
    const YeeComponent& f=component(fields,a);
    const double q[3]={i+.5-SHIFT[a][0],j+.5-SHIFT[a][1],globalZ-SHIFT[a][2]};
    const std::ptrdiff_t base[3]={static_cast<std::ptrdiff_t>(std::floor(q[0])),
      static_cast<std::ptrdiff_t>(std::floor(q[1])),static_cast<std::ptrdiff_t>(std::floor(q[2]))};
    const double fraction[3]={q[0]-base[0],q[1]-base[1],q[2]-base[2]};
    double result=0.;
    for(unsigned dk=0;dk<2;++dk)
      for(unsigned dj=0;dj<2;++dj)
        for(unsigned di=0;di<2;++di) {
          const double weight=(di?fraction[0]:1-fraction[0])*
            (dj?fraction[1]:1-fraction[1])*(dk?fraction[2]:1-fraction[2]);
          if(weight==0) continue;
          const std::ptrdiff_t x=base[0]+di,y=base[1]+dj,z=base[2]+dk-offset;
          require(x>=0 && y>=0 && z>=0 && x<static_cast<std::ptrdiff_t>(f.nx()) &&
                  y<static_cast<std::ptrdiff_t>(f.ny()) && z<static_cast<std::ptrdiff_t>(f.nz()),
                  "Diagnostic sample left the lattice");
          result+=weight*f(x,y,z);
        }
    return result;
  }

  struct Point { double value[6]; };

  Point sampleLab(const EBFieldGrid& fields,std::size_t offset,std::size_t i,
                   std::size_t j,double globalZ,const aprl::BoostFrameTransform& frame)
  {
    aprl::FieldVector<Double> e(0.),b(0.),labE(0.),labB(0.);
    for(unsigned a=0;a<3;++a) {
      e[a]=sampleComponent(fields,a,offset,i,j,globalZ);
      b[a]=sampleComponent(fields,a+3,offset,i,j,globalZ);
    }
    frame.fieldsBoxToLab(e,b,labE,labB);
    Point out;
    for(unsigned a=0;a<3;++a) {out.value[a]=labE[a];out.value[a+3]=aprl::SI::c*labB[a];}
    return out;
  }

  struct Stats
  {
    unsigned steps=0;double beta=0.;std::size_t fullNz=0;
    std::uint64_t initialBytes=0,boundaryBytes=0,movedBytes=0;
    std::uint64_t fieldValuesChecked=0,fieldUnequal=0;
    double maxStateDifference=0.,maxOracleError=0.,maxZeroError=0.;
    double maxWeakReference=0.,maxWeakError=0.;
    double fullKernelSeconds=0.,windowKernelSeconds=0.,negativeKernelSeconds=0.,elapsedSeconds=0.;
    std::size_t fullMemory=0,windowMemory=0;
    std::string csv;
  };

  void checkState(const EBFieldGrid& window,const EBFieldGrid& reference,
                   std::size_t offset,Stats& result)
  {
    for(unsigned a=0;a<6;++a) {
      const YeeComponent& f=component(window,a),&r=component(reference,a);
      for(std::size_t k=2;k+2<f.nz();++k)
        for(std::size_t j=0;j<f.ny();++j)
          for(std::size_t i=0;i<f.nx();++i) {
            const double value=f(i,j,k),expected=r(i,j,k+offset);
            require(std::isfinite(value) && std::isfinite(expected),"Nonfinite oracle state");
            ++result.fieldValuesChecked;
            if(std::memcmp(&value,&expected,sizeof(double))!=0) ++result.fieldUnequal;
            result.maxStateDifference=std::max(result.maxStateDifference,
              std::abs(value-expected)*(a<3?1.:aprl::SI::c));
          }
    }
  }

  double seconds(std::chrono::steady_clock::time_point start)
  { return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(); }

  void timedMagnetic(EBFieldGrid& fields,const aprl::EBElectricHaloView& halo,double& elapsed)
  { const auto start=std::chrono::steady_clock::now();fields.advanceMagnetic(halo);elapsed+=seconds(start); }

  void timedElectric(EBFieldGrid& fields,double& elapsed)
  { const auto start=std::chrono::steady_clock::now();fields.advanceElectric();elapsed+=seconds(start); }

  Stats runOracle(const std::string& directory,unsigned steps,double beta,unsigned id)
  {
    Stats result;result.steps=steps;result.beta=beta;result.fullNz=steps+WIDTH+32;
    const auto begin=std::chrono::steady_clock::now();
    EBFieldGrid reference(geometry(result.fullNz)),background(geometry(result.fullNz));
    EBFieldGrid window(geometry(WIDTH)),windowBackground(geometry(WIDTH)),negative(geometry(WIDTH));
    window.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(new aprl::PerfectElectricConductorBoundary(false,false)));
    windowBackground.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(new aprl::PerfectElectricConductorBoundary(false,false)));
    initialize(reference,true);initialize(background,false);
    const std::size_t initialOffset=steps+16;
    std::size_t offset=initialOffset;
    copyInitial(window,reference,offset,result.initialBytes);
    copyInitial(windowBackground,background,offset,result.initialBytes);
    std::uint64_t negativeInitial=0,negativeMoved=0,ignoredBytes=0;
    copyInitial(negative,reference,offset,negativeInitial);
    result.fullMemory=reference.memoryFootprint().totalBytes();
    result.windowMemory=window.memoryFootprint().totalBytes();
    aprl::BoostFrameTransform frame;
    frame.setFromGamma(1/std::sqrt(1-beta*beta),aprl::SI::c,0.);
    std::ostringstream name;name<<"oracle-"<<id<<"-steps-"<<steps<<".csv";
    result.csv=name.str();std::ofstream csv=freshFile(directory+"/"+result.csv);
    csv<<"step,i,j,time_box_s,z_box_m";
    const char* groups[]={"reference","oracle","zero_boundary","reference_background","oracle_background"};
    const char* axes[]={"Ex","Ey","Ez","cBx","cBy","cBz"};
    for(const char* group:groups) for(const char* axis:axes) csv<<','<<group<<'_'<<axis;
    csv<<'\n';
    for(unsigned step=0;step<=steps;++step) {
      if(step) {
        // E halos are from the reference at n; updating reference B leaves E intact.
        const aprl::EBElectricHaloView h=referenceHalos(reference,window,offset,result.boundaryBytes);
        const aprl::EBElectricHaloView hb=referenceHalos(background,windowBackground,offset,result.boundaryBytes);
        timedMagnetic(reference,aprl::EBElectricHaloView(),result.fullKernelSeconds);
        timedMagnetic(background,aprl::EBElectricHaloView(),result.fullKernelSeconds);
        timedMagnetic(window,h,result.windowKernelSeconds);
        timedMagnetic(windowBackground,hb,result.windowKernelSeconds);
        timedMagnetic(negative,aprl::EBElectricHaloView(),result.negativeKernelSeconds);
        prescribedCurrent(reference,0,step,true);prescribedCurrent(background,0,step,false);
        prescribedCurrent(window,offset,step,true);prescribedCurrent(windowBackground,offset,step,false);
        prescribedCurrent(negative,offset,step,true);
        timedElectric(reference,result.fullKernelSeconds);timedElectric(background,result.fullKernelSeconds);
        timedElectric(window,result.windowKernelSeconds);timedElectric(windowBackground,result.windowKernelSeconds);
        timedElectric(negative,result.negativeKernelSeconds);
        electricFaces(window,reference,offset,result.boundaryBytes);
        electricFaces(windowBackground,background,offset,result.boundaryBytes);
        const std::size_t next=initialOffset-static_cast<std::size_t>(std::floor(beta*step));
        shiftWindow(window,reference,offset,next,true,result.boundaryBytes,result.movedBytes);
        shiftWindow(windowBackground,background,offset,next,true,result.boundaryBytes,result.movedBytes);
        shiftWindow(negative,reference,offset,next,false,ignoredBytes,negativeMoved);
        offset=next;
      }
      checkState(window,reference,offset,result);checkState(windowBackground,background,offset,result);
      const double q=initialOffset+.5*WIDTH-beta*step;
      for(std::size_t j=0;j<NY;++j) for(std::size_t i=0;i<NX;++i) {
        const Point r=sampleLab(reference,0,i,j,q,frame),w=sampleLab(window,offset,i,j,q,frame),
          z=sampleLab(negative,offset,i,j,q,frame),b=sampleLab(background,0,i,j,q,frame),
          wb=sampleLab(windowBackground,offset,i,j,q,frame);
        const Point values[]={r,w,z,b,wb};
        csv<<step<<','<<i<<','<<j<<','<<step*geometry(WIDTH).dt<<','<<q*DZ;
        for(const Point& p:values) for(double v:p.value) {require(std::isfinite(v),"Nonfinite plane value");csv<<','<<v;}
        csv<<'\n';
        for(unsigned a=0;a<6;++a) {
          result.maxOracleError=std::max(result.maxOracleError,std::abs(w.value[a]-r.value[a]));
          result.maxZeroError=std::max(result.maxZeroError,std::abs(z.value[a]-r.value[a]));
          result.maxWeakReference=std::max(result.maxWeakReference,std::abs(r.value[a]-b.value[a]));
          result.maxWeakError=std::max(result.maxWeakError,
            std::abs((w.value[a]-wb.value[a])-(r.value[a]-b.value[a])));
        }
      }
    }
    csv.close();require(static_cast<bool>(csv),"Failed to write plane CSV");
    result.elapsedSeconds=seconds(begin);
    require(result.fieldUnequal==0 && result.maxOracleError==0 && result.maxWeakError==0,
            "Oracle-boundary window did not reproduce the reference bitwise");
    require(result.maxWeakReference>0 && result.maxZeroError>100*result.maxWeakReference,
            "Missing-boundary negative control was not distinguishable from the weak signal");
    return result;
  }

  void jsonStats(std::ostream& out,const Stats& s)
  {
    out<<"{\"steps\":"<<s.steps<<",\"beta\":"<<s.beta<<",\"full_nz\":"<<s.fullNz
       <<",\"window_nz\":"<<WIDTH<<",\"plane_pixels\":"<<NX*NY<<",\"csv\":\""<<s.csv<<"\""
       <<",\"field_values_checked\":"<<s.fieldValuesChecked<<",\"field_bitwise_mismatches\":"<<s.fieldUnequal
       <<",\"max_state_error_E_units\":"<<s.maxStateDifference
       <<",\"max_lab_plane_error_E_units\":"<<s.maxOracleError
       <<",\"weak_lab_plane_max_E_units\":"<<s.maxWeakReference
       <<",\"weak_lab_plane_error_E_units\":"<<s.maxWeakError
       <<",\"zero_boundary_error_E_units\":"<<s.maxZeroError
       <<",\"zero_boundary_error_over_weak\":"<<s.maxZeroError/s.maxWeakReference
       <<",\"reference_cell_updates_per_grid\":"<<static_cast<std::uint64_t>(NX*NY*s.fullNz)*s.steps
       <<",\"window_cell_updates_per_grid\":"<<static_cast<std::uint64_t>(NX*NY*WIDTH)*s.steps
       <<",\"reference_grid_count\":2,\"oracle_window_grid_count\":2,\"negative_window_grid_count\":1"
       <<",\"reference_bytes_per_grid\":"<<s.fullMemory<<",\"window_bytes_per_grid\":"<<s.windowMemory
       <<",\"oracle_initial_bytes_two_windows\":"<<s.initialBytes
       <<",\"oracle_supplied_boundary_bytes_two_windows\":"<<s.boundaryBytes
       <<",\"window_shift_bytes_two_windows\":"<<s.movedBytes
       <<",\"reference_kernel_seconds_two_grids\":"<<s.fullKernelSeconds
       <<",\"window_kernel_seconds_two_grids\":"<<s.windowKernelSeconds
       <<",\"negative_kernel_seconds\":"<<s.negativeKernelSeconds
       <<",\"total_fixture_seconds_including_reference_and_csv\":"<<s.elapsedSeconds<<'}';
  }
}

int main(int argc,char** argv)
{
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int size=0;MPI_Comm_size(MPI_COMM_WORLD,&size);
  int result=0;
  try {
    require(size==1,"This bounded field-kernel fixture uses exactly one MPI rank");
    require(argc==2,"Usage: aprl_window_streaming_probe OUTPUT_DIRECTORY");
    aprl::configureComputeThreads(provided);
    makeDirectory(argv[1]);
    std::vector<Stats> outcomes;
    const double beta[]={.9,aprl::BoostFrameTransform::betaFromGamma(929.1305492858543)};
    const unsigned steps[]={48,96,192};
    for(unsigned b=0;b<2;++b) for(unsigned n:steps) {
      outcomes.push_back(runOracle(argv[1],n,beta[b],b));
      std::cout<<"PASS oracle beta="<<std::setprecision(17)<<beta[b]<<" steps="<<n
               <<" bitwise_mismatches="<<outcomes.back().fieldUnequal
               <<" negative_over_weak="<<outcomes.back().maxZeroError/outcomes.back().maxWeakReference<<std::endl;
    }
    std::ofstream json=freshFile(std::string(argv[1])+"/oracle-results.json");
    json<<"{\"format\":\"aprl-window-oracle-probe-v1\",\"pass\":true,\"nx\":"<<NX<<",\"ny\":"<<NY
        <<",\"dz_m\":"<<DZ<<",\"dt_s\":"<<geometry(WIDTH).dt<<",\"weak_initial_scale\":"<<WEAK
        <<",\"solver\":\"unmodified Cowan-z EBFieldGrid\",\"cases\":[";
    for(std::size_t i=0;i<outcomes.size();++i) {if(i)json<<',';jsonStats(json,outcomes[i]);}
    json<<"],\"limitations\":[\"Oracle E faces/halos and entering layers require the full reference simulation; not an independent exterior-field solver.\","
        <<"\"Manufactured prescribed current, no self-consistent particles, held activation or CPML.\","
        <<"\"Single MPI rank; no distributed ownership or checkpoint implementation.\","
        <<"\"Native leapfrog B half-time is preserved; E/B time alignment is not changed.\","
        <<"\"Reported cell-update counts show local work scaling only; total fixture includes reference generation and is not a production speedup measurement.\"]}\n";
    json.close();require(static_cast<bool>(json),"Failed to write result JSON");
  } catch(const std::exception& error) {std::cerr<<"FAIL window streaming probe: "<<error.what()<<std::endl;result=1;}
  MPI_Finalize();return result;
}
