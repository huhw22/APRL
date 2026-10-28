// Native particle interpolation oracle. No PIC/production acceptance implied.
#include "particle_field_halo.h"
#include "eb_particles.h"
#include "compute_threads.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using namespace aprl;
const std::size_t NX=10,NY=12,NZ=25;
const Double DZ=1.2543262415359034e-6;
unsigned long long samples=0,negativeSamples=0,unchangedSamples=0,physicalSamples=0,borisSamples=0;
Double maximumE=0.,maximumB=0.,maximumU=0.,maximumPosition=0.,negativeE=0.,negativeB=0.;
void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
bool equal(const FieldVector<Double>& a,const FieldVector<Double>& b) {
  for(unsigned k=0;k<3;++k)if(std::memcmp(&a[k],&b[k],sizeof(Double)))return false;
  return true;
}
Double error(const FieldVector<Double>& a,const FieldVector<Double>& b) {
  long double n=0,d=0;
  for(unsigned k=0;k<3;++k) {const long double v=a[k],delta=static_cast<long double>(b[k])-v;n+=v*v;d+=delta*delta;}
  return static_cast<Double>(std::sqrt(d/std::max(n,1e-90L)));
}
EBGridGeometry geometry(std::size_t nz) {
  return EBGridGeometry(NX,NY,nz,5e-6,5e-6,DZ,DZ/SI::c,EBMaxwellSolver::CowanZ);
}
YeeComponent& component(EBFieldGrid& f,unsigned a) {
  switch(a){case 0:return f.ex();case 1:return f.ey();case 2:return f.ez();
    case 3:return f.bx();case 4:return f.by();case 5:return f.bz();}
  throw std::out_of_range("Check component");
}
void fill(EBFieldGrid& f,std::size_t offset,unsigned pattern,unsigned refresh) {
  for(unsigned a=0;a<6;++a) {
    auto& c=component(f,a);
    for(std::size_t k=0;k<c.nz();++k)for(std::size_t j=0;j<c.ny();++j)
      for(std::size_t i=0;i<c.nx();++i) {
        const Double z=static_cast<Double>(offset+k);
        Double v=1.+.3*a+.013*i+.021*j+.071*z+.19*refresh;
        if(pattern) v+=.21*std::sin(.3*i+.7*j+.4*z+.1*a+.23*refresh)+.001*i*j*z;
        c(i,j,k)=v*(a<3?1e6:1e-3);
      }
  }
}
struct Point {
  FieldVector<Double> position,legacyE,legacyB;
  bool haloNeeded,physicalEnd;
};
void check(int rank,int ranks) {
  const std::size_t localNz=NZ/ranks+(static_cast<std::size_t>(rank)<NZ%ranks?1:0);
  const std::size_t offset=NZ/ranks*rank+std::min(static_cast<std::size_t>(rank),NZ%ranks);
  EBFieldGrid local(geometry(localNz)),whole(geometry(NZ));
  EBParticleFieldHalo halo(local.geometry(),MPI_COMM_WORLD);
  FieldVector<Double> origin(0.0);origin[0]=-4.25*5e-6;origin[1]=-5.75*5e-6;origin[2]=-31.25*DZ;
  FieldVector<Double> localOrigin(origin);localOrigin[2]+=offset*DZ;
  const Double xy[6]={0.,.1,.49,3.375,7.51,9.9};
  const Double edge[5]={0.,.1,.49,.5,.51};
  std::vector<Double> z{1.25,static_cast<Double>(localNz)*.5};
  for(Double d:edge){z.push_back(d);z.push_back(localNz-d);}
  for(unsigned pattern=0;pattern<2;++pattern)for(unsigned refresh=0;refresh<2;++refresh) {
    fill(local,offset,pattern,refresh);fill(whole,0,pattern,refresh);
    std::vector<Point> points;
    for(Double x:xy)for(Double y:xy)for(Double q:z) {
      Point p;p.position=origin;p.position[0]+=x*5e-6;p.position[1]+=y*5e-6;
      p.position[2]+=(offset+q)*DZ;
      p.haloNeeded=(rank>0&&q<.5)||(rank+1<ranks&&q>localNz-.5);
      p.physicalEnd=(rank==0&&q<.5)||(rank+1==ranks&&q>localNz-.5);
      local.sampleFieldsPosition(p.position,localOrigin,p.legacyE,p.legacyB);
      FieldVector<Double> e(0.0),b(0.0);whole.sampleFieldsPosition(p.position,origin,e,b);
      if(p.haloNeeded) {
        const Double de=error(e,p.legacyE),db=error(b,p.legacyB);
        negativeE=std::max(negativeE,de);negativeB=std::max(negativeB,db);
        if(de>1e-8||db>1e-8)++negativeSamples;
      }
      points.push_back(p);
    }
    bool rejected=false;
    try {halo.clear(local);}catch(const std::logic_error&){rejected=true;}
    require(rejected,"Clear without install was accepted");
    EBParticleInterpolationHaloView partial;Double dummy=0.;partial.lowerEz=&dummy;
    rejected=false;try{local.setParticleInterpolationHalo(partial);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected&&!local.hasParticleInterpolationHalo(),"Incomplete native halo accepted");
    halo.exchangeAndInstall(local);
    require(local.hasParticleInterpolationHalo()==(ranks>1),"Physical halo mask differs from rank layout");
    rejected=false;try{halo.exchangeAndInstall(local);}catch(const std::logic_error&){rejected=true;}
    require(rejected,"Repeated exchange without clear was accepted");
    for(const Point& p:points) {
      FieldVector<Double> expectedE(0.0),expectedB(0.0),e(0.0),b(0.0);
      whole.sampleFieldsPosition(p.position,origin,expectedE,expectedB);
      local.sampleFieldsPosition(p.position,localOrigin,e,b);
      const Double de=error(expectedE,e),db=error(expectedB,b);
      maximumE=std::max(maximumE,de);maximumB=std::max(maximumB,db);++samples;
      require(de<=2e-12&&db<=2e-12,"MPI particle interpolation differs from global native oracle");
      if(!p.haloNeeded) {
        // .5-cell endpoints can cross a floating-point branch by a few ULPs;
        // the exact old local arithmetic is required away from those endpoints.
        const Double lattice=(p.position[2]-localOrigin[2])/DZ-.5;
        if((lattice>=0.&&lattice<=localNz-1.)||p.physicalEnd) {
          require(equal(e,p.legacyE)&&equal(b,p.legacyB),"Interior/physical sampling changed bitwise");
          ++unchangedSamples;if(p.physicalEnd)++physicalSamples;
        }
      }
      RelativisticParticleSI a,c;
      a.position=p.position;a.properVelocity[0]=.04;a.properVelocity[1]=-.03;a.properVelocity[2]=.23;
      a.charge=-SI::elementaryCharge;a.mass=SI::electronMass;c=a;
      RelativisticBorisPusher::pushFromGrid(a,whole,origin,whole.geometry().dt);
      RelativisticBorisPusher::pushFromGrid(c,local,localOrigin,local.geometry().dt);
      const Double du=error(a.properVelocity,c.properVelocity),dp=error(a.position,c.position);
      maximumU=std::max(maximumU,du);maximumPosition=std::max(maximumPosition,dp);++borisSamples;
      require(du<=2e-12&&dp<=2e-12,"Original Boris push differs with complete neighbour sampling");
    }
    halo.clear(local);require(!local.hasParticleInterpolationHalo(),"Clear retained a stale field view");
    for(const Point& p:points) {
      FieldVector<Double> e(0.0),b(0.0);local.sampleFieldsPosition(p.position,localOrigin,e,b);
      require(equal(e,p.legacyE)&&equal(b,p.legacyB),"Clear did not restore the legacy path bitwise");
    }
  }
  require(halo.exchanges()==4,"Unexpected exchange lifetime count");
  if(ranks==1)require(halo.sentBytes()==0&&halo.receivedBytes()==0&&halo.memoryBytes()==0,"Single rank allocated/communicated particle halo planes");
  else require(halo.sentBytes()>0&&halo.sentBytes()==halo.receivedBytes(),"Missing native plane communication");
}
}

int main(int argc,char** argv) {
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int rank=0,ranks=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  try {
    configureComputeThreads(provided);require(ranks==1||ranks==2||ranks==4,"Use 1, 2 or 4 MPI ranks");
    check(rank,ranks);
    unsigned long long count[5]={samples,negativeSamples,unchangedSamples,physicalSamples,borisSamples},totals[5]={};
    MPI_Allreduce(count,totals,5,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    Double values[6]={maximumE,maximumB,maximumU,maximumPosition,negativeE,negativeB},all[6]={};
    MPI_Allreduce(values,all,6,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(ranks>1)require(totals[1]>0&&all[4]>1e-8&&all[5]>1e-8,"No-halo negative control did not expose truncation");
    if(rank==0) {
      std::ofstream file;if(argc>1)file.open(argv[1]);if(argc>1)require(file.good(),"Cannot open output");
      std::ostream& out=argc>1?static_cast<std::ostream&>(file):std::cout;
      out<<std::setprecision(17)<<"{\"status\":\"PASS\",\"ranks\":"<<ranks
        <<",\"global_nz\":25,\"uneven_slabs\":"<<(ranks>1?"true":"false")
        <<",\"field_samples\":"<<totals[0]<<",\"negative_control_samples\":"<<totals[1]
        <<",\"legacy_unchanged_samples_bitwise\":"<<totals[2]<<",\"physical_end_samples_bitwise\":"<<totals[3]
        <<",\"boris_push_samples\":"<<totals[4]<<",\"maximum_E_relative_error\":"<<all[0]
        <<",\"maximum_B_relative_error\":"<<all[1]<<",\"maximum_Boris_u_relative_error\":"<<all[2]
        <<",\"maximum_Boris_position_relative_error\":"<<all[3]
        <<",\"legacy_missing_plane_E_error\":"<<all[4]<<",\"legacy_missing_plane_B_error\":"<<all[5]
        <<",\"relative_tolerance\":2e-12,\"refreshes_per_pattern\":2,\"patterns\":[\"affine\",\"mixed\"]"
        <<",\"scope\":\"Native neighbour data, original trilinear interpolation and one original Boris push; no full PIC or production acceptance\"}\n";
    }
    MPI_Finalize();return 0;
  }catch(const std::exception& e){std::cerr<<"rank="<<rank<<": "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;}
}
