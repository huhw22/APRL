// Production communicates three layers: two nominal support layers plus one
// IEEE/access guard. No particle/field equation or nonzero weight changes.
#include "source_coupling.h"
#include "source_coupling_math.h"
#include "longitudinal_shape.h"
#include "z_partition.h"
#include "eb_mpi.h"
#include "runtime_util.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using namespace aprl;
namespace {
using Vec=FieldVector<Double>;using Shape=longitudinal_shape::Shape;
unsigned long long currentPaths=0,magneticPoints=0,zeroCorners=0,ulpRejections=0,multicellRejections=0,pushes=0;
unsigned long long binaryExactSamples=0,binaryExactZeroCorners=0,productionSamples=0,roundoffTwoRejections=0,pairedInvalidSeedRejections=0;
double roundoffThirdNodeMagnitude=0.,roundoffThirdEdgeMagnitude=0.;
double maximumCurrentError=0.,maximumWorkPairError=0.,maximumMagneticEnergyError=0.;
void require(bool yes,const char* message){if(!yes)throw std::runtime_error(message);}
bool same(double a,double b){return std::memcmp(&a,&b,sizeof(double))==0;}
EBGridGeometry geometry(std::size_t n){const double dz=1.2543262415359034e-6;return EBGridGeometry(12,14,n,5e-6,5e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);}
std::vector<std::size_t> cutsFor(int size){const std::size_t widths[6]={2,2,2,2,3,8};std::vector<std::size_t> cuts(size+1);for(int r=0;r<size;++r)cuts[r+1]=cuts[r]+widths[r%6];return cuts;}
YeeComponent& component(EBFieldGrid& f,unsigned a){switch(a){case 0:return f.ex();case 1:return f.ey();case 2:return f.ez();case 3:return f.bx();case 4:return f.by();case 5:return f.bz();case 6:return f.jx();case 7:return f.jy();case 8:return f.jz();}throw std::out_of_range("Component");}
void fill(EBFieldGrid& f,std::size_t offset,bool electric=true,bool magnetic=true,unsigned phase=0){
 for(unsigned a=0;a<6;++a){auto& v=component(f,a);for(std::size_t k=0;k<v.nz();++k)for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i){
  const double integer=512*(a+1)+17*i+23*j+31*(k+offset)+37*phase;
  v(i,j,k)=a<3?(electric?std::ldexp(integer,6):0.):(magnetic?std::ldexp(integer,-20):0.);
 }}
}
Vec position(const EBGridGeometry& g,const Vec& origin,double z){Vec p(origin);p[0]+=4.25*g.dx;p[1]+=5.375*g.dy;p[2]+=z*g.dz;return p;}
template<class Visit> void paths(const EBGridGeometry& g,const Vec& origin,Visit visit){
 const double n=g.nz;const double start[7]={0.,std::nextafter(0.,1.),.01,.5*n,n-.01,std::nextafter(n,0.),n};
 const double displacements[7]={-1.,-.99,-.01,0.,.01,.99,1.};unsigned count=0;
 for(double z:start)for(double delta:displacements){Vec a=position(g,origin,z),b(a);b[0]+=.17*g.dx;b[1]-=.23*g.dy;b[2]+=delta*g.dz;visit(a,b,-SI::elementaryCharge*(1+count++%3));}
}
double compareCurrent(EBFieldGrid& a,EBFieldGrid& b){double error=0.,scale=0.;for(unsigned d=6;d<9;++d){auto& x=component(a,d);auto& y=component(b,d);for(std::size_t k=0;k<x.size();++k){error=std::max(error,std::abs(x.data()[k]-y.data()[k]));scale=std::max(scale,std::abs(y.data()[k]));}}
 const double relative=error/std::max(scale,1e-300);maximumCurrentError=std::max(maximumCurrentError,relative);require(relative<2e-12,"Compact current differs from depth-eight current");return relative;}
void currentCheck(int rank,int size){const auto cuts=cutsFor(size);const auto off=cuts[rank],n=cuts[rank+1]-off;EBFieldGrid small(geometry(n)),full(geometry(n));Vec origin(0.);origin[2]=off*small.geometry().dz;
 longitudinal_shape::Current a(small,MPI_COMM_WORLD,Shape::LinearBinomial3,3),b(full,MPI_COMM_WORLD,Shape::LinearBinomial3,8);a.clear();b.clear();
 paths(small.geometry(),origin,[&](const Vec& begin,const Vec& end,double charge){a.deposit(origin,begin,end,charge);b.deposit(origin,begin,end,charge);++currentPaths;});
 a.complete();b.complete();EBZSlabHaloExchange exchange(MPI_COMM_WORLD);exchange.sumSharedTangentialCurrent(small);exchange.sumSharedTangentialCurrent(full);compareCurrent(small,full);
 EBFieldGrid negative(geometry(8));longitudinal_shape::Current guard(negative,MPI_COMM_WORLD,Shape::LinearBinomial3,3);guard.clear();
 Vec negativeOrigin(0.);negativeOrigin[2]=8*rank*negative.geometry().dz;Vec first=position(negative.geometry(),negativeOrigin,2.25),last(first);last[2]+=3*negative.geometry().dz;
 bool failed=false;try{guard.deposit(negativeOrigin,first,last,-SI::elementaryCharge);}catch(const std::exception&){failed=true;}require(failed,"Compact current accepted a multicell orbit");++multicellRejections;
}

// Reproduce the r4 rank-2 IEEE endpoint in a domain which retains its
// third-layer contribution. At four/eight ranks the source is still rank 2;
// at two ranks a common world origin of 4*dz gives the identical SI values.
void roundoffRegression(int rank,int size){
 if(size==1)return;
 const int selected=size>=4?2:0;std::vector<std::size_t> cuts(size+1,0);
 for(int r=0;r<size;++r)cuts[r+1]=cuts[r]+(r+1==size?8:2);
 const auto partition=ZPartition::fromCuts(cuts);const auto off=partition.offset(rank),n=partition.cells(rank);
 EBFieldGrid nominal(geometry(n)),guarded(geometry(n)),full(geometry(n));const double dz=nominal.geometry().dz;
 const double base=selected==2?0.:4*dz;Vec origin(0.);origin[2]=base+off*dz;
 longitudinal_shape::Current two(nominal,MPI_COMM_WORLD,Shape::LinearBinomial3,2),three(guarded,MPI_COMM_WORLD,Shape::LinearBinomial3,3),eight(full,MPI_COMM_WORLD,Shape::LinearBinomial3,8);
 two.clear();three.clear();eight.clear();int rejected=0;
 if(rank==selected){Vec a=position(nominal.geometry(),origin,2.),b(a);b[0]+=.17*nominal.geometry().dx;b[1]-=.23*nominal.geometry().dy;b[2]+=dz;
  const double first=(a[2]-origin[2])/dz,delta=(b[2]-a[2])/dz;
  const double cut=(3.-first)/delta,width=1.-cut,u=.5*(1.+cut+width*.57735026918962576451),q=first+u*delta;
  require(first>2.&&delta>1.&&width>0.&&q>3.&&std::floor(q)==3.,"Required production-grid IEEE third-support case disappeared");
  const double weighted=.5*width*.25*(q-std::floor(q));require(weighted>0.,"Third node weight was rounded or thresholded to zero");
  try{two.deposit(origin,a,b,-SI::elementaryCharge);}catch(const std::out_of_range&){rejected=1;++roundoffTwoRejections;}
  three.deposit(origin,a,b,-SI::elementaryCharge);eight.deposit(origin,a,b,-SI::elementaryCharge);
 }
 int totalRejected=0;MPI_Allreduce(&rejected,&totalRejected,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);require(totalRejected==1,"Nominal two-layer IEEE negative control did not reject exactly its source rank");
 three.complete();eight.complete();EBZSlabHaloExchange exchange(MPI_COMM_WORLD);exchange.sumSharedTangentialCurrent(guarded);exchange.sumSharedTangentialCurrent(full);
 for(unsigned a=6;a<9;++a){auto& x=component(guarded,a);auto& y=component(full,a);require(x.size()==y.size(),"Roundoff regression current shape changed");for(std::size_t i=0;i<x.size();++i)require(same(x.data()[i],y.data()[i]),"Three layers did not preserve depth-eight IEEE current bitwise");}
 const std::size_t node=cuts[selected]+5,edge=cuts[selected]+4;
 require(node<=partition.totalCells()&&edge<partition.totalCells(),"IEEE regression must retain its nonzero third layer in the global domain");
 if(partition.ownerNode(node)==rank){const auto k=node-off;const auto& j=guarded.jx();for(std::size_t y=0;y<j.ny();++y)for(std::size_t x=0;x<j.nx();++x)roundoffThirdNodeMagnitude=std::max(roundoffThirdNodeMagnitude,std::abs(j(x,y,k)));}
 if(partition.ownerCell(edge)==rank){const auto k=edge-off;const auto& j=guarded.jz();for(std::size_t y=0;y<j.ny();++y)for(std::size_t x=0;x<j.nx();++x)roundoffThirdEdgeMagnitude=std::max(roundoffThirdEdgeMagnitude,std::abs(j(x,y,k)));}
 double local[2]={roundoffThirdNodeMagnitude,roundoffThirdEdgeMagnitude},global[2]={};MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
 require(global[0]>0.&&global[1]>0.,"Three-layer IEEE regression lost a real retained current contribution");
}

// Reconstruct actual native accesses independently after IEEE normalization.
// Nominal-two tests skip only exact zero corners. Production-three tests keep
// the original sampler, including zero-weight corner accesses; no truncation.
void support(const EBFieldGrid& f,const Vec& origin,const Vec& point,int rank,int size,unsigned depth,bool skipZero,bool& nonzeroOutside,bool& zeroOutside){
 nonzeroOutside=zeroOutside=false;const double z=(point[2]-origin[2])/f.geometry().dz;
 for(unsigned a=0;a<3;++a){const auto& v=a==0?f.bx():(a==1?f.by():f.bz());const double coordinate=z-(a==2?0.:.5),base=std::floor(coordinate),fraction=coordinate-base;
  for(unsigned corner=0;corner<2;++corner){const long long k=static_cast<long long>(base)+corner;const double weight=corner?fraction:1.-fraction;
   const bool physical=(k<0&&rank==0)||(k>=static_cast<long long>(v.nz())&&rank+1==size);
   const bool outside=!physical&&(k< -static_cast<long long>(depth)||k>=static_cast<long long>(v.nz()+depth));
   if(outside){if(skipZero&&weight==0.)zeroOutside=true;else nonzeroOutside=true;}
  }
 }
}
void magneticCheck(int rank,int size,bool binaryExact,unsigned depth,bool skipZero){const auto cuts=cutsFor(size);const auto off=cuts[rank],n=cuts[rank+1]-off;
 const auto beforeSamples=magneticPoints,beforeZeros=zeroCorners,beforeRejections=ulpRejections;
 EBGridGeometry g=geometry(n);if(binaryExact){const double dz=std::ldexp(1.,-20);g=EBGridGeometry(12,14,n,4*dz,4*dz,dz,dz/SI::c,EBMaxwellSolver::CowanZ);}
 EBFieldGrid f(g);fill(f,off);
 source_coupling_math::NativeHalo small,full;small.exchange(f,MPI_COMM_WORLD,true,depth);full.exchange(f,MPI_COMM_WORLD,true,8);Vec origin(0.);origin[2]=off*f.geometry().dz;
 const double middle[12]={-.5,std::nextafter(-.5,-1.),std::nextafter(-.5,0.),0.,.25,.5,.5*n,n-.5,n+.5,std::nextafter(n+.5,0.),std::nextafter(n+.5,std::numeric_limits<double>::infinity()),n};
 for(double m:middle)for(int shift=-1;shift<=1;++shift){Vec p=position(f.geometry(),origin,m+shift);bool outside=false,zero=false;support(f,origin,p,rank,size,depth,skipZero,outside,zero);
  const auto expected=source_coupling_math::magneticSample(f,origin,p,false,&full,true,false);
  bool failed=false;source_coupling_math::MagneticSample actual;try{actual=source_coupling_math::magneticSample(f,origin,p,false,&small,true,skipZero);}catch(const std::out_of_range&){failed=true;}
  if(outside){require(failed,"Compact magnetic sample silently discarded a nonzero out-of-support corner");++ulpRejections;}
  else{require(!failed,"Compact magnetic sample read an unnecessary zero corner outside its collar");for(unsigned a=0;a<3;++a)require(same(actual.value[a],expected.value[a]),"Compact magnetic gather changed a nonzero arithmetic result");++magneticPoints;if(zero)++zeroCorners;}
 }
 unsigned long long local[2]={zeroCorners-beforeZeros,ulpRejections-beforeRejections},global[2]={};MPI_Allreduce(local,global,2,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
 if(binaryExact&&depth==2&&size>1)require(global[0]>0,"Binary-exact zero third-corner regression was not exercised");
 if(binaryExact){binaryExactSamples+=magneticPoints-beforeSamples;binaryExactZeroCorners+=zeroCorners-beforeZeros;}else productionSamples+=magneticPoints-beforeSamples;
 // Depending on coordinate rounding a particular SI nextafter point can
 // collapse onto the endpoint. Explicitly place a genuinely nonzero third
 // corner in the interior-rank collar so rejection is always tested.
 if(size>1){Vec p=position(f.geometry(),origin,rank+1<size?double(n)+double(depth)-.25:-double(depth)+.25);bool outside=false,zero=false;support(f,origin,p,rank,size,depth,skipZero,outside,zero);require(outside,"Nonzero outside-collar negative control is not outside");bool failed=false;try{source_coupling_math::magneticSample(f,origin,p,false,&small,true,skipZero);}catch(const std::out_of_range&){failed=true;}require(failed,"Compact magnetic gather accepted actual outside-collar support");++ulpRejections;}
}

void constraintCheck(int rank,int size,const std::string& directory){const auto cuts=cutsFor(size);const auto off=cuts[rank],n=cuts[rank+1]-off;EBFieldGrid small(geometry(n)),full(geometry(n));fill(small,off);fill(full,off);Vec origin(0.);origin[2]=off*small.geometry().dz;
 createDirectories(directory+"/small",MPI_COMM_WORLD);createDirectories(directory+"/full",MPI_COMM_WORLD);
 ConservativeSourceCoupling a(MPI_COMM_WORLD,directory+"/small",false,false,true,true),b(MPI_COMM_WORLD,directory+"/full",false,false,true,false);
 a.begin(small,origin,{});b.begin(full,origin,{});fill(small,off,true,true,1);fill(full,off,true,true,1);
 paths(small.geometry(),origin,[&](const Vec& first,const Vec& last,double charge){a.depositFrozenSourceSupport(first,last,charge,true);b.depositFrozenSourceSupport(first,last,charge,true);});
 a.completeCurrent();b.completeCurrent();EBZSlabHaloExchange exchange(MPI_COMM_WORLD);exchange.sumSharedTangentialCurrent(small);exchange.sumSharedTangentialCurrent(full);compareCurrent(small,full);
 const auto wa=a.finalizeConstraints(),wb=b.finalizeConstraints();for(unsigned i=0;i<wa.size();++i)require(wa[i]==wb[i],"Compact constraint electric gather/work changed");
 long double work=0.,absolute=0.;for(unsigned d=0;d<3;++d){const auto& j=component(small,d+6);const auto& e=component(small,d);const auto& g=small.geometry();for(std::size_t k=0;k<j.nz();++k)for(std::size_t y=0;y<j.ny();++y)for(std::size_t x=0;x<j.nx();++x){
  // Old and new electric fields differ by one known dyadic constant.
  long double weight=1.;if(d<2&&((rank>0&&k==0)||(rank+1<size&&k==g.nz)))weight=.5;
  const long double value=weight*((long double)e(x,y,k)-std::ldexp(37.,5))*j(x,y,k)*g.dt*g.dx*g.dy*g.dz;work+=value;absolute+=std::abs(value);
 }}
 const long double local[3]={work,wa[0],absolute};long double global[3]={};MPI_Allreduce(local,global,3,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
 maximumWorkPairError=double(std::abs(global[0]+global[1])/std::max(global[2],1e-100L));require(maximumWorkPairError<2e-12,"Compact current and electric feedback lost adjoint work pairing");a.close(0);b.close(0);
}
void particleEqual(const RelativisticParticleSI& a,const RelativisticParticleSI& b){require(a.id==b.id&&a.sourceId==b.sourceId&&same(a.charge,b.charge)&&same(a.mass,b.mass)&&same(a.weight,b.weight),"Compact coupling changed particle identity");for(unsigned d=0;d<3;++d)require(same(a.position[d],b.position[d])&&same(a.properVelocity[d],b.properVelocity[d]),"Compact coupling changed particle state");}
void workEqual(const BorisWork& a,const BorisWork& b){require(a.gridWork==b.gridWork&&a.prescribedWork==b.prescribedWork&&a.totalDeltaK==b.totalDeltaK&&a.magneticRotationDeltaK==b.magneticRotationDeltaK&&a.roundoffResidual==b.roundoffResidual&&a.absoluteWork==b.absoluteWork,"Compact coupling changed replay work");for(unsigned k=0;k<6;++k)require(a.compensation[k]==b.compensation[k],"Compact coupling changed work accumulation state");}
void pushCheck(int rank,int size,const std::string& directory){const auto cuts=cutsFor(size);const auto off=cuts[rank],n=cuts[rank+1]-off;const double z[5]={0.,.01,.5*double(n),double(n)-.01,double(n)};
 for(unsigned mode=0;mode<3;++mode)for(unsigned point=0;point<5;++point)for(int direction:{-1,1}){
  const std::string name="case-"+std::to_string(mode)+"-"+std::to_string(point)+"-"+std::to_string(direction);
  createDirectories(directory+"/"+name+"-small",MPI_COMM_WORLD);createDirectories(directory+"/"+name+"-full",MPI_COMM_WORLD);
  EBFieldGrid small(geometry(n)),full(geometry(n));fill(small,off,mode!=1,mode!=0);fill(full,off,mode!=1,mode!=0);Vec origin(0.);origin[2]=off*small.geometry().dz;
  RelativisticParticleSI initial;initial.id=rank*100+point+1;initial.sourceId=initial.id;initial.position=position(small.geometry(),origin,z[point]);initial.charge=-SI::elementaryCharge;initial.mass=SI::electronMass;initial.properVelocity=Vec(0.);initial.properVelocity[0]=.23;initial.properVelocity[1]=-.17;initial.properVelocity[2]=direction*(point%2?10.:.7);
  ConservativeSourceCoupling a(MPI_COMM_WORLD,directory+"/"+name+"-small",false,false,true,true),b(MPI_COMM_WORLD,directory+"/"+name+"-full",false,false,true,false);
  a.begin(small,origin,{initial});b.begin(full,origin,{initial});a.exchangeCandidate();b.exchangeCandidate();SIFieldSourceSet sources;BoostFrameTransform frame;frame.setFromGamma(1.,SI::c,0.);RelativisticParticleSI pa,pb;
  // Independent geometric interval check for a seed at the true global
  // sacrificial face. Normalization can put the intended endpoint just
  // outside; if its whole initial orbit is outside, both modes must reject.
  Vec initialGuess(initial.position);initialGuess.pmv(small.geometry().dt*SI::c/BoostFrameTransform::gammaFromProperVelocity(initial.properVelocity),initial.properVelocity);
  const double z0=(initial.position[2]-origin[2])/small.geometry().dz;
  const double z1=z0+(initialGuess[2]-initial.position[2])/small.geometry().dz;
  const bool noGlobalDomain=(rank==0&&std::max(z0,z1)<0.)||(rank+1==size&&std::min(z0,z1)>n);
  bool invalidSeed=false;
  for(unsigned trial=0;trial<3;++trial){a.beginTrial();b.beginTrial();pa=initial;pb=initial;bool rejectedSmall=false,rejectedFull=false;std::string smallError,fullError;
   try{a.push(pa,sources,frame,0.,small.geometry().dt,7,nullptr,nullptr,nullptr);}catch(const std::out_of_range& error){rejectedSmall=true;smallError=error.what();}
   try{b.push(pb,sources,frame,0.,full.geometry().dt,7,nullptr,nullptr,nullptr);}catch(const std::out_of_range& error){rejectedFull=true;fullError=error.what();}
   if(rejectedSmall||rejectedFull){
    if(!(rejectedSmall&&rejectedFull&&noGlobalDomain&&trial==0&&smallError=="Quintic active orbit outside domain"&&fullError==smallError)){
     std::cerr<<std::setprecision(18)<<"Push rejection mismatch rank="<<rank<<" case="<<name<<" trial="<<trial<<" origin_z="<<origin[2]<<" begin_z="<<initial.position[2]<<" guess_z="<<initialGuess[2]<<" normalized_begin="<<z0<<" normalized_end="<<z1<<" compact_error="<<smallError<<" full_error="<<fullError<<'\n';
     throw std::runtime_error("Production three-layer and depth-eight rejection are not equivalent");
    }
    ++pairedInvalidSeedRejections;invalidSeed=true;break;
   }
   require(!noGlobalDomain,"Both couplers accepted a seed whose complete initial orbit misses the true global domain");particleEqual(pa,pb);
   const auto ra=a.localResidual(),rb=b.localResidual();for(unsigned i=0;i<ra.size();++i)require(same(ra[i],rb[i]),"Compact initial/residual electric or magnetic force changed");++pushes;
  }
  if(invalidSeed){a.close(0);b.close(0);continue;}
  if(mode==1){const auto delta=source_coupling_math::deltaGamma(initial.properVelocity,pa.properVelocity);const double error=double(std::abs(delta)/source_coupling_math::gammaLD(initial.properVelocity));maximumMagneticEnergyError=std::max(maximumMagneticEnergyError,error);require(error<2e-13,"Magnetic-only coupling performed spurious electric work");}
  a.setReplay(true);b.setReplay(true);pa=initial;pb=initial;BorisWork wa,wb;a.push(pa,sources,frame,0.,small.geometry().dt,7,nullptr,&wa,nullptr);b.push(pb,sources,frame,0.,full.geometry().dt,7,nullptr,&wb,nullptr);particleEqual(pa,pb);workEqual(wa,wb);
  a.close(0);b.close(0);
 }
 // The compact coupler must reject a super-CFL gather even if it would stay
 // inside a large local slab. Default depth-eight multi-cell tests stay valid.
 EBFieldGrid f(geometry(8));fill(f,8*rank);Vec origin(0.);origin[2]=8*rank*f.geometry().dz;RelativisticParticleSI p;p.id=p.sourceId=1;p.charge=-SI::elementaryCharge;p.mass=SI::electronMass;p.position=position(f.geometry(),origin,3.);p.properVelocity=Vec(0.);p.properVelocity[2]=10.;
 createDirectories(directory+"/multicell-guard",MPI_COMM_WORLD);ConservativeSourceCoupling guard(MPI_COMM_WORLD,directory+"/multicell-guard",false,false,true,true);guard.begin(f,origin,{p});guard.exchangeCandidate();SIFieldSourceSet sources;BoostFrameTransform frame;frame.setFromGamma(1.,SI::c,0.);bool failed=false;try{guard.push(p,sources,frame,0.,3*f.geometry().dt,7,nullptr,nullptr,nullptr);}catch(const std::out_of_range&){failed=true;}require(failed,"Compact coupling accepted a super-CFL multicell gather");++multicellRejections;guard.close(0);
}
}
int main(int argc,char** argv){MPI_Init(&argc,&argv);int rank=0,size=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);try{require(argc==2,"Need independent output directory");createDirectories(argv[1],MPI_COMM_WORLD);currentCheck(rank,size);roundoffRegression(rank,size);magneticCheck(rank,size,true,2,true);magneticCheck(rank,size,false,3,false);constraintCheck(rank,size,std::string(argv[1])+"/constraints");pushCheck(rank,size,std::string(argv[1])+"/push");
 unsigned long long local[11]={currentPaths,magneticPoints,zeroCorners,ulpRejections,multicellRejections,pushes,binaryExactSamples,binaryExactZeroCorners,productionSamples,roundoffTwoRejections,pairedInvalidSeedRejections},total[11]={};MPI_Allreduce(local,total,11,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);if(size==8)require(total[10]>=3,"Expected global-face same-rejection controls were not exercised");double values[5]={maximumCurrentError,maximumWorkPairError,maximumMagneticEnergyError,roundoffThirdNodeMagnitude,roundoffThirdEdgeMagnitude},maximum[5]={};MPI_Allreduce(values,maximum,5,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
 if(rank==0)std::cout<<std::setprecision(18)<<"{\"passed\":true,\"physical_acceptance\":false,\"ranks\":"<<size<<",\"one_cell_current_paths\":"<<total[0]<<",\"magnetic_samples_bitwise\":"<<total[1]<<",\"zero_outside_corner_cases\":"<<total[2]<<",\"nonzero_outside_support_rejections\":"<<total[3]<<",\"multicell_rejections\":"<<total[4]<<",\"initial_and_residual_pushes_bitwise\":"<<total[5]<<",\"binary_exact_magnetic_samples_bitwise\":"<<total[6]<<",\"binary_exact_zero_outside_corner_cases\":"<<total[7]<<",\"production_grid_magnetic_samples_bitwise\":"<<total[8]<<",\"production_halo_depth\":3,\"roundoff_nominal_two_rejections\":"<<total[9]<<",\"roundoff_retained_third_node_current\":"<<maximum[3]<<",\"roundoff_retained_third_edge_current\":"<<maximum[4]<<",\"paired_invalid_global_face_seed_rejections\":"<<total[10]<<",\"maximum_current_relative_difference\":"<<maximum[0]<<",\"maximum_adjoint_work_relative_error\":"<<maximum[1]<<",\"maximum_magnetic_only_relative_energy_error\":"<<maximum[2]<<"}\n";
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}MPI_Finalize();return 0;}
