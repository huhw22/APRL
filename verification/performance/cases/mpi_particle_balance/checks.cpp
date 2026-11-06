// Independent partition oracle, thin-slab routing, and native state transport.
// This is program acceptance, not a new FEL/space-charge physical benchmark.
#include "particle_balance.h"
#include "z_partition.h"
#include "source_coupling_math.h"
#include "longitudinal_shape.h"
#include "moving_window_fields.h"
#include "eb_mpi.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

using namespace aprl;
namespace {
using Vec=FieldVector<Double>;
using Shape=longitudinal_shape::Shape;
unsigned long long plannerCases=0,rejections=0,haloValues=0,stateValues=0;
double plannerError=0.,currentError=0.,continuationError=0.;
void require(bool good,const char* message){if(!good)throw std::runtime_error(message);}
bool same(double a,double b){return std::memcmp(&a,&b,sizeof(double))==0;}
void rejected(const std::function<void()>& action){bool failed=false;try{action();}catch(const std::exception&){failed=true;}require(failed,"Invalid planner input accepted");++rejections;}

long double oracle(const std::vector<double>& w,std::size_t ranks,std::size_t minimum,std::size_t maximum){
 const auto n=w.size();std::vector<std::vector<long double> > interval(n+1,std::vector<long double>(n+1,0.));
 for(std::size_t first=0;first<n;++first)for(std::size_t last=first+1;last<=n;++last)interval[first][last]=interval[first][last-1]+w[last-1];
 const long double inf=std::numeric_limits<long double>::infinity();
 std::vector<long double> previous(n+1,inf),next(n+1,inf);previous[0]=0.;
 for(std::size_t p=0;p<ranks;++p){std::fill(next.begin(),next.end(),inf);
  for(std::size_t end=minimum;end<=n;++end)for(std::size_t length=minimum;length<=maximum&&length<=end;++length)
   if(std::isfinite(previous[end-length]))next[end]=std::min(next[end],std::max(previous[end-length],interval[end-length][end]));
  previous.swap(next);
 }return previous[n];
}
void onePlanner(const std::vector<double>& work,std::size_t ranks,std::size_t minimum,std::size_t maximum=0){
 const auto original=work;const auto cuts=particle_balance::plan(work,ranks,minimum,maximum);
 require(work==original,"Planner modified its input cost array");
 require(cuts.size()==ranks+1&&cuts.front()==0&&cuts.back()==work.size(),"Planner coverage invalid");
 long double observed=0.,sum=0.;
 for(std::size_t r=0;r<ranks;++r){const auto length=cuts[r+1]-cuts[r];require(length>=minimum&&(!maximum||length<=maximum),"Planner cell limits violated");
  long double value=0.;for(std::size_t k=cuts[r];k<cuts[r+1];++k)value+=work[k];observed=std::max(observed,value);sum+=value;}
 const auto expected=oracle(work,ranks,minimum,maximum?std::min(maximum,work.size()):work.size());
 const double error=double(std::abs(observed-expected)/std::max(std::abs(expected),std::numeric_limits<long double>::min()));
 plannerError=std::max(plannerError,error);require(error<2e-13,"Integer planner differs from independent exact DP optimum");
 const auto measured=particle_balance::measure(work,cuts);
 require(std::abs(measured.total-sum)<2e-13L*std::max(sum,1e-300L),"Reported partition total differs from independent sum");
 require(std::abs(measured.maximum-observed)<2e-13L*std::max(observed,1e-300L),"Reported maximum differs from independent sum");
 require(std::abs(measured.imbalance-(sum>0?observed/(sum/ranks):1.L))<2e-13L,"Reported true imbalance differs");
 require(cuts==particle_balance::plan(work,ranks,minimum,maximum),"Planner tie handling is nondeterministic");++plannerCases;
}
void planner(){
 onePlanner({0,0,10,0,10,10,0,0},3,2); // Greedy longest-prefix fails at limit 10.
 onePlanner({0,0,10,0,10,10,0,0},3,2,3);
 onePlanner(std::vector<double>(16,0.),4,2);onePlanner(std::vector<double>(17,1.),4,2,5);
 onePlanner(std::vector<double>(16,1.),8,2);onePlanner({1e-100,2e-100,9e-100,1e-100,2e-100,3e-100},2,2);
 onePlanner({1e200,2e200,9e200,1e200,2e200,3e200},2,2);
 std::mt19937 random(20261106);
 for(unsigned n=0;n<512;++n){const std::size_t cells=2+random()%23,minimum=1+random()%std::min<std::size_t>(3,cells),ranks=1+random()%(cells/minimum);
  const std::size_t leastMax=cells/ranks+(cells%ranks?1:0),maximum=leastMax+random()%(cells-leastMax+1);
  std::vector<double> work(cells);for(auto& x:work)x=(random()%23)*.125+(n%3?.03125:0.);
  onePlanner(work,ranks,minimum,n%2?maximum:0);
 }
 const std::vector<double> w={4,4,4,4,1,1,1,1};const std::vector<std::size_t> old={0,4,8},tied={0,2,8};
 require(particle_balance::plan(w,2,2,0,old,.3)==old,"Hysteresis did not retain insufficient improvement");
 require(particle_balance::plan(w,2,2,0,old,.2)!=old,"Hysteresis rejected sufficient improvement");
 require(particle_balance::plan(w,2,2,0,tied,0.)==tied,"Planner did not retain an equally optimal current partition");
 rejected([]{particle_balance::plan({},1);});rejected([]{particle_balance::plan({1,1},0);});
 rejected([]{particle_balance::plan({1,1},1,0);});rejected([]{particle_balance::plan({1,1,1},2,2);});
 rejected([]{particle_balance::plan(std::vector<double>(9,1.),4,2,2);});
 rejected([]{particle_balance::plan({1,-1,1},1);});rejected([]{particle_balance::plan({1,std::numeric_limits<double>::infinity()},1);});
 rejected([]{particle_balance::plan({1,std::numeric_limits<double>::quiet_NaN()},1);});
 rejected([]{particle_balance::plan({1,1},1,2,0,{},-.1);});
 rejected([]{particle_balance::plan({1,1},1,2,0,{},std::numeric_limits<double>::quiet_NaN());});
 rejected([]{particle_balance::plan({1,1,1,1},2,2,0,{0,1,4});});
 rejected([]{particle_balance::measure({1,1},{0,0,2});});
}

YeeComponent& component(EBFieldGrid& f,unsigned a){switch(a){case 0:return f.ex();case 1:return f.ey();case 2:return f.ez();case 3:return f.bx();case 4:return f.by();case 5:return f.bz();case 6:return f.jx();case 7:return f.jy();case 8:return f.jz();}throw std::out_of_range("Test component");}
EBGridGeometry geometry(std::size_t nz){const double dz=1.2543262415359034e-6;return EBGridGeometry(12,14,nz,5e-6,5e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);}
double mark(unsigned a,std::size_t i,std::size_t j,long long k){return std::ldexp(double(1+a),22)+std::ldexp(double(k+1024),10)+32.*j+i;}
std::vector<std::size_t> mixed(int size,bool tiny){std::vector<std::size_t> cuts(size+1,0);const std::size_t values[6]={2,2,2,2,3,8};for(int r=0;r<size;++r)cuts[r+1]=cuts[r]+(tiny?values[r%6]:12);return cuts;}
void fill(EBFieldGrid& f,std::size_t offset){for(unsigned a=0;a<9;++a){auto& v=component(f,a);for(std::size_t k=0;k<v.nz();++k)for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i)v(i,j,k)=mark(a,i,j,offset+k);}}
void nativeHalo(int rank,int size){
 // Pure indexing can represent a one-cell slab; the unchanged EBFieldGrid
 // deliberately requires at least two cells. Do not weaken that kernel guard.
 const auto oneCell=ZPartition::fromCuts({0,1,3,6,14});
 require(oneCell.minimumCells()==1&&oneCell.ownerCell(0)==0&&oneCell.ownerNode(1)==1&&oneCell.ownerNode(14)==3,"One-cell integer ownership convention differs");
 const auto cuts=mixed(size,true);const auto partition=ZPartition::fromCuts(cuts);const auto gathered=ZPartition::gather(partition.cells(rank),MPI_COMM_WORLD);require(gathered.cuts==cuts,"Gathered integer partition differs");
 for(std::size_t k=0;k<partition.totalCells();++k){const int owner=partition.ownerCell(k);require(cuts[owner]<=k&&k<cuts[owner+1],"Integer cell owner wrong");require(partition.ownerNode(k)==owner,"Shared node does not use right owner");}
 require(partition.ownerNode(partition.totalCells())==size-1,"Final node owner wrong");
 EBFieldGrid f(geometry(partition.cells(rank)));fill(f,partition.offset(rank));
 for(bool magnetic:{false,true}){source_coupling_math::NativeHalo halo;halo.exchange(f,MPI_COMM_WORLD,magnetic,8);
  for(unsigned a=0;a<3;++a){auto& v=component(f,a+(magnetic?3:0));const bool node=magnetic?a==2:a<2;const long long last=partition.totalCells()-(node?0:1);
   for(long long k=-8;k<static_cast<long long>(v.nz())+8;++k){const long long global=static_cast<long long>(partition.offset(rank))+k;
    for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i){
     if((k>=0&&k<static_cast<long long>(v.nz()))||halo.contains(k,v.nz())){const double want=magnetic?mark(a+3,i,j,std::max(0LL,std::min(global,last))):(global<0||global>last?0.:mark(a,i,j,global));require(same(halo.value(v,a,i,j,k),want),"Thin multi-rank native halo differs bitwise");++haloValues;}
     if(magnetic){const double want=mark(a+3,i,j,std::max(0LL,std::min(global,last)));require(same(halo.openMagneticValue(v,a,i,j,k),want),"Open magnetic halo lost real global endpoint");++haloValues;}
    }
   }
  }
 }
}
void routedCurrent(int rank,int size,Shape shape){
 auto cuts=mixed(size,true);if(size==1)cuts[1]=2;const auto partition=ZPartition::fromCuts(cuts);const auto off=partition.offset(rank),n=partition.totalCells();
 EBGridGeometry g=geometry(partition.cells(rank)),global=geometry(n);EBFieldGrid f(g),whole(global);
 longitudinal_shape::Current distributed(f,MPI_COMM_WORLD,shape),serial(whole,MPI_COMM_SELF,shape);distributed.clear();serial.clear();
 Vec origin(0.),local(0.);local[2]=off*g.dz;const unsigned particles=static_cast<unsigned>(2*n+17);
 for(unsigned p=0;p<particles;++p){Vec begin(0.);begin[0]=(4.13+.031*(p%17))*g.dx;begin[1]=(5.21+.037*(p%13))*g.dy;
  begin[2]=(.25+.5*((p*37)%particles)/double(particles))*n*g.dz;Vec end(begin);end[0]+=.13*g.dx;end[1]-=.17*g.dy;end[2]+=(p%2?.071:-.083)*g.dz;
  const double charge=-SI::elementaryCharge*(1+p%3);serial.deposit(origin,begin,end,charge);
  const double a=begin[2]/g.dz,b=end[2]/g.dz,lower=off,upper=off+g.nz;double lo=0.,hi=1.;
  if(b>a){lo=std::max(lo,(lower-a)/(b-a));hi=std::min(hi,(upper-a)/(b-a));}
  else{lo=std::max(lo,(upper-a)/(b-a));hi=std::min(hi,(lower-a)/(b-a));}
  if(hi>lo){Vec first(begin),last(begin);for(unsigned d=0;d<3;++d){first[d]+=lo*(end[d]-begin[d]);last[d]+=hi*(end[d]-begin[d]);}distributed.deposit(local,first,last,charge);}
 }
 distributed.complete();EBZSlabHaloExchange exchange(MPI_COMM_WORLD);exchange.sumSharedTangentialCurrent(f);
 double error=0.,scale=0.;for(unsigned a=6;a<9;++a){auto& v=component(f,a);auto& expected=component(whole,a);for(std::size_t k=0;k<v.nz();++k)for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i){error=std::max(error,std::abs(v(i,j,k)-expected(i,j,k+off)));scale=std::max(scale,std::abs(expected(i,j,k+off)));}}
 currentError=std::max(currentError,error/std::max(scale,1e-300));require(error<2e-11*std::max(scale,1e-300),"Thin multi-rank source differs from serial complete source");
}

EBCPMLParameters parameters(){EBCPMLParameters p;p.cells[0]=p.cells[1]=2;p.cells[2]=0;p.polynomialOrder=3.;p.kappaMax=8.;p.alphaFraction=0.;p.targetReflection=1e-8;return p;}
struct History {std::size_t row,layers,slots;};
History history(unsigned h,std::size_t nz){switch(h){case 0:return {13,nz,4};case 1:case 2:return {0,nz,0};case 3:return {15,nz,4};case 4:return {14,nz+1,4};case 5:case 6:return {12,nz+1,4};case 7:case 8:return {0,nz+1,0};case 9:return {14,nz+1,4};case 10:return {15,nz,4};case 11:return {13,nz,4};}throw std::out_of_range("Test history");}
double historyMark(unsigned h,std::size_t slot,std::size_t q,std::size_t k){return mark(h+9,slot,q,k);}
void histories(EBFieldGrid& f,std::size_t offset,bool write){auto& cpml=dynamic_cast<EBConvolutionalPML&>(f.boundary());require(cpml.checkpointHistoryCount()==12,"Wrong CPML history count");for(unsigned h=0;h<12;++h){auto& v=cpml.checkpointHistory(h);const auto s=history(h,f.geometry().nz);require(v.size()==s.row*s.layers*s.slots,"Independent CPML history layout differs");for(std::size_t slot=0;slot<s.slots;++slot)for(std::size_t k=0;k<s.layers;++k)for(std::size_t q=0;q<s.row;++q){const auto index=(slot*s.layers+k)*s.row+q;const double want=historyMark(h,slot,q,offset+k);if(write)v[index]=want;else{require(same(v[index],want),"Repartition changed a CPML history bit");++stateValues;}}}}
void state(EBFieldGrid& f,std::size_t offset){for(unsigned a=0;a<9;++a){auto& v=component(f,a);for(std::size_t k=0;k<v.nz();++k)for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i){require(same(v(i,j,k),mark(a,i,j,offset+k)),"Repartition changed a native field/current bit");++stateValues;}}histories(f,offset,false);}
double reference(void*,unsigned a,std::size_t i,std::size_t j,std::ptrdiff_t k,std::size_t){return mark(a,i,j,k);}
void observe(void* context,const char*,const EBFieldGrid&){++*static_cast<unsigned*>(context);}
void ledgerEqual(const MovingWindowFieldEnergyLedger& a,const MovingWindowFieldEnergyLedger& b){require(a.enabled==b.enabled&&a.inletOverwriteEnergyJ==b.inletOverwriteEnergyJ&&a.insertedEnergyJ==b.insertedEnergyJ&&a.croppedEnergyJ==b.croppedEnergyJ&&a.remapEnergyJ==b.remapEnergyJ&&a.remapInterfaceCorrectionJ==b.remapInterfaceCorrectionJ&&a.outwardEnergyFace==b.outwardEnergyFace&&a.overwriteEvents==b.overwriteEvents&&a.fluxSteps==b.fluxSteps&&a.remapSteps==b.remapSteps&&a.diagnosticCellSamples==b.diagnosticCellSamples&&a.diagnosticSeconds==b.diagnosticSeconds,"Repartition changed cumulative ledger state");}
void step(MovingWindowFields& window,EBFieldGrid& f,std::size_t number,std::ptrdiff_t offset){f.clearCurrent();EBZSlabHaloExchange exchange(MPI_COMM_WORLD);window.beginStep(f,number,offset);exchange.advanceMagnetic(f,window.physicalElectricHalos(f,number));exchange.advanceElectric(f);window.finishEnergyLedgerStep(f);window.shift(f,number+1,offset-1);}
void repartition(int rank,int size){
 const std::size_t total=12*size;std::vector<std::size_t> cuts(size+1,0);const std::size_t lengths[3]={2,3,8};for(int r=0;r+1<size;++r)cuts[r+1]=cuts[r]+lengths[r%3];cuts[size]=total;
 const auto next=ZPartition::fromCuts(cuts);const std::size_t oldOffset=12*rank;
 EBFieldGrid old(geometry(12)),moved(geometry(next.cells(rank))),back(geometry(12));MovingWindowReference ref;ref.value=reference;
 MovingWindowFields first(old.geometry(),total,oldOffset,MPI_COMM_WORLD,parameters(),ref);first.installBoundary(old);first.enableEnergyLedger();unsigned callbacks=0;first.setStageObserver(observe,&callbacks);
 fill(old,oldOffset);histories(old,oldOffset,true);const auto before=first.energyLedger();
 auto second=first.repartition(old,moved,next.offset(rank));ledgerEqual(before,second->energyLedger());require(callbacks==0,"Repartition invoked a physical observer");
 state(old,oldOffset);state(moved,next.offset(rank));auto third=second->repartition(moved,back,oldOffset);state(back,oldOffset);ledgerEqual(before,third->energyLedger());
 // Propagate the transported nonzero CPML histories on both partitions.
 step(first,old,0,0);step(*second,moved,0,0);require(callbacks>0,"Repartition dropped the stage observer");
 EBFieldGrid expected(geometry(next.cells(rank)));auto after=first.repartition(old,expected,next.offset(rank));ledgerEqual(first.energyLedger(),after->energyLedger());
 double error=0.,scale=0.;for(unsigned a=0;a<9;++a){auto& v=component(moved,a);auto& w=component(expected,a);for(std::size_t i=0;i<v.size();++i){error=std::max(error,std::abs(v.data()[i]-w.data()[i]));scale=std::max(scale,std::abs(w.data()[i]));}}
 auto& cpml=dynamic_cast<EBConvolutionalPML&>(moved.boundary());auto& expectedPml=dynamic_cast<EBConvolutionalPML&>(expected.boundary());
 for(unsigned h=0;h<12;++h){const auto& a=cpml.checkpointHistory(h);const auto& b=expectedPml.checkpointHistory(h);require(a.size()==b.size(),"Continuation history size differs");for(std::size_t i=0;i<a.size();++i){error=std::max(error,std::abs(a[i]-b[i]));scale=std::max(scale,std::abs(b[i]));}}
 continuationError=error/std::max(scale,1e-300);require(continuationError<2e-12,"Repartition changed continued Cowan/CPML evolution");
 // A transported established state must accept its next exact step/offset,
 // and collective repartition during an unfinished field step must fail.
 after->beginStep(expected,1,-1);EBFieldGrid activeDestination(geometry(next.cells(rank)));
 bool activeRejected=false;try{auto wrong=after->repartition(expected,activeDestination,next.offset(rank));}catch(const std::exception&){activeRejected=true;}
 require(activeRejected,"Repartition during an unfinished field step accepted");
 if(size>1){ // Deliberately inconsistent old shared node must never be averaged.
  EBFieldGrid bad(geometry(12)),destination(geometry(next.cells(rank)));MovingWindowFields guard(bad.geometry(),total,oldOffset,MPI_COMM_WORLD,parameters(),ref);guard.installBoundary(bad);fill(bad,oldOffset);histories(bad,oldOffset,true);
  if(rank==0)bad.ex()(0,0,bad.geometry().nz)+=1.;bool failed=false;try{auto wrong=guard.repartition(bad,destination,next.offset(rank));}catch(const std::exception&){failed=true;}require(failed,"Inconsistent duplicated field node accepted");
  fill(bad,oldOffset);if(rank==0){auto& v=dynamic_cast<EBConvolutionalPML&>(bad.boundary()).checkpointHistory(4);const auto s=history(4,12);v[(s.layers-1)*s.row]+=1.;}
  failed=false;try{auto wrong=guard.repartition(bad,destination,next.offset(rank));}catch(const std::exception&){failed=true;}require(failed,"Inconsistent duplicated CPML node accepted");
 }
}
}
int main(int argc,char** argv){MPI_Init(&argc,&argv);int rank=0,size=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
 try{planner();nativeHalo(rank,size);routedCurrent(rank,size,Shape::LinearBinomial3);routedCurrent(rank,size,Shape::LinearNotch);routedCurrent(rank,size,Shape::QuinticCompensated);repartition(rank,size);
  double local[3]={plannerError,currentError,continuationError},global[3]={};MPI_Allreduce(local,global,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  unsigned long long counts[2]={haloValues,stateValues},totals[2]={};MPI_Allreduce(counts,totals,2,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  if(rank==0)std::cout<<std::setprecision(18)<<"{\"passed\":true,\"physical_acceptance\":false,\"ranks\":"<<size<<",\"planner_oracle_cases_per_rank\":"<<plannerCases<<",\"planner_invalid_controls_per_rank\":"<<rejections<<",\"planner_max_relative_optimality_error\":"<<global[0]<<",\"native_halo_values_bitwise\":"<<totals[0]<<",\"native_field_and_CPML_values_bitwise\":"<<totals[1]<<",\"thin_slab_current_max_relative_error\":"<<global[1]<<",\"continued_field_CPML_max_relative_error\":"<<global[2]<<",\"three_wide_shapes\":true,\"inconsistent_shared_node_negative_controls\":"<<(size>1?2:0)<<"}\n";
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}MPI_Finalize();return 0;
}
