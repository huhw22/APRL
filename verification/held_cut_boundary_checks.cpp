// Operator-only indexing/causality checks. The scientific probe retains the
// full production transverse mesh; this small array fixture tests the adapter.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <mpi.h>

#include "held_cut_boundary.h"

namespace {
using namespace aprl;
using namespace aprl::verification;
constexpr std::size_t NX=24,NY=26,W=32,FULL=112,START=64;
constexpr Double DZ=1.2543262415359034e-6;
std::uint64_t checks=0,fieldValues=0,historyValues=0,negativeDifferences=0;
std::uint64_t skippedShiftDifferences=0,missingHistoryDifferences=0;
Double maximumFieldDifference=0.,maximumHistoryDifference=0.;

void check(bool ok,const std::string& message) {
  ++checks;if(!ok)throw std::runtime_error(message);
}
bool equal(Double a,Double b) {return std::memcmp(&a,&b,sizeof(Double))==0;}
YeeComponent& c(EBFieldGrid& f,unsigned int a) {
  switch(a) {
    case 0:return f.ex();case 1:return f.ey();case 2:return f.ez();
    case 3:return f.bx();case 4:return f.by();case 5:return f.bz();
  }
  throw std::out_of_range("Bad field");
}
const YeeComponent& c(const EBFieldGrid& f,unsigned int a) {
  return c(const_cast<EBFieldGrid&>(f),a);
}
EBGridGeometry geometry(std::size_t nz) {
  return EBGridGeometry(NX,NY,nz,5e-6,5e-6,DZ,DZ/SI::c,EBMaxwellSolver::CowanZ);
}
EBCPMLParameters parameters() {
  EBCPMLParameters p;p.cells[0]=p.cells[1]=10;p.cells[2]=0;
  p.polynomialOrder=3.;p.targetReflection=1e-8;p.kappaMax=8.;p.alphaFraction=0.;
  return p;
}
void seed(EBFieldGrid& f) {
  for(unsigned int a=0;a<6;++a) {
    auto& v=c(f,a);
    for(std::size_t k=0;k<v.nz();++k)for(std::size_t j=0;j<v.ny();++j)
      for(std::size_t i=0;i<v.nx();++i) {
        const int tag=static_cast<int>((13*i+17*j+19*k+23*a)%101)-50;
        v(i,j,k)=std::ldexp(static_cast<Double>(tag),a<3?-4:-24);
      }
  }
  PerfectElectricConductorBoundary transverse(false,false);
  transverse.afterElectricUpdate(f);
  auto& p=dynamic_cast<EBConvolutionalPML&>(f.boundary());
  for(std::size_t h=0;h<12;++h) {
    auto& a=p.checkpointHistory(h);
    for(std::size_t i=0;i<a.size();++i)
      a[i]=std::ldexp(static_cast<Double>(static_cast<int>((i+3*h)%31)-15),
                     h<6?-8:-28);
  }
}
struct Donor {EBFieldGrid* fields;unsigned int state;};
Double query(void* raw,unsigned int component,std::size_t i,std::size_t j,
             std::ptrdiff_t k,unsigned int state) {
  auto& d=*static_cast<Donor*>(raw);
  if(state!=d.state)throw std::logic_error("Reference native time mismatch");
  const auto& f=c(*d.fields,component);
  if(k<0||static_cast<std::size_t>(k)>=f.nz())
    throw std::out_of_range("Reference native position mismatch");
  return f(i,j,static_cast<std::size_t>(k));
}
struct Fill {Donor* donor;std::size_t* offset;};
Double incoming(void* raw,unsigned int component,std::size_t i,std::size_t j,
                std::size_t k,Double) {
  auto& f=*static_cast<Fill*>(raw);
  if(component>=6)return 0.;
  return query(f.donor,component,i,j,static_cast<std::ptrdiff_t>(*f.offset+k),
               f.donor->state);
}
void poison(EBFieldGrid& f,unsigned int n) {
  // The integer-z terminal planes are not a usable physical boundary. The
  // half-z components end one half cell earlier and remain trusted initially.
  for(unsigned int a:{0u,1u,5u}) {
    auto& v=c(f,a);
    for(std::size_t j=0;j<v.ny();++j)for(std::size_t i=0;i<v.nx();++i)
      v(i,j,W)=static_cast<Double>((a+1)*(n+1))*(a<3?1e7:1e-1)+i+3*j;
  }
}
std::uint64_t compare(const EBFieldGrid& shortGrid,const MovingWindowRemap& shortLayout,
    const EBFieldGrid& fullGrid,const MovingWindowRemap& fullLayout,
    std::size_t offset,bool requireEqual) {
  std::uint64_t differences=0;
  for(unsigned int a=0;a<6;++a) {
    const auto& s=c(shortGrid,a);const auto& f=c(fullGrid,a);
    for(std::size_t k=0;k<W;++k)for(std::size_t j=0;j<s.ny();++j)
      for(std::size_t i=0;i<s.nx();++i) {
        const Double x=s(i,j,k),y=f(i,j,k+offset);
        ++fieldValues;
        if(!equal(x,y))++differences;
        if(requireEqual)maximumFieldDifference=std::max(maximumFieldDifference,std::abs(x-y));
      }
  }
  const auto& p=dynamic_cast<const EBConvolutionalPML&>(shortGrid.boundary());
  const auto& q=dynamic_cast<const EBConvolutionalPML&>(fullGrid.boundary());
  for(std::size_t h=0;h<12;++h) {
    const auto& s=shortLayout.historyLayout(h);const auto& f=fullLayout.historyLayout(h);
    if(s.derivativeAxis==2) {
      check(p.checkpointHistory(h).empty()&&q.checkpointHistory(h).empty(),"Longitudinal CPML is not empty");
      continue;
    }
    const std::size_t width=s.dimensions[s.derivativeAxis==0?1:0];
    const auto& ps=p.checkpointHistory(h);const auto& pf=q.checkpointHistory(h);
    for(std::size_t slot=0;slot<s.activeAxisIndices.size();++slot)
      for(std::size_t k=0;k<W;++k)for(std::size_t ij=0;ij<width;++ij) {
        const Double x=ps[(slot*s.dimensions[2]+k)*width+ij];
        const Double y=pf[(slot*f.dimensions[2]+offset+k)*width+ij];
        ++historyValues;
        if(!equal(x,y))++differences;
        if(requireEqual)maximumHistoryDifference=std::max(maximumHistoryDifference,std::abs(x-y));
      }
  }
  if(requireEqual)check(differences==0,"Right cut failed bitwise full-domain causality");
  return differences;
}

void run(bool skipShift,bool importHistory) {
  auto p=parameters();EBFieldGrid full(geometry(FULL)),shortGrid(geometry(W));
  full.setBoundary(std::unique_ptr<EBBoundaryOperator>(new EBConvolutionalPML(
      full.geometry(),FULL,0,p,false,false)));
  seed(full);Donor donor{&full,0};std::size_t offset=START;
  HeldCutReference reference;reference.context=&donor;reference.field=query;
  auto* boundary=new HeldCutBoundary(shortGrid.geometry(),p,offset,reference,1e9);
  shortGrid.setBoundary(std::unique_ptr<EBBoundaryOperator>(boundary));
  MovingWindowRemap fullLayout(full.geometry(),p,{-1});
  MovingWindowRemap shortLayout(shortGrid.geometry(),p,{-1});
  copyHeldCutFields(shortGrid,full,offset);
  copyHeldCutHistories(shortGrid,shortLayout,full,fullLayout,offset);
  compare(shortGrid,shortLayout,full,fullLayout,offset,true);
  Fill fill{&donor,&offset};MovingWindowIncoming in;in.context=&fill;in.value=incoming;
  for(unsigned int n=0;n<12;++n) {
    poison(shortGrid,n);boundary->setRightHaloValue((n+1)*1e9);
    boundary->beginStep(shortGrid,n);
    const auto halo=boundary->electricHalos(shortGrid,n);
    full.advanceMagnetic();full.advanceElectric();++donor.state;
    shortGrid.advanceMagnetic(halo);shortGrid.advanceElectric();
    const bool shift=!(skipShift&&n==2);
    if(shift) {
      --offset;
      shortLayout.shift(shortGrid,-1,offset*DZ,in,MovingWindowCurrentPolicy::Reset);
      if(importHistory)
        copyHeldCutHistories(shortGrid,shortLayout,full,fullLayout,offset,0,1);
    }
    const bool pass=!skipShift&&importHistory;
    const std::uint64_t d=compare(shortGrid,shortLayout,full,fullLayout,offset,pass);
    if(!pass&&d) {
      negativeDifferences+=d;
      if(skipShift)skippedShiftDifferences=d;
      else missingHistoryDifferences=d;
      return;
    }
  }
  if(skipShift||!importHistory)
    throw std::runtime_error("Negative control did not reveal missing causality/state");
  check(boundary->referenceValues()>0,"No held field queries recorded");
}

void rejects() {
  auto p=parameters();EBFieldGrid full(geometry(FULL)),small(geometry(W));
  full.setBoundary(std::unique_ptr<EBBoundaryOperator>(new EBConvolutionalPML(
      full.geometry(),FULL,0,p,false,false)));
  Donor donor{&full,0};HeldCutReference r;r.context=&donor;r.field=query;
  std::size_t offset=START;HeldCutBoundary b(small.geometry(),p,offset,r);
  bool rejected=false;try{b.electricHalos(small,0);}catch(const std::logic_error&){rejected=true;}
  check(rejected,"Unprepared halo request accepted");
  b.beginStep(small,0);rejected=false;
  try{b.beginStep(small,1);}catch(const std::logic_error&){rejected=true;}
  check(rejected,"Overlapping step accepted");
  ++offset;rejected=false;
  try{b.electricHalos(small,0);}catch(const std::logic_error&){rejected=true;}
  check(rejected,"Mid-step remap accepted");
  p.cells[2]=2;rejected=false;
  try{HeldCutBoundary invalid(small.geometry(),p,offset,r);}
  catch(const std::invalid_argument&){rejected=true;}
  check(rejected,"Longitudinal CPML accepted");
  rejected=false;
  try{copyHeldCutFields(small,full,FULL);}
  catch(const std::out_of_range&){rejected=true;}
  check(rejected,"Out-of-range donor copy accepted");
}
}

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);int rank=0,size=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  try {
    check(size==1,"Held cut operator checks require one rank");
    run(false,true);run(true,true);run(false,false);rejects();
    std::ostream* out=&std::cout;std::ofstream file;
    if(argc==2) {file.open(argv[1]);check(static_cast<bool>(file),"Cannot open report");out=&file;}
    *out<<"{\"format\":\"aprl-held-cut-operator-checks-v1\",\"pass\":true,"
        <<"\"checks\":"<<checks<<",\"field_values_checked\":"<<fieldValues
        <<",\"history_values_checked\":"<<historyValues
        <<",\"negative_control_differences\":"<<negativeDifferences
        <<",\"skipped_shift_differences\":"<<skippedShiftDifferences
        <<",\"missing_incoming_history_differences\":"<<missingHistoryDifferences
        <<",\"right_halo_and_terminal_faces_poisoned\":true"
        <<",\"maximum_valid_field_difference\":"<<maximumFieldDifference
        <<",\"maximum_valid_history_difference\":"<<maximumHistoryDifference
        <<",\"scientific_geometry\":false,\"scope\":\"native Cowan+transverse CPML field dependency, poisoned right halo and terminal integer faces, skipped-shift and missing-side-history negative controls; no particle or production acceptance\"}\n";
    MPI_Finalize();return 0;
  }catch(const std::exception& e) {
    if(rank==0)std::cerr<<"held cut checks FAILED: "<<e.what()<<"\n";
    MPI_Finalize();return 1;
  }
}
