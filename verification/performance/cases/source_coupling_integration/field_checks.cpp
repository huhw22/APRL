#include "source_coupling.h"
#include "moving_window_fields.h"
#include "discrete_field_ledger.h"
#include "eb_mpi.h"
#include "runtime_util.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
using namespace aprl;
Double reference(void*,unsigned a,std::size_t i,std::size_t j,std::ptrdiff_t k,std::size_t n){return(a<3?1e5:3e-4)*std::sin(.17*i+.23*j+.11*k+.07*a+.03*n);}
YeeComponent& component(EBFieldGrid& f,unsigned a){switch(a){case 0:return f.ex();case 1:return f.ey();case 2:return f.ez();case 3:return f.bx();case 4:return f.by();case 5:return f.bz();case 6:return f.jx();case 7:return f.jy();default:return f.jz();}}
int main(int argc,char** argv){MPI_Init(&argc,&argv);int rank=0,size=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
 try{if(argc!=2)throw std::runtime_error("Need test directory");const std::string out=argv[1];createDirectories(out,MPI_COMM_WORLD);
  const Double dz=1.2543262415359034e-6;EBGridGeometry g(16,16,12,5e-6,5e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);EBFieldGrid a(g),b(g);EBZSlabHaloExchange halo(MPI_COMM_WORLD);EBCPMLParameters cp;cp.cells[0]=cp.cells[1]=2;cp.cells[2]=0;MovingWindowReference ref;ref.value=reference;
  MovingWindowFields wa(g,g.nz*size,g.nz*rank,MPI_COMM_WORLD,cp,ref),wb(g,g.nz*size,g.nz*rank,MPI_COMM_WORLD,cp,ref);wa.installBoundary(a);wb.installBoundary(b);
  for(unsigned q=0;q<6;++q){auto& x=component(a,q);auto& y=component(b,q);for(std::size_t k=0;k<x.nz();++k)for(std::size_t j=0;j<x.ny();++j)for(std::size_t i=0;i<x.nx();++i)x(i,j,k)=y(i,j,k)=reference(nullptr,q,i,j,k+g.nz*rank,0);}
  ConservativeSourceCoupling coupling(MPI_COMM_WORLD,out);DiscreteFieldLedger ledger(out,g,g.nz*size,g.nz*rank,2,2,MPI_COMM_WORLD);wb.setDiscreteObserver(&DiscreteFieldLedger::observe,&ledger);Double maximum=0.;
  for(unsigned n=0;n<12;++n){const std::ptrdiff_t offset=-std::ptrdiff_t(n);ledger.begin(b,n);a.clearCurrent();b.clearCurrent();wa.beginStep(a,n,offset);wb.beginStep(b,n,offset);halo.advanceMagnetic(a,wa.physicalElectricHalos(a,n));halo.advanceMagnetic(b,wb.physicalElectricHalos(b,n));FieldVector<Double> origin(0.);origin[2]=(offset+static_cast<std::ptrdiff_t>(g.nz*rank))*g.dz;coupling.begin(b,origin,{});wb.previewElectricBoundary(true);halo.advanceElectric(b);wb.previewElectricBoundary(false);
   for(unsigned q=6;q<9;++q){auto& x=component(a,q);auto& y=component(b,q);for(std::size_t k=0;k<x.nz();++k)for(std::size_t j=0;j<x.ny();++j)for(std::size_t i=0;i<x.nx();++i)x(i,j,k)=y(i,j,k)=1e3*std::cos(.19*i+.13*j+.07*(offset+static_cast<std::ptrdiff_t>(k+g.nz*rank))+.3*q+.17*n);}
   halo.advanceElectric(a);wb.undoPreviewElectricBoundary(b);halo.sumSharedTangentialCurrent(b);coupling.applyCurrentSource(-1.);wb.commitPreparedElectricBoundary(b);ledger.finishElectric(b);
   for(unsigned q=0;q<6;++q){const auto& x=component(a,q);const auto& y=component(b,q);Double scale=0.,error=0.;for(std::size_t k=0;k<x.size();++k){scale=std::max(scale,std::abs(x.data()[k]));error=std::max(error,std::abs(x.data()[k]-y.data()[k]));}maximum=std::max(maximum,error/std::max(scale,1e-300));}
   wa.shift(a,n+1,offset-1);wb.shift(b,n+1,offset-1);ledger.afterRemap(b);
  }
  coupling.close(12);ledger.close(true);Double global=0.;MPI_Allreduce(&maximum,&global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);if(global>1024*std::numeric_limits<Double>::epsilon())throw std::runtime_error("Prepared electric update changed original Cowan/CPML fields");if(rank==0)std::cout<<std::setprecision(18)<<"{\"passed\":true,\"ranks\":"<<size<<",\"steps\":12,\"maximum_scaled_field_difference\":"<<global<<"}\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}MPI_Finalize();return 0;
}
