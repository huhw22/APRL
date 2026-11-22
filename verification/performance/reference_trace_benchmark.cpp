// Bounded initialization-only timing tool; never produces simulation fields.
#include "open_transverse_trace.h"
#include "compute_threads.h"
#include <mpi.h>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace aprl;
using aprl::verification::OpenTransverseTrace;
int main(int argc,char** argv) {
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int rank=0,ranks=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  try {
    if(argc!=5)throw std::runtime_error("Usage: trace_benchmark SOURCES OUTPUT_PREFIX NZ GAMMA");
    configureComputeThreads(provided);
    const EBGridGeometry g(144,96,std::stoul(argv[3]),4.e-5,4.e-5,4.e-6,4.e-6/SI::c,EBMaxwellSolver::CowanZ);
    const double gamma=std::stod(argv[4]);
    OpenTransverseTrace::Sources source;
    if(rank==0) {
      std::ifstream in(argv[1],std::ios::binary);std::uint64_t n=0;
      in.read(reinterpret_cast<char*>(&n),sizeof(n));
      if(!in || n>100000000)throw std::runtime_error("Bad source file");
      for(auto* v:{&source.x,&source.y,&source.z,&source.q}) {
        v->resize(n);in.read(reinterpret_cast<char*>(v->data()),n*sizeof(double));
      }
      if(!in || in.peek()!=EOF)throw std::runtime_error("Wrong source file length");
    }
    OpenTransverseTrace trace;
    if(rank==0)trace.configure(g);
    MPI_Barrier(MPI_COMM_WORLD);
    const auto begin=std::chrono::steady_clock::now();
#ifdef APRL_BENCHMARK_LEGACY
    if(ranks!=1)throw std::runtime_error("Legacy benchmark requires one rank");
    trace.build(source,g,gamma);
#else
    trace.build(source,g,gamma,MPI_COMM_WORLD);
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    if(rank==0) {
      std::ofstream binary(std::string(argv[2])+".bin",std::ios::binary);
      binary.write(reinterpret_cast<const char*>(trace.values.data()),trace.values.size()*sizeof(double));
      if(!binary)throw std::runtime_error("Cannot write trace");
      std::ofstream report(std::string(argv[2])+".json");
      report<<std::setprecision(17)<<"{\"ranks\":"<<ranks<<",\"threads_per_rank\":"<<computeThreadCount()
            <<",\"source_nodes\":"<<source.q.size()<<",\"target_cells\":[144,96,"<<g.nz
            <<"],\"wall_seconds\":"<<elapsed<<",\"direct_checks\":"<<trace.directChecks
            <<",\"max_relative_interpolation_error\":"<<trace.maxRelativeError;
#ifndef APRL_BENCHMARK_LEGACY
      report<<",\"coarse_points\":"<<trace.coarsePoints<<",\"fallback_tiles\":"<<trace.fallbackTiles;
#endif
      report<<"}\n";
      std::cout<<"wall_seconds="<<std::setprecision(10)<<elapsed<<'\n';
    }
    MPI_Finalize();return 0;
  }catch(const std::exception& e){std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
  return 1;
}
