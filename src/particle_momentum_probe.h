#ifndef APRL_PARTICLE_MOMENTUM_PROBE_H
#define APRL_PARTICLE_MOMENTUM_PROBE_H

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <unistd.h>
#ifdef APRL_HAS_OPENMP
#include <omp.h>
#endif
#include "eb_particles.h"

namespace aprl {
// Diagnostic-only binary trace of already sampled fields and actual states.
// Each OpenMP thread owns its fd; there are no MPI/HDF5 calls or buffered data.
// Enabling this does not resample fields or modify a state-producing expression.
class ParticleMomentumProbe {
public:
  struct Record {
    std::uint64_t step, sourceId;
    Double time, dt, charge;
    Double position[3], before[3], after[3], electric[3], magnetic[3];
  };
  static_assert(sizeof(Record)==160,"Momentum probe record layout changed");
  static bool selected(Double time) {
    const Settings& s=settings();
    if(s.directory.empty())return false;
    if(s.firstOnly)return true;
    const auto step=static_cast<std::uint64_t>(std::floor(time/s.dt+1.e-7));
    return step%s.stride==0;
  }
  static void record(const RelativisticParticleSI& before,
      const RelativisticParticleSI& after,const FieldVector<Double>& electric,
      const FieldVector<Double>& magnetic,Double time,Double dt) {
    const Settings& s=settings();
    if(s.firstOnly) {
      static thread_local std::unordered_set<std::uint64_t> seen;
      if(!seen.insert(before.sourceId).second)return;
    }
    Record r{};r.step=static_cast<std::uint64_t>(std::floor(time/s.dt+1.e-7));
    r.sourceId=before.sourceId;r.time=time;r.dt=dt;r.charge=before.charge;
    for(unsigned a=0;a<3;++a) {
      r.position[a]=before.position[a];r.before[a]=before.properVelocity[a];
      r.after[a]=after.properVelocity[a];r.electric[a]=electric[a];r.magnetic[a]=magnetic[a];
    }
    static thread_local int fd=openFile(s.directory);
    const char* data=reinterpret_cast<const char*>(&r);std::size_t left=sizeof(r);
    while(left) {
      const ssize_t n=::write(fd,data,left);
      if(n<0&&errno==EINTR)continue;
      if(n<=0)throw std::runtime_error("Momentum probe write failed");
      data+=n;left-=static_cast<std::size_t>(n);
    }
  }
private:
  struct Settings {
    std::string directory;Double dt=0.;std::uint64_t stride=20;bool firstOnly=false;
    Settings() {
      const char* d=std::getenv("APRL_PARTICLE_MOMENTUM_PROBE");
      if(!d||!*d)return;
      directory=d;const char* t=std::getenv("APRL_PARTICLE_MOMENTUM_PROBE_DT");
      const char* first=std::getenv("APRL_PARTICLE_MOMENTUM_PROBE_FIRST");
      firstOnly=first&&std::string(first)=="1";
      if(!t)throw std::runtime_error("Momentum probe requires nominal field dt");
      char* end=nullptr;dt=std::strtod(t,&end);
      if(!end||*end||!(dt>0.)||!std::isfinite(dt))throw std::runtime_error("Invalid momentum probe dt");
    }
  };
  static const Settings& settings() {static const Settings s;return s;}
  static int openFile(const std::string& directory) {
    const char* rank=std::getenv("SLURM_PROCID");
    if(!rank)rank=std::getenv("OMPI_COMM_WORLD_RANK");
    if(!rank)throw std::runtime_error("Momentum probe requires rank identity");
    int thread=0;
#ifdef APRL_HAS_OPENMP
    thread=omp_get_thread_num();
#endif
    char suffix[96];std::snprintf(suffix,sizeof(suffix),"/momentum-rank-%05d-thread-%03d.bin",std::atoi(rank),thread);
    const int fd=::open((directory+suffix).c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);
    if(fd<0)throw std::runtime_error("Cannot create momentum probe trace");
    return fd;
  }
};
}
#endif
