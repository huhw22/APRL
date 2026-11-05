#ifndef APRL_PARTICLE_FORCE_PROBE_H
#define APRL_PARTICLE_FORCE_PROBE_H

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>
#include <unistd.h>
#include <vector>
#ifdef APRL_HAS_OPENMP
#include <omp.h>
#endif
#include "eb_particles.h"

namespace aprl {

// Read-only selected-particle trace of the actual converged active intervals.
// No grid is sampled here. Analytic device fields are reevaluated at the same
// midpoint events as the already completed pusher, then checked against its
// independently saved device work. No callback, MPI, HDF5 or state update occurs.
class ParticleForceProbe {
public:
  struct Header {
    char magic[8];
    std::uint32_t version, recordBytes, headerBytes, endianMarker;
    std::uint64_t rank, thread, selectedIds, reserved0, reserved1;
  };
  struct Record {
    std::uint64_t step, sourceId, runtimeId;
    std::uint32_t substeps, flags;
    Double time, dt, charge, mass, weight;
    Double positionBefore[3], positionAfter[3];
    Double properBefore[3], properAfter[3];
    Double gridElectric[3], gridMagnetic[3];
    Double gammaBefore, gammaAfter;
    Double prescribedImpulseElectric[3], prescribedImpulseMagnetic[3];
    Double prescribedDeltaGamma;
  };
  static_assert(sizeof(Double)==8,"Force probe requires float64 Double");
  static_assert(sizeof(Header)==64,"Force probe header layout changed");
  static_assert(offsetof(Record,time)==32,"Force probe double payload moved");
  static_assert(sizeof(Record)==288,"Force probe record layout changed");

  static bool enabled() {return !settings().directory.empty();}
  static bool selected(std::uint64_t sourceId) {
    const Settings& s=settings();
    return !s.directory.empty() && s.ids.find(sourceId)!=s.ids.end();
  }

  // Pure record construction is public for numerical diagnostic tests. All
  // endpoints are actual on-shell pusher states, not interpolated plane data.
  static Record makeRecord(std::uint64_t step,
      const RelativisticParticleSI& before,const RelativisticParticleSI& after,
      const FieldVector<Double>& gridElectric,const FieldVector<Double>& gridMagnetic,
      Double time,Double dt,unsigned substeps,
      const std::vector<RelativisticParticleSI>& states,
      const std::vector<BorisWork>& works,const SIFieldSourceSet& sources,
      const BoostFrameTransform& frame) {
    if(!substeps || states.size()!=static_cast<std::size_t>(substeps)+1 ||
       works.size()!=substeps || !(dt>0.) || !std::isfinite(dt) ||
       !(before.mass>0.) || before.sourceId!=after.sourceId)
      throw std::invalid_argument("Invalid force-probe interval history");
    for(unsigned a=0;a<3;++a)
      if(states.front().position[a]!=before.position[a] ||
         states.front().properVelocity[a]!=before.properVelocity[a] ||
         states.back().position[a]!=after.position[a] ||
         states.back().properVelocity[a]!=after.properVelocity[a])
        throw std::invalid_argument("Force-probe endpoints differ from converged history");
    Record r{};r.step=step;r.sourceId=before.sourceId;r.runtimeId=before.id;
    r.substeps=substeps;r.flags=1; // Bit 0: independent device-work check passed.
    r.time=time;r.dt=dt;r.charge=before.charge;r.mass=before.mass;r.weight=before.weight;
    r.gammaBefore=BoostFrameTransform::gammaFromProperVelocity(before.properVelocity);
    r.gammaAfter=BoostFrameTransform::gammaFromProperVelocity(after.properVelocity);
    for(unsigned a=0;a<3;++a) {
      r.positionBefore[a]=before.position[a];r.positionAfter[a]=after.position[a];
      r.properBefore[a]=before.properVelocity[a];r.properAfter[a]=after.properVelocity[a];
      r.gridElectric[a]=gridElectric[a];r.gridMagnetic[a]=gridMagnetic[a];
    }
    long double electric[3]={},magnetic[3]={},work=0.,savedWork=0.,absoluteWork=0.;
    const Double h=dt/substeps;
    const long double rest=static_cast<long double>(before.mass)*SI::c*SI::c;
    for(unsigned n=0;n<substeps;++n) {
      const auto& first=states[n];const auto& last=states[n+1];
      if(first.sourceId!=before.sourceId || last.sourceId!=before.sourceId ||
         first.charge!=before.charge || first.mass!=before.mass ||
         last.charge!=before.charge || last.mass!=before.mass)
        throw std::invalid_argument("Force-probe particle identity changed inside history");
      const Double g0=BoostFrameTransform::gammaFromProperVelocity(first.properVelocity);
      const Double g1=BoostFrameTransform::gammaFromProperVelocity(last.properVelocity);
      FieldVector<Double> point(first.position);
      point.pmv(.5*h*SI::c/g0,first.properVelocity);
      FieldVector<Double> ee(0.),bb(0.);
      sources.addPrescribedBox(point,time+(n+.5)*h,frame,ee,bb);
      long double velocity[3];
      for(unsigned a=0;a<3;++a)
        velocity[a]=(static_cast<long double>(first.properVelocity[a])+last.properVelocity[a])/
          (static_cast<long double>(g0)+g1);
      for(unsigned a=0;a<3;++a) {
        const unsigned b=(a+1)%3,c=(a+2)%3;
        electric[a]+=static_cast<long double>(before.charge)*h*ee[a]/(before.mass*SI::c);
        magnetic[a]+=static_cast<long double>(before.charge)*h/before.mass*
          (velocity[b]*bb[c]-velocity[c]*bb[b]);
        work+=static_cast<long double>(before.charge)*ee[a]*
          (static_cast<long double>(last.position[a])-first.position[a]);
      }
      savedWork+=works[n].prescribedWork;
      absoluteWork+=std::abs(works[n].prescribedWork);
    }
    const long double scale=absoluteWork+std::abs(work)+std::abs(savedWork);
    if(std::abs(work-savedWork)>256*std::numeric_limits<Double>::epsilon()*
       std::max(scale,std::numeric_limits<long double>::min()))
      throw std::runtime_error("Force-probe analytic device resampling differs from actual work");
    for(unsigned a=0;a<3;++a) {
      r.prescribedImpulseElectric[a]=static_cast<Double>(electric[a]);
      r.prescribedImpulseMagnetic[a]=static_cast<Double>(magnetic[a]);
    }
    r.prescribedDeltaGamma=static_cast<Double>(work/rest);
    Double values[32];std::memcpy(values,reinterpret_cast<const char*>(&r)+32,sizeof(values));
    for(Double value:values)if(!std::isfinite(value))
      throw std::runtime_error("Nonfinite force-probe record");
    return r;
  }

  static void record(std::uint64_t step,int rank,
      const RelativisticParticleSI& before,const RelativisticParticleSI& after,
      const FieldVector<Double>& gridElectric,const FieldVector<Double>& gridMagnetic,
      Double time,Double dt,unsigned substeps,
      const std::vector<RelativisticParticleSI>& states,
      const std::vector<BorisWork>& works,const SIFieldSourceSet& sources,
      const BoostFrameTransform& frame) {
    if(!selected(before.sourceId))return;
    const Record r=makeRecord(step,before,after,gridElectric,gridMagnetic,
      time,dt,substeps,states,works,sources,frame);
    Registry& registry=writers();
    if(registry.closed)throw std::logic_error("Force probe already closed");
    static thread_local Writer* writer=nullptr;
    if(!writer) {
      int thread=0;
#ifdef APRL_HAS_OPENMP
      thread=omp_get_thread_num();
#endif
      std::lock_guard<std::mutex> lock(registry.mutex);
      registry.files.emplace_back(new Writer(rank,thread,settings()));
      writer=registry.files.back().get();
    }
    if(writer->rank!=rank)throw std::logic_error("Force-probe writer rank changed");
    writeAll(writer->fd,reinterpret_cast<const char*>(&r),sizeof(r));
    const auto* data=reinterpret_cast<const unsigned char*>(&r);
    for(std::size_t i=0;i<sizeof(r);++i){writer->checksum^=data[i];writer->checksum*=1099511628211ULL;}
    if(!writer->count)writer->firstStep=r.step;
    writer->lastStep=r.step;++writer->count;
  }

  // Call once on each rank, on its main thread after all preparation workers
  // have joined. This closes every local fd; there is deliberately no MPI.
  static void close(int rank,std::uint64_t step,bool complete) {
    if(!enabled())return;
    Registry& registry=writers();
    if(registry.closed)throw std::logic_error("Force probe close repeated");
    std::ostringstream report;report<<std::setprecision(17)
      <<"{\n  \"schema\": \"APRL_particle_force_interval_v1\",\n"
      <<"  \"frame\": \"boost\",\n  \"rank\": "<<rank<<",\n"
      <<"  \"complete\": "<<(complete?"true":"false")<<",\n  \"step\": "<<step<<",\n"
      <<"  \"record_bytes\": "<<sizeof(Record)<<",\n  \"header_bytes\": "<<sizeof(Header)<<",\n"
      <<"  \"endianness_marker\": 16909060,\n  \"ids_file\": "<<quoted(settings().idsFile)<<",\n"
      <<"  \"selected_ids\": [";
    for(std::size_t i=0;i<settings().orderedIds.size();++i)
      report<<(i?", ":"")<<settings().orderedIds[i];
    report<<"],\n  \"files\": [";std::uint64_t count=0;
    for(std::size_t i=0;i<registry.files.size();++i) {
      Writer& f=*registry.files[i];
      if(f.rank!=rank)throw std::logic_error("Force-probe close rank differs");
      struct stat stat{};
      if(::fstat(f.fd,&stat)!=0 || stat.st_size!=static_cast<off_t>(sizeof(Header)+f.count*sizeof(Record)))
        throw std::runtime_error("Force-probe trace length differs");
      if(::close(f.fd)!=0)throw std::runtime_error("Cannot close force-probe trace");
      f.fd=-1;count+=f.count;
      report<<(i?",\n    ":"\n    ")<<"{\"path\": "<<quoted(f.path)
        <<", \"thread\": "<<f.thread<<", \"records\": "<<f.count
        <<", \"first_step\": "<<f.firstStep<<", \"last_step\": "<<f.lastStep
        <<", \"records_fnv1a64\": \""<<std::hex<<f.checksum<<std::dec<<"\"}";
    }
    report<<"\n  ],\n  \"records\": "<<count<<",\n"
      <<"  \"scope\": \"Full converged active pusher intervals; physical terminal clipping remains external.\",\n"
      <<"  \"units\": \"SI x/time/E/B/q/m; u=gamma*v/c; impulses are delta-u; prescribedDeltaGamma is dimensionless; mass and charge include macro weight.\"\n}\n";
    const std::string stem=rankStem(rank);
    writeText(stem+".json",report.str());
    std::ostringstream marker;marker<<(complete?1:0)<<' '<<step<<' '<<count<<'\n';
    writeText(stem+".complete",marker.str());registry.closed=true;
  }

private:
  struct Settings {
    std::string directory,idsFile;
    std::unordered_set<std::uint64_t> ids;
    std::vector<std::uint64_t> orderedIds;
    Settings() {
      const char* d=std::getenv("APRL_PARTICLE_FORCE_PROBE");
      if(!d||!*d)return;directory=d;
      const char* p=std::getenv("APRL_PARTICLE_FORCE_PROBE_IDS");
      if(!p||!*p)throw std::invalid_argument("Force probe requires an explicit source-ID file");
      idsFile=p;std::ifstream input(idsFile);
      if(!input)throw std::runtime_error("Cannot open force-probe source-ID file");
      std::string line;
      while(std::getline(input,line)) {
        const auto comment=line.find('#');if(comment!=std::string::npos)line.resize(comment);
        std::istringstream words(line);std::string word;
        while(words>>word) {
          char* end=nullptr;errno=0;const auto id=std::strtoull(word.c_str(),&end,10);
          if(word.empty()||word[0]<'0'||word[0]>'9'||!end||*end||errno==ERANGE||!id||!ids.insert(id).second)
            throw std::invalid_argument("Invalid, zero or duplicate force-probe source ID");
          orderedIds.push_back(id);
        }
      }
      if(!input.eof()||orderedIds.empty())throw std::invalid_argument("Empty or unreadable force-probe source-ID file");
      std::sort(orderedIds.begin(),orderedIds.end());
      std::unordered_map<std::uint64_t,std::uint32_t> families;
      for(auto id:orderedIds)families[(id-1)/16]|=1u<<((id-1)%16);
      for(const auto& family:families)if(family.second!=65535u)
        throw std::invalid_argument("Force-probe selection must contain complete 16-particle phase families");
    }
  };
  struct Writer {
    int fd=-1,rank,thread;std::string path;
    std::uint64_t count=0,firstStep=0,lastStep=0,checksum=14695981039346656037ULL;
    Writer(int r,int t,const Settings& s):rank(r),thread(t) {
      char name[96];std::snprintf(name,sizeof(name),"/force-rank-%05d-thread-%03d.bin",rank,thread);path=s.directory+name;
      fd=::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);
      if(fd<0)throw std::runtime_error("Cannot create force-probe trace");
      Header h{};const char magic[8]={'A','P','R','L','F','P','1','\0'};std::memcpy(h.magic,magic,8);
      h.version=1;h.recordBytes=sizeof(Record);h.headerBytes=sizeof(Header);h.endianMarker=0x01020304;
      h.rank=rank;h.thread=thread;h.selectedIds=s.orderedIds.size();
      try{writeAll(fd,reinterpret_cast<const char*>(&h),sizeof(h));}
      catch(...){::close(fd);fd=-1;throw;}
    }
    ~Writer(){if(fd>=0)::close(fd);}
  };
  struct Registry {std::mutex mutex;std::vector<std::unique_ptr<Writer>> files;bool closed=false;};
  static const Settings& settings(){static const Settings value;return value;}
  static Registry& writers(){static Registry value;return value;}
  static void writeAll(int fd,const char* data,std::size_t bytes) {
    while(bytes) {const ssize_t n=::write(fd,data,bytes);if(n<0&&errno==EINTR)continue;
      if(n<=0)throw std::runtime_error("Cannot write force-probe output");data+=n;bytes-=static_cast<std::size_t>(n);}
  }
  static std::string rankStem(int rank) {
    char name[64];std::snprintf(name,sizeof(name),"/force-rank-%05d",rank);return settings().directory+name;
  }
  static void writeText(const std::string& path,const std::string& text) {
    const int fd=::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);
    if(fd<0)throw std::runtime_error("Cannot create force-probe completion output");
    try{writeAll(fd,text.data(),text.size());}catch(...){::close(fd);throw;}
    if(::close(fd)!=0)throw std::runtime_error("Cannot close force-probe completion output");
  }
  static std::string quoted(const std::string& value) {
    std::ostringstream text;text<<'"';
    for(unsigned char c:value) {
      if(c=='"'||c=='\\')text<<'\\'<<c;
      else if(c<32)text<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(c)<<std::dec;
      else text<<c;
    }text<<'"';return text.str();
  }
};
}
#endif
