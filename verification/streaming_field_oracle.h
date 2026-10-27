#ifndef APRL_VERIFICATION_STREAMING_FIELD_ORACLE_H
#define APRL_VERIFICATION_STREAMING_FIELD_ORACLE_H

// Diagnostic-only native-state replay. Version 2 distinguishes physical PIC
// donors from independent prescribed-uniform-source donors. Both currently
// require a full-grid donor evolution and provide no production speed claim.
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>

#include "eb_field.h"

namespace aprl { namespace verification {

struct NativeOracleConfig {
  EBGridGeometry geometry;
  unsigned steps,guard,seed;
  unsigned sourceKind; // 0: physical PIC history; 1: prescribed uniform source.
  Double beta,charge,originZ,releaseZ,detectorZ,windowBeta;
  NativeOracleConfig(const EBGridGeometry& g,unsigned n,unsigned p,unsigned s,
      Double b,Double q,Double z,Double release,Double detector,Double window,
      unsigned sourceKindValue=0)
    :geometry(g),steps(n),guard(p),seed(s),sourceKind(sourceKindValue),beta(b),charge(q),originZ(z),
     releaseZ(release),detectorZ(detector),windowBeta(window) {}
};

namespace native_oracle_detail {
inline void require(bool ok,const std::string& message) {
  if(!ok)throw std::runtime_error("Native oracle: "+message);
}
inline void platform() {
  static_assert(sizeof(Double)==8 && std::numeric_limits<Double>::is_iec559,
                "Oracle requires IEEE binary64 fields");
  const std::uint16_t one=1;
  require(*reinterpret_cast<const unsigned char*>(&one)==1,"requires little-endian storage");
}
inline std::ofstream create(const std::string& path) {
  require(!std::ifstream(path).good(),"refusing overwrite: "+path);
  std::ofstream stream(path,std::ios::out|std::ios::binary);
  require(stream.good(),"cannot create "+path);return stream;
}
inline void read(std::ifstream& stream,void* destination,std::size_t bytes) {
  stream.read(static_cast<char*>(destination),static_cast<std::streamsize>(bytes));
  require(stream.good(),"truncated binary input");
}
inline void eof(std::ifstream& stream) {
  char extra=0;stream.read(&extra,1);
  require(stream.eof() && stream.gcount()==0,"unexpected trailing bytes");
}
inline void magic(std::ifstream& stream,const char* expected) {
  char bytes[8];read(stream,bytes,8);
  require(std::string(bytes,8)==std::string(expected,8),"wrong binary magic/version");
}
inline void finish(std::ofstream& stream) {
  stream.flush();require(stream.good(),"write/flush failed");
  stream.close();require(!stream.fail(),"close failed");
}
inline std::string framePath(const std::string& directory,unsigned step) {
  std::ostringstream path;path<<directory<<"/state-"<<std::setw(6)<<std::setfill('0')<<step<<".bin";
  return path.str();
}
inline const YeeComponent& component(const EBFieldGrid& f,unsigned a) {
  if(a==0)return f.ex();
  if(a==1)return f.ey();
  if(a==2)return f.ez();
  if(a==3)return f.bx();
  if(a==4)return f.by();
  return f.bz();
}
inline std::array<std::uint64_t,8> integers(const NativeOracleConfig& c) {
  return {{c.geometry.nx,c.geometry.ny,c.geometry.nz,c.steps,c.guard,c.seed,6,c.sourceKind}};
}
inline const char* sourceKindName(const NativeOracleConfig& c) {
  return c.sourceKind==1?"prescribed-uniform-initial-source":"physical-pic-history";
}
inline std::array<Double,10> reals(const NativeOracleConfig& c) {
  return {{c.geometry.dx,c.geometry.dy,c.geometry.dz,c.geometry.dt,c.beta,
            c.charge,c.originZ,c.releaseZ,c.detectorZ,c.windowBeta}};
}
inline bool vertexZ(unsigned component) {return component==0||component==1||component==5;}
inline std::size_t slabCells(const NativeOracleConfig& c) {return 16+2*c.guard+40+4;}
inline std::size_t firstK(const NativeOracleConfig& c,unsigned state) {
  // The bounded fixture anchors its release and origin on integer cells.
  // Round those anchors once, matching the probe's exact 218 / -84 constants.
  const long long release=std::llround(c.releaseZ/c.geometry.dz);
  const long long origin=std::llround(c.originZ/c.geometry.dz);
  const long long first=static_cast<long long>(std::floor(
    static_cast<Double>(release)-c.windowBeta*state))-origin-c.guard-20-2;
  require(first>=0,"negative recording slab origin");return static_cast<std::size_t>(first);
}
inline std::string jsonString(const std::string& value) {
  std::ostringstream out;out<<'"';
  for(unsigned char c:value) {
    if(c=='"'||c=='\\')out<<'\\'<<c;
    else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<unsigned(c)<<std::dec;
    else out<<c;
  }
  out<<'"';return out.str();
}
} // namespace native_oracle_detail

class NativeFieldOracleWriter {
public:
  NativeFieldOracleWriter(const std::string& directory,const NativeOracleConfig& config)
    :directory_(directory),config_(config) {
    using namespace native_oracle_detail;
    platform();require(config_.sourceKind<=1,"unknown donor source kind");
    require(mkdir(directory.c_str(),0755)==0,"recording directory must be new: "+directory);
    auto stream=create(directory+"/header.bin");stream.write("APRLOH02",8);
    const auto u=integers(config_);const auto d=reals(config_);
    stream.write(reinterpret_cast<const char*>(u.data()),sizeof(u));
    stream.write(reinterpret_cast<const char*>(d.data()),sizeof(d));native_oracle_detail::finish(stream);
  }

  void record(const EBFieldGrid& fields,unsigned stateStep,std::size_t firstK) {
    using namespace native_oracle_detail;
    require(stateStep==frames_ && stateStep<=config_.steps,"record state order changed");
    const auto& g=fields.geometry();
    require(g.nx==config_.geometry.nx&&g.ny==config_.geometry.ny&&g.nz==config_.geometry.nz,
            "recording is not the full mesh");
    const std::size_t cells=slabCells(config_);
    require(firstK==native_oracle_detail::firstK(config_,stateStep)&&
            firstK<=g.nz && cells<=g.nz-firstK,"record slab outside prescribed full-mesh location");
    auto stream=create(framePath(directory_,stateStep));stream.write("APRLOF02",8);
    const std::uint64_t header[4]={stateStep,firstK,cells,6};
    stream.write(reinterpret_cast<const char*>(header),sizeof(header));
    std::uint64_t bytes=8+sizeof(header);
    for(unsigned a=0;a<6;++a) {
      const auto& source=component(fields,a);
      const std::size_t planes=cells+(vertexZ(a)?1:0),plane=source.nx()*source.ny();
      require(firstK+planes<=source.nz(),"native component slab outside mesh");
      const std::uint64_t description[4]={a,source.nx(),source.ny(),planes};
      stream.write(reinterpret_cast<const char*>(description),sizeof(description));
      const Double* values=source.data()+firstK*plane;
      for(std::size_t i=0;i<plane*planes;++i)
        require(std::isfinite(values[i]),"nonfinite recorded field");
      const std::size_t payload=plane*planes*sizeof(Double);
      stream.write(reinterpret_cast<const char*>(values),payload);
      bytes+=sizeof(description)+payload;
    }
    native_oracle_detail::finish(stream);++frames_;bytes_+=bytes;
  }

  void finish() {
    using namespace native_oracle_detail;
    require(frames_==config_.steps+1,"incomplete recording");
    auto manifest=create(directory_+"/manifest.json");
    manifest<<"{\"format\":\"aprl-native-field-oracle-v2\",\"completed\":true,\"frames\":"<<frames_
      <<",\"frame_bytes\":"<<bytes_<<",\"byte_order\":\"little\",\"field_type\":\"IEEE binary64\""
      <<",\"source_kind\":"<<config_.sourceKind<<",\"source_provenance\":\""<<sourceKindName(config_)<<"\""
      <<",\"prescribed_source_contract\":\"source kind 1 uses only initial positions/momenta and uniform prescribed trajectories; no physical PIC future states\""
      <<",\"state_contract\":\"completed state n: E(n), B(n-1/2); native staggered components\""
      <<",\"header_format\":\"APRLOH02; uint64[8] nx,ny,full_nz,steps,guard,seed,6,source_kind; double[10] dx,dy,dz,dt,beta,charge,origin_z,release_z,detector_z,window_beta\""
      <<",\"frame_format\":\"APRLOF02; uint64[4] state,first_global_k,slab_cells,6; six components each uint64[4] id,nx,ny,nz then double[nz][ny][nx]\""
      <<",\"scope\":\"diagnostic donor requiring complete full-grid field evolution; not a production reference supply or speed optimization\"}\n";
    native_oracle_detail::finish(manifest);
    auto marker=create(directory_+"/complete.bin");marker.write("APRLOC02",8);
    marker.write(reinterpret_cast<const char*>(&frames_),sizeof(frames_));native_oracle_detail::finish(marker);
  }
  std::uint64_t frames() const{return frames_;}
  std::uint64_t bytes() const{return bytes_;}
private:
  std::string directory_;NativeOracleConfig config_;
  std::uint64_t frames_=0,bytes_=0;
};

class NativeFieldOracleReader {
  struct Component {std::size_t nx=0,ny=0,nz=0;std::vector<Double> values;};
  struct Frame {
    bool valid=false;unsigned state=0;std::size_t first=0;
    std::uint64_t touched=0;std::array<Component,6> components;
  };
public:
  NativeFieldOracleReader(const std::string& directory,const NativeOracleConfig& expected)
    :directory_(directory),expected_(expected) {
    using namespace native_oracle_detail;
    platform();require(expected_.sourceKind<=1,"unknown expected donor source kind");
    std::ifstream marker(directory+"/complete.bin",std::ios::binary);
    require(marker.good(),"no completed full-domain recording at "+directory);
    magic(marker,"APRLOC02");std::uint64_t frames=0;read(marker,&frames,sizeof(frames));eof(marker);
    require(frames==expected_.steps+1,"completed recording has wrong length");
    std::ifstream header(directory+"/header.bin",std::ios::binary);
    require(header.good(),"missing recording header");magic(header,"APRLOH02");
    std::array<std::uint64_t,8> u;std::array<Double,10> d;
    read(header,u.data(),sizeof(u));read(header,d.data(),sizeof(d));eof(header);
    require(u[7]==expected_.sourceKind,"recording donor source kind mismatch");
    require(u==integers(expected_)&&d==reals(expected_),"recording configuration mismatch");
  }

  static Double query(void* context,unsigned component,std::size_t i,std::size_t j,
                      Double globalK,Double stateStep) {
    native_oracle_detail::require(context!=nullptr,"null callback context");
    return static_cast<NativeFieldOracleReader*>(context)->value(component,i,j,globalK,stateStep);
  }
  Double value(unsigned a,std::size_t i,std::size_t j,Double globalK,Double stateStep) {
    using namespace native_oracle_detail;
    require(a<6&&std::isfinite(globalK)&&std::isfinite(stateStep)&&globalK>=0&&
            stateStep>=0&&stateStep<=expected_.steps&&globalK<=expected_.geometry.nz&&
            globalK==std::floor(globalK)&&stateStep==std::floor(stateStep),
            "query must use valid integer native state and world lattice");
    Frame& frame=load(static_cast<unsigned>(stateStep));
    const auto& c=frame.components[a];const std::size_t k=static_cast<std::size_t>(globalK);
    require(i<c.nx&&j<c.ny&&k>=frame.first&&k-frame.first<c.nz,
            "requested native value outside recorded slab");
    const Double result=c.values[((k-frame.first)*c.ny+j)*c.nx+i];
    require(std::isfinite(result),"nonfinite replay value");++queries_;return result;
  }
  std::size_t memoryBytes() const {
    std::size_t n=0;for(const auto& f:frames_)for(const auto& c:f.components)n+=c.values.capacity()*sizeof(Double);
    return n;
  }
  std::uint64_t queries() const{return queries_;}
  std::uint64_t loads() const{return loads_;}
  std::uint64_t bytesRead() const{return bytesRead_;}
private:
  Frame& load(unsigned state) {
    using namespace native_oracle_detail;
    for(auto& f:frames_)if(f.valid&&f.state==state){f.touched=++clock_;return f;}
    Frame& f=!frames_[0].valid?frames_[0]:!frames_[1].valid?frames_[1]:
      (frames_[0].touched<=frames_[1].touched?frames_[0]:frames_[1]);
    f.valid=false;
    std::ifstream stream(framePath(directory_,state),std::ios::binary);
    require(stream.good(),"missing state "+std::to_string(state));magic(stream,"APRLOF02");
    std::uint64_t header[4];read(stream,header,sizeof(header));
    require(header[0]==state&&header[1]==firstK(expected_,state)&&
            header[2]==slabCells(expected_)&&header[3]==6&&
            header[1]<=expected_.geometry.nz&&header[2]<=expected_.geometry.nz-header[1],
            "invalid frame geometry or time");
    f.state=state;f.first=header[1];std::uint64_t bytes=8+sizeof(header);
    const std::size_t nx=expected_.geometry.nx,ny=expected_.geometry.ny;
    const std::size_t x[6]={nx,nx+1,nx+1,nx+1,nx,nx};
    const std::size_t y[6]={ny+1,ny,ny+1,ny,ny+1,ny};
    for(unsigned a=0;a<6;++a) {
      std::uint64_t d[4];read(stream,d,sizeof(d));auto& c=f.components[a];
      const std::size_t z=slabCells(expected_)+(vertexZ(a)?1:0);
      require(d[0]==a&&d[1]==x[a]&&d[2]==y[a]&&d[3]==z,"component geometry mismatch");
      c.nx=x[a];c.ny=y[a];c.nz=z;c.values.resize(c.nx*c.ny*c.nz);
      read(stream,c.values.data(),c.values.size()*sizeof(Double));
      for(Double v:c.values)require(std::isfinite(v),"nonfinite cached frame");
      bytes+=sizeof(d)+c.values.size()*sizeof(Double);
    }
    eof(stream);f.valid=true;f.touched=++clock_;++loads_;bytesRead_+=bytes;return f;
  }
  std::string directory_;NativeOracleConfig expected_;std::array<Frame,2> frames_;
  std::uint64_t queries_=0,loads_=0,bytesRead_=0,clock_=0;
};

} } // namespace aprl::verification
#endif
