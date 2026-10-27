// Diagnostic-file integrity checks only. No field or particle time advance.
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <mpi.h>

#include "streaming_field_oracle.h"

namespace {
using aprl::Double;
using aprl::EBFieldGrid;
using aprl::EBGridGeometry;
using aprl::YeeComponent;
using aprl::verification::NativeOracleConfig;
using aprl::verification::NativeFieldOracleWriter;
using aprl::verification::NativeFieldOracleReader;

std::uint64_t assertions=0,errors=0,bitwiseFailures=0,replayedValues=0;
std::uint64_t rejections=0,frameLoads=0,readBytes=0,recordedBytes=0;

void check(bool value,const std::string& message) {
  ++assertions;
  if(!value){++errors;throw std::runtime_error(message);}
}

bool exists(const std::string& path) {
  struct stat info;
  return stat(path.c_str(),&info)==0;
}

void newDirectory(const std::string& path) {
  check(mkdir(path.c_str(),0755)==0,"Directory must be new: "+path);
}

template<typename Operation>
void rejects(Operation operation,const std::string& expected,const std::string& label) {
  bool rejected=false;
  try {operation();}
  catch(const std::exception& error) {
    rejected=std::string(error.what()).find(expected)!=std::string::npos;
    check(rejected,label+": wrong rejection: "+error.what());
  }
  check(rejected,label+": invalid input was accepted");
  ++rejections;
}

YeeComponent& component(EBFieldGrid& fields,unsigned a) {
  switch(a) {
    case 0:return fields.ex();case 1:return fields.ey();case 2:return fields.ez();
    case 3:return fields.bx();case 4:return fields.by();case 5:return fields.bz();
  }
  throw std::logic_error("Invalid test component");
}

Double tag(unsigned state,unsigned a,std::size_t i,std::size_t j,std::size_t k) {
  // Dyadic, exactly representable tags, including both signs of zero.
  const std::uint64_t parity=state+a+i+3*j+5*k;
  if(parity%113==0)return std::copysign(0.,parity%2?-1.:1.);
  const Double value=static_cast<Double>(1+state*1000000+a*100000+k*64+j*8+i)/8.;
  return parity%2?-value:value;
}

void fill(EBFieldGrid& fields,unsigned state) {
  for(unsigned a=0;a<6;++a) {
    auto& f=component(fields,a);
    for(std::size_t k=0;k<f.nz();++k)
      for(std::size_t j=0;j<f.ny();++j)
        for(std::size_t i=0;i<f.nx();++i)f(i,j,k)=tag(state,a,i,j,k);
  }
}

std::size_t fileSize(const std::string& path) {
  struct stat info;
  check(stat(path.c_str(),&info)==0&&info.st_size>=0,"Cannot stat "+path);
  return static_cast<std::size_t>(info.st_size);
}

// Broken fixtures are independent copies; valid donor files are never edited.
void copyFile(const std::string& from,const std::string& to,std::size_t trim=0) {
  check(!exists(to),"Refusing to overwrite "+to);
  const std::size_t size=fileSize(from);
  check(trim<=size,"Invalid truncation size");
  std::vector<char> data(size);
  std::ifstream in(from,std::ios::binary);
  in.read(data.data(),static_cast<std::streamsize>(data.size()));
  check(in.good(),"Cannot read fixture "+from);
  std::ofstream out(to,std::ios::binary);
  out.write(data.data(),static_cast<std::streamsize>(size-trim));out.close();
  check(!out.fail(),"Cannot write fixture "+to);
}

void copyHeaderAndMarker(const std::string& valid,const std::string& target) {
  newDirectory(target);
  copyFile(valid+"/header.bin",target+"/header.bin");
  copyFile(valid+"/complete.bin",target+"/complete.bin");
}

void runChecks(const std::string& directory) {
  const Double dz=1e-6;
  const EBGridGeometry geometry(4,4,424,2*dz,2*dz,dz,dz/aprl::SI::c,
                                aprl::EBMaxwellSolver::CowanZ);
  const NativeOracleConfig config(geometry,1,24,710271,.23,-1e-13,
                                  -84e-6,218e-6,234e-6,.8,1);
  EBFieldGrid fields(geometry);
  const std::string valid=directory+"/valid";
  NativeFieldOracleWriter writer(valid,config);
  // Independent integer geometry: floor(218-.8*n)+84-24-20-2.
  const std::size_t first[2]={256,255},cells=108;
  for(unsigned state=0;state<2;++state) {
    fill(fields,state);writer.record(fields,state,first[state]);
  }
  writer.finish();
  check(writer.frames()==2,"Wrong recorded state count");
  recordedBytes=writer.bytes();
  check(recordedBytes==fileSize(valid+"/state-000000.bin")+
                       fileSize(valid+"/state-000001.bin"),"Recorded byte count differs");

  NativeFieldOracleReader reader(valid,config);
  check(reader.loads()==0&&reader.queries()==0,"Reader eagerly loaded states");
  // Revisit both states in reverse order to verify the two-frame cache.
  const unsigned order[4]={0,1,1,0};
  for(unsigned state:order)for(unsigned a=0;a<6;++a) {
    const auto& f=component(fields,a);
    const bool vertex=a==0||a==1||a==5;
    const std::size_t planes=cells+(vertex?1:0);
    for(std::size_t k=0;k<planes;++k)
      for(std::size_t j=0;j<f.ny();++j)
        for(std::size_t i=0;i<f.nx();++i) {
          const Double actual=reader.value(a,i,j,first[state]+k,state);
          const Double expected=tag(state,a,i,j,first[state]+k);
          const bool same=std::memcmp(&actual,&expected,sizeof(Double))==0;
          if(!same)++bitwiseFailures;
          check(same,"Native value did not replay bitwise");++replayedValues;
        }
  }
  check(reader.loads()==2,"Two-frame revisits caused additional loads");
  check(reader.queries()==replayedValues,"Successful replay query count differs");
  check(reader.bytesRead()==recordedBytes,"Reader byte count differs");
  const Double callback=NativeFieldOracleReader::query(&reader,4,2,2,first[1]+3,1);
  const Double expected=tag(1,4,2,2,first[1]+3);
  check(std::memcmp(&callback,&expected,sizeof(Double))==0,"Callback replay mismatch");

  NativeOracleConfig wrongKind=config;wrongKind.sourceKind=0;
  rejects([&](){NativeFieldOracleReader wrong(valid,wrongKind);},
          "donor source kind mismatch","Physical-history/independent-source confusion");
  rejects([&](){reader.value(0,0,0,first[0]+.25,0);},"integer native","Fractional coordinate");
  rejects([&](){reader.value(0,0,0,first[0],.5);},"integer native","Fractional state");
  rejects([&](){reader.value(0,0,0,first[0],2);},"integer native","State beyond recording");
  rejects([&](){reader.value(0,0,0,first[0],-1);},"integer native","Negative state");
  rejects([&](){reader.value(0,0,0,-1,0);},"integer native","Negative global coordinate");
  rejects([&](){reader.value(0,0,0,425,0);},"integer native","Beyond global mesh");
  rejects([&](){reader.value(0,0,0,first[0]-1,0);},"outside recorded slab","Below recorded slab");
  rejects([&](){reader.value(0,0,0,first[0]+cells+1,0);},"outside recorded slab","Above vertex slab");
  rejects([&](){reader.value(3,0,0,first[0]+cells,0);},"outside recorded slab","Past half-cell component");
  rejects([&](){reader.value(0,4,0,first[0],0);},"outside recorded slab","Invalid transverse x");
  rejects([&](){reader.value(0,0,5,first[0],0);},"outside recorded slab","Invalid transverse y");
  rejects([&](){reader.value(6,0,0,first[0],0);},"integer native","Invalid component");
  rejects([&](){reader.value(0,0,0,std::numeric_limits<Double>::quiet_NaN(),0);},
          "integer native","Nonfinite coordinate");
  rejects([&](){NativeFieldOracleReader::query(nullptr,0,0,0,first[0],0);},
          "null callback context","Null callback");
  check(reader.loads()==2,"Rejected coordinates unexpectedly reloaded cached frames");
  frameLoads=reader.loads();readBytes=reader.bytesRead();

  const std::string missing=directory+"/missing-record";newDirectory(missing);
  rejects([&](){NativeFieldOracleReader broken(missing,config);},
          "no completed full-domain recording","Missing recording");
  const std::string incomplete=directory+"/incomplete-record";
  NativeFieldOracleWriter unfinished(incomplete,config);
  fill(fields,0);unfinished.record(fields,0,first[0]);
  rejects([&](){unfinished.finish();},"incomplete recording","Premature writer completion");
  rejects([&](){NativeFieldOracleReader broken(incomplete,config);},
          "no completed full-domain recording","Incomplete recording");

  const std::string missingHeader=directory+"/missing-header";newDirectory(missingHeader);
  copyFile(valid+"/complete.bin",missingHeader+"/complete.bin");
  rejects([&](){NativeFieldOracleReader broken(missingHeader,config);},
          "missing recording header","Missing header");
  const std::string truncatedHeader=directory+"/truncated-header";newDirectory(truncatedHeader);
  copyFile(valid+"/complete.bin",truncatedHeader+"/complete.bin");
  copyFile(valid+"/header.bin",truncatedHeader+"/header.bin",1);
  rejects([&](){NativeFieldOracleReader broken(truncatedHeader,config);},
          "truncated binary input","Truncated header");
  const std::string missingFrame=directory+"/missing-frame";
  copyHeaderAndMarker(valid,missingFrame);
  rejects([&](){NativeFieldOracleReader broken(missingFrame,config);
                broken.value(0,0,0,first[0],0);},"missing state","Missing frame");
  const std::string truncatedFrame=directory+"/truncated-frame";
  copyHeaderAndMarker(valid,truncatedFrame);
  copyFile(valid+"/state-000000.bin",truncatedFrame+"/state-000000.bin",1);
  rejects([&](){NativeFieldOracleReader broken(truncatedFrame,config);
                broken.value(0,0,0,first[0],0);},"truncated binary input","Truncated frame");
  // Existing source/valid files remain immutable, including on rejected writes.
  rejects([&](){NativeFieldOracleWriter duplicate(valid,config);},
          "recording directory must be new","Refuse recording overwrite");
}

void writeMetrics(const std::string& directory,bool completed) {
  const std::string path=directory+"/metrics.json";
  if(exists(path))throw std::runtime_error("Refusing to overwrite metrics");
  std::ofstream out(path);
  if(!out.good())throw std::runtime_error("Cannot write oracle-check metrics");
  out<<"{\"format\":\"aprl-streaming-oracle-checks-v1\",\"completed\":"
     <<(completed?"true":"false")<<",\"assertions\":"<<assertions<<",\"checks\":"<<assertions
     <<",\"errors\":"<<errors<<",\"bitwise_failures\":"<<bitwiseFailures
     <<",\"replayed_values\":"<<replayedValues<<",\"expected_rejections\":"<<rejections
     <<",\"reader_frame_loads\":"<<frameLoads<<",\"reader_bytes\":"<<readBytes
     <<",\"recorded_bytes\":"<<recordedBytes
     <<",\"recording_format\":\"aprl-native-field-oracle-v2\",\"source_kind\":1"
       ",\"scope\":\"Single-rank synthetic diagnostic-file integrity; no field evolution, physical acceptance or production-speed claim\"}\n";
  out.close();if(out.fail())throw std::runtime_error("Cannot close oracle-check metrics");
}
} // namespace

int main(int argc,char** argv) {
  int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int size=0;MPI_Comm_size(MPI_COMM_WORLD,&size);
  std::string directory;bool ownOutput=false;
  try {
    check(size==1,"Oracle checks require one MPI rank");
    check(argc==2,"Usage: aprl_streaming_oracle_checks OUT_DIRECTORY");
    directory=argv[1];
    check(mkdir(directory.c_str(),0755)==0||errno==EEXIST,"Cannot create checks directory");
    check(!exists(directory+"/metrics.json"),"Refusing to overwrite metrics");
    ownOutput=true;
    runChecks(directory);writeMetrics(directory,true);
    std::cout<<"Native field oracle checks PASS: "<<assertions<<" assertions\n";
    MPI_Finalize();return 0;
  } catch(const std::exception& error) {
    if(errors==0)++errors;
    std::cerr<<"Native field oracle checks FAILED: "<<error.what()<<'\n';
    if(ownOutput&&!exists(directory+"/metrics.json"))
      try {writeMetrics(directory,false);}catch(...) {}
    MPI_Finalize();return 1;
  }
}
