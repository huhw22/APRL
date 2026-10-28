// Directory visibility check, deliberately warming missing-path lookups before
// the real collective helper. No directory/file retry is allowed in new mode.
// Run negative modes separately so a broken collective cannot hide later cases.
#include "runtime_util.h"

#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <mpi.h>

namespace {
void require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}
std::string failure(const char* operation,const std::string& path,int error) {
  std::ostringstream out;out<<operation<<" '"<<path<<"': errno="<<error<<" ("<<std::strerror(error)<<')';
  return out.str();
}
void absent(const std::string& path) {
  struct stat value;errno=0;const int result=::stat(path.c_str(),&value);const int error=errno;
  require(result==-1&&error==ENOENT,
    failure("Expected fresh path with stat returning ENOENT",path,error));
}
struct Descriptor {
  int value;
  explicit Descriptor(int fd):value(fd) {}
  ~Descriptor(){if(value>=0)::close(value);}
  void close() {
    const int fd=value;value=-1;
    require(::close(fd)==0,"Cannot close output-directory check file");
  }
};
void writeRead(const std::string& path,const std::string& payload) {
  // This must be the first operation after createDirectories in the new case:
  // no stat, access, sleep, extra collective or transient-error retry.
  Descriptor out(::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0664));
  const int openError=errno;
  require(out.value>=0,failure("Immediate exclusive open",path,openError));
  std::size_t offset=0;
  while(offset<payload.size()) {
    const ssize_t count=::write(out.value,payload.data()+offset,payload.size()-offset);
    if(count<0&&errno==EINTR)continue;
    const int error=errno;require(count>0,failure("write",path,error));offset+=count;
  }
  const int syncStatus=::fsync(out.value),syncError=errno;
  require(syncStatus==0,failure("fsync",path,syncError));out.close();
  Descriptor in(::open(path.c_str(),O_RDONLY));
  const int readOpenError=errno;
  require(in.value>=0,failure("Read-back open",path,readOpenError));
  std::string result(payload.size()+1,'\0');offset=0;
  while(offset<result.size()) {
    const ssize_t count=::read(in.value,&result[offset],result.size()-offset);
    if(count<0&&errno==EINTR)continue;
    const int error=errno;require(count>=0,failure("read",path,error));
    if(count==0)break;
    offset+=count;
  }
  result.resize(offset);require(result==payload,"Read-back bytes/length differ: "+path);in.close();
}
bool safePrefix(const std::string& value) {
  if(value.empty()||value=="."||value=="..")return false;
  for(char c:value)if(!((c>='a'&&c<='z')||(c>='A'&&c<='Z')||
      (c>='0'&&c<='9')||c=='-'||c=='_'||c=='.'))return false;
  return true;
}
}

int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);int rank=0,ranks=1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  try {
    require((ranks==1||ranks==2)&&argc>=3,
      "Usage: CHECK new EXISTING_PARENT UNIQUE_PREFIX | CHECK file EXISTING_FILE | CHECK parent-file EXISTING_FILE; use 1 or 2 ranks");
    const std::string mode=argv[1],input=argv[2];
    require((mode=="new"&&argc==4)||(mode=="file"&&argc==3)||(mode=="parent-file"&&argc==3),
      "Invalid output-directory check mode/arguments");
    char host[MPI_MAX_PROCESSOR_NAME]={};int hostLength=0;
    MPI_Get_processor_name(host,&hostLength);
    std::vector<char> hosts(static_cast<std::size_t>(ranks)*MPI_MAX_PROCESSOR_NAME,0);
    MPI_Gather(host,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,hosts.data(),MPI_MAX_PROCESSOR_NAME,MPI_CHAR,0,MPI_COMM_WORLD);
    struct stat inputStatus;
    require(::stat(input.c_str(),&inputStatus)==0,"Input must already exist on every rank: "+input);
    unsigned rounds=0,rejections=0;
    if(mode=="new") {
      require(S_ISDIR(inputStatus.st_mode),"new mode requires an existing parent directory");
      const std::string prefix=argv[3];require(safePrefix(prefix),"Use a fresh simple prefix of letters/digits/._-");
      for(unsigned round=0;round<16;++round) {
        std::ostringstream name;name<<prefix<<"-r"<<std::setw(2)<<std::setfill('0')<<round;
        const std::string directory=aprl::joinPath(input,name.str()+"/inner/leaf");
        const std::string file=aprl::joinPath(directory,"rank-"+std::to_string(rank)+".dat");
        const std::string rootFile=aprl::joinPath(directory,"root.dat");
        const std::string payload="directory-check round="+std::to_string(round)+" rank="+std::to_string(rank)+"\n";
        absent(directory);absent(file);absent(rootFile);
        MPI_Barrier(MPI_COMM_WORLD);
        aprl::createDirectories(directory,MPI_COMM_WORLD);
        writeRead(file,payload);
        if(rank==0)writeRead(rootFile,"root-output round="+std::to_string(round)+"\n");
        // Only after the opens and read-backs, coordinate the next fresh case.
        MPI_Barrier(MPI_COMM_WORLD);++rounds;
      }
    } else {
      require(S_ISREG(inputStatus.st_mode),"Negative modes require an existing regular file");
      const std::string requested=mode=="file"?input:aprl::joinPath(input,"must-not-be-created");
      MPI_Barrier(MPI_COMM_WORLD);unsigned rejected=0;
      try{aprl::createDirectories(requested,MPI_COMM_WORLD);}
      catch(const std::exception&){rejected=1;}
      MPI_Allreduce(&rejected,&rejections,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
      require(rejections==static_cast<unsigned>(ranks),
        "createDirectories must reject the invalid path on every rank");
      struct stat after;require(::stat(input.c_str(),&after)==0&&S_ISREG(after.st_mode)&&
        after.st_size==inputStatus.st_size&&after.st_ino==inputStatus.st_ino,
        "Negative check changed the original regular file");
    }
    if(rank==0) {
      std::cout<<"{\"status\":\"PASS\",\"mode\":\""<<mode<<"\",\"ranks\":"<<ranks
        <<",\"rounds\":"<<rounds<<",\"rank_files\":"<<rounds*ranks<<",\"root_files\":"<<rounds
        <<",\"collective_rejections\":"<<rejections
        <<",\"stat_before_create\":\"ENOENT in new mode\",\"leaf_precreated\":false"
        <<",\"extra_barrier_before_open\":false,\"retrying_opens\":false,\"hosts\":[";
      for(int r=0;r<ranks;++r){if(r)std::cout<<',';std::cout<<'"'<<&hosts[static_cast<std::size_t>(r)*MPI_MAX_PROCESSOR_NAME]<<'"';}
      std::cout<<"],\"scope\":\"shared-directory immediate visibility and collective failure only; not HDF5 or simulation acceptance\"}"<<std::endl;
    }
    MPI_Finalize();return 0;
  }catch(const std::exception& error) {
    std::cerr<<"output-directory-check rank="<<rank<<": "<<error.what()<<std::endl;
    MPI_Abort(MPI_COMM_WORLD,2);return 2;
  }
}
