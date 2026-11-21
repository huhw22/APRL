#ifndef APRL_BEAM_CPML_START_H
#define APRL_BEAM_CPML_START_H

#include "config.h"
#include "eb_cpml.h"
#include <yaml-cpp/yaml.h>
#include <mpi.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace aprl { namespace beam_cpml_start {

inline std::size_t product(std::size_t a,std::size_t b) {
  if(b && a>std::numeric_limits<std::size_t>::max()/b) throw std::overflow_error("CPML startup array size overflow");
  return a*b;
}
inline std::vector<Double> readSlab(const std::string& path,std::size_t globalPlanes,
    std::size_t first,std::size_t planes,std::size_t width) {
  static_assert(sizeof(Double)==8,"CPML startup arrays require binary64");
  const std::uint16_t one=1;
  if(*reinterpret_cast<const unsigned char*>(&one)!=1) throw std::runtime_error("CPML startup arrays require little endian doubles");
  if(!width || first>globalPlanes || planes>globalPlanes-first) throw std::runtime_error("Invalid CPML startup slab "+path);
  const auto row=product(width,sizeof(Double)), size=product(row,globalPlanes);
  if(size>std::size_t(std::numeric_limits<std::streamoff>::max())) throw std::overflow_error("CPML startup file exceeds stream offset");
  std::ifstream in(path,std::ios::binary|std::ios::ate);
  if(!in || in.tellg()!=static_cast<std::streamoff>(size)) throw std::runtime_error("Incorrect CPML startup file size "+path);
  in.seekg(static_cast<std::streamoff>(first*row));
  std::vector<Double> data(product(planes,width));
  if(!in.read(reinterpret_cast<char*>(data.data()),static_cast<std::streamsize>(planes*row)))
    throw std::runtime_error("Short CPML startup slab "+path);
  for(const auto v:data) if(!std::isfinite(v)) throw std::runtime_error("Nonfinite CPML startup sample "+path);
  return data;
}
inline std::string fingerprint(const std::string& path) {
  std::ifstream in(path,std::ios::binary);
  if(!in)throw std::runtime_error("Cannot open CPML startup input "+path);
  std::array<unsigned char,65536> buf;
  std::uint64_t hash=14695981039346656037ULL;
  for(;;) {
    in.read(reinterpret_cast<char*>(buf.data()),buf.size());
    const auto n=in.gcount();
    for(std::streamsize i=0;i<n;++i){hash^=buf[i];hash*=1099511628211ULL;}
    if(in.bad() || (!in && !in.eof())) throw std::runtime_error("Cannot read CPML startup input "+path);
    if(in.eof()) break;
  }
  std::ostringstream out;
  out<<std::setfill('0')<<std::setw(16)<<std::hex<<hash;
  return out.str();
}
inline std::vector<unsigned> historyIds(bool left) {
  return left ? std::vector<unsigned>{3,4,9,10} : std::vector<unsigned>{0,3,4,5,6,9,10,11};
}
inline std::vector<std::string> names(bool left) {
  std::vector<std::string> list={"Ex","Ey","Ez","Bx","By","Bz"};
  for(auto h:historyIds(left)) list.push_back("psi"+std::to_string(h));
  return list;
}
inline void require(bool cond,const std::string& key) {
  if(!cond) throw std::runtime_error("CPML startup package does not match "+key);
}
inline void checkFloat(const YAML::Node& node,const char* key,double actual,double floor=0.) {
  const double expected=node.as<double>();
  const double scale=std::max({std::abs(actual),std::abs(expected),floor});
  require(std::isfinite(expected)&&std::isfinite(actual)&&std::abs(expected-actual)<=2e-13*scale,key);
}
inline void verifyPackage(const SimulationConfig& c,const EBGridGeometry& ref,const FieldVector<Double>& origin,
    double beta,std::size_t windowNz,std::ptrdiff_t offset) {
  const auto& mw=c.mesh.movingWindow;
  const auto& dir=mw.cpmlStartupDirectory;
  const auto p=YAML::LoadFile(dir+"/startup.yaml");
  require(p["version"].as<int>()==1 && p["method"].as<std::string>()==mw.cpmlStartup,"method/version");
  require(mw.cpmlStartup=="osiris" || mw.cpmlStartup=="left-20056","method");
  require(mw.openTransverseReference && mw.referenceThroughCpml && !mw.referenceSurfaceCurrent,"open transverse Cowan reference");
  require(c.beam.type==BeamInputType::Hdf5,"HDF5 beam");
  const std::vector<std::size_t> dims={ref.nx,ref.ny,ref.nz},window={ref.nx,ref.ny,windowNz};
  require(p["window_cells"].as<std::vector<std::size_t>>()==window,"window cells");
  require(p["reference_cells"].as<std::vector<std::size_t>>()==dims,"reference cells");
  require(p["reference_offset"].as<std::ptrdiff_t>()==offset,"window offset");
  const auto& cp=c.boundary.cpml;
  require(p["cpml_cells"].as<std::vector<std::size_t>>()==std::vector<std::size_t>(cp.cells,cp.cells+3),"CPML cells");
  require(cp.cells[0]>0 && cp.cells[1]>0 && cp.cells[2]==0,"four transverse sides");
  checkFloat(p["target_reflection"],"target reflection",cp.targetReflection);
  checkFloat(p["polynomial_order"],"polynomial order",cp.polynomialOrder);
  checkFloat(p["kappa_max"],"kappa",cp.kappaMax);
  checkFloat(p["alpha_fraction"],"alpha",cp.alphaFraction);
  const double h[]={ref.dx,ref.dy,ref.dz};
  for(unsigned i=0;i<3;++i) {
    checkFloat(p["reference_origin_m"][i],"reference origin",origin[i],h[i]);
    checkFloat(p["spacing_m"][i],"spacing",h[i]);
    checkFloat(p["beam_position_offset_m"][i],"beam offset",c.beam.positionOffset[i],h[i]);
  }
  checkFloat(p["dt_s"],"timestep",ref.dt);
  checkFloat(p["reference_beta"],"reference beta",beta);
  checkFloat(p["frame_gamma"],"boost gamma",c.mesh.boostGamma);
  checkFloat(p["beam_electrons"],"beam electrons",c.beam.electrons);
  const auto b=p["beam_reference"];
  checkFloat(b["initial_center_z"],"initial beam center",c.reference.initialCenterZ,ref.dz);
  require(bool(b["input_plane_z"])==c.reference.inputPlaneZSet,"input plane presence");
  if(c.reference.inputPlaneZSet) checkFloat(b["input_plane_z"],"input plane",c.reference.inputPlaneZ,ref.dz);
  require(bool(b["activation_plane_z"])==c.reference.activationPlaneZSet,"activation plane presence");
  if(c.reference.activationPlaneZSet) checkFloat(b["activation_plane_z"],"activation plane",c.reference.activationPlaneZ,ref.dz);
  require(b["frozen_lab_plane_injection"].as<bool>(false)==c.reference.frozenLabPlaneInjection,"frozen injection");
  require(b["hold_transverse_until_activation"].as<bool>(false)==c.reference.holdTransverseUntilActivation,"transverse hold");
  require(p["transverse_x_binomial3"].as<bool>()==c.mesh.transverseXBinomial3,"transverse shape");
  const std::string shape=c.mesh.longitudinalBinomial3?"linear-binomial3":c.mesh.longitudinalNotch?"linear-notch":
      c.mesh.longitudinalQuintic?"quintic-compensated":"linear";
  require(p["longitudinal_shape"].as<std::string>()==shape,"longitudinal shape");
  require(p["beam_fnv1a64"].as<std::string>()==fingerprint(c.beam.file),"beam fingerprint");
  const auto entries=names(mw.cpmlStartup=="left-20056");
  const auto hashes=p["arrays_fnv1a64"];
  require(hashes.IsMap() && hashes.size()==entries.size(),"array inventory");
  for(const auto& name:entries) require(hashes[name].as<std::string>()==fingerprint(dir+"/"+name+".bin"),name+" fingerprint");
}
inline std::vector<std::size_t> exteriorIndices(const YeeComponent& u,std::size_t lx,std::size_t ly) {
  require(lx && ly && lx<(u.nx()+1)/2 && ly<(u.ny()+1)/2,"transverse collar");
  std::vector<std::size_t> indices;
  for(std::size_t j=0;j<u.ny();++j)for(std::size_t i=0;i<u.nx();++i)
    if(i<lx || i>=u.nx()-lx || j<ly || j>=u.ny()-ly) indices.push_back(j*u.nx()+i);
  return indices;
}
inline void load(const std::string& dir,bool left,EBFieldGrid& grid,std::size_t localOffset,
    std::size_t globalNz,std::size_t lx,std::size_t ly) {
  auto& pml=dynamic_cast<EBConvolutionalPML&>(grid.boundary());
  std::array<YeeComponent*,6> fields={{&grid.ex(),&grid.ey(),&grid.ez(),&grid.bx(),&grid.by(),&grid.bz()}};
  const auto list=names(left); const auto ids=historyIds(left);
  for(std::size_t c=0;c<list.size();++c) {
    const bool memory=c>=fields.size();
    const unsigned term=memory?ids[c-fields.size()]:0;
    const unsigned target=memory?(term<6?3:0)+(term%6)/2:static_cast<unsigned>(c);
    auto& u=*fields[target];
    const auto path=dir+"/"+list[c]+".bin";
    if(left) {
      const auto rows=grid.geometry().ny+1,width=product(rows,lx);
      const auto data=readSlab(path,globalNz+1,localOffset,u.nz(),width);
      std::vector<Double>* hist=memory?&pml.checkpointHistory(term):nullptr;
      if(hist)require(hist->size()==product(2*lx,product(u.nz(),u.ny())),"left history layout");
      for(std::size_t k=0;k<u.nz();++k)for(std::size_t j=0;j<u.ny();++j)for(std::size_t i=0;i<lx;++i) {
        const double v=data[(k*rows+j)*lx+i];
        if(hist)(*hist)[(i*u.nz()+k)*u.ny()+j]=v;
        else u(i,j,k)+=v;
      }
    } else if(memory) {
      const auto layout=pml.nativeZHistoryLayout(term);
      const bool alongX=term%6==3 || term%6==4;
      require(layout.layers==u.nz() && layout.valuesPerPlane==2*(alongX?lx*u.ny():ly*u.nx()),"history layout");
      const auto data=readSlab(path,globalNz+layout.nodeCentered,localOffset,layout.layers,layout.valuesPerPlane);
      pml.restoreNativeZHistory(term,data);
    } else {
      const auto ix=exteriorIndices(u,lx,ly);
      const auto data=readSlab(path,globalNz+(u.nz()==grid.geometry().nz+1),localOffset,u.nz(),ix.size());
      const auto area=u.nx()*u.ny();
      for(std::size_t k=0;k<u.nz();++k)for(std::size_t n=0;n<ix.size();++n)
        u.data()[k*area+ix[n]]=data[k*ix.size()+n];
    }
  }
}
// Propagate the first reader error before continuing with any collective or update.
inline void checkRanks(const std::string& error,int rank,MPI_Comm comm) {
  int size=0;MPI_Comm_size(comm,&size);
  int failing=error.empty()?size:rank;
  MPI_Allreduce(MPI_IN_PLACE,&failing,1,MPI_INT,MPI_MIN,comm);
  if(failing==size)return;
  int length=rank==failing?static_cast<int>(error.size()):0;
  MPI_Bcast(&length,1,MPI_INT,failing,comm);
  std::string message=rank==failing?error:std::string(length,' ');
  if(length) MPI_Bcast(&message[0],length,MPI_CHAR,failing,comm);
  throw std::runtime_error(message);
}
inline void initialize(const SimulationConfig& c,const EBGridGeometry& ref,const FieldVector<Double>& origin,
    double beta,EBFieldGrid& grid,std::size_t localOffset,std::size_t windowNz,
    std::ptrdiff_t offset,int rank,MPI_Comm comm) {
  if(c.mesh.movingWindow.cpmlStartup=="standard") return;
  std::string error;
  if(rank==0) try {verifyPackage(c,ref,origin,beta,windowNz,offset);}
  catch(const std::exception& e){error=e.what();}
  checkRanks(error,rank,comm);
  try {load(c.mesh.movingWindow.cpmlStartupDirectory,c.mesh.movingWindow.cpmlStartup=="left-20056",
      grid,localOffset,windowNz,c.boundary.cpml.cells[0],c.boundary.cpml.cells[1]);}
  catch(const std::exception& e){error=e.what();}
  checkRanks(error,rank,comm);
}
}}
#endif
