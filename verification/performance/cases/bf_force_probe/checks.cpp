#include "particle_force_probe.h"
#include "source_coupling_math.h"
#include "runtime_util.h"
#include <yaml-cpp/yaml.h>
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace aprl;
using Vec=FieldVector<Double>;
namespace {
void need(bool value,const char* text){if(!value)throw std::runtime_error(text);}
struct Fixture {
  RelativisticParticleSI first,last;
  Vec electric=Vec(0.),magnetic=Vec(0.);
  std::vector<RelativisticParticleSI> states;
  std::vector<BorisWork> works;
  SIFieldSourceSet sources;
  BoostFrameTransform frame;
  double time=0.,dt=0.;unsigned substeps=7;
};
Fixture fixture(unsigned n,bool device,bool self=true) {
  Fixture f;f.frame.setFromGamma(2.,SI::c,0.);
  f.first.id=100+n;f.first.sourceId=n%16+1;f.first.weight=100.+n;
  f.first.charge=-SI::elementaryCharge*f.first.weight;
  f.first.mass=SI::electronMass*f.first.weight;
  f.first.position[0]=2e-6;f.first.position[1]=-3e-6;f.first.position[2]=.02;
  f.first.properVelocity[0]=.04+.001*n;f.first.properVelocity[1]=-.03;
  f.first.properVelocity[2]=.6+.003*n;
  f.time=1e-12+n*1e-14;f.dt=1e-13*(1.+.031*n);
  if(self)for(unsigned a=0;a<3;++a) {
    f.electric[a]=2e5*std::sin(.71*n+.37*a);
    f.magnetic[a]=.04*std::cos(.43*n+.51*a);
  }
  if(device) {
    SIMagneticElement magnet;magnet.type=SIMagnetType::UniformDipole;
    magnet.center=Vec(0.);magnet.length=1.;magnet.peakMagneticField=.12;
    magnet.polarizationAngle=.37;magnet.prepare();
    f.sources.addPrescribedLabMagnet(magnet);
  }
  auto p=f.first;f.states.push_back(p);const double h=f.dt/f.substeps;
  for(unsigned k=0;k<f.substeps;++k) {
    const auto before=p;const double g0=BoostFrameTransform::gammaFromProperVelocity(p.properVelocity);
    Vec point(p.position);point.pmv(.5*h*SI::c/g0,p.properVelocity);
    Vec ee(0.),bb(0.);f.sources.addPrescribedBox(point,f.time+(k+.5)*h,f.frame,ee,bb);
    Vec E=f.electric;E+=ee;Vec B=f.magnetic;B+=bb;
    source_coupling_math::momentum(p.properVelocity,E,B,p.charge,p.mass,h,g0);
    const double g1=BoostFrameTransform::gammaFromProperVelocity(p.properVelocity);
    Vec mean(before.properVelocity);mean+=p.properVelocity;
    p.position.pmv(h*SI::c/(g0+g1),mean);
    BorisWork work;work.gridWork=source_coupling_math::chargeWork(p.charge,f.electric,before.position,p.position);
    work.prescribedWork=source_coupling_math::chargeWork(p.charge,ee,before.position,p.position);
    work.totalDeltaK=static_cast<long double>(p.mass)*SI::c*SI::c*
      source_coupling_math::deltaGamma(before.properVelocity,p.properVelocity);
    f.works.push_back(work);f.states.push_back(p);
  }
  f.last=p;return f;
}
ParticleForceProbe::Record make(const Fixture& f,unsigned step) {
  return ParticleForceProbe::makeRecord(step,f.first,f.last,f.electric,f.magnetic,
    f.time,f.dt,f.substeps,f.states,f.works,f.sources,f.frame);
}
struct Errors {double momentum=0.,energy=0.,deviceLab=0.;};
Errors budgets() {
  Errors errors;
  for(unsigned n=0;n<32;++n)for(bool device:{false,true}) {
    const auto f=fixture(n,device);const auto r=make(f,n+1);
    long double displacement[3];for(unsigned a=0;a<3;++a)
      displacement[a]=static_cast<long double>(r.positionAfter[a])-r.positionBefore[a];
    long double gammaGrid=0.;
    for(unsigned a=0;a<3;++a) {
      const unsigned b=(a+1)%3,c=(a+2)%3;
      const long double electric=static_cast<long double>(r.charge)*r.dt/(r.mass*SI::c)*r.gridElectric[a];
      const long double magnetic=static_cast<long double>(r.charge)/(r.mass*SI::c)*
        (displacement[b]*r.gridMagnetic[c]-displacement[c]*r.gridMagnetic[b]);
      const long double total=electric+magnetic+r.prescribedImpulseElectric[a]+r.prescribedImpulseMagnetic[a];
      errors.momentum=std::max(errors.momentum,static_cast<double>(std::abs(
        (static_cast<long double>(r.properAfter[a])-r.properBefore[a])-total)));
      gammaGrid+=static_cast<long double>(r.charge)/(r.mass*SI::c*SI::c)*r.gridElectric[a]*displacement[a];
    }
    const long double actual=source_coupling_math::deltaGamma(f.first.properVelocity,f.last.properVelocity);
    errors.energy=std::max(errors.energy,static_cast<double>(std::abs(actual-gammaGrid-r.prescribedDeltaGamma)));
  }
  for(unsigned n=0;n<16;++n) {
    const auto f=fixture(n,true,false);const auto r=make(f,n+1);
    const long double lab=f.frame.gamma()*(static_cast<long double>(r.prescribedDeltaGamma)+f.frame.beta()*
      (static_cast<long double>(r.prescribedImpulseElectric[2])+r.prescribedImpulseMagnetic[2]));
    errors.deviceLab=std::max(errors.deviceLab,static_cast<double>(std::abs(lab)));
  }
  need(errors.momentum<2e-12,"Interval self/device momentum budget failed");
  need(errors.energy<2e-12,"Interval self/device energy budget failed");
  need(errors.deviceLab<2e-12,"Static lab device unexpectedly supplies lab work");
  // An altered saved device work must be rejected; no false matching through
  // redefining a missing contribution as an arithmetic residual.
  auto bad=fixture(7,true);bad.works[0].prescribedWork*=1.1;
  bool rejected=false;try{make(bad,1);}catch(const std::runtime_error&){rejected=true;}
  need(rejected,"Altered actual device work was not rejected");
  return errors;
}
unsigned readFiles(const std::string& output,int rank,unsigned expected) {
  char stem[64];std::snprintf(stem,sizeof(stem),"/force-rank-%05d",rank);
  const YAML::Node manifest=YAML::LoadFile(output+stem+".json");
  need(manifest["complete"].as<bool>(),"Probe completion false");
  need(manifest["record_bytes"].as<unsigned>()==288,"Reader record schema differs");
  need(manifest["records"].as<unsigned>()==expected,"Missing selected intervals");
  unsigned total=0;
  for(const auto& node:manifest["files"]) {
    std::ifstream input(node["path"].as<std::string>(),std::ios::binary);
    ParticleForceProbe::Header header;input.read(reinterpret_cast<char*>(&header),sizeof(header));
    need(bool(input)&&std::string(header.magic)=="APRLFP1"&&header.recordBytes==288&&header.headerBytes==64,
      "Invalid binary trace framing");
    need(header.rank==static_cast<unsigned>(rank)&&header.endianMarker==0x01020304,"Header rank/endian differs");
    ParticleForceProbe::Record record;std::uint64_t hash=14695981039346656037ULL;unsigned count=0;
    while(input.read(reinterpret_cast<char*>(&record),sizeof(record))) {
      need(record.sourceId>=1&&record.sourceId<=16&&record.flags==1&&record.substeps==7,
        "Unexpected trace identity or validity flags");
      const auto* bytes=reinterpret_cast<const unsigned char*>(&record);
      for(std::size_t i=0;i<sizeof(record);++i){hash^=bytes[i];hash*=1099511628211ULL;}
      ++count;
    }
    need(input.eof()&&input.gcount()==0,"Truncated binary record accepted");
    need(count==node["records"].as<unsigned>(),"File and manifest counts differ");
    need(hash==std::stoull(node["records_fnv1a64"].as<std::string>(),nullptr,16),"Trace checksum differs");
    total+=count;
  }
  need(total==expected,"Trace total differs");return total;
}
}
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);int rank=0,size=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  try {
    need(argc>=2,"Pass a new output directory");const std::string output=argv[1];
    createDirectories(output,MPI_COMM_WORLD);const bool invalid=argc>2&&std::string(argv[2])=="--invalid-ids";
    const bool disabled=argc>2&&std::string(argv[2])=="--disabled";
    const std::string ids=output+"/selected-ids.txt";
    if(rank==0){std::ofstream stream(ids);for(unsigned i=1;i<=(invalid?15:16);++i)stream<<i<<'\n';}
    MPI_Barrier(MPI_COMM_WORLD);
    if(disabled)unsetenv("APRL_PARTICLE_FORCE_PROBE");
    else{setenv("APRL_PARTICLE_FORCE_PROBE",output.c_str(),1);setenv("APRL_PARTICLE_FORCE_PROBE_IDS",ids.c_str(),1);}
    if(invalid) {
      bool rejected=false;try{ParticleForceProbe::enabled();}catch(const std::invalid_argument&){rejected=true;}
      need(rejected,"Incomplete phase family accepted");
      if(rank==0)std::cout<<"{\"passed\":true,\"incomplete_family_rejected\":true}\n";
      MPI_Finalize();return 0;
    }
    const Errors error=budgets();unsigned count=0;
    for(unsigned i=0;i<16;++i)if(i%size==static_cast<unsigned>(rank))++count;
#pragma omp parallel for num_threads(2) schedule(static)
    for(int i=0;i<16;++i)if(i%size==rank) {
      const auto f=fixture(i,true);ParticleForceProbe::record(i+1,rank,f.first,f.last,
        f.electric,f.magnetic,f.time,f.dt,f.substeps,f.states,f.works,f.sources,f.frame);
    }
    ParticleForceProbe::close(rank,16,true);
    if(!disabled)readFiles(output,rank,count);
    else need(!ParticleForceProbe::enabled(),"Default-off probe became active");
    unsigned global=0;MPI_Reduce(&count,&global,1,MPI_UNSIGNED,MPI_SUM,0,MPI_COMM_WORLD);
    if(rank==0)std::cout<<std::setprecision(17)<<"{\"passed\":true,\"momentum_abs_error\":"<<error.momentum
      <<",\"energy_abs_error\":"<<error.energy<<",\"static_lab_device_work_abs\":"<<error.deviceLab
      <<",\"records\":"<<(disabled?0:global)<<",\"disabled\":"<<(disabled?"true":"false")<<"}\n";
    MPI_Finalize();return 0;
  }catch(const std::exception& e){std::cerr<<"force probe check failed: "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}
}
