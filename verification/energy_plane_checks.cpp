/* Diagnostic fixtures only: affine E/B sampling and fixed-energy endpoint
 * rotations. These are not a Maxwell or full PIC conservation solution. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <mpi.h>
#include <hdf5.h>
#include "lab_detectors.h"

namespace
{
  void require(bool ok, const std::string& text)
  { if (!ok) throw std::runtime_error(text); }
  void closeEnough(long double value, long double reference, const char* text)
  {
    require(std::isfinite(value) && std::abs(value-reference) <=
      3.e-12L*std::max(std::abs(reference),1.e-29L), text);
  }
  std::size_t offset(int rank, int size)
  { return (32/size)*rank+std::min(rank,32%size); }

  template<class T> std::vector<T> column(hid_t dataset, const char* name,
                                         hid_t native, std::size_t count)
  {
    hid_t type=H5Tcreate(H5T_COMPOUND,sizeof(T));
    require(type>=0 && H5Tinsert(type,name,0,native)>=0,"Cannot define read field");
    std::vector<T> result(count);
    require(H5Dread(dataset,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,result.data())>=0,"Cannot read energy field");
    H5Tclose(type); return result;
  }

  long double particleEnergy(hid_t d,const char* name,std::size_t count)
  {
    const std::string low=std::string(name)+"_low";
    return static_cast<long double>(column<double>(d,name,H5T_NATIVE_DOUBLE,count).back())+
      column<double>(d,low.c_str(),H5T_NATIVE_DOUBLE,count).back();
  }

  void fill(aprl::EBFieldGrid& fields,std::size_t first,unsigned step)
  {
    aprl::YeeComponent* arrays[]={&fields.ex(),&fields.ey(),&fields.ez(),
      &fields.bx(),&fields.by(),&fields.bz()};
    const double shifts[]={0,0,.5,.5,.5,0};
    for(unsigned a=0;a<6;++a)
      for(std::size_t k=0;k<arrays[a]->nz();++k)
        for(std::size_t j=0;j<arrays[a]->ny();++j)
          for(std::size_t i=0;i<arrays[a]->nx();++i)
            {
              const double q=first+k+shifts[a];
              (*arrays[a])(i,j,k)=a==0?1.+.1*q+.2*step:
                (a==4?(2.+.15*q-.07*step)/aprl::SI::c:0.);
            }
  }

  long double analyticPower(const aprl::BoostFrameTransform& frame,double plane,
                            double dt,unsigned n,long double area)
  {
    const double q=frame.boxZFromLabZAndBoxT(plane,n*dt)/.05;
    aprl::FieldVector<aprl::Double> e(0),b(0),el(0),bl(0);
    e[0]=1.+.1*q+.2*n;b[1]=(2.+.15*q-.07*n)/aprl::SI::c;
    frame.fieldsBoxToLab(e,b,el,bl);
    return area*el[0]*bl[1]/aprl::SI::mu0;
  }

  void run(const std::string& directory,int rank,int size)
  {
    using namespace aprl;
    const double dt=.01/SI::c;
    const EBGridGeometry global(8,6,32,.1,.2,.05,dt);
    const EBGridGeometry local(8,6,offset(rank+1,size)-offset(rank,size),.1,.2,.05,dt);
    FieldVector<Double> origin(0),localOrigin(0);localOrigin[2]=offset(rank,size)*.05;
    BoostFrameTransform frame;frame.setFromGamma(2,SI::c,0);
    DetectorConfig cfg;cfg.directory=directory;
    cfg.physicalGuardCells[0]=cfg.physicalGuardCells[1]=1;
    for(unsigned i=0;i<2;++i)
      {
        EnergyDetectorPlaneConfig p;p.name=i?"sparse":"dense";
        p.z=frame.labZFromBoxZT(16.125*.05,0);
        p.rhythm=(i?100:.25)*dt/frame.gamma();p.bufferSamples=2;p.compression=i;
        cfg.energyPlanes.push_back(p);
      }
    require(cfg.enabled()&&cfg.fieldPlanes.empty()&&cfg.particlePlanes.empty(),
            "Energy planes are not independent");
    LabDetectorManager manager(cfg,global,origin,localOrigin,frame,MPI_COMM_WORLD);
    require(manager.hasEnergyPlanes(),"Energy manager did not enable");
    EBFieldGrid fields(local);
    fill(fields,offset(rank,size),0);manager.sampleFieldPlanes(fields,0);
    RelativisticParticleSI before,after;
    before.mass=after.mass=SI::electronMass*1.e7;
    before.charge=after.charge=-SI::elementaryCharge*1.e7;
    before.weight=after.weight=1.e7;
    before.id=after.id=rank;
    before.position[0]=after.position[0]=.3;
    before.position[1]=after.position[1]=.4;
    const double plane=cfg.energyPlanes[0].z;
    before.position[2]=frame.boxZFromLabZAndBoxT(plane-.02,0);
    after.position[2]=frame.boxZFromLabZAndBoxT(plane+.02,dt);
    FieldVector<Double> u0(0),u1(0);u0[0]=5;u1[0]=-5;u0[2]=u1[2]=100;
    frame.properVelocityLabToBox(u0,before.properVelocity);
    frame.properVelocityLabToBox(u1,after.properVelocity);
    std::vector<EnergyPlaneCrossing> events;
    manager.previewEnergyParticleStep(before,after,0,dt,events);
    require(events.size()==2&&events[0].direction==1,"Missing forward substep crossings");
    const long double mc2=static_cast<long double>(before.mass)*SI::c*SI::c;
    const long double gamma=std::sqrt(10026.L);
    const long double expectedK=mc2*10025.L/(gamma+1.L);
    closeEnough(events[0].kineticLab,expectedK,"Affine lab energy changed in pure magnetic endpoint rotation");
    closeEnough(events[0].totalLab,frame.gamma()*events[0].totalBox+
      frame.gammaBeta()*events[0].cpzBox,"Crossing four-energy Lorentz identity failed");
    closeEnough(events[0].cpzLab,frame.gamma()*events[0].cpzBox+
      frame.gammaBeta()*events[0].totalBox,"Crossing momentum Lorentz identity failed");
    FieldVector<Double> midpoint(0),midLab(0);
    for(unsigned k=0;k<3;++k)midpoint[k]=.5*(before.properVelocity[k]+after.properVelocity[k]);
    frame.properVelocityBoxToLab(midpoint,midLab);
    const long double wrongK=mc2*(BoostFrameTransform::gammaFromProperVelocity(midLab)-1.L);
    require(std::abs(wrongK-expectedK)>1.e-4L*expectedK,
            "Legacy three-momentum interpolation negative control is insensitive");
    manager.commitEnergyParticleCrossings(events,-dt); // rejected after true terminal time
    manager.commitEnergyParticleCrossings(events,dt);
    // One downstream and one upstream crossing outside the physical aperture
    // cancel in signed transport but remain separately countable.
    before.position[0]=after.position[0]=.01;
    manager.captureEnergyParticleStep(before,after,0,dt);
    std::swap(before,after);
    before.position[2]=frame.boxZFromLabZAndBoxT(plane+.02,0);
    after.position[2]=frame.boxZFromLabZAndBoxT(plane-.02,dt);
    manager.captureEnergyParticleStep(before,after,0,dt);
    // Adjacent substeps landing exactly on the plane must count only once.
    before.position[2]=frame.boxZFromLabZAndBoxT(plane-.02,0);
    after.position[2]=frame.boxZFromLabZAndBoxT(plane,dt);
    events.clear();manager.previewEnergyParticleStep(before,after,0,dt,events);
    require(events.size()==2,"End-on-plane crossing absent");
    before=after;after.position[2]=frame.boxZFromLabZAndBoxT(plane+.02,2*dt);
    events.clear();manager.previewEnergyParticleStep(before,after,dt,2*dt,events);
    require(events.empty(),"Start-on-plane crossing duplicated");
    manager.collectParticleCrossings();
    for(unsigned n=1;n<=5;++n)
      {
        FieldVector<Double> translated=origin;
        if(n==2)translated[2]=100.; // one deliberately unavailable field step
        FieldVector<Double> translatedLocal=translated;translatedLocal[2]+=offset(rank,size)*.05;
        manager.updateGridOrigins(translated,translatedLocal);
        fill(fields,offset(rank,size),n);manager.sampleFieldPlanes(fields,n*dt);
      }
    manager.close(true);MPI_Barrier(MPI_COMM_WORLD);
    if(rank!=0)return;
    const long double area=8*.1L*6*.2L,physicalArea=6*.1L*4*.2L;
    long double expectedIntegral=0;
    const unsigned left[]={0,3,4};
    for(unsigned i=0;i<3;++i)
      expectedIntegral+=.5L*dt/frame.gamma()*(analyticPower(frame,plane,dt,left[i],area)+
                                           analyticPower(frame,plane,dt,left[i]+1,area));
    for(unsigned p=0;p<2;++p)
      {
        const std::string path=directory+"/"+(p?"sparse":"dense")+"-energy.h5";
        hid_t f=H5Fopen(path.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
        require(f>=0,"Cannot open energy output");
        hid_t g=H5Gopen2(f,"/energy_plane",H5P_DEFAULT);
        hid_t d=H5Dopen2(g,"committed_records",H5P_DEFAULT);
        std::uint64_t count=0;unsigned char complete=0;
        require(H5Dread(d,H5T_NATIVE_UINT64,H5S_ALL,H5S_ALL,H5P_DEFAULT,&count)>=0,"No commit marker");H5Dclose(d);
        d=H5Dopen2(g,"complete",H5P_DEFAULT);
        require(H5Dread(d,H5T_NATIVE_UCHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,&complete)>=0&&complete==1,"Not complete");H5Dclose(d);
        require(count==(p?2:6),"Output cadence/final flush mismatch");
        d=H5Dopen2(g,"records",H5P_DEFAULT);
        const auto samples=column<std::uint64_t>(d,"field_sample_count",H5T_NATIVE_UINT64,count);
        const auto outside=column<std::uint64_t>(d,"outside_field_window_steps",H5T_NATIVE_UINT64,count);
        const auto gaps=column<std::uint64_t>(d,"field_gap_count",H5T_NATIVE_UINT64,count);
        require(samples.back()==5&&outside.back()==1&&gaps.back()==1,"Field coverage counters differ");
        const auto energy=column<double>(d,"field_energy_lab_J",H5T_NATIVE_DOUBLE,count);
        closeEnough(energy.back(),expectedIntegral,"Field integral cadence/gap error");
        closeEnough(column<double>(d,"physical_field_energy_lab_J",H5T_NATIVE_DOUBLE,count).back(),
                    expectedIntegral*physicalArea/area,"Physical field aperture error");
        closeEnough(column<double>(d,"field_power_lab_W",H5T_NATIVE_DOUBLE,count).back(),
                    analyticPower(frame,plane,dt,5,area),"MPI plane sampling error");
        if(!p)
          {
            require(column<unsigned char>(d,"field_available",H5T_NATIVE_UCHAR,count)[2]==0,
                    "Missing unavailable field flag");
            require(std::isnan(column<double>(d,"field_power_lab_W",H5T_NATIVE_DOUBLE,count)[2]),
                    "Unavailable field incorrectly reported zero");
          }
        require(column<std::uint64_t>(d,"particle_downstream_count",H5T_NATIVE_UINT64,count).back()==2U*size&&
                column<std::uint64_t>(d,"particle_upstream_count",H5T_NATIVE_UINT64,count).back()==static_cast<unsigned>(size)&&
                column<std::uint64_t>(d,"physical_particle_downstream_count",H5T_NATIVE_UINT64,count).back()==static_cast<unsigned>(size)&&
                column<std::uint64_t>(d,"physical_particle_upstream_count",H5T_NATIVE_UINT64,count).back()==0,
                "Signed counts, clipping or aperture mismatch");
        closeEnough(particleEnergy(d,"particle_kinetic_lab_J",count),size*expectedK,"Wrong particle K or double macro weighting");
        closeEnough(particleEnergy(d,"physical_particle_kinetic_lab_J",count),size*expectedK,"Wrong physical particle K");
        closeEnough(particleEnergy(d,"particle_cpz_lab_J",count),size*mc2*100,"Wrong lab particle cPz");
        const long double totalLab=particleEnergy(d,"particle_total_lab_J",count);
        const long double totalBox=particleEnergy(d,"particle_total_box_J",count);
        const long double cpzBox=particleEnergy(d,"particle_cpz_box_J",count);
        closeEnough(totalLab,frame.gamma()*totalBox+frame.gammaBeta()*cpzBox,"Stored four-energy Lorentz identity mismatch");
        H5Dclose(d);H5Gclose(g);H5Fclose(f);
      }
    std::ofstream report((directory+"/summary.json").c_str());
    report<<std::setprecision(17)<<"{\"status\":\"PASS\",\"ranks\":"<<size
      <<",\"energy_planes_independent\":true,\"pure_B_endpoint_energy_invariant\":true"
      <<",\"legacy_momentum_negative_relative_bias\":"<<static_cast<double>((wrongK-expectedK)/expectedK)
      <<",\"terminal_clip_and_half_open_crossings\":true,\"macro_weight_once\":true"
      <<",\"signed_transport_and_physical_aperture\":true,\"field_step_integration_independent_of_write_rhythm\":true"
      <<",\"no_unavailable_gap_bridge\":true,\"scope\":\"analytic diagnostic fixtures, not full PIC conservation\"}\n";
    require(report.good(),"Cannot write summary");
    std::cout<<"energy planes PASS ranks="<<size<<'\n';
  }
}
int main(int argc,char** argv)
{
  MPI_Init(&argc,&argv);int rank=0,size=1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  try
    {
      require(argc==2&&size>=1&&size<=8,"Usage: aprl_energy_plane_checks NEW_DIRECTORY (1..8 ranks)");
      if(rank==0)require(::mkdir(argv[1],0775)==0,"Use a fresh output directory");
      MPI_Barrier(MPI_COMM_WORLD);run(argv[1],rank,size);MPI_Finalize();return 0;
    }
  catch(const std::exception& e)
    { std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,2);return 2; }
}
