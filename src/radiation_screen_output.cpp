#include "radiation_screen_output.h"
#include "runtime_util.h"
#include "eb_field.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace aprl
{
  namespace
  {
    const std::size_t batchSize=4096;

    RadiationScreenSpec timeScreen(RadiationScreenSpec spec)
    {
      if (!spec.nt)
        throw std::invalid_argument("radiation_screen requires an arrival-time window");
      // The optional frequency list is metadata for offline processing only.
      spec.frequenciesHz.clear();
      return spec;
    }

    void check(herr_t status, const char* operation)
    {
      if (status<0) throw std::runtime_error(std::string("radiation_screen: ")+operation);
    }

    struct Handle
    {
      hid_t id;
      herr_t (*release)(hid_t);
      Handle(hid_t value, herr_t (*closer)(hid_t)) : id(value), release(closer)
      {
        if (id<0) throw std::runtime_error("radiation_screen: cannot open HDF5 object");
      }
      ~Handle() { release(id); }
      Handle(const Handle&)=delete;
      Handle& operator=(const Handle&)=delete;
    };

    void scalar(hid_t file, const char* name, hid_t type, const void* value)
    {
      Handle space(H5Screate(H5S_SCALAR),H5Sclose);
      Handle attribute(H5Acreate2(file,name,type,space.id,H5P_DEFAULT,H5P_DEFAULT),H5Aclose);
      check(H5Awrite(attribute.id,type,value),name);
    }

    void label(hid_t file, const char* name, const std::string& text)
    {
      Handle type(H5Tcopy(H5T_C_S1),H5Tclose);
      check(H5Tset_size(type.id,text.size()+1),name);
      scalar(file,name,type.id,text.c_str());
    }

    void axis(hid_t file, const char* name, const std::vector<double>& values)
    {
      hsize_t count=values.size();
      Handle space(H5Screate_simple(1,&count,NULL),H5Sclose);
      Handle data(H5Dcreate2(file,name,H5T_IEEE_F64LE,space.id,H5P_DEFAULT,
                             H5P_DEFAULT,H5P_DEFAULT),H5Dclose);
      check(H5Dwrite(data.id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,values.data()),name);
    }

    void linearAxis(hid_t file, const char* name, std::size_t count, double first, double spacing)
    {
      std::vector<double> values(count);
      for (std::size_t i=0;i<count;++i) values[i]=std::fma(double(i),spacing,first);
      axis(file,name,values);
    }

    void field(hid_t file, const RadiationScreenSpec& spec, const char* name,
               const std::vector<double>& samples, bool magnetic)
    {
      hsize_t dims[]={spec.nt,spec.ny,spec.nx,3};
      Handle space(H5Screate_simple(4,dims,NULL),H5Sclose);
      Handle data(H5Dcreate2(file,name,H5T_IEEE_F64LE,space.id,H5P_DEFAULT,
                             H5P_DEFAULT,H5P_DEFAULT),H5Dclose);
      label(data.id,"units",magnetic ? "T" : "V/m");
      const double conversion=magnetic ? 1./SI::c : 1.;
      const std::size_t componentOffset=magnetic ? 3 : 0;
      hsize_t slab[]={1,spec.ny,spec.nx,3};
      Handle memory(H5Screate_simple(4,slab,NULL),H5Sclose);
      std::vector<double> plane(spec.nx*spec.ny*3);
      for (std::size_t k=0;k<spec.nt;++k)
        {
          for (std::size_t p=0;p<spec.nx*spec.ny;++p)
            {
              const std::size_t input=(p*spec.nt+k)*6+componentOffset;
              for (std::size_t c=0;c<3;++c)
                plane[p*3+c]=samples[input+c]*conversion;
            }
          hsize_t offset[]={k,0,0,0};
          check(H5Sselect_hyperslab(space.id,H5S_SELECT_SET,offset,NULL,slab,NULL),name);
          check(H5Dwrite(data.id,H5T_NATIVE_DOUBLE,memory.id,space.id,H5P_DEFAULT,plane.data()),name);
        }
    }

    void reduce(std::vector<double>& samples, MPI_Comm comm, int rank)
    {
      for (std::size_t i=0;i<samples.size();)
        {
          const int count=static_cast<int>(std::min<std::size_t>(samples.size()-i,1048576));
          MPI_Reduce(rank==0 ? MPI_IN_PLACE : samples.data()+i,
                     rank==0 ? samples.data()+i : NULL,count,MPI_DOUBLE,MPI_SUM,0,comm);
          i+=count;
        }
    }

    // All ranks must leave collective initialization/finalization together.
    void rootStatus(std::string& error, int rank, MPI_Comm comm)
    {
      int length=static_cast<int>(error.size());
      MPI_Bcast(&length,1,MPI_INT,0,comm);
      if (!length) return;
      if (rank!=0) error.resize(length);
      MPI_Bcast(&error[0],length,MPI_CHAR,0,comm);
      throw std::runtime_error(error);
    }
  }

  RadiationScreenOutput::RadiationScreenOutput(const RadiationScreenSpec& spec,
      const std::string& path, bool overwrite, MPI_Comm comm):
    screen_(timeScreen(spec)), communicator_(comm), rank_(0), file_(-1)
  {
    MPI_Comm_rank(comm,&rank_);
    const auto slash=path.find_last_of('/');
    if (slash!=std::string::npos)
      createDirectories(slash==0 ? "/" : path.substr(0,slash),comm);
    std::string error;
    if (rank_==0)
      try
        {
          file_=H5Fcreate(path.c_str(),overwrite ? H5F_ACC_TRUNC : H5F_ACC_EXCL,
                          H5P_DEFAULT,H5P_DEFAULT);
          if (file_<0) throw std::runtime_error("Cannot create radiation screen output: "+path);
          const int incomplete=0;
          scalar(file_,"complete",H5T_NATIVE_INT,&incomplete);
          if (!spec.frequenciesHz.empty())
            axis(file_,"requested_frequency_hz",spec.frequenciesHz);
          check(H5Fflush(file_,H5F_SCOPE_LOCAL),"write initial status");
        }
      catch (const std::exception& e)
        {
          if (file_>=0) H5Fclose(file_);
          file_=-1;
          error=e.what();
        }
    rootStatus(error,rank_,comm);
  }

  RadiationScreenOutput::~RadiationScreenOutput()
  {
    if (file_>=0) H5Fclose(file_);
  }

  void RadiationScreenOutput::append(const std::vector<RadiationSegment>& boxSubsteps,
      double boxStop, const BoostFrameTransform& frame)
  {
    for (const auto& box : boxSubsteps)
      {
        if (box.t0>=boxStop) break;
        ++committed_;
        const bool partial=boxStop<box.t1;
        if (partial) ++clipped_;
        if (box.u0==box.u1 || box.charge==0.) continue;
        RadiationSegment lab=box;
        FieldVector<double> u(0.),uLab(0.);
        if (partial)
          {
            // The particle manager locates a terminal hit within the field step.
            // Use only the affected Boris substep's linear endpoint prefix.
            const double fraction=(boxStop-box.t0)/(box.t1-box.t0);
            for (unsigned a=0;a<3;++a)
              {
                lab.x1[a]=std::fma(fraction,box.x1[a]-box.x0[a],box.x0[a]);
                lab.u1[a]=std::fma(fraction,box.u1[a]-box.u0[a],box.u0[a]);
              }
            lab.t1=boxStop;
          }
        for (int end=0;end<2;++end)
          {
            auto& x=end ? lab.x1 : lab.x0;
            auto& proper=end ? lab.u1 : lab.u0;
            double& t=end ? lab.t1 : lab.t0;
            const double tBox=t,zBox=x[2];
            frame.boxToLab(tBox,zBox,t,x[2]);
            for (unsigned a=0;a<3;++a) u[a]=proper[a];
            frame.properVelocityBoxToLab(u,uLab);
            for (unsigned a=0;a<3;++a) proper[a]=uLab[a];
          }
        if (lab.u0==lab.u1) continue;
        ++accelerated_;
        batch_.push_back(lab);
        if (batch_.size()==batchSize) flush();
      }
  }

  void RadiationScreenOutput::flush()
  {
    if (batch_.empty()) return;
    const double started=MPI_Wtime();
    screen_.accumulate(batch_);
    accumulationSeconds_+=MPI_Wtime()-started;
    batch_.clear();
  }

  void RadiationScreenOutput::close(bool complete)
  {
    flush();
    const unsigned long long storage=screen_.storageBytes();
    const unsigned long long local[]={committed_,accelerated_,clipped_};
    unsigned long long totals[3]={};
    MPI_Reduce(local,totals,3,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,communicator_);
    double maxSeconds=0.;
    MPI_Reduce(&accumulationSeconds_,&maxSeconds,1,MPI_DOUBLE,MPI_MAX,0,communicator_);
    auto time=screen_.releaseTimeField();
    reduce(time,communicator_,rank_);
    std::string error;
    if (rank_==0)
      try
        {
          const auto& spec=screen_.spec();
          linearAxis(file_,"x_m",spec.nx,spec.x0,spec.dx);
          linearAxis(file_,"y_m",spec.ny,spec.y0,spec.dy);
          scalar(file_,"z_m",H5T_NATIVE_DOUBLE,&spec.z);
          scalar(file_,"time_origin_s",H5T_NATIVE_DOUBLE,&spec.phaseOrigin);
          label(file_,"model","finite-distance vacuum Lienard-Wiechert acceleration term");
          label(file_,"frame","laboratory; xyz Cartesian components, pixel centres");
          linearAxis(file_,"time_s",spec.nt,spec.timeStart+.5*spec.timeStep,spec.timeStep);
          scalar(file_,"time_start_s",H5T_NATIVE_DOUBLE,&spec.timeStart);
          scalar(file_,"time_step_s",H5T_NATIVE_DOUBLE,&spec.timeStep);
          label(file_,"time_order","arrival time bin centre,y,x,component; mean over bin");
          field(file_,spec,"E_time",time,false);
          field(file_,spec,"B_time",time,true);
          scalar(file_,"committed_substeps",H5T_NATIVE_ULLONG,&totals[0]);
          scalar(file_,"accelerated_segments",H5T_NATIVE_ULLONG,&totals[1]);
          scalar(file_,"clipped_substeps",H5T_NATIVE_ULLONG,&totals[2]);
          scalar(file_,"max_accumulation_s",H5T_NATIVE_DOUBLE,&maxSeconds);
          scalar(file_,"storage_bytes_per_rank",H5T_NATIVE_ULLONG,&storage);
          const int done=complete ? 1 : 0;
          {
            Handle status(H5Aopen(file_,"complete",H5P_DEFAULT),H5Aclose);
            check(H5Awrite(status.id,H5T_NATIVE_INT,&done),"complete");
          }
          check(H5Fclose(file_),"close");
          file_=-1;
        }
      catch (const std::exception& e) {error=e.what();}
    rootStatus(error,rank_,communicator_);
  }
}
