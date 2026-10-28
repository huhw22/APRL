/* Fixed laboratory diagnostics on translated MPI slabs. These affine native
 * E/B data are an interpolation fixture, not a Maxwell/PIC simulation. */
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
  using aprl::Double;
  const std::size_t NX=8, NY=6, NZ=95, STEPS=320;
  const double DZ=1.2543262415359034e-6, GAMMA=929.1305492858543;
  const double TOLERANCE=2.e-10;
  const double STAGGER[6][3]={{.5,0,0},{0,.5,0},{0,0,.5},
    {0,.5,.5},{.5,0,.5},{.5,.5,0}};
  struct Point { double e[3],b[3]; };
  struct Output { std::vector<double> time; std::vector<Point> fields; };
  static_assert(sizeof(Point)==6*sizeof(double),"Unexpected field record padding");

  void require(bool ok,const std::string& message)
  { if(!ok) throw std::runtime_error(message); }

  std::size_t offset(int rank,int size)
  { return (NZ/size)*rank+std::min<std::size_t>(rank,NZ%size); }

  aprl::FieldVector<Double> initialOrigin()
  {
    aprl::FieldVector<Double> out(0.);
    out[0]=-2.7e-5;out[1]=1.9e-5;out[2]=.0314159;
    return out;
  }

  aprl::FieldVector<Double> translatedOrigin(std::size_t step)
  {
    auto out=initialOrigin();out[2]-=static_cast<double>(step/16)*DZ;
    return out;
  }

  aprl::FieldVector<Double> slabOrigin(const aprl::FieldVector<Double>& global,
                                     int rank,int size)
  { auto out=global;out[2]+=offset(rank,size)*DZ;return out; }

  double affine(unsigned component,double x,double y,double z,double step)
  {
    const double value=1.+.19*component+.013*x+.021*y+.017*z+
      .0007*(1.+.2*component)*step;
    return component<3?value:value/aprl::SI::c;
  }

  void fill(aprl::EBFieldGrid& fields,std::size_t first,
            const aprl::FieldVector<Double>& origin,std::size_t step)
  {
    aprl::YeeComponent* components[]={&fields.ex(),&fields.ey(),&fields.ez(),
      &fields.bx(),&fields.by(),&fields.bz()};
    const double translation=(origin[2]-initialOrigin()[2])/DZ;
    for(unsigned a=0;a<6;++a)
      for(std::size_t k=0;k<components[a]->nz();++k)
        for(std::size_t j=0;j<components[a]->ny();++j)
          for(std::size_t i=0;i<components[a]->nx();++i)
            (*components[a])(i,j,k)=affine(a,i+STAGGER[a][0],j+STAGGER[a][1],
              translation+first+k+STAGGER[a][2],step-(a>=3?.5:0.));
  }

  std::vector<Point> analytic(const aprl::BoostFrameTransform& frame,
      const aprl::FieldVector<Double>& gridOrigin,double boxZ,
      std::size_t step,bool wrongBTime=false)
  {
    std::vector<Point> out(NX*NY);
    const double coordinate=(boxZ-gridOrigin[2])/DZ;
    const double translation=(gridOrigin[2]-initialOrigin()[2])/DZ;
    for(std::size_t j=0;j<NY;++j)
      for(std::size_t i=0;i<NX;++i)
        {
          aprl::FieldVector<Double> e(0.),b(0.),el(0.),bl(0.);
          for(unsigned a=0;a<6;++a)
            {
              // Match only the physical global clamp, never a local MPI edge.
              const double low=STAGGER[a][2],high=NZ-low;
              const double z=translation+std::max(low,std::min(high,coordinate));
              const double value=affine(a,i+.5,j+.5,z,
                step-(a>=3&&!wrongBTime?.5:0.));
              (a<3?e[a]:b[a-3])=value;
            }
          frame.fieldsBoxToLab(e,b,el,bl);
          for(unsigned a=0;a<3;++a)
            { out[j*NX+i].e[a]=el[a];out[j*NX+i].b[a]=bl[a]; }
        }
    return out;
  }

  double difference(const Point* a,const Point* b,std::size_t count)
  {
    double maximum=0.;
    for(std::size_t p=0;p<count;++p)for(unsigned k=0;k<6;++k)
      {
        const double x=k<3?a[p].e[k]:aprl::SI::c*a[p].b[k-3];
        const double y=k<3?b[p].e[k]:aprl::SI::c*b[p].b[k-3];
        require(std::isfinite(x)&&std::isfinite(y),"Nonfinite detector field");
        maximum=std::max(maximum,std::abs(x-y)/std::max(1.,std::abs(y)));
      }
    return maximum;
  }

  double attribute(hid_t group,const char* name)
  {
    const hid_t handle=H5Aopen(group,name,H5P_DEFAULT);double value=0.;
    require(handle>=0&&H5Aread(handle,H5T_NATIVE_DOUBLE,&value)>=0,
      "Cannot read detector coordinate attribute");H5Aclose(handle);return value;
  }

  Output read(const std::string& filename,const aprl::EBGridGeometry& geometry)
  {
    const hid_t file=H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
    require(file>=0,"Cannot open "+filename);
    const hid_t group=H5Gopen2(file,"/field_plane",H5P_DEFAULT);
    require(group>=0,"Missing detector group");
    require(attribute(group,"x_first_m")==initialOrigin()[0]+.5*geometry.dx&&
            attribute(group,"y_first_m")==initialOrigin()[1]+.5*geometry.dy,
      "Laboratory aperture changed while moving the grid");
    std::uint64_t committed=0;unsigned char complete=0;
    hid_t d=H5Dopen2(group,"committed_samples",H5P_DEFAULT);
    require(d>=0&&H5Dread(d,H5T_NATIVE_UINT64,H5S_ALL,H5S_ALL,H5P_DEFAULT,&committed)>=0,
      "Cannot read detector commit count");H5Dclose(d);
    d=H5Dopen2(group,"complete",H5P_DEFAULT);
    require(d>=0&&H5Dread(d,H5T_NATIVE_UCHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,&complete)>=0,
      "Cannot read detector completion");H5Dclose(d);
    require(complete==1,"Detector did not complete");
    Output out;out.time.resize(committed);out.fields.resize(committed*NX*NY);
    d=H5Dopen2(group,"time_s",H5P_DEFAULT);require(d>=0,"Missing times");
    hid_t space=H5Dget_space(d);hsize_t shape[3]={0,0,0};
    require(H5Sget_simple_extent_ndims(space)==1&&
      H5Sget_simple_extent_dims(space,shape,NULL)>=0&&shape[0]==committed,"Time dimensions differ");
    H5Sclose(space);
    if(committed)require(H5Dread(d,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,out.time.data())>=0,
      "Cannot read times");H5Dclose(d);
    d=H5Dopen2(group,"fields",H5P_DEFAULT);require(d>=0,"Missing fields");space=H5Dget_space(d);
    require(H5Sget_simple_extent_ndims(space)==3&&H5Sget_simple_extent_dims(space,shape,NULL)>=0&&
      shape[0]==committed&&shape[1]==NY&&shape[2]==NX,"Field dimensions differ");H5Sclose(space);
    const hid_t type=H5Tcreate(H5T_COMPOUND,sizeof(Point));const hsize_t three=3;
    const hid_t vector=H5Tarray_create2(H5T_NATIVE_DOUBLE,1,&three);
    require(type>=0&&vector>=0&&H5Tinsert(type,"electric_V_per_m",HOFFSET(Point,e),vector)>=0&&
      H5Tinsert(type,"magnetic_T",HOFFSET(Point,b),vector)>=0,"Cannot create field read type");
    if(committed)require(H5Dread(d,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,out.fields.data())>=0,
      "Cannot read fields");
    H5Tclose(vector);H5Tclose(type);H5Dclose(d);H5Gclose(group);H5Fclose(file);return out;
  }

  aprl::DetectorConfig configuration(const std::string& directory,
      const aprl::BoostFrameTransform& frame,double dt)
  {
    aprl::DetectorConfig out;out.directory=directory;
    const char* names[]={"sweep","sparse","upper","lower"};
    const double positions[]={83.375,83.375,NZ+.375,3.375};
    for(unsigned p=0;p<4;++p)
      {
        aprl::FieldDetectorPlaneConfig plane;plane.name=names[p];
        plane.z=frame.labZFromBoxZT(initialOrigin()[2]+positions[p]*DZ,0.);
        plane.rhythm=(p==1?2.25:.5)*dt/GAMMA;plane.bufferSamples=7;plane.compression=0;
        out.fieldPlanes.push_back(plane);
      }
    return out;
  }

  int owner(double coordinate,int size)
  {
    for(int rank=0;rank<size;++rank)
      if(coordinate<offset(rank+1,size))return rank;
    return size-1;
  }

  std::vector<std::size_t> expectedSteps(const aprl::FieldDetectorPlaneConfig& plane,
      const aprl::BoostFrameTransform& frame,double dt,bool moving,unsigned index)
  {
    std::vector<std::size_t> out;
    for(std::size_t n=0;n<STEPS;++n)
      {
        if(index==1)
          {
            // Independent exact cadence construction; no copy of manager's
            // accumulated-next-time scheduling or its origin update logic.
            bool selected=false;
            for(std::size_t k=0;k<STEPS&&std::ceil(2.25*k)<=n;++k)
              if(std::ceil(2.25*k)==n){selected=true;break;}
            if(!selected)continue;
          }
        const auto origin=moving?translatedOrigin(n):initialOrigin();
        const double q=(frame.boxZFromLabZAndBoxT(plane.z,n*dt)-origin[2])/DZ;
        if(q>=0.&&q<=NZ)out.push_back(n);
      }
    return out;
  }

  void run(const std::string& directory,int rank,int size)
  {
    const auto initial=initialOrigin();
    aprl::BoostFrameTransform frame;frame.setFromGamma(GAMMA,aprl::SI::c,initial[2]);
    const double dt=.25*DZ/(frame.beta()*aprl::SI::c);
    aprl::EBGridGeometry global(NX,NY,NZ,5e-6,5e-6,DZ,dt),slab(global);
    slab.nz=offset(rank+1,size)-offset(rank,size);
    aprl::EBFieldGrid fixedFields(slab),movingFields(slab);
    const auto fixedConfig=configuration(directory+"/fixed",frame,dt);
    const auto movingConfig=configuration(directory+"/moving",frame,dt);
    const auto staleConfig=configuration(directory+"/stale-origin-control",frame,dt);
    aprl::LabDetectorManager fixed(fixedConfig,global,initial,slabOrigin(initial,rank,size),frame,MPI_COMM_WORLD);
    aprl::LabDetectorManager moving(movingConfig,global,initial,slabOrigin(initial,rank,size),frame,MPI_COMM_WORLD);
    aprl::LabDetectorManager stale(staleConfig,global,initial,slabOrigin(initial,rank,size),frame,MPI_COMM_WORLD);
    unsigned rejected=0;
    for(unsigned invalid=0;invalid<3;++invalid)
      {
        auto g=initial,l=slabOrigin(initial,rank,size);
        if(invalid==0)g[2]=std::numeric_limits<double>::quiet_NaN();
        if(invalid==1){g[0]+=global.dx;l[0]+=global.dx;}
        if(invalid==2)l[2]+=DZ;
        try{moving.updateGridOrigins(g,l);}
        catch(const std::invalid_argument&){++rejected;}
      }
    require(rejected==3,"Invalid translated grid was accepted");
    unsigned fixedChanges=0,movingChanges=0;int previousFixed=-1,previousMoving=-1;
    for(std::size_t n=0;n<STEPS;++n)
      {
        const auto shifted=translatedOrigin(n);
        fill(fixedFields,offset(rank,size),initial,n);
        fill(movingFields,offset(rank,size),shifted,n);
        moving.updateGridOrigins(shifted,slabOrigin(shifted,rank,size));
        moving.updateGridOrigins(shifted,slabOrigin(shifted,rank,size));
        fixed.sampleFieldPlanes(fixedFields,n*dt);
        moving.sampleFieldPlanes(movingFields,n*dt);
        moving.sampleFieldPlanes(movingFields,n*dt); // A move must not reset cadence.
        stale.sampleFieldPlanes(movingFields,n*dt); // Deliberately omit the update.
        const double z=frame.boxZFromLabZAndBoxT(fixedConfig.fieldPlanes[0].z,n*dt);
        const int f=owner((z-initial[2])/DZ,size),m=owner((z-shifted[2])/DZ,size);
        if(previousFixed>=0&&f!=previousFixed)++fixedChanges;
        if(previousMoving>=0&&m!=previousMoving)++movingChanges;
        previousFixed=f;previousMoving=m;
      }
    fixed.close(true);moving.close(true);stale.close(true);
    bool closedRejected=false;
    try{moving.updateGridOrigins(initial,slabOrigin(initial,rank,size));}
    catch(const std::logic_error&){closedRejected=true;}
    require(closedRejected,"Closed detector accepted an origin update");
    if(rank!=0)return;
    require(size==1||(fixedChanges>0&&movingChanges>0),"Fixture did not exercise MPI ownership changes");
    double fixedError=0.,movingError=0.,pairError=0.,staleError=0.,wrongBError=0.;
    std::size_t compared=0,lowerMoving=0,lowerFixed=0;
    for(unsigned p=0;p<fixedConfig.fieldPlanes.size();++p)
      {
        const auto& plane=fixedConfig.fieldPlanes[p];
        const Output f=read(fixedConfig.directory+"/"+plane.name+".h5",global);
        const Output m=read(movingConfig.directory+"/"+plane.name+".h5",global);
        for(unsigned mode=0;mode<2;++mode)
          {
            const auto steps=expectedSteps(plane,frame,dt,mode!=0,p);const Output& output=mode?m:f;
            require(output.time.size()==steps.size(),"Cadence or window coverage changed: "+plane.name);
            for(std::size_t row=0;row<steps.size();++row)
              {
                const std::size_t n=steps[row];
                require(output.time[row]==frame.labTFromLabZBoxT(plane.z,n*dt),"Laboratory event time changed");
                const auto origin=mode?translatedOrigin(n):initial;
                const double z=frame.boxZFromLabZAndBoxT(plane.z,n*dt);
                const auto expected=analytic(frame,origin,z,n);
                const double error=difference(output.fields.data()+row*NX*NY,expected.data(),NX*NY);
                if(mode)movingError=std::max(movingError,error);else fixedError=std::max(fixedError,error);
                compared+=NX*NY*6;
                if(mode&&p==0)
                  {
                    const auto wrong=analytic(frame,origin,z,n,true);
                    wrongBError=std::max(wrongBError,difference(expected.data(),wrong.data(),NX*NY));
                  }
              }
          }
        if(p<2)
          {
            require(f.time==m.time,"Moving grid changed common-plane cadence");
            pairError=std::max(pairError,difference(m.fields.data(),f.fields.data(),f.fields.size()));
          }
        if(p==0)
          {
            const auto s=read(staleConfig.directory+"/"+plane.name+".h5",global);
            require(s.time==f.time,"Stale-origin control changed cadence");
            staleError=difference(s.fields.data(),f.fields.data(),f.fields.size());
          }
        if(p==3){lowerMoving=m.time.size();lowerFixed=f.time.size();}
      }
    require(lowerMoving>lowerFixed,"Lower-plane re-entry control was not exercised");
    require(fixedError<TOLERANCE&&movingError<TOLERANCE&&pairError<TOLERANCE,
      "Moving detector field interpolation differs from fixed-grid/analytic reference");
    require(staleError>1.e-4&&wrongBError>1.e-6,"Origin/time negative controls are insensitive");
    std::ofstream report((directory+"/summary.json").c_str());report<<std::setprecision(17)
      <<"{\"status\":\"PASS\",\"ranks\":"<<size<<",\"global_nz\":"<<NZ
      <<",\"samples\":"<<STEPS<<",\"components_compared\":"<<compared
      <<",\"relative_tolerance\":"<<TOLERANCE<<",\"fixed_analytic_error\":"<<fixedError
      <<",\"moving_analytic_error\":"<<movingError<<",\"moving_fixed_error\":"<<pairError
      <<",\"stale_origin_negative_error\":"<<staleError<<",\"wrong_B_time_negative_error\":"<<wrongBError
      <<",\"fixed_owner_changes\":"<<fixedChanges<<",\"moving_owner_changes\":"<<movingChanges
      <<",\"lower_fixed_samples\":"<<lowerFixed<<",\"lower_moving_samples\":"<<lowerMoving
      <<",\"invalid_origin_controls\":3,\"closed_update_rejected\":true,\"cadence_and_lab_times_identical\":true"
      <<",\"scope\":\"Analytic native staggered E/B on translated grids; original global clamp and MPI neighbour sampling; no Maxwell/PIC or checkpoint acceptance.\"}\n";
    require(report.good(),"Cannot write moving detector report");
    std::cout<<"moving detectors PASS ranks="<<size<<" fixed="<<fixedError
      <<" moving="<<movingError<<" paired="<<pairError<<'\n';
  }
}

int main(int argc,char** argv)
{
  MPI_Init(&argc,&argv);int rank=0,size=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
  try
    {
      require(argc==2&&size>=1&&size<=8,
        "Usage: aprl_moving_detector_checks NEW_OUTPUT_DIRECTORY (1..8 ranks)");
      if(rank==0)require(::mkdir(argv[1],0775)==0,"Use a fresh output directory with an existing parent");
      MPI_Barrier(MPI_COMM_WORLD);run(argv[1],rank,size);MPI_Finalize();return 0;
    }
  catch(const std::exception& e)
    { std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,2);return 2; }
}
