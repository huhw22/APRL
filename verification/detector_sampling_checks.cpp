/* Diagnostic interpolation fixture, not a solution of Maxwell's equations.
 * The unsplit Yee grid is the reference; no filtering is applied to outputs. */
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <mpi.h>
#include <hdf5.h>
#include "lab_detectors.h"
#include "radiation_screen_output.h"

namespace
{
  using aprl::Double;
  const std::size_t NX = 6, NY = 5, SAMPLES = 1016;
  const unsigned TONE = 37;
  const double GAMMA = 929.1305492858543;
  const double DZ = 1.2543262415359034e-6;
  const double PI = std::acos(-1.0);
  struct Point { double e[3]; double b[3]; };
  static_assert(sizeof(Point) == 6 * sizeof(double), "MPI fixture point must be packed doubles");
  struct Output { std::vector<double> time; std::vector<Point> field; };
  struct Errors { double boundary = 0.0, sweep = 0.0, sparse = 0.0;
                  double legacy = 0.0; std::size_t boundaryPlanes = 0; };

  void require(bool ok, const std::string& message)
  { if (!ok) throw std::runtime_error(message); }

  void mkdirChecked(const std::string& path)
  { require(::mkdir(path.c_str(), 0775) == 0 || errno == EEXIST,
      "Cannot create fixture output directory: " + path); }

  void checkRadiationScreenOutput(const std::string& directory, int rank)
  {
    aprl::RadiationScreenSpec spec;
    spec.nx=spec.ny=1; spec.dx=spec.dy=1.; spec.z=.5;
    spec.nt=4; spec.timeStep=1.e-15; spec.frequenciesHz={0.,1.e12};
    const std::string path=directory+"/radiation-waveform.h5";
    aprl::RadiationScreenOutput writer(spec,path,true,MPI_COMM_WORLD);
    writer.close(true);
    if (rank!=0) return;
    hid_t file=H5Fopen(path.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
    require(file>=0,"Cannot read the radiation waveform");
    for (const char* name:{"time_s","E_time","B_time","requested_frequency_hz"})
      require(H5Lexists(file,name,H5P_DEFAULT)>0,"Missing waveform item: "+std::string(name));
    for (const char* name:{"frequency_hz","E_spectrum","B_spectrum"})
      require(H5Lexists(file,name,H5P_DEFAULT)==0,"Unexpected online spectrum: "+std::string(name));
    hid_t hints=H5Dopen2(file,"requested_frequency_hz",H5P_DEFAULT);
    hid_t space=H5Dget_space(hints);
    require(H5Sget_simple_extent_ndims(space)==1 && H5Sget_simple_extent_npoints(space)==2,
            "Offline frequency list changed shape");
    double frequencies[2];
    require(H5Dread(hints,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,frequencies)>=0 &&
            frequencies[0]==0. && frequencies[1]==1.e12,"Offline frequency list lost its order or units");
    H5Sclose(space); H5Dclose(hints); H5Fclose(file);
  }

  std::size_t offset(std::size_t nz, int rank, int size)
  { return (nz / size) * rank + std::min<std::size_t>(rank, nz % size); }

  std::size_t count(std::size_t nz, int rank, int size)
  { return nz / size + (static_cast<std::size_t>(rank) < nz % size ? 1 : 0); }

  int owner(double q, std::size_t nz, int size)
  {
    for (int rank = 0; rank < size; ++rank)
      if (q < offset(nz, rank, size) + count(nz, rank, size)) return rank;
    return size - 1;
  }

  aprl::EBGridGeometry geometry(std::size_t nz)
  { return aprl::EBGridGeometry(NX, NY, nz, 4.0 * DZ, 5.0 * DZ,
      DZ, 0.2 * DZ / aprl::SI::c); }

  aprl::FieldVector<Double> origin()
  {
    aprl::FieldVector<Double> value(0.0);
    value[0] = -2.7e-5; value[1] = 1.9e-5; value[2] = 0.0314159;
    return value;
  }

  void fill(aprl::EBFieldGrid& fields, std::size_t first,
            bool sweep, double q = 0.0, std::size_t sample = 0)
  {
    aprl::YeeComponent* arrays[] = { &fields.ex(), &fields.ey(), &fields.ez(),
      &fields.bx(), &fields.by(), &fields.bz() };
    const double shifts[6][3] = {{.5,0,0}, {0,.5,0}, {0,0,.5},
      {0,.5,.5}, {.5,0,.5}, {.5,.5,0}};
    for (unsigned a = 0; a < 6; ++a)
      for (std::size_t k = 0; k < arrays[a]->nz(); ++k)
        for (std::size_t j = 0; j < arrays[a]->ny(); ++j)
          for (std::size_t i = 0; i < arrays[a]->nx(); ++i)
            {
              const double x = i + shifts[a][0], y = j + shifts[a][1];
              const double z = k + first + shifts[a][2];
              double value = 1.0 + .17 * a + .013 * x + .021 * y;
              if (sweep)
                value += .05 * (1.0 + .1 * a) * (z - q) +
                  1.0e-5 * (1.0 + .2 * a) * std::cos(
                    2.0 * PI * TONE * sample / SAMPLES + .12 * x + .17 * y + .2 * a);
              else
                value += .019 * z + .03 * std::sin(.4 * x - .2 * y + .15 * z + .3 * a);
              (*arrays[a])(i,j,k) = value / (a < 3 ? 1.0 : aprl::SI::c);
            }
  }

  std::vector<Point> sample(const aprl::EBFieldGrid& fields,
      const aprl::FieldVector<Double>& localOrigin,
      const aprl::BoostFrameTransform& frame, double z)
  {
    std::vector<Point> result(NX * NY);
    aprl::FieldVector<Double> p(origin()); p[2] = z;
    for (std::size_t j = 0; j < NY; ++j)
      for (std::size_t i = 0; i < NX; ++i)
        {
          p[0] = origin()[0] + (i + .5) * geometry(2).dx;
          p[1] = origin()[1] + (j + .5) * geometry(2).dy;
          aprl::FieldVector<Double> e(0.0), b(0.0);
          fields.sampleFieldsPositionLab(p, localOrigin, frame, e, b);
          for (unsigned a = 0; a < 3; ++a)
            { result[j * NX + i].e[a] = e[a]; result[j * NX + i].b[a] = b[a]; }
        }
    return result;
  }

  double scaled(const Point& value, unsigned a)
  { return a < 3 ? value.e[a] : aprl::SI::c * value.b[a - 3]; }

  double difference(const std::vector<Point>& found, const std::vector<Point>& expected)
  {
    require(found.size() == expected.size(), "Field record count mismatch");
    double error = 0.0;
    for (std::size_t p = 0; p < found.size(); ++p)
      for (unsigned a = 0; a < 6; ++a)
        {
          const double value = scaled(found[p], a), reference = scaled(expected[p], a);
          require(std::isfinite(value) && std::isfinite(reference), "Nonfinite field");
          error = std::max(error, std::abs(value-reference) / std::max(1.0, std::abs(reference)));
        }
    return error;
  }

  Output read(const std::string& path)
  {
    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    require(file >= 0, "Cannot read detector fixture: " + path);
    hid_t group = H5Gopen2(file, "/field_plane", H5P_DEFAULT);
    require(group >= 0, "Missing field_plane group");
    std::uint64_t committed = 0;
    unsigned char complete = 0;
    hid_t d = H5Dopen2(group, "committed_samples", H5P_DEFAULT);
    require(d >= 0 && H5Dread(d,H5T_NATIVE_UINT64,H5S_ALL,H5S_ALL,H5P_DEFAULT,&committed)>=0,
      "Cannot read committed_samples"); H5Dclose(d);
    d = H5Dopen2(group, "complete", H5P_DEFAULT);
    require(d >= 0 && H5Dread(d,H5T_NATIVE_UCHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,&complete)>=0,
      "Cannot read complete"); H5Dclose(d);
    require(complete == 1, "Detector output not complete");
    Output result; result.time.resize(committed); result.field.resize(committed * NX * NY);
    d = H5Dopen2(group, "time_s", H5P_DEFAULT);
    require(d >= 0, "Cannot open time dataset");
    if (committed) require(H5Dread(d,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,
      result.time.data())>=0, "Cannot read time dataset");
    H5Dclose(d);
    d = H5Dopen2(group, "fields", H5P_DEFAULT);
    require(d >= 0, "Cannot open fields dataset");
    hid_t space = H5Dget_space(d); hsize_t shape[3] = {0,0,0};
    require(H5Sget_simple_extent_ndims(space)==3 &&
      H5Sget_simple_extent_dims(space,shape,NULL)>=0 && shape[0]==committed &&
      shape[1]==NY && shape[2]==NX, "Unexpected field dataset dimensions");
    H5Sclose(space);
    hid_t type = H5Tcreate(H5T_COMPOUND,sizeof(Point)); hsize_t n=3;
    hid_t vector = H5Tarray_create2(H5T_NATIVE_DOUBLE,1,&n);
    require(type>=0 && vector>=0 &&
      H5Tinsert(type,"electric_V_per_m",HOFFSET(Point,e),vector)>=0 &&
      H5Tinsert(type,"magnetic_T",HOFFSET(Point,b),vector)>=0, "Cannot create field memory type");
    if (committed) require(H5Dread(d,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,
      result.field.data())>=0, "Cannot read field dataset");
    H5Tclose(vector); H5Tclose(type); H5Dclose(d); H5Gclose(group); H5Fclose(file);
    return result;
  }

  aprl::FieldDetectorPlaneConfig plane(const std::string& name, double z, double rhythm)
  {
    aprl::FieldDetectorPlaneConfig result;
    result.name=name; result.z=z; result.rhythm=rhythm;
    result.bufferSamples=8; result.compression=0; return result;
  }

  void boundaries(const std::string& directory, int rank, int size, Errors& errors)
  {
    const std::size_t nz=257;
    aprl::EBGridGeometry global=geometry(nz), slab(global); slab.nz=count(nz,rank,size);
    aprl::EBFieldGrid local(slab), whole(global);
    fill(local,offset(nz,rank,size),false); if(rank==0) fill(whole,0,false);
    aprl::FieldVector<Double> localOrigin(origin()); localOrigin[2]+=offset(nz,rank,size)*DZ;
    std::set<double> positions={0.0,.125,.499,.5, nz-.5,nz-.125,static_cast<double>(nz)};
    for (int edge : std::set<int>{1,size/2,size-1})
      if(edge>0 && edge<size)
        for(double delta : {-.5,-.25,0.0,.25,.5}) positions.insert(offset(nz,edge,size)+delta);
    for(unsigned mode=0;mode<3;++mode)
      {
        aprl::BoostFrameTransform frame(mode==0?1.0:GAMMA,
          mode==0?0.0:(mode==1?1.0:-1.0)*aprl::BoostFrameTransform::betaFromGamma(GAMMA),aprl::SI::c,origin()[2]);
        aprl::DetectorConfig config; config.directory=directory;
        std::vector<std::vector<Point> > expected;
        unsigned index=0;
        for(double q:positions)
          {
            const double zLab=frame.labZFromBoxZT(origin()[2]+q*DZ,0.0);
            config.fieldPlanes.push_back(plane("boundary-"+std::to_string(mode)+"-"+
              std::to_string(index++),zLab,1.e-15));
            if(rank==0) expected.push_back(sample(whole,origin(),frame,
              frame.boxZFromLabZAndBoxT(zLab,0.0)));
          }
        aprl::LabDetectorManager detector(config,global,origin(),localOrigin,frame,MPI_COMM_WORLD,NULL,true);
        detector.sampleFieldPlanes(local,0.0);
        detector.sampleFieldPlanes(local,0.0); // All planes not due: no duplicate records.
        detector.close(true);
        if(rank==0)
          for(std::size_t p=0;p<config.fieldPlanes.size();++p)
            {
              const Output output=read(directory+"/"+config.fieldPlanes[p].name+".h5");
              require(output.time.size()==1,"Repeated time produced duplicate detector sample");
              errors.boundary=std::max(errors.boundary,difference(output.field,expected[p]));
              ++errors.boundaryPlanes;
            }
      }
  }

  void writeScaled(std::ofstream& stream, const std::vector<Point>& points)
  {
    for(const Point& point:points)
      for(unsigned a=0;a<6;++a)
        { const double value=scaled(point,a); stream.write(reinterpret_cast<const char*>(&value),sizeof(value)); }
    require(stream.good(),"Cannot write reference binary");
  }

  void sweep(const std::string& directory, int rank, int size, Errors& errors)
  {
    const std::size_t nz=256;
    aprl::EBGridGeometry global=geometry(nz), slab(global); slab.nz=count(nz,rank,size);
    aprl::EBFieldGrid local(slab), whole(global);
    aprl::FieldVector<Double> localOrigin(origin()); localOrigin[2]+=offset(nz,rank,size)*DZ;
    aprl::BoostFrameTransform frame; frame.setFromGamma(GAMMA,aprl::SI::c,origin()[2]);
    const double dtBox=.25*DZ/(frame.beta()*aprl::SI::c), dtLab=dtBox/GAMMA;
    // Avoid exact vertex/half-cell collisions in this time trace. The legacy
    // sampler derives its coordinate from a separately rounded local origin;
    // at an exact handoff that can throw before its known clamp is recorded.
    // Dedicated boundary cases still test exact interfaces and global edges.
    const double startCell=nz-.7;
    const double zLab=frame.labZFromBoxZT(origin()[2]+startCell*DZ,0.0);
    aprl::DetectorConfig config; config.directory=directory;
    config.fieldPlanes.push_back(plane("sweep",zLab,.5*dtLab));
    config.fieldPlanes.push_back(plane("sweep-sparse",zLab,2.25*dtLab));
    aprl::LabDetectorManager detector(config,global,origin(),localOrigin,frame,MPI_COMM_WORLD,NULL,true);
    std::vector<Point> expected, legacy;
    std::vector<double> expectedTime, coordinates; std::vector<int> owners;
    std::ofstream referenceFile, legacyFile;
    if(rank==0)
      {
        referenceFile.open((directory+"/sweep-reference.bin").c_str(),std::ios::binary);
        legacyFile.open((directory+"/sweep-legacy.bin").c_str(),std::ios::binary);
      }
    for(std::size_t n=0;n<SAMPLES;++n)
      {
        const double t=n*dtBox, z=frame.boxZFromLabZAndBoxT(zLab,t), q=(z-origin()[2])/DZ;
        fill(local,offset(nz,rank,size),true,q,n);
        detector.sampleFieldPlanes(local,t);
        std::vector<Point> contribution(NX*NY), old(NX*NY);
        const int selected=owner(q,nz,size);
        if(rank==selected) contribution=sample(local,localOrigin,frame,z);
        MPI_Reduce(reinterpret_cast<double*>(contribution.data()),reinterpret_cast<double*>(old.data()),
          NX*NY*6,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        if(rank==0)
          {
            fill(whole,0,true,q,n);
            const std::vector<Point> target=sample(whole,origin(),frame,z);
            expected.insert(expected.end(),target.begin(),target.end());
            legacy.insert(legacy.end(),old.begin(),old.end());
            writeScaled(referenceFile,target); writeScaled(legacyFile,old);
            expectedTime.push_back(frame.labTFromLabZBoxT(zLab,t));
            coordinates.push_back(q); owners.push_back(selected);
          }
      }
    detector.close(true);
    if(rank!=0) return;
    referenceFile.close(); legacyFile.close();
    const Output output=read(directory+"/sweep.h5"), sparse=read(directory+"/sweep-sparse.h5");
    require(output.time.size()==SAMPLES,"Sweep cadence lost samples");
    errors.sweep=difference(output.field,expected); errors.legacy=difference(legacy,expected);
    for(std::size_t n=0;n<SAMPLES;++n)
      require(output.time[n]==expectedTime[n],"Detector time differs from transform time");
    require(sparse.time.size()==1+static_cast<std::size_t>((SAMPLES-1)/2.25),"Sparse cadence count mismatch");
    std::vector<Point> sparseExpected;
    for(std::size_t n=0;n<sparse.time.size();++n)
      {
        const std::size_t index=static_cast<std::size_t>(std::ceil(n*2.25));
        require(sparse.time[n]==expectedTime[index],"Sparse cadence selected wrong sample");
        sparseExpected.insert(sparseExpected.end(),expected.begin()+index*NX*NY,
          expected.begin()+(index+1)*NX*NY);
      }
    errors.sparse=difference(sparse.field,sparseExpected);
    std::ofstream trace((directory+"/sweep-trace.csv").c_str()); trace<<std::setprecision(17);
    trace<<"sample,time_s,z_cell,owner";
    for(const char* prefix:{"actual","reference","legacy"})
      for(const char* name:{"Ex","Ey","Ez","cBx","cBy","cBz"}) trace<<','<<prefix<<'_'<<name;
    trace<<'\n'; const std::size_t pixel=(NY/2)*NX+NX/2;
    for(std::size_t n=0;n<SAMPLES;++n)
      {
        trace<<n<<','<<output.time[n]<<','<<coordinates[n]<<','<<owners[n];
        const std::vector<Point>* traces[]={&output.field,&expected,&legacy};
        for(const std::vector<Point>* fields:traces)
          for(unsigned a=0;a<6;++a) trace<<','<<scaled((*fields)[n*NX*NY+pixel],a);
        trace<<'\n';
      }
    std::ofstream metadata((directory+"/sweep-metadata.json").c_str()); metadata<<std::setprecision(17);
    metadata<<"{\"fixture\":\"analytic diagnostic interpolation, not Maxwell evolution\","
      "\"binary_dtype\":\"native float64\",\"component_order\":[\"Ex\",\"Ey\",\"Ez\",\"cBx\",\"cBy\",\"cBz\"],"
      "\"shape\":["<<SAMPLES<<','<<NY<<','<<NX<<",6],\"ranks\":"<<size<<",\"global_nz\":"<<nz<<
      ",\"tone_bin\":"<<TONE<<",\"box_tone_amplitude\":1e-5,\"gamma\":"<<GAMMA<<
      ",\"dt_lab_s\":"<<dtLab<<",\"step_cells\":0.25,\"start_cell\":"<<startCell<<
      ",\"legacy\":\"same generic local-only interpolation as before fix\"}\n";
  }
}

int main(int argc,char** argv)
{
  MPI_Init(&argc,&argv); int rank=0,size=1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
  try
    {
      require(argc==2 || (argc==3 && std::string(argv[2])=="--sweep-only"),
        "Usage: aprl_detector_sampling_checks OUTPUT_DIRECTORY [--sweep-only]");
      const bool sweepOnly=argc==3;
      require(size<=128,"Fixture requires at least two local cells");
      if(rank==0) mkdirChecked(argv[1]);
      MPI_Barrier(MPI_COMM_WORLD);
      Errors errors;
      if(!sweepOnly)
        {
          checkRadiationScreenOutput(argv[1],rank);
          boundaries(argv[1],rank,size,errors);
        }
      sweep(argv[1],rank,size,errors);
      int failed=0;
      if(rank==0)
        {
          const double tolerance=2.e-10;
          failed=errors.boundary>tolerance || errors.sweep>tolerance || errors.sparse>tolerance;
          std::ofstream report((std::string(argv[1])+"/summary.json").c_str()); report<<std::setprecision(17);
          report<<"{\"pass\":"<<(failed?"false":"true")<<",\"ranks\":"<<size<<
            ",\"sweep_only\":"<<(sweepOnly?"true":"false")<<
            ",\"relative_tolerance\":"<<tolerance<<",\"boundary_planes\":"<<errors.boundaryPlanes<<
            ",\"boundary_max_relative_error\":"<<errors.boundary<<",\"sweep_max_relative_error\":"<<errors.sweep<<
            ",\"sparse_max_relative_error\":"<<errors.sparse<<",\"legacy_max_relative_error\":"<<errors.legacy<<"}\n";
          std::cout<<std::setprecision(17)<<"sampling ranks="<<size<<" boundary="<<errors.boundary<<
            " sweep="<<errors.sweep<<" sparse="<<errors.sparse<<" legacy="<<errors.legacy<<
            (failed?" FAIL":" PASS")<<'\n';
        }
      MPI_Bcast(&failed,1,MPI_INT,0,MPI_COMM_WORLD); MPI_Finalize(); return failed;
    }
  catch(const std::exception& error)
    { std::cerr<<"rank "<<rank<<": "<<error.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,2); return 2; }
}
