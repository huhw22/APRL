#include "particle_file.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <stdexcept>

#include <hdf5.h>

namespace aprl
{
  namespace
  {
    void requireHandle(hid_t handle, const std::string& message)
    {
      if (handle < 0) throw std::runtime_error(message);
    }

    void requireStatus(herr_t status, const std::string& message)
    {
      if (status < 0) throw std::runtime_error(message);
    }

    hid_t memoryRecordType(int version)
    {
      hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(ParticleInputRecord));
      requireHandle(type, "Cannot create particle input memory datatype");
      const bool timedPlane =
        version >= ParticleHdf5File::formatVersion;
      hsize_t positionDimensions[1] = {timedPlane ? 2U : 3U};
      hsize_t velocityDimensions[1] = {3};
      hid_t positionType = H5Tarray_create2(
        H5T_NATIVE_DOUBLE, 1, positionDimensions);
      hid_t velocityType = H5Tarray_create2(
        H5T_NATIVE_DOUBLE, 1, velocityDimensions);
      if (positionType < 0 || velocityType < 0)
        {
          if (positionType >= 0) H5Tclose(positionType);
          if (velocityType >= 0) H5Tclose(velocityType);
          H5Tclose(type);
          throw std::runtime_error("Cannot create particle input vector datatype");
        }
      const herr_t first = H5Tinsert(type,
        timedPlane ? "plane_position_m" : "position_m",
        HOFFSET(ParticleInputRecord, position), positionType);
      const herr_t time = timedPlane ? H5Tinsert(type,
        "arrival_time_offset_s",
        HOFFSET(ParticleInputRecord, position) + 2 * sizeof(Double),
        H5T_NATIVE_DOUBLE) : 0;
      const herr_t second = H5Tinsert(type, "proper_velocity",
        HOFFSET(ParticleInputRecord, properVelocity), velocityType);
      const herr_t third = H5Tinsert(type, "source_id",
        HOFFSET(ParticleInputRecord, sourceId), H5T_NATIVE_UINT64);
      const herr_t fourth =
        version >= ParticleHdf5File::snapshotFormatVersion ?
        H5Tinsert(type, "macro_weight",
          HOFFSET(ParticleInputRecord, macroWeight), H5T_NATIVE_DOUBLE) : 0;
      H5Tclose(positionType);
      H5Tclose(velocityType);
      if (first < 0 || time < 0 || second < 0 || third < 0 || fourth < 0)
        {
          H5Tclose(type);
          throw std::runtime_error("Cannot define particle input memory datatype");
        }
      return type;
    }

    int readFormatVersion(hid_t group)
    {
      hid_t attribute = H5Aopen(group, "format_version", H5P_DEFAULT);
      requireHandle(attribute, "Particle file is missing format_version");
      int version = 0;
      const herr_t status = H5Aread(attribute, H5T_NATIVE_INT, &version);
      H5Aclose(attribute);
      requireStatus(status, "Cannot read particle format_version");
      return version;
    }
  }

  ParticleInputRecord::ParticleInputRecord()
    : sourceId(0), macroWeight(1.0)
  {
    for (unsigned int i = 0; i < 3; ++i)
      position[i] = properVelocity[i] = 0.0;
  }

  std::vector<RelativisticParticleSI> ParticleHdf5File::readDistributed(
      const std::string& filename, Double totalElectrons,
      const FieldVector<Double>& positionOffsetSI,
      MPI_Comm communicator, unsigned long long& globalRecords,
      int& inputFormatVersion)
  {
    if (communicator == MPI_COMM_NULL)
      throw std::invalid_argument("Particle input communicator is null");
    if (filename.empty())
      throw std::invalid_argument("Particle HDF5 filename cannot be empty");
    if (!(totalElectrons > 0.0) || !std::isfinite(totalElectrons))
      throw std::invalid_argument("Beam electron count must be positive");

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

#ifndef H5_HAVE_PARALLEL
    if (size > 1)
      throw std::runtime_error(
        "MPI particle input requires an HDF5 build with parallel I/O support");
#endif

    hid_t access = H5Pcreate(H5P_FILE_ACCESS);
    requireHandle(access, "Cannot create HDF5 particle file access properties");
#ifdef H5_HAVE_PARALLEL
    requireStatus(H5Pset_fapl_mpio(access, communicator, MPI_INFO_NULL),
                  "Cannot enable collective MPI particle input");
#endif
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, access);
    H5Pclose(access);
    requireHandle(file, "Cannot open particle HDF5 file: " + filename);

    hid_t group = H5Gopen2(file, "/particles", H5P_DEFAULT);
    if (group < 0)
      {
        H5Fclose(file);
        throw std::runtime_error("Particle file is missing /particles group");
      }
    const int version = readFormatVersion(group);
    if (version != legacyFormatVersion &&
        version != snapshotFormatVersion &&
        version != labPlaneOffsetFormatVersion &&
        version != formatVersion)
      {
        H5Gclose(group);
        H5Fclose(file);
        throw std::runtime_error(
          "Unsupported particle HDF5 format version; supported versions are 1, 2, 3 and 4");
      }
    inputFormatVersion = version;
    hid_t dataset = H5Dopen2(group, "records", H5P_DEFAULT);
    if (dataset < 0)
      {
        H5Gclose(group);
        H5Fclose(file);
        throw std::runtime_error("Particle file is missing /particles/records");
      }
    hid_t fileSpace = H5Dget_space(dataset);
    requireHandle(fileSpace, "Cannot inspect particle record dataset");
    if (H5Sget_simple_extent_ndims(fileSpace) != 1)
      {
        H5Sclose(fileSpace);
        H5Dclose(dataset);
        H5Gclose(group);
        H5Fclose(file);
        throw std::runtime_error("Particle records must be one-dimensional");
      }
    hsize_t dimensions[1] = {0};
    H5Sget_simple_extent_dims(fileSpace, dimensions, NULL);
    if (dimensions[0] == 0 || dimensions[0] >
        static_cast<hsize_t>(ULLONG_MAX))
      {
        H5Sclose(fileSpace);
        H5Dclose(dataset);
        H5Gclose(group);
        H5Fclose(file);
        throw std::runtime_error("Particle record dataset has invalid size");
      }
    globalRecords = static_cast<unsigned long long>(dimensions[0]);

    const unsigned long long base = globalRecords /
      static_cast<unsigned long long>(size);
    const unsigned long long remainder = globalRecords %
      static_cast<unsigned long long>(size);
    const unsigned long long localCount = base +
      (static_cast<unsigned long long>(rank) < remainder ? 1 : 0);
    const unsigned long long localOffset =
      static_cast<unsigned long long>(rank) * base +
      std::min(static_cast<unsigned long long>(rank), remainder);
    std::vector<ParticleInputRecord> input(
      static_cast<std::size_t>(localCount));

    hid_t memorySpace = -1;
    if (localCount > 0)
      {
        const hsize_t start[1] = {static_cast<hsize_t>(localOffset)};
        const hsize_t count[1] = {static_cast<hsize_t>(localCount)};
        requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
                                          start, NULL, count, NULL),
                      "Cannot select particle input hyperslab");
        memorySpace = H5Screate_simple(1, count, NULL);
      }
    else
      {
        requireStatus(H5Sselect_none(fileSpace),
                      "Cannot select empty particle hyperslab");
        memorySpace = H5Screate(H5S_NULL);
      }
    requireHandle(memorySpace, "Cannot create particle input memory space");
    hid_t transfer = H5Pcreate(H5P_DATASET_XFER);
    requireHandle(transfer, "Cannot create particle HDF5 transfer properties");
#ifdef H5_HAVE_PARALLEL
    requireStatus(H5Pset_dxpl_mpio(transfer, H5FD_MPIO_COLLECTIVE),
                  "Cannot enable collective particle dataset read");
#endif
    hid_t memoryType = memoryRecordType(version);
    const herr_t readStatus = H5Dread(dataset, memoryType, memorySpace,
      fileSpace, transfer, input.empty() ? NULL : &input[0]);
    H5Tclose(memoryType);
    H5Pclose(transfer);
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
    H5Dclose(dataset);
    H5Gclose(group);
    H5Fclose(file);
    requireStatus(readStatus, "Cannot read particle HDF5 records");

    Double localWeightSum = 0.0;
    for (std::size_t i = 0; i < input.size(); ++i)
      {
        if (!(input[i].macroWeight > 0.0) ||
            !std::isfinite(input[i].macroWeight))
          throw std::runtime_error(
            "Particle HDF5 macro_weight must be positive and finite");
        localWeightSum += input[i].macroWeight;
      }
    Double globalWeightSum = 0.0;
    MPI_Allreduce(&localWeightSum, &globalWeightSum, 1, MPI_DOUBLE,
                  MPI_SUM, communicator);
    if (!(globalWeightSum > 0.0) || !std::isfinite(globalWeightSum))
      throw std::runtime_error(
        "Particle HDF5 macro weights have an invalid global sum");
    std::vector<RelativisticParticleSI> result;
    result.reserve(input.size());
    for (std::size_t i = 0; i < input.size(); ++i)
      {
        RelativisticParticleSI particle;
        for (unsigned int component = 0; component < 3; ++component)
          {
            if (!std::isfinite(input[i].position[component]) ||
                !std::isfinite(input[i].properVelocity[component]))
              throw std::runtime_error("Particle HDF5 record contains non-finite values");
            particle.position[component] = input[i].position[component];
            if (version < formatVersion || component < 2)
              particle.position[component] += positionOffsetSI[component];
            particle.properVelocity[component] =
              input[i].properVelocity[component];
          }
        const Double macroElectrons = totalElectrons *
          input[i].macroWeight / globalWeightSum;
        particle.charge = -SI::elementaryCharge * macroElectrons;
        particle.mass = SI::electronMass * macroElectrons;
        particle.weight = input[i].macroWeight;
        particle.sourceId = input[i].sourceId;
        result.push_back(particle);
      }
    return result;
  }
}
