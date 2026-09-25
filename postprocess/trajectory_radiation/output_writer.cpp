#include "output_writer.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <sys/stat.h>

#include "radiation_kernel.h"

namespace radiation
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

    std::string parentDirectory(const std::string& path)
    {
      const std::string::size_type separator = path.find_last_of('/');
      if (separator == std::string::npos) return std::string();
      if (separator == 0) return "/";
      return path.substr(0, separator);
    }

    void createDirectories(const std::string& directory)
    {
      if (directory.empty() || directory == "." || directory == "/") return;
      std::string current;
      std::size_t start = 0;
      if (directory[0] == '/')
        {
          current = "/";
          start = 1;
        }
      while (start <= directory.size())
        {
          const std::size_t separator = directory.find('/', start);
          const std::string part = directory.substr(start,
            separator == std::string::npos ? std::string::npos :
            separator - start);
          if (!part.empty())
            {
              if (!current.empty() && current[current.size() - 1] != '/')
                current += "/";
              current += part;
              if (::mkdir(current.c_str(), 0775) != 0 && errno != EEXIST)
                throw std::runtime_error(
                  "Cannot create output directory: " + current +
                  ": " + std::strerror(errno));
            }
          if (separator == std::string::npos) break;
          start = separator + 1;
        }
    }

    std::vector<hsize_t> fullFieldDimensions(const RadiationConfig& config)
    {
      return std::vector<hsize_t>{
        static_cast<hsize_t>(config.shots.size()),
        static_cast<hsize_t>(config.photonEnergyEV.count),
        static_cast<hsize_t>(config.thetaY.count),
        static_cast<hsize_t>(config.thetaX.count), 2
      };
    }

    std::vector<hsize_t> summaryFieldDimensions(
        const RadiationConfig& config)
    {
      return std::vector<hsize_t>{
        static_cast<hsize_t>(config.photonEnergyEV.count),
        static_cast<hsize_t>(config.thetaY.count),
        static_cast<hsize_t>(config.thetaX.count), 2
      };
    }

    std::vector<hsize_t> fullScalarDimensions(const RadiationConfig& config)
    {
      return std::vector<hsize_t>{
        static_cast<hsize_t>(config.shots.size()),
        static_cast<hsize_t>(config.photonEnergyEV.count),
        static_cast<hsize_t>(config.thetaY.count),
        static_cast<hsize_t>(config.thetaX.count)
      };
    }

    std::vector<hsize_t> summaryScalarDimensions(
        const RadiationConfig& config)
    {
      return std::vector<hsize_t>{
        static_cast<hsize_t>(config.photonEnergyEV.count),
        static_cast<hsize_t>(config.thetaY.count),
        static_cast<hsize_t>(config.thetaX.count)
      };
    }

    std::vector<hsize_t> fieldChunks(const RadiationConfig& config,
                                     bool withShot)
    {
      std::vector<hsize_t> result;
      if (withShot) result.push_back(1);
      result.push_back(static_cast<hsize_t>(std::min(
        config.frequencyBlock, config.photonEnergyEV.count)));
      result.push_back(static_cast<hsize_t>(std::min(
        config.thetaYBlock, config.thetaY.count)));
      result.push_back(static_cast<hsize_t>(std::min<std::size_t>(
        64, config.thetaX.count)));
      result.push_back(2);
      return result;
    }

    std::vector<hsize_t> scalarChunks(const RadiationConfig& config,
                                      bool withShot,
                                      std::size_t trailing = 0)
    {
      std::vector<hsize_t> result;
      if (withShot) result.push_back(1);
      result.push_back(static_cast<hsize_t>(std::min(
        config.frequencyBlock, config.photonEnergyEV.count)));
      result.push_back(static_cast<hsize_t>(std::min(
        config.thetaYBlock, config.thetaY.count)));
      result.push_back(static_cast<hsize_t>(std::min<std::size_t>(
        64, config.thetaX.count)));
      if (trailing > 0) result.push_back(static_cast<hsize_t>(trailing));
      return result;
    }

    ComplexValue multiplyConjugate(const ComplexValue& left,
                                   const ComplexValue& right)
    {
      ComplexValue value;
      value.real = left.real * right.real + left.imag * right.imag;
      value.imag = left.real * right.imag - left.imag * right.real;
      return value;
    }
  }

  FarFieldWriter::FarFieldWriter(const RadiationConfig& config)
    : config_(config), photonEnergyEV_(config.photonEnergyEV.values()),
      thetaX_(config.thetaX.values()), thetaY_(config.thetaY.values()),
      omega_(photonEnergyToOmega(photonEnergyEV_)), file_(-1), group_(-1),
      complexType_(-1), shotField_(-1), shotIntensity_(-1),
      shotStokes_(-1), meanField_(-1), meanIntensity_(-1),
      coherentIntensity_(-1), incoherentIntensity_(-1),
      meanStokes_(-1), completeDataset_(-1), open_(false)
  {
    createDirectories(parentDirectory(config_.outputFile));
    file_ = H5Fcreate(config_.outputFile.c_str(),
      config_.overwrite ? H5F_ACC_TRUNC : H5F_ACC_EXCL,
                      H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(file_, "Cannot create far-field output: " +
      config_.outputFile);
    try
      {
        group_ = H5Gcreate2(file_, "/far_field", H5P_DEFAULT,
                            H5P_DEFAULT, H5P_DEFAULT);
        requireHandle(group_, "Cannot create /far_field output group");
        complexType_ = createComplexType();

        writeAxis("photon_energy_eV", photonEnergyEV_, "eV");
        std::vector<double> omegaDouble(omega_.size());
        std::vector<double> frequency(omega_.size());
        std::vector<double> wavelength(omega_.size());
        for (std::size_t index = 0; index < omega_.size(); ++index)
          {
            omegaDouble[index] = static_cast<double>(omega_[index]);
            frequency[index] = static_cast<double>(
              omega_[index] / (2.0L * constants::pi));
            wavelength[index] = static_cast<double>(
              2.0L * constants::pi * constants::c / omega_[index]);
          }
        writeAxis("omega_rad_per_s", omegaDouble, "rad/s");
        writeAxis("frequency_Hz", frequency, "Hz");
        writeAxis("wavelength_m", wavelength, "m");
        writeAxis("theta_x_rad", thetaX_, "rad");
        writeAxis("theta_y_rad", thetaY_, "rad");

        const double basis[9] = {
          static_cast<double>(config_.observationAxis[0]),
          static_cast<double>(config_.observationAxis[1]),
          static_cast<double>(config_.observationAxis[2]),
          static_cast<double>(config_.horizontalAxis[0]),
          static_cast<double>(config_.horizontalAxis[1]),
          static_cast<double>(config_.horizontalAxis[2]),
          static_cast<double>(config_.verticalAxis[0]),
          static_cast<double>(config_.verticalAxis[1]),
          static_cast<double>(config_.verticalAxis[2])
        };
        hsize_t basisDimensions[2] = {3, 3};
        hid_t basisSpace = H5Screate_simple(2, basisDimensions, NULL);
        hid_t basisDataset = H5Dcreate2(group_, "observation_basis",
          H5T_IEEE_F64LE, basisSpace, H5P_DEFAULT, H5P_DEFAULT,
          H5P_DEFAULT);
        H5Sclose(basisSpace);
        requireHandle(basisDataset, "Cannot create observation basis");
        requireStatus(H5Dwrite(basisDataset, H5T_NATIVE_DOUBLE,
          H5S_ALL, H5S_ALL, H5P_DEFAULT, basis),
          "Cannot write observation basis");
        writeStringAttribute(basisDataset, "row_order",
          "central_axis,horizontal,vertical");
        H5Dclose(basisDataset);

        std::size_t maximumShotName = 1;
        for (std::size_t shot = 0; shot < config_.shots.size(); ++shot)
          maximumShotName = std::max(maximumShotName,
            config_.shots[shot].name.size());
        hid_t shotNameType = H5Tcopy(H5T_C_S1);
        requireHandle(shotNameType, "Cannot create shot-name datatype");
        requireStatus(H5Tset_size(shotNameType, maximumShotName + 1),
                      "Cannot size shot-name datatype");
        hsize_t shotNameDimensions[1] = {
          static_cast<hsize_t>(config_.shots.size())};
        hid_t shotNameSpace = H5Screate_simple(1, shotNameDimensions, NULL);
        hid_t shotNameDataset = H5Dcreate2(group_, "shot_names",
          shotNameType, shotNameSpace, H5P_DEFAULT, H5P_DEFAULT,
          H5P_DEFAULT);
        H5Sclose(shotNameSpace);
        requireHandle(shotNameDataset, "Cannot create shot-name dataset");
        std::vector<char> shotNames(config_.shots.size() *
          (maximumShotName + 1), '\0');
        for (std::size_t shot = 0; shot < config_.shots.size(); ++shot)
          std::copy(config_.shots[shot].name.begin(),
            config_.shots[shot].name.end(),
            shotNames.begin() + shot * (maximumShotName + 1));
        requireStatus(H5Dwrite(shotNameDataset, shotNameType,
          H5S_ALL, H5S_ALL, H5P_DEFAULT, &shotNames[0]),
          "Cannot write shot names");
        H5Dclose(shotNameDataset);
        H5Tclose(shotNameType);

        hid_t scalarSpace = H5Screate(H5S_SCALAR);
        completeDataset_ = H5Dcreate2(group_, "complete", H5T_STD_U8LE,
          scalarSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(scalarSpace);
        requireHandle(completeDataset_, "Cannot create completion marker");
        const unsigned char incomplete = 0;
        requireStatus(H5Dwrite(completeDataset_, H5T_NATIVE_UCHAR,
          H5S_ALL, H5S_ALL, H5P_DEFAULT, &incomplete),
          "Cannot initialize completion marker");

        const std::vector<hsize_t> fieldDimensions =
          fullFieldDimensions(config_);
        const std::vector<hsize_t> scalarDimensions =
          fullScalarDimensions(config_);
        shotField_ = createDataset("electric_field_spectral",
          complexType_, fieldDimensions, fieldChunks(config_, true));
        shotIntensity_ = createDataset("spectral_energy_density",
          H5T_IEEE_F64LE, scalarDimensions,
          scalarChunks(config_, true));
        std::vector<hsize_t> stokesDimensions = scalarDimensions;
        stokesDimensions.push_back(4);
        shotStokes_ = createDataset("stokes_spectral_energy_density",
          H5T_IEEE_F64LE, stokesDimensions,
          scalarChunks(config_, true, 4));

        meanField_ = createDataset("mean_electric_field_spectral",
          complexType_, summaryFieldDimensions(config_),
          fieldChunks(config_, false));
        meanIntensity_ = createDataset("mean_spectral_energy_density",
          H5T_IEEE_F64LE, summaryScalarDimensions(config_),
          scalarChunks(config_, false));
        coherentIntensity_ = createDataset(
          "coherent_spectral_energy_density", H5T_IEEE_F64LE,
          summaryScalarDimensions(config_), scalarChunks(config_, false));
        incoherentIntensity_ = createDataset(
          "fluctuation_spectral_energy_density", H5T_IEEE_F64LE,
          summaryScalarDimensions(config_), scalarChunks(config_, false));
        std::vector<hsize_t> meanStokesDimensions =
          summaryScalarDimensions(config_);
        meanStokesDimensions.push_back(4);
        meanStokes_ = createDataset("mean_stokes_spectral_energy_density",
          H5T_IEEE_F64LE, meanStokesDimensions,
          scalarChunks(config_, false, 4));

        writeStringAttribute(shotField_, "unit", "V*s/m");
        writeStringAttribute(meanField_, "unit", "V*s/m");
        writeStringAttribute(shotIntensity_, "unit",
          "J*s/sr = d^2W/(domega*dOmega)");
        writeStringAttribute(meanIntensity_, "unit",
          "J*s/sr = d^2W/(domega*dOmega)");
        writeStringAttribute(coherentIntensity_, "unit",
          "J*s/sr = d^2W/(domega*dOmega)");
        writeStringAttribute(incoherentIntensity_, "unit",
          "J*s/sr = d^2W/(domega*dOmega)");
        writeStringAttribute(shotStokes_, "unit", "J*s/sr");
        writeStringAttribute(meanStokes_, "unit", "J*s/sr");

        writeUnsignedAttribute(group_, "format_version", 1);
        writeUnsignedAttribute(group_, "ensemble_shots",
          config_.shots.size());
        writeDoubleAttribute(group_, "observer_distance_m",
          config_.distanceM);
        writeStringAttribute(group_, "fourier_convention",
          "E_tilde(omega)=integral E(t) exp(+i omega t) dt");
        writeStringAttribute(group_, "trajectory_model",
          "piecewise-linear path; internal velocity-kink endpoint sum; artificial first/last endpoints suppressed");
        writeStringAttribute(group_, "angle_parameterization",
          "n=normalize(axis+tan(theta_x)*horizontal+tan(theta_y)*vertical)");
        writeStringAttribute(group_, "polarization_order",
          "0=projected horizontal,1=n_cross_horizontal");
        writeStringAttribute(group_, "field_unit", "V*s/m");
        writeStringAttribute(group_, "spectral_energy_density_unit",
          "J*s/sr = d^2W/(domega*dOmega)");
        writeStringAttribute(group_, "stokes_order",
          "I,Q,U,V with V=-2 Im(Eh*conj(Ev)) under the stored Fourier convention");
        writeStringAttribute(group_, "macroparticle_interpretation",
          "charge_C is summed coherently; fluctuation terms are physical shot noise only when the input ensemble and macroparticle model represent that noise");
        open_ = true;
      }
    catch (...)
      {
        close();
        throw;
      }
  }

  FarFieldWriter::~FarFieldWriter()
  {
    try { close(); } catch (...) {}
  }

  hid_t FarFieldWriter::fileHandle() const
  {
    if (!open_) throw std::logic_error("Far-field writer is closed");
    return file_;
  }

  hid_t FarFieldWriter::createComplexType() const
  {
    hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(ComplexValue));
    requireHandle(type, "Cannot create complex HDF5 datatype");
    try
      {
        requireStatus(H5Tinsert(type, "real", HOFFSET(ComplexValue, real),
          H5T_NATIVE_DOUBLE), "Cannot add complex real field");
        requireStatus(H5Tinsert(type, "imag", HOFFSET(ComplexValue, imag),
          H5T_NATIVE_DOUBLE), "Cannot add complex imaginary field");
      }
    catch (...)
      {
        H5Tclose(type);
        throw;
      }
    return type;
  }

  hid_t FarFieldWriter::createDataset(
      const char* name, hid_t type,
      const std::vector<hsize_t>& dimensions,
      const std::vector<hsize_t>& chunks) const
  {
    hid_t space = H5Screate_simple(static_cast<int>(dimensions.size()),
      &dimensions[0], NULL);
    requireHandle(space, std::string("Cannot create dataspace: ") + name);
    hid_t creation = H5Pcreate(H5P_DATASET_CREATE);
    if (creation < 0)
      {
        H5Sclose(space);
        throw std::runtime_error(std::string("Cannot create properties: ") + name);
      }
    try
      {
        requireStatus(H5Pset_chunk(creation,
          static_cast<int>(chunks.size()), &chunks[0]),
          std::string("Cannot set chunks: ") + name);
        if (config_.compression > 0)
          requireStatus(H5Pset_deflate(creation, config_.compression),
            std::string("Cannot set compression: ") + name);
        hid_t dataset = H5Dcreate2(group_, name, type, space, H5P_DEFAULT,
                                   creation, H5P_DEFAULT);
        H5Pclose(creation);
        H5Sclose(space);
        requireHandle(dataset, std::string("Cannot create dataset: ") + name);
        return dataset;
      }
    catch (...)
      {
        H5Pclose(creation);
        H5Sclose(space);
        throw;
      }
  }

  void FarFieldWriter::writeAxis(const char* name,
                                 const std::vector<double>& values,
                                 const char* unit) const
  {
    hsize_t count[1] = {static_cast<hsize_t>(values.size())};
    hid_t space = H5Screate_simple(1, count, NULL);
    requireHandle(space, std::string("Cannot create axis: ") + name);
    hid_t dataset = H5Dcreate2(group_, name, H5T_IEEE_F64LE, space,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(dataset, std::string("Cannot create axis dataset: ") + name);
    requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, values.empty() ? NULL : &values[0]),
      std::string("Cannot write axis: ") + name);
    writeStringAttribute(dataset, "unit", unit);
    H5Dclose(dataset);
  }

  void FarFieldWriter::writeStringAttribute(
      hid_t object, const char* name, const std::string& value) const
  {
    hid_t type = H5Tcopy(H5T_C_S1);
    requireHandle(type, "Cannot copy string datatype");
    requireStatus(H5Tset_size(type, value.size() + 1),
                  "Cannot size string attribute");
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create string attribute space");
    hid_t attribute = H5Acreate2(object, name, type, space,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute: ") + name);
    requireStatus(H5Awrite(attribute, type, value.c_str()),
                  std::string("Cannot write attribute: ") + name);
    H5Aclose(attribute);
    H5Tclose(type);
  }

  void FarFieldWriter::writeDoubleAttribute(
      hid_t object, const char* name, double value) const
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create scalar attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_IEEE_F64LE, space,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute: ") + name);
    requireStatus(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &value),
                  std::string("Cannot write attribute: ") + name);
    H5Aclose(attribute);
  }

  void FarFieldWriter::writeUnsignedAttribute(
      hid_t object, const char* name, unsigned long long value) const
  {
    hid_t space = H5Screate(H5S_SCALAR);
    requireHandle(space, "Cannot create scalar attribute space");
    hid_t attribute = H5Acreate2(object, name, H5T_STD_U64LE, space,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(attribute, std::string("Cannot create attribute: ") + name);
    requireStatus(H5Awrite(attribute, H5T_NATIVE_ULLONG, &value),
                  std::string("Cannot write attribute: ") + name);
    H5Aclose(attribute);
  }

  void FarFieldWriter::writeComplexBlock(
      hid_t dataset, std::size_t shot, const RadiationBlock& block,
      const std::vector<ComplexValue>& values) const
  {
    hsize_t start[5] = {static_cast<hsize_t>(shot),
      static_cast<hsize_t>(block.frequencyOffset),
      static_cast<hsize_t>(block.thetaYOffset), 0, 0};
    hsize_t count[5] = {1,
      static_cast<hsize_t>(block.frequencyCount),
      static_cast<hsize_t>(block.thetaYCount),
      static_cast<hsize_t>(block.thetaXCount), 2};
    hid_t fileSpace = H5Dget_space(dataset);
    requireHandle(fileSpace, "Cannot get complex output space");
    requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
      start, NULL, count, NULL), "Cannot select complex output block");
    hid_t memorySpace = H5Screate_simple(5, count, NULL);
    requireHandle(memorySpace, "Cannot create complex memory block");
    const herr_t status = H5Dwrite(dataset, complexType_, memorySpace,
      fileSpace, H5P_DEFAULT, values.empty() ? NULL : &values[0]);
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
    requireStatus(status, "Cannot write complex output block");
  }

  void FarFieldWriter::writeScalarBlock(
      hid_t dataset, std::size_t shot, const RadiationBlock& block,
      const std::vector<double>& values, std::size_t trailing) const
  {
    const int dimensions = trailing == 0 ? 4 : 5;
    hsize_t start[5] = {static_cast<hsize_t>(shot),
      static_cast<hsize_t>(block.frequencyOffset),
      static_cast<hsize_t>(block.thetaYOffset), 0, 0};
    hsize_t count[5] = {1,
      static_cast<hsize_t>(block.frequencyCount),
      static_cast<hsize_t>(block.thetaYCount),
      static_cast<hsize_t>(block.thetaXCount),
      static_cast<hsize_t>(trailing)};
    hid_t fileSpace = H5Dget_space(dataset);
    requireHandle(fileSpace, "Cannot get scalar output space");
    requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
      start, NULL, count, NULL), "Cannot select scalar output block");
    hid_t memorySpace = H5Screate_simple(dimensions, count, NULL);
    requireHandle(memorySpace, "Cannot create scalar memory block");
    const herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memorySpace,
      fileSpace, H5P_DEFAULT, values.empty() ? NULL : &values[0]);
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
    requireStatus(status, "Cannot write scalar output block");
  }

  void FarFieldWriter::readComplexBlock(
      hid_t dataset, std::size_t shot, const RadiationBlock& block,
      std::vector<ComplexValue>& values) const
  {
    values.resize(block.frequencyCount * block.thetaYCount *
                  block.thetaXCount * 2);
    hsize_t start[5] = {static_cast<hsize_t>(shot),
      static_cast<hsize_t>(block.frequencyOffset),
      static_cast<hsize_t>(block.thetaYOffset), 0, 0};
    hsize_t count[5] = {1,
      static_cast<hsize_t>(block.frequencyCount),
      static_cast<hsize_t>(block.thetaYCount),
      static_cast<hsize_t>(block.thetaXCount), 2};
    hid_t fileSpace = H5Dget_space(dataset);
    requireHandle(fileSpace, "Cannot get complex input space");
    requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
      start, NULL, count, NULL), "Cannot select complex input block");
    hid_t memorySpace = H5Screate_simple(5, count, NULL);
    requireHandle(memorySpace, "Cannot create complex input block");
    const herr_t status = H5Dread(dataset, complexType_, memorySpace,
      fileSpace, H5P_DEFAULT, values.empty() ? NULL : &values[0]);
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
    requireStatus(status, "Cannot read complex input block");
  }

  void FarFieldWriter::readScalarBlock(
      hid_t dataset, std::size_t shot, const RadiationBlock& block,
      std::vector<double>& values) const
  {
    values.resize(block.frequencyCount * block.thetaYCount *
                  block.thetaXCount);
    hsize_t start[4] = {static_cast<hsize_t>(shot),
      static_cast<hsize_t>(block.frequencyOffset),
      static_cast<hsize_t>(block.thetaYOffset), 0};
    hsize_t count[4] = {1,
      static_cast<hsize_t>(block.frequencyCount),
      static_cast<hsize_t>(block.thetaYCount),
      static_cast<hsize_t>(block.thetaXCount)};
    hid_t fileSpace = H5Dget_space(dataset);
    requireHandle(fileSpace, "Cannot get scalar input space");
    requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
      start, NULL, count, NULL), "Cannot select scalar input block");
    hid_t memorySpace = H5Screate_simple(4, count, NULL);
    requireHandle(memorySpace, "Cannot create scalar input block");
    const herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memorySpace,
      fileSpace, H5P_DEFAULT, values.empty() ? NULL : &values[0]);
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
    requireStatus(status, "Cannot read scalar input block");
  }

  void FarFieldWriter::writeShotBlock(
      std::size_t shot, const RadiationBlock& block,
      const std::vector<long double>& omega)
  {
    if (!open_ || shot >= config_.shots.size() ||
        omega.size() != photonEnergyEV_.size())
      throw std::logic_error("Invalid far-field block write");
    const std::size_t points = block.frequencyCount *
      block.thetaYCount * block.thetaXCount;
    std::vector<ComplexValue> field(points * 2);
    std::vector<double> intensity(points);
    std::vector<double> stokes(points * 4);
    const long double fieldScale = 1.0L /
      (4.0L * constants::pi * constants::epsilon0 * constants::c *
       static_cast<long double>(config_.distanceM));
    const long double energyScale = 1.0L /
      (16.0L * constants::pi * constants::pi * constants::pi *
       constants::epsilon0 * constants::c);
    for (std::size_t frequency = 0;
         frequency < block.frequencyCount; ++frequency)
      {
        const long double propagationPhase = std::remainder(
          omega[block.frequencyOffset + frequency] *
          static_cast<long double>(config_.distanceM) / constants::c,
          2.0L * constants::pi);
        const long double phaseReal = std::cos(propagationPhase);
        const long double phaseImag = std::sin(propagationPhase);
        for (std::size_t y = 0; y < block.thetaYCount; ++y)
          for (std::size_t x = 0; x < block.thetaXCount; ++x)
            {
              const std::size_t point =
                (frequency * block.thetaYCount + y) *
                block.thetaXCount + x;
              long double amplitudeReal[2] = {};
              long double amplitudeImag[2] = {};
              for (std::size_t polarization = 0;
                   polarization < 2; ++polarization)
                {
                  amplitudeReal[polarization] = block.amplitude[
                    block.scalarIndex(frequency, y, x, polarization, 0)];
                  amplitudeImag[polarization] = block.amplitude[
                    block.scalarIndex(frequency, y, x, polarization, 1)];
                  field[point * 2 + polarization].real =
                    static_cast<double>(fieldScale *
                      (amplitudeReal[polarization] * phaseReal -
                       amplitudeImag[polarization] * phaseImag));
                  field[point * 2 + polarization].imag =
                    static_cast<double>(fieldScale *
                      (amplitudeReal[polarization] * phaseImag +
                       amplitudeImag[polarization] * phaseReal));
                }
              const long double h2 = amplitudeReal[0] * amplitudeReal[0] +
                amplitudeImag[0] * amplitudeImag[0];
              const long double v2 = amplitudeReal[1] * amplitudeReal[1] +
                amplitudeImag[1] * amplitudeImag[1];
              const long double crossReal = amplitudeReal[0] *
                amplitudeReal[1] + amplitudeImag[0] * amplitudeImag[1];
              const long double crossImag = amplitudeImag[0] *
                amplitudeReal[1] - amplitudeReal[0] * amplitudeImag[1];
              intensity[point] = static_cast<double>(energyScale * (h2 + v2));
              stokes[point * 4] = intensity[point];
              stokes[point * 4 + 1] =
                static_cast<double>(energyScale * (h2 - v2));
              stokes[point * 4 + 2] =
                static_cast<double>(energyScale * 2.0L * crossReal);
              stokes[point * 4 + 3] =
                static_cast<double>(energyScale * -2.0L * crossImag);
            }
      }
    writeComplexBlock(shotField_, shot, block, field);
    writeScalarBlock(shotIntensity_, shot, block, intensity);
    writeScalarBlock(shotStokes_, shot, block, stokes, 4);
  }

  void FarFieldWriter::finalizeEnsemble()
  {
    if (!open_) throw std::logic_error("Far-field writer is closed");
    const long double fieldToEnergy =
      static_cast<long double>(config_.distanceM) * config_.distanceM /
      (constants::pi * constants::mu0 * constants::c);
    for (std::size_t frequencyOffset = 0;
         frequencyOffset < config_.photonEnergyEV.count;
         frequencyOffset += config_.frequencyBlock)
      for (std::size_t thetaYOffset = 0;
           thetaYOffset < config_.thetaY.count;
           thetaYOffset += config_.thetaYBlock)
        {
          RadiationBlock block;
          block.frequencyOffset = frequencyOffset;
          block.frequencyCount = std::min(config_.frequencyBlock,
            config_.photonEnergyEV.count - frequencyOffset);
          block.thetaYOffset = thetaYOffset;
          block.thetaYCount = std::min(config_.thetaYBlock,
            config_.thetaY.count - thetaYOffset);
          block.thetaXCount = config_.thetaX.count;
          const std::size_t points = block.frequencyCount *
            block.thetaYCount * block.thetaXCount;
          std::vector<long double> fieldSum(points * 4, 0.0L);
          std::vector<long double> intensitySum(points, 0.0L);
          std::vector<long double> stokesSum(points * 4, 0.0L);
          std::vector<ComplexValue> shotField;
          std::vector<double> shotIntensity;
          for (std::size_t shot = 0; shot < config_.shots.size(); ++shot)
            {
              readComplexBlock(shotField_, shot, block, shotField);
              readScalarBlock(shotIntensity_, shot, block, shotIntensity);
              for (std::size_t point = 0; point < points; ++point)
                {
                  intensitySum[point] += shotIntensity[point];
                  const ComplexValue& horizontal = shotField[point * 2];
                  const ComplexValue& vertical = shotField[point * 2 + 1];
                  fieldSum[point * 4] += horizontal.real;
                  fieldSum[point * 4 + 1] += horizontal.imag;
                  fieldSum[point * 4 + 2] += vertical.real;
                  fieldSum[point * 4 + 3] += vertical.imag;
                  const long double h2 = horizontal.real * horizontal.real +
                    horizontal.imag * horizontal.imag;
                  const long double v2 = vertical.real * vertical.real +
                    vertical.imag * vertical.imag;
                  const long double crossReal = horizontal.real * vertical.real +
                    horizontal.imag * vertical.imag;
                  const long double crossImag = horizontal.imag * vertical.real -
                    horizontal.real * vertical.imag;
                  stokesSum[point * 4] += fieldToEnergy * (h2 + v2);
                  stokesSum[point * 4 + 1] += fieldToEnergy * (h2 - v2);
                  stokesSum[point * 4 + 2] +=
                    fieldToEnergy * 2.0L * crossReal;
                  stokesSum[point * 4 + 3] +=
                    fieldToEnergy * -2.0L * crossImag;
                }
            }
          const long double inverseShots = 1.0L /
            static_cast<long double>(config_.shots.size());
          std::vector<ComplexValue> meanField(points * 2);
          std::vector<double> meanIntensity(points);
          std::vector<double> coherent(points);
          std::vector<double> fluctuation(points);
          std::vector<double> meanStokes(points * 4);
          for (std::size_t point = 0; point < points; ++point)
            {
              for (std::size_t polarization = 0;
                   polarization < 2; ++polarization)
                {
                  meanField[point * 2 + polarization].real =
                    static_cast<double>(inverseShots *
                      fieldSum[point * 4 + 2 * polarization]);
                  meanField[point * 2 + polarization].imag =
                    static_cast<double>(inverseShots *
                      fieldSum[point * 4 + 2 * polarization + 1]);
                }
              meanIntensity[point] = static_cast<double>(
                inverseShots * intensitySum[point]);
              long double fieldMagnitude2 = 0.0L;
              for (std::size_t polarization = 0;
                   polarization < 2; ++polarization)
                {
                  const ComplexValue& value =
                    meanField[point * 2 + polarization];
                  fieldMagnitude2 += value.real * value.real +
                    value.imag * value.imag;
                }
              coherent[point] = static_cast<double>(
                fieldToEnergy * fieldMagnitude2);
              fluctuation[point] = std::max(0.0,
                meanIntensity[point] - coherent[point]);
              for (std::size_t component = 0; component < 4; ++component)
                meanStokes[point * 4 + component] = static_cast<double>(
                  inverseShots * stokesSum[point * 4 + component]);
            }

          /* Summary datasets omit the shot dimension, so use an in-memory
           * view whose first dimension starts at frequency. */
          hsize_t complexStart[4] = {
            static_cast<hsize_t>(block.frequencyOffset),
            static_cast<hsize_t>(block.thetaYOffset), 0, 0};
          hsize_t complexCount[4] = {
            static_cast<hsize_t>(block.frequencyCount),
            static_cast<hsize_t>(block.thetaYCount),
            static_cast<hsize_t>(block.thetaXCount), 2};
          hid_t fieldSpace = H5Dget_space(meanField_);
          requireStatus(H5Sselect_hyperslab(fieldSpace, H5S_SELECT_SET,
            complexStart, NULL, complexCount, NULL),
            "Cannot select mean field output block");
          hid_t fieldMemory = H5Screate_simple(4, complexCount, NULL);
          requireStatus(H5Dwrite(meanField_, complexType_, fieldMemory,
            fieldSpace, H5P_DEFAULT, &meanField[0]),
            "Cannot write mean field output block");
          H5Sclose(fieldMemory);
          H5Sclose(fieldSpace);

          hsize_t scalarStart[4] = {
            static_cast<hsize_t>(block.frequencyOffset),
            static_cast<hsize_t>(block.thetaYOffset), 0, 0};
          hsize_t scalarCount[4] = {
            static_cast<hsize_t>(block.frequencyCount),
            static_cast<hsize_t>(block.thetaYCount),
            static_cast<hsize_t>(block.thetaXCount), 4};
          hid_t scalarMemory3 = H5Screate_simple(3, scalarCount, NULL);
          const hid_t scalarDatasets[3] = {
            meanIntensity_, coherentIntensity_, incoherentIntensity_
          };
          const std::vector<double>* scalarValues[3] = {
            &meanIntensity, &coherent, &fluctuation
          };
          for (unsigned int datasetIndex = 0; datasetIndex < 3;
               ++datasetIndex)
            {
              hid_t scalarSpace = H5Dget_space(scalarDatasets[datasetIndex]);
              requireStatus(H5Sselect_hyperslab(scalarSpace, H5S_SELECT_SET,
                scalarStart, NULL, scalarCount, NULL),
                "Cannot select summary intensity block");
              requireStatus(H5Dwrite(scalarDatasets[datasetIndex],
                H5T_NATIVE_DOUBLE, scalarMemory3, scalarSpace,
                H5P_DEFAULT, &(*scalarValues[datasetIndex])[0]),
                "Cannot write summary intensity block");
              H5Sclose(scalarSpace);
            }
          H5Sclose(scalarMemory3);
          hid_t stokesSpace = H5Dget_space(meanStokes_);
          requireStatus(H5Sselect_hyperslab(stokesSpace, H5S_SELECT_SET,
            scalarStart, NULL, scalarCount, NULL),
            "Cannot select mean Stokes block");
          hid_t stokesMemory = H5Screate_simple(4, scalarCount, NULL);
          requireStatus(H5Dwrite(meanStokes_, H5T_NATIVE_DOUBLE,
            stokesMemory, stokesSpace, H5P_DEFAULT, &meanStokes[0]),
            "Cannot write mean Stokes block");
          H5Sclose(stokesMemory);
          H5Sclose(stokesSpace);
        }
    finalizeSpectra();
    finalizeCoherence();
    const unsigned char complete = 1;
    requireStatus(H5Dwrite(completeDataset_, H5T_NATIVE_UCHAR,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &complete),
      "Cannot commit completion marker");
    requireStatus(H5Fflush(file_, H5F_SCOPE_GLOBAL),
                  "Cannot flush far-field output");
  }

  double FarFieldWriter::solidAngleWeight(
      std::size_t y, std::size_t x) const
  {
    if (thetaX_.size() < 2 || thetaY_.size() < 2) return 0.0;
    const long double dx = static_cast<long double>(
      (thetaX_.back() - thetaX_.front()) / (thetaX_.size() - 1));
    const long double dy = static_cast<long double>(
      (thetaY_.back() - thetaY_.front()) / (thetaY_.size() - 1));
    const long double tx = std::tan(static_cast<long double>(thetaX_[x]));
    const long double ty = std::tan(static_cast<long double>(thetaY_[y]));
    const long double secX2 = 1.0L + tx * tx;
    const long double secY2 = 1.0L + ty * ty;
    const long double jacobian = secX2 * secY2 /
      std::pow(1.0L + tx * tx + ty * ty, 1.5L);
    const long double endpointX =
      (x == 0 || x + 1 == thetaX_.size()) ? 0.5L : 1.0L;
    const long double endpointY =
      (y == 0 || y + 1 == thetaY_.size()) ? 0.5L : 1.0L;
    return static_cast<double>(dx * dy * endpointX * endpointY * jacobian);
  }

  void FarFieldWriter::finalizeSpectra()
  {
    const std::size_t shots = config_.shots.size();
    const std::size_t frequencies = photonEnergyEV_.size();
    std::vector<double> perShot(shots * frequencies, 0.0);
    RadiationBlock block;
    block.frequencyOffset = 0;
    block.frequencyCount = std::min(config_.frequencyBlock, frequencies);
    block.thetaXCount = thetaX_.size();
    for (std::size_t frequencyOffset = 0;
         frequencyOffset < frequencies;
         frequencyOffset += config_.frequencyBlock)
      {
        block.frequencyOffset = frequencyOffset;
        block.frequencyCount = std::min(config_.frequencyBlock,
          frequencies - frequencyOffset);
        for (std::size_t yOffset = 0; yOffset < thetaY_.size();
             yOffset += config_.thetaYBlock)
          {
            block.thetaYOffset = yOffset;
            block.thetaYCount = std::min(config_.thetaYBlock,
              thetaY_.size() - yOffset);
            std::vector<double> values;
            for (std::size_t shot = 0; shot < shots; ++shot)
              {
                readScalarBlock(shotIntensity_, shot, block, values);
                for (std::size_t f = 0; f < block.frequencyCount; ++f)
                  for (std::size_t y = 0; y < block.thetaYCount; ++y)
                    for (std::size_t x = 0; x < block.thetaXCount; ++x)
                      {
                        const std::size_t point =
                          (f * block.thetaYCount + y) *
                          block.thetaXCount + x;
                        perShot[shot * frequencies + frequencyOffset + f] +=
                          values[point] * solidAngleWeight(yOffset + y, x);
                      }
              }
          }
      }
    hsize_t perShotDimensions[2] = {
      static_cast<hsize_t>(shots), static_cast<hsize_t>(frequencies)};
    hid_t space = H5Screate_simple(2, perShotDimensions, NULL);
    requireHandle(space, "Cannot create integrated spectrum space");
    hid_t dataset = H5Dcreate2(group_, "energy_spectrum",
      H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(dataset, "Cannot create integrated energy spectrum");
    requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &perShot[0]), "Cannot write integrated energy spectrum");
    writeStringAttribute(dataset, "unit", "J*s = dW/domega");
    H5Dclose(dataset);

    std::vector<double> mean(frequencies, 0.0);
    std::vector<double> logarithmicEnergy(frequencies, 0.0);
    for (std::size_t f = 0; f < frequencies; ++f)
      {
        for (std::size_t shot = 0; shot < shots; ++shot)
          mean[f] += perShot[shot * frequencies + f];
        mean[f] /= static_cast<double>(shots);
        logarithmicEnergy[f] = mean[f] * static_cast<double>(omega_[f]);
      }
    hsize_t meanDimensions[1] = {static_cast<hsize_t>(frequencies)};
    space = H5Screate_simple(1, meanDimensions, NULL);
    dataset = H5Dcreate2(group_, "mean_energy_spectrum",
      H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(dataset, "Cannot create mean energy spectrum");
    requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &mean[0]), "Cannot write mean energy spectrum");
    writeStringAttribute(dataset, "unit", "J*s = dW/domega");
    H5Dclose(dataset);
    dataset = H5Dcreate2(group_, "mean_energy_per_log_frequency",
      H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(dataset, "Cannot create logarithmic energy spectrum");
    requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &logarithmicEnergy[0]),
      "Cannot write logarithmic energy spectrum");
    writeStringAttribute(dataset, "unit", "J = dW/dln(omega)");
    H5Dclose(dataset);

    std::vector<double> bandEnergy(shots, 0.0);
    if (frequencies >= 2)
      for (std::size_t shot = 0; shot < shots; ++shot)
        for (std::size_t f = 1; f < frequencies; ++f)
          bandEnergy[shot] += 0.5 *
            (perShot[shot * frequencies + f - 1] +
             perShot[shot * frequencies + f]) *
            static_cast<double>(omega_[f] - omega_[f - 1]);
    hsize_t energyDimensions[1] = {static_cast<hsize_t>(shots)};
    space = H5Screate_simple(1, energyDimensions, NULL);
    dataset = H5Dcreate2(group_, "band_energy_J", H5T_IEEE_F64LE,
      space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    requireHandle(dataset, "Cannot create band-energy dataset");
    requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &bandEnergy[0]), "Cannot write band energy");
    writeStringAttribute(dataset, "unit", "J");
    writeStringAttribute(dataset, "definition",
      "trapezoidal integral of energy_spectrum over configured omega axis and angular grid");
    H5Dclose(dataset);
    double meanBandEnergy = 0.0;
    for (std::size_t shot = 0; shot < shots; ++shot)
      meanBandEnergy += bandEnergy[shot];
    meanBandEnergy /= static_cast<double>(shots);
    writeDoubleAttribute(group_, "mean_band_energy_J", meanBandEnergy);
    writeDoubleAttribute(group_, "band_min_photon_energy_eV",
      photonEnergyEV_.front());
    writeDoubleAttribute(group_, "band_max_photon_energy_eV",
      photonEnergyEV_.back());
    writeUnsignedAttribute(group_, "angular_integral_available",
      thetaX_.size() >= 2 && thetaY_.size() >= 2 ? 1 : 0);
  }

  std::size_t FarFieldWriter::nearestIndex(
      const std::vector<double>& values, double requested) const
  {
    if (values.empty()) throw std::logic_error("Cannot search an empty axis");
    std::size_t best = 0;
    double distance = std::abs(values[0] - requested);
    for (std::size_t index = 1; index < values.size(); ++index)
      if (std::abs(values[index] - requested) < distance)
        {
          best = index;
          distance = std::abs(values[index] - requested);
        }
    return best;
  }

  ComplexValue FarFieldWriter::readFieldPoint(
      std::size_t shot, std::size_t frequency, std::size_t thetaY,
      std::size_t thetaX, std::size_t polarization) const
  {
    hsize_t start[5] = {static_cast<hsize_t>(shot),
      static_cast<hsize_t>(frequency), static_cast<hsize_t>(thetaY),
      static_cast<hsize_t>(thetaX), static_cast<hsize_t>(polarization)};
    hsize_t count[5] = {1, 1, 1, 1, 1};
    hid_t fileSpace = H5Dget_space(shotField_);
    requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
      start, NULL, count, NULL), "Cannot select field point");
    hid_t memorySpace = H5Screate_simple(5, count, NULL);
    ComplexValue value = {0.0, 0.0};
    requireStatus(H5Dread(shotField_, complexType_, memorySpace,
      fileSpace, H5P_DEFAULT, &value), "Cannot read field point");
    H5Sclose(memorySpace);
    H5Sclose(fileSpace);
    return value;
  }

  void FarFieldWriter::finalizeCoherence()
  {
    const CoherenceConfig& coherence = config_.coherence;
    if (coherence.spatialPhotonEnergyEV.empty() &&
        coherence.temporalPhotonEnergyEV.empty()) return;
    hid_t coherenceGroup = H5Gcreate2(file_, "/coherence", H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(coherenceGroup, "Cannot create /coherence group");
    writeUnsignedAttribute(coherenceGroup, "ensemble_shots",
      config_.shots.size());
    writeStringAttribute(coherenceGroup, "definition",
      "W_ab(p,q)=mean(conj(E_a(p))*E_b(q)); a is reference polarization");
    if (config_.shots.size() == 1)
      writeStringAttribute(coherenceGroup, "single_shot_warning",
        "A single deterministic field produces a rank-one outer product and cannot establish statistical partial coherence");

    if (!coherence.spatialPhotonEnergyEV.empty())
      {
        const std::size_t references =
          coherence.spatialReferenceAngles.size();
        const std::size_t selectedFrequencies =
          coherence.spatialPhotonEnergyEV.size();
        std::vector<std::size_t> frequencyIndex(selectedFrequencies);
        std::vector<double> actualEnergy(selectedFrequencies);
        for (std::size_t f = 0; f < selectedFrequencies; ++f)
          {
            frequencyIndex[f] = nearestIndex(photonEnergyEV_,
              coherence.spatialPhotonEnergyEV[f]);
            actualEnergy[f] = photonEnergyEV_[frequencyIndex[f]];
          }
        std::vector<std::size_t> referenceX(references);
        std::vector<std::size_t> referenceY(references);
        std::vector<double> actualAngles(references * 2);
        for (std::size_t reference = 0; reference < references; ++reference)
          {
            referenceX[reference] = nearestIndex(thetaX_,
              static_cast<double>(
                coherence.spatialReferenceAngles[reference][0]));
            referenceY[reference] = nearestIndex(thetaY_,
              static_cast<double>(
                coherence.spatialReferenceAngles[reference][1]));
            actualAngles[2 * reference] = thetaX_[referenceX[reference]];
            actualAngles[2 * reference + 1] = thetaY_[referenceY[reference]];
          }
        hsize_t dimensions[6] = {
          static_cast<hsize_t>(references),
          static_cast<hsize_t>(selectedFrequencies),
          static_cast<hsize_t>(thetaY_.size()),
          static_cast<hsize_t>(thetaX_.size()), 2, 2};
        hid_t space = H5Screate_simple(6, dimensions, NULL);
        hid_t dataset = H5Dcreate2(coherenceGroup,
          "spatial_cross_spectral_density", complexType_, space,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(space);
        requireHandle(dataset, "Cannot create spatial CSD dataset");
        for (std::size_t reference = 0; reference < references; ++reference)
          for (std::size_t selected = 0;
               selected < selectedFrequencies; ++selected)
            for (std::size_t yOffset = 0; yOffset < thetaY_.size();
                 yOffset += config_.thetaYBlock)
              {
                const std::size_t yCount = std::min(config_.thetaYBlock,
                  thetaY_.size() - yOffset);
                std::vector<ComplexValue> result(
                  yCount * thetaX_.size() * 4, ComplexValue{0.0, 0.0});
                RadiationBlock block;
                block.frequencyOffset = frequencyIndex[selected];
                block.frequencyCount = 1;
                block.thetaYOffset = yOffset;
                block.thetaYCount = yCount;
                block.thetaXCount = thetaX_.size();
                std::vector<ComplexValue> target;
                for (std::size_t shot = 0;
                     shot < config_.shots.size(); ++shot)
                  {
                    ComplexValue referenceField[2] = {
                      readFieldPoint(shot, frequencyIndex[selected],
                        referenceY[reference], referenceX[reference], 0),
                      readFieldPoint(shot, frequencyIndex[selected],
                        referenceY[reference], referenceX[reference], 1)
                    };
                    readComplexBlock(shotField_, shot, block, target);
                    for (std::size_t point = 0;
                         point < yCount * thetaX_.size(); ++point)
                      for (std::size_t a = 0; a < 2; ++a)
                        for (std::size_t b = 0; b < 2; ++b)
                          {
                            const ComplexValue product = multiplyConjugate(
                              referenceField[a], target[point * 2 + b]);
                            const std::size_t index = point * 4 + a * 2 + b;
                            result[index].real += product.real /
                              static_cast<double>(config_.shots.size());
                            result[index].imag += product.imag /
                              static_cast<double>(config_.shots.size());
                          }
                  }
                hsize_t start[6] = {
                  static_cast<hsize_t>(reference),
                  static_cast<hsize_t>(selected),
                  static_cast<hsize_t>(yOffset), 0, 0, 0};
                hsize_t count[6] = {1, 1,
                  static_cast<hsize_t>(yCount),
                  static_cast<hsize_t>(thetaX_.size()), 2, 2};
                hid_t fileSpace = H5Dget_space(dataset);
                requireStatus(H5Sselect_hyperslab(fileSpace,
                  H5S_SELECT_SET, start, NULL, count, NULL),
                  "Cannot select spatial CSD block");
                hid_t memorySpace = H5Screate_simple(6, count, NULL);
                requireStatus(H5Dwrite(dataset, complexType_, memorySpace,
                  fileSpace, H5P_DEFAULT, &result[0]),
                  "Cannot write spatial CSD block");
                H5Sclose(memorySpace);
                H5Sclose(fileSpace);
              }
        writeStringAttribute(dataset, "unit", "(V*s/m)^2");
        H5Dclose(dataset);

        hsize_t energySize[1] = {
          static_cast<hsize_t>(actualEnergy.size())};
        space = H5Screate_simple(1, energySize, NULL);
        dataset = H5Dcreate2(coherenceGroup,
          "spatial_photon_energy_eV", H5T_IEEE_F64LE, space,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(space);
        requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
          H5S_ALL, H5P_DEFAULT, &actualEnergy[0]),
          "Cannot write spatial CSD energies");
        H5Dclose(dataset);
        hsize_t angleSize[2] = {
          static_cast<hsize_t>(references), 2};
        space = H5Screate_simple(2, angleSize, NULL);
        dataset = H5Dcreate2(coherenceGroup,
          "spatial_reference_angles_rad", H5T_IEEE_F64LE, space,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(space);
        requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
          H5S_ALL, H5P_DEFAULT, &actualAngles[0]),
          "Cannot write spatial CSD reference angles");
        H5Dclose(dataset);
      }

    if (!coherence.temporalPhotonEnergyEV.empty())
      {
        const std::size_t references =
          coherence.temporalReferenceAngles.size();
        const std::size_t selectedFrequencies =
          coherence.temporalPhotonEnergyEV.size();
        std::vector<std::size_t> frequencyIndex(selectedFrequencies);
        std::vector<double> actualEnergy(selectedFrequencies);
        for (std::size_t f = 0; f < selectedFrequencies; ++f)
          {
            frequencyIndex[f] = nearestIndex(photonEnergyEV_,
              coherence.temporalPhotonEnergyEV[f]);
            actualEnergy[f] = photonEnergyEV_[frequencyIndex[f]];
          }
        hsize_t dimensions[5] = {
          static_cast<hsize_t>(references),
          static_cast<hsize_t>(selectedFrequencies),
          static_cast<hsize_t>(selectedFrequencies), 2, 2};
        hid_t space = H5Screate_simple(5, dimensions, NULL);
        hid_t dataset = H5Dcreate2(coherenceGroup,
          "temporal_cross_spectral_density", complexType_, space,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(space);
        requireHandle(dataset, "Cannot create temporal CSD dataset");
        std::vector<double> actualAngles(references * 2);
        for (std::size_t reference = 0; reference < references; ++reference)
          {
            const std::size_t x = nearestIndex(thetaX_,
              static_cast<double>(
                coherence.temporalReferenceAngles[reference][0]));
            const std::size_t y = nearestIndex(thetaY_,
              static_cast<double>(
                coherence.temporalReferenceAngles[reference][1]));
            actualAngles[2 * reference] = thetaX_[x];
            actualAngles[2 * reference + 1] = thetaY_[y];
            std::vector<ComplexValue> fields(config_.shots.size() *
              selectedFrequencies * 2);
            for (std::size_t shot = 0; shot < config_.shots.size(); ++shot)
              for (std::size_t f = 0; f < selectedFrequencies; ++f)
                for (std::size_t polarization = 0;
                     polarization < 2; ++polarization)
                  fields[(shot * selectedFrequencies + f) * 2 +
                    polarization] = readFieldPoint(shot,
                      frequencyIndex[f], y, x, polarization);
            std::vector<ComplexValue> result(
              selectedFrequencies * selectedFrequencies * 4,
              ComplexValue{0.0, 0.0});
            for (std::size_t first = 0;
                 first < selectedFrequencies; ++first)
              for (std::size_t second = 0;
                   second < selectedFrequencies; ++second)
                for (std::size_t a = 0; a < 2; ++a)
                  for (std::size_t b = 0; b < 2; ++b)
                    {
                      ComplexValue& value = result[
                        ((first * selectedFrequencies + second) * 2 + a) *
                        2 + b];
                      for (std::size_t shot = 0;
                           shot < config_.shots.size(); ++shot)
                        {
                          const ComplexValue product = multiplyConjugate(
                            fields[(shot * selectedFrequencies + first) *
                              2 + a],
                            fields[(shot * selectedFrequencies + second) *
                              2 + b]);
                          value.real += product.real /
                            static_cast<double>(config_.shots.size());
                          value.imag += product.imag /
                            static_cast<double>(config_.shots.size());
                        }
                    }
            hsize_t start[5] = {
              static_cast<hsize_t>(reference), 0, 0, 0, 0};
            hsize_t count[5] = {1,
              static_cast<hsize_t>(selectedFrequencies),
              static_cast<hsize_t>(selectedFrequencies), 2, 2};
            hid_t fileSpace = H5Dget_space(dataset);
            requireStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET,
              start, NULL, count, NULL),
              "Cannot select temporal CSD block");
            hid_t memorySpace = H5Screate_simple(5, count, NULL);
            requireStatus(H5Dwrite(dataset, complexType_, memorySpace,
              fileSpace, H5P_DEFAULT, &result[0]),
              "Cannot write temporal CSD block");
            H5Sclose(memorySpace);
            H5Sclose(fileSpace);
          }
        writeStringAttribute(dataset, "unit", "(V*s/m)^2");
        H5Dclose(dataset);

        hsize_t energySize[1] = {
          static_cast<hsize_t>(actualEnergy.size())};
        space = H5Screate_simple(1, energySize, NULL);
        dataset = H5Dcreate2(coherenceGroup,
          "temporal_photon_energy_eV", H5T_IEEE_F64LE, space,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(space);
        requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
          H5S_ALL, H5P_DEFAULT, &actualEnergy[0]),
          "Cannot write temporal CSD energies");
        H5Dclose(dataset);
        hsize_t angleSize[2] = {
          static_cast<hsize_t>(references), 2};
        space = H5Screate_simple(2, angleSize, NULL);
        dataset = H5Dcreate2(coherenceGroup,
          "temporal_reference_angles_rad", H5T_IEEE_F64LE, space,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(space);
        requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
          H5S_ALL, H5P_DEFAULT, &actualAngles[0]),
          "Cannot write temporal CSD reference angles");
        H5Dclose(dataset);
      }
    H5Gclose(coherenceGroup);
  }

  void FarFieldWriter::close()
  {
    if (shotField_ >= 0) H5Dclose(shotField_);
    if (shotIntensity_ >= 0) H5Dclose(shotIntensity_);
    if (shotStokes_ >= 0) H5Dclose(shotStokes_);
    if (meanField_ >= 0) H5Dclose(meanField_);
    if (meanIntensity_ >= 0) H5Dclose(meanIntensity_);
    if (coherentIntensity_ >= 0) H5Dclose(coherentIntensity_);
    if (incoherentIntensity_ >= 0) H5Dclose(incoherentIntensity_);
    if (meanStokes_ >= 0) H5Dclose(meanStokes_);
    if (completeDataset_ >= 0) H5Dclose(completeDataset_);
    if (complexType_ >= 0) H5Tclose(complexType_);
    if (group_ >= 0) H5Gclose(group_);
    if (file_ >= 0) H5Fclose(file_);
    shotField_ = shotIntensity_ = shotStokes_ = -1;
    meanField_ = meanIntensity_ = coherentIntensity_ = -1;
    incoherentIntensity_ = meanStokes_ = completeDataset_ = -1;
    complexType_ = group_ = file_ = -1;
    open_ = false;
  }
}
