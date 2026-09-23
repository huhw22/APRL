#include "time_average.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "radiation_kernel.h"

namespace radiation
{
  namespace
  {
    struct ComplexValue
    {
      double real;
      double imag;
    };

    void requireHandle(hid_t handle, const std::string& message)
    {
      if (handle < 0) throw std::runtime_error(message);
    }

    void requireStatus(herr_t status, const std::string& message)
    {
      if (status < 0) throw std::runtime_error(message);
    }

    std::size_t checkedProduct(std::size_t left, std::size_t right,
                               const char* description)
    {
      if (right != 0 && left >
          std::numeric_limits<std::size_t>::max() / right)
        throw std::overflow_error(std::string("Size overflow in ") +
          description);
      return left * right;
    }

    std::size_t checkedSum(std::size_t left, std::size_t right,
                           const char* description)
    {
      if (left > std::numeric_limits<std::size_t>::max() - right)
        throw std::overflow_error(std::string("Size overflow in ") +
          description);
      return left + right;
    }

    void addProduct(long double leftReal, long double leftImag,
                    long double rightReal, long double rightImag,
                    long double& resultReal, long double& resultImag)
    {
      resultReal += leftReal * rightReal + leftImag * rightImag;
      resultImag += leftReal * rightImag - leftImag * rightReal;
    }

    hid_t createComplexFileType()
    {
      hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(ComplexValue));
      requireHandle(type, "Cannot create time-average complex file type");
      requireStatus(H5Tinsert(type, "real", HOFFSET(ComplexValue, real),
        H5T_IEEE_F64LE), "Cannot create complex real member");
      requireStatus(H5Tinsert(type, "imag", HOFFSET(ComplexValue, imag),
        H5T_IEEE_F64LE), "Cannot create complex imaginary member");
      return type;
    }

    hid_t createComplexMemoryType()
    {
      hid_t type = H5Tcreate(H5T_COMPOUND, sizeof(ComplexValue));
      requireHandle(type, "Cannot create time-average complex memory type");
      requireStatus(H5Tinsert(type, "real", HOFFSET(ComplexValue, real),
        H5T_NATIVE_DOUBLE), "Cannot map complex real member");
      requireStatus(H5Tinsert(type, "imag", HOFFSET(ComplexValue, imag),
        H5T_NATIVE_DOUBLE), "Cannot map complex imaginary member");
      return type;
    }

    hid_t createDataset(hid_t group, const char* name, hid_t type,
                        const std::vector<hsize_t>& dimensions,
                        unsigned int compression)
    {
      hid_t space = H5Screate_simple(static_cast<int>(dimensions.size()),
        &dimensions[0], NULL);
      requireHandle(space, std::string("Cannot create time-average space: ") +
        name);
      hid_t creation = H5P_DEFAULT;
      if (compression > 0)
        {
          creation = H5Pcreate(H5P_DATASET_CREATE);
          requireHandle(creation, "Cannot create time-average properties");
          std::vector<hsize_t> chunks(dimensions);
          for (std::size_t index = 0; index < chunks.size(); ++index)
            chunks[index] = std::min<hsize_t>(chunks[index], 16);
          requireStatus(H5Pset_chunk(creation,
            static_cast<int>(chunks.size()), &chunks[0]),
            "Cannot set time-average chunks");
          requireStatus(H5Pset_deflate(creation, compression),
            "Cannot set time-average compression");
        }
      hid_t dataset = H5Dcreate2(group, name, type, space, H5P_DEFAULT,
        creation, H5P_DEFAULT);
      if (creation != H5P_DEFAULT) H5Pclose(creation);
      H5Sclose(space);
      requireHandle(dataset, std::string("Cannot create time-average dataset: ") +
        name);
      return dataset;
    }

    void writeStringAttribute(hid_t object, const char* name,
                              const std::string& value)
    {
      hid_t type = H5Tcopy(H5T_C_S1);
      requireHandle(type, "Cannot copy string type");
      requireStatus(H5Tset_size(type, value.size() + 1),
        "Cannot size string attribute");
      hid_t space = H5Screate(H5S_SCALAR);
      hid_t attribute = H5Acreate2(object, name, type, space,
        H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(space);
      requireHandle(attribute, std::string("Cannot create attribute: ") + name);
      requireStatus(H5Awrite(attribute, type, value.c_str()),
        std::string("Cannot write attribute: ") + name);
      H5Aclose(attribute);
      H5Tclose(type);
    }

    void writeDoubleAttribute(hid_t object, const char* name, double value)
    {
      hid_t space = H5Screate(H5S_SCALAR);
      hid_t attribute = H5Acreate2(object, name, H5T_IEEE_F64LE, space,
        H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(space);
      requireHandle(attribute, std::string("Cannot create attribute: ") + name);
      requireStatus(H5Awrite(attribute, H5T_NATIVE_DOUBLE, &value),
        std::string("Cannot write attribute: ") + name);
      H5Aclose(attribute);
    }

    void writeUnsignedAttribute(hid_t object, const char* name,
                                unsigned long long value)
    {
      hid_t space = H5Screate(H5S_SCALAR);
      hid_t attribute = H5Acreate2(object, name, H5T_STD_U64LE, space,
        H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(space);
      requireHandle(attribute, std::string("Cannot create attribute: ") + name);
      requireStatus(H5Awrite(attribute, H5T_NATIVE_ULLONG, &value),
        std::string("Cannot write attribute: ") + name);
      H5Aclose(attribute);
    }

    void writeRealDataset(hid_t group, const char* name,
                          const std::vector<hsize_t>& dimensions,
                          const std::vector<double>& values,
                          unsigned int compression,
                          const char* unit)
    {
      hid_t dataset = createDataset(group, name, H5T_IEEE_F64LE,
        dimensions, compression);
      requireStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, values.empty() ? NULL : &values[0]),
        std::string("Cannot write time-average dataset: ") + name);
      if (unit != NULL) writeStringAttribute(dataset, "unit", unit);
      H5Dclose(dataset);
    }

    void writeComplexDataset(hid_t group, const char* name,
                             const std::vector<hsize_t>& dimensions,
                             const std::vector<ComplexValue>& values,
                             unsigned int compression,
                             const char* unit)
    {
      hid_t fileType = createComplexFileType();
      hid_t memoryType = createComplexMemoryType();
      hid_t dataset = createDataset(group, name, fileType, dimensions,
        compression);
      requireStatus(H5Dwrite(dataset, memoryType, H5S_ALL, H5S_ALL,
        H5P_DEFAULT, values.empty() ? NULL : &values[0]),
        std::string("Cannot write time-average dataset: ") + name);
      if (unit != NULL) writeStringAttribute(dataset, "unit", unit);
      H5Dclose(dataset);
      H5Tclose(memoryType);
      H5Tclose(fileType);
    }
  }

  TimeAverageAccumulator::TimeAverageAccumulator(
      const RadiationConfig& config)
    : config_(config), photonEnergyEV_(config.timeAverage.photonEnergyEV),
      omega_(photonEnergyToOmega(photonEnergyEV_)),
      thetaX_(config.thetaX.values()), thetaY_(config.thetaY.values()),
      windowCenters_(config.timeAverage.windowCenters()),
      referenceDirections_(), currentReferenceField_(), spatialCsd_(),
      targetIntensity_(), referenceIntensity_(), sampleReferenceIntensity_(),
      temporalCsd_(), expectedSamples_(0), estimatedBytes_(0),
      sampleCount_(0), sampleOpen_(false)
  {
    if (!config_.timeAverage.enabled)
      throw std::invalid_argument(
        "Cannot construct a disabled time-average accumulator");
    for (std::size_t reference = 0;
         reference < config_.timeAverage.referenceAngles.size(); ++reference)
      referenceDirections_.push_back(makeObservationDirection(config_,
        static_cast<double>(config_.timeAverage.referenceAngles[reference][0]),
        static_cast<double>(config_.timeAverage.referenceAngles[reference][1])));

    const std::size_t references = referenceDirections_.size();
    const std::size_t frequencies = omega_.size();
    const std::size_t angles = checkedProduct(thetaX_.size(), thetaY_.size(),
      "time-average angular grid");
    std::size_t spatialValues = checkedProduct(references, frequencies,
      "time-average spatial CSD");
    spatialValues = checkedProduct(spatialValues, angles,
      "time-average spatial CSD");
    spatialValues = checkedProduct(spatialValues, 8,
      "time-average spatial CSD");
    const std::size_t targetValues = checkedProduct(frequencies, angles,
      "time-average angular intensity");
    const std::size_t referenceValues = checkedProduct(references, frequencies,
      "time-average reference intensity");
    expectedSamples_ = checkedProduct(config_.shots.size(),
      windowCenters_.size(), "time-average sample count");
    const std::size_t sampleReferenceValues = checkedProduct(
      expectedSamples_, referenceValues,
      "time-average per-window reference intensity");
    std::size_t temporalValues = checkedProduct(referenceValues, frequencies,
      "time-average two-frequency CSD");
    temporalValues = checkedProduct(temporalValues, 8,
      "time-average two-frequency CSD");
    const std::size_t currentValues = checkedProduct(referenceValues, 4,
      "time-average reference field");
    std::size_t allValues = checkedSum(spatialValues, targetValues,
      "time-average accumulator");
    allValues = checkedSum(allValues, referenceValues,
      "time-average accumulator");
    allValues = checkedSum(allValues, sampleReferenceValues,
      "time-average accumulator");
    allValues = checkedSum(allValues, temporalValues,
      "time-average accumulator");
    allValues = checkedSum(allValues, currentValues,
      "time-average accumulator");
    estimatedBytes_ = checkedProduct(allValues, sizeof(long double),
      "time-average accumulator bytes");
    const std::size_t maximumBytes = checkedProduct(
      config_.timeAverage.maximumAccumulatorMiB, 1024U * 1024U,
      "time-average memory limit");
    if (estimatedBytes_ > maximumBytes)
      throw std::invalid_argument(
        "time_average accumulator exceeds maximum_accumulator_mib; "
        "reduce references, frequencies, or angular grid");

    currentReferenceField_.assign(currentValues, 0.0L);
    spatialCsd_.assign(spatialValues, 0.0L);
    targetIntensity_.assign(targetValues, 0.0L);
    referenceIntensity_.assign(referenceValues, 0.0L);
    sampleReferenceIntensity_.assign(sampleReferenceValues, 0.0L);
    temporalCsd_.assign(temporalValues, 0.0L);
  }

  const std::vector<long double>& TimeAverageAccumulator::omega() const
  {
    return omega_;
  }

  const std::vector<double>& TimeAverageAccumulator::windowCenters() const
  {
    return windowCenters_;
  }

  const std::vector<ObservationDirection>&
  TimeAverageAccumulator::referenceDirections() const
  {
    return referenceDirections_;
  }

  std::size_t TimeAverageAccumulator::estimatedBytes() const
  {
    return estimatedBytes_;
  }

  unsigned long long TimeAverageAccumulator::sampleCount() const
  {
    return sampleCount_;
  }

  std::size_t TimeAverageAccumulator::spatialIndex(
      std::size_t reference, std::size_t frequency,
      std::size_t thetaY, std::size_t thetaX,
      std::size_t referencePolarization, std::size_t targetPolarization,
      std::size_t realOrImag) const
  {
    return ((((((reference * omega_.size() + frequency) * thetaY_.size() +
      thetaY) * thetaX_.size() + thetaX) * 2 + referencePolarization) *
      2 + targetPolarization) * 2 + realOrImag);
  }

  std::size_t TimeAverageAccumulator::temporalIndex(
      std::size_t reference, std::size_t firstFrequency,
      std::size_t secondFrequency, std::size_t firstPolarization,
      std::size_t secondPolarization, std::size_t realOrImag) const
  {
    return ((((((reference * omega_.size() + firstFrequency) * omega_.size() +
      secondFrequency) * 2 + firstPolarization) * 2 +
      secondPolarization) * 2 + realOrImag));
  }

  void TimeAverageAccumulator::beginSample(
      const std::vector<long double>& referenceAmplitude)
  {
    if (sampleOpen_) throw std::logic_error("A time-average sample is open");
    if (sampleCount_ >= expectedSamples_)
      throw std::logic_error("Too many time-average samples");
    const std::size_t expected = omega_.size() *
      referenceDirections_.size() * 4;
    if (referenceAmplitude.size() != expected)
      throw std::invalid_argument("Invalid reference radiation sample size");
    const long double fieldScale = 1.0L /
      (4.0L * constants::pi * constants::epsilon0 * constants::c *
       static_cast<long double>(config_.distanceM));
    for (std::size_t frequency = 0; frequency < omega_.size(); ++frequency)
      {
        const long double phase = std::remainder(
          omega_[frequency] * static_cast<long double>(config_.distanceM) /
          constants::c, 2.0L * constants::pi);
        const long double phaseReal = std::cos(phase);
        const long double phaseImag = std::sin(phase);
        for (std::size_t reference = 0;
             reference < referenceDirections_.size(); ++reference)
          for (std::size_t polarization = 0; polarization < 2;
               ++polarization)
            {
              const std::size_t source =
                ((frequency * referenceDirections_.size() + reference) *
                  2 + polarization) * 2;
              const std::size_t target =
                ((reference * omega_.size() + frequency) * 2 +
                  polarization) * 2;
              currentReferenceField_[target] = fieldScale *
                (referenceAmplitude[source] * phaseReal -
                 referenceAmplitude[source + 1] * phaseImag);
              currentReferenceField_[target + 1] = fieldScale *
                (referenceAmplitude[source] * phaseImag +
                 referenceAmplitude[source + 1] * phaseReal);
            }
      }

    for (std::size_t reference = 0;
         reference < referenceDirections_.size(); ++reference)
      {
        for (std::size_t first = 0; first < omega_.size(); ++first)
          {
            long double intensity = 0.0L;
            for (std::size_t polarization = 0; polarization < 2;
                 ++polarization)
              {
                const std::size_t index =
                  ((reference * omega_.size() + first) * 2 +
                    polarization) * 2;
                intensity += currentReferenceField_[index] *
                  currentReferenceField_[index] +
                  currentReferenceField_[index + 1] *
                  currentReferenceField_[index + 1];
              }
            referenceIntensity_[reference * omega_.size() + first] +=
              intensity;
            sampleReferenceIntensity_[
              (static_cast<std::size_t>(sampleCount_) *
                referenceDirections_.size() + reference) * omega_.size() +
                first] = intensity;
            for (std::size_t second = 0; second < omega_.size(); ++second)
              for (std::size_t a = 0; a < 2; ++a)
                for (std::size_t b = 0; b < 2; ++b)
                  {
                    const std::size_t left =
                      ((reference * omega_.size() + first) * 2 + a) * 2;
                    const std::size_t right =
                      ((reference * omega_.size() + second) * 2 + b) * 2;
                    const std::size_t result = temporalIndex(reference,
                      first, second, a, b, 0);
                    addProduct(currentReferenceField_[left],
                      currentReferenceField_[left + 1],
                      currentReferenceField_[right],
                      currentReferenceField_[right + 1],
                      temporalCsd_[result], temporalCsd_[result + 1]);
                  }
          }
      }
    sampleOpen_ = true;
  }

  void TimeAverageAccumulator::accumulateAngularBlock(
      const std::vector<long double>& targetAmplitude,
      std::size_t thetaYOffset, std::size_t thetaYCount)
  {
    if (!sampleOpen_) throw std::logic_error("No time-average sample is open");
    if (thetaYOffset + thetaYCount > thetaY_.size())
      throw std::out_of_range("Time-average angular block is outside grid");
    const std::size_t directions = thetaYCount * thetaX_.size();
    if (targetAmplitude.size() != omega_.size() * directions * 4)
      throw std::invalid_argument("Invalid angular radiation sample size");
    const long double fieldScale = 1.0L /
      (4.0L * constants::pi * constants::epsilon0 * constants::c *
       static_cast<long double>(config_.distanceM));
    for (std::size_t frequency = 0; frequency < omega_.size(); ++frequency)
      {
        const long double phase = std::remainder(
          omega_[frequency] * static_cast<long double>(config_.distanceM) /
          constants::c, 2.0L * constants::pi);
        const long double phaseReal = std::cos(phase);
        const long double phaseImag = std::sin(phase);
        for (std::size_t localY = 0; localY < thetaYCount; ++localY)
          for (std::size_t x = 0; x < thetaX_.size(); ++x)
            {
              const std::size_t direction = localY * thetaX_.size() + x;
              long double targetField[4] = {};
              for (std::size_t polarization = 0; polarization < 2;
                   ++polarization)
                {
                  const std::size_t source =
                    ((frequency * directions + direction) * 2 +
                      polarization) * 2;
                  targetField[2 * polarization] = fieldScale *
                    (targetAmplitude[source] * phaseReal -
                     targetAmplitude[source + 1] * phaseImag);
                  targetField[2 * polarization + 1] = fieldScale *
                    (targetAmplitude[source] * phaseImag +
                     targetAmplitude[source + 1] * phaseReal);
                }
              const std::size_t angle =
                (frequency * thetaY_.size() + thetaYOffset + localY) *
                thetaX_.size() + x;
              targetIntensity_[angle] +=
                targetField[0] * targetField[0] +
                targetField[1] * targetField[1] +
                targetField[2] * targetField[2] +
                targetField[3] * targetField[3];
              for (std::size_t reference = 0;
                   reference < referenceDirections_.size(); ++reference)
                for (std::size_t a = 0; a < 2; ++a)
                  for (std::size_t b = 0; b < 2; ++b)
                    {
                      const std::size_t left =
                        ((reference * omega_.size() + frequency) * 2 + a) * 2;
                      const std::size_t result = spatialIndex(reference,
                        frequency, thetaYOffset + localY, x, a, b, 0);
                      addProduct(currentReferenceField_[left],
                        currentReferenceField_[left + 1],
                        targetField[2 * b], targetField[2 * b + 1],
                        spatialCsd_[result], spatialCsd_[result + 1]);
                    }
            }
      }
  }

  void TimeAverageAccumulator::endSample()
  {
    if (!sampleOpen_) throw std::logic_error("No time-average sample is open");
    ++sampleCount_;
    sampleOpen_ = false;
  }

  void TimeAverageAccumulator::write(hid_t file) const
  {
    if (sampleOpen_ || sampleCount_ == 0 ||
        sampleCount_ != expectedSamples_)
      throw std::logic_error("Time-average samples are incomplete");
    hid_t group = H5Gcreate2(file, "/time_average", H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);
    requireHandle(group, "Cannot create /time_average group");
    const long double inverseSamples = 1.0L /
      static_cast<long double>(sampleCount_);
    const long double fieldToEnergy =
      static_cast<long double>(config_.distanceM) * config_.distanceM /
      (constants::pi * constants::mu0 * constants::c);

    std::vector<ComplexValue> spatial(spatialCsd_.size() / 2);
    for (std::size_t index = 0; index < spatial.size(); ++index)
      {
        spatial[index].real = static_cast<double>(
          spatialCsd_[2 * index] * inverseSamples);
        spatial[index].imag = static_cast<double>(
          spatialCsd_[2 * index + 1] * inverseSamples);
      }
    std::vector<double> meanEnergy(targetIntensity_.size());
    for (std::size_t index = 0; index < meanEnergy.size(); ++index)
      meanEnergy[index] = static_cast<double>(targetIntensity_[index] *
        inverseSamples * fieldToEnergy);
    std::vector<double> referenceEnergy(referenceIntensity_.size());
    for (std::size_t index = 0; index < referenceEnergy.size(); ++index)
      referenceEnergy[index] = static_cast<double>(referenceIntensity_[index] *
        inverseSamples * fieldToEnergy);
    std::vector<double> sampleReferenceEnergy(
      sampleReferenceIntensity_.size());
    for (std::size_t index = 0; index < sampleReferenceEnergy.size(); ++index)
      sampleReferenceEnergy[index] = static_cast<double>(
        sampleReferenceIntensity_[index] * fieldToEnergy);
    std::vector<ComplexValue> temporal(temporalCsd_.size() / 2);
    for (std::size_t index = 0; index < temporal.size(); ++index)
      {
        temporal[index].real = static_cast<double>(
          temporalCsd_[2 * index] * inverseSamples);
        temporal[index].imag = static_cast<double>(
          temporalCsd_[2 * index + 1] * inverseSamples);
      }
    std::vector<double> coherenceSquared(
      referenceDirections_.size() * omega_.size() * thetaY_.size() *
      thetaX_.size(), 0.0);
    for (std::size_t reference = 0;
         reference < referenceDirections_.size(); ++reference)
      for (std::size_t frequency = 0; frequency < omega_.size(); ++frequency)
        for (std::size_t y = 0; y < thetaY_.size(); ++y)
          for (std::size_t x = 0; x < thetaX_.size(); ++x)
            {
              long double numerator = 0.0L;
              for (std::size_t a = 0; a < 2; ++a)
                for (std::size_t b = 0; b < 2; ++b)
                  {
                    const std::size_t index = spatialIndex(reference,
                      frequency, y, x, a, b, 0);
                    const long double real = spatialCsd_[index] *
                      inverseSamples;
                    const long double imag = spatialCsd_[index + 1] *
                      inverseSamples;
                    numerator += real * real + imag * imag;
                  }
              const long double referenceValue = referenceIntensity_[
                reference * omega_.size() + frequency] * inverseSamples;
              const long double targetValue = targetIntensity_[
                (frequency * thetaY_.size() + y) * thetaX_.size() + x] *
                inverseSamples;
              const long double denominator = referenceValue * targetValue;
              const std::size_t output =
                ((reference * omega_.size() + frequency) * thetaY_.size() +
                  y) * thetaX_.size() + x;
              if (denominator > 0.0L)
                coherenceSquared[output] = static_cast<double>(
                  std::min(1.0L, std::max(0.0L, numerator / denominator)));
            }
    std::vector<double> solidAngleWeight(
      thetaY_.size() * thetaX_.size(), 0.0);
    if (thetaX_.size() >= 2 && thetaY_.size() >= 2)
      {
        const long double dx = static_cast<long double>(
          (thetaX_.back() - thetaX_.front()) / (thetaX_.size() - 1));
        const long double dy = static_cast<long double>(
          (thetaY_.back() - thetaY_.front()) / (thetaY_.size() - 1));
        for (std::size_t y = 0; y < thetaY_.size(); ++y)
          for (std::size_t x = 0; x < thetaX_.size(); ++x)
            {
              const long double tx = std::tan(
                static_cast<long double>(thetaX_[x]));
              const long double ty = std::tan(
                static_cast<long double>(thetaY_[y]));
              const long double jacobian = (1.0L + tx * tx) *
                (1.0L + ty * ty) /
                std::pow(1.0L + tx * tx + ty * ty, 1.5L);
              const long double endpointX =
                (x == 0 || x + 1 == thetaX_.size()) ? 0.5L : 1.0L;
              const long double endpointY =
                (y == 0 || y + 1 == thetaY_.size()) ? 0.5L : 1.0L;
              solidAngleWeight[y * thetaX_.size() + x] =
                static_cast<double>(dx * dy * jacobian *
                  endpointX * endpointY);
            }
      }

    const std::vector<hsize_t> spatialDimensions = {
      static_cast<hsize_t>(referenceDirections_.size()),
      static_cast<hsize_t>(omega_.size()),
      static_cast<hsize_t>(thetaY_.size()),
      static_cast<hsize_t>(thetaX_.size()), 2, 2};
    writeComplexDataset(group, "spatial_cross_spectral_density",
      spatialDimensions, spatial, config_.compression, "(V*s/m)^2");
    const std::vector<hsize_t> angularDimensions = {
      static_cast<hsize_t>(omega_.size()),
      static_cast<hsize_t>(thetaY_.size()),
      static_cast<hsize_t>(thetaX_.size())};
    writeRealDataset(group, "mean_window_spectral_energy_density",
      angularDimensions, meanEnergy, config_.compression,
      "J*s/sr, arithmetic mean per Hann window");
    const std::vector<hsize_t> referenceDimensions = {
      static_cast<hsize_t>(referenceDirections_.size()),
      static_cast<hsize_t>(omega_.size())};
    writeRealDataset(group, "reference_mean_window_spectral_energy_density",
      referenceDimensions, referenceEnergy, config_.compression,
      "J*s/sr, arithmetic mean per Hann window");
    const std::vector<hsize_t> sampleReferenceDimensions = {
      static_cast<hsize_t>(config_.shots.size()),
      static_cast<hsize_t>(windowCenters_.size()),
      static_cast<hsize_t>(referenceDirections_.size()),
      static_cast<hsize_t>(omega_.size())};
    writeRealDataset(group, "reference_window_spectral_energy_density",
      sampleReferenceDimensions, sampleReferenceEnergy,
      config_.compression, "J*s/sr per Hann window");
    const std::vector<hsize_t> coherenceDimensions = {
      static_cast<hsize_t>(referenceDirections_.size()),
      static_cast<hsize_t>(omega_.size()),
      static_cast<hsize_t>(thetaY_.size()),
      static_cast<hsize_t>(thetaX_.size())};
    writeRealDataset(group, "spectral_degree_of_coherence_squared",
      coherenceDimensions, coherenceSquared, config_.compression, "1");
    writeRealDataset(group, "solid_angle_quadrature_weight",
      std::vector<hsize_t>{static_cast<hsize_t>(thetaY_.size()),
        static_cast<hsize_t>(thetaX_.size())},
      solidAngleWeight, 0, "sr");
    const std::vector<hsize_t> temporalDimensions = {
      static_cast<hsize_t>(referenceDirections_.size()),
      static_cast<hsize_t>(omega_.size()),
      static_cast<hsize_t>(omega_.size()), 2, 2};
    writeComplexDataset(group, "two_frequency_cross_spectral_density",
      temporalDimensions, temporal, config_.compression, "(V*s/m)^2");

    std::vector<double> omegaDouble(omega_.size());
    for (std::size_t index = 0; index < omega_.size(); ++index)
      omegaDouble[index] = static_cast<double>(omega_[index]);
    std::vector<double> referenceAngles(referenceDirections_.size() * 2);
    for (std::size_t index = 0; index < referenceDirections_.size(); ++index)
      {
        referenceAngles[2 * index] = referenceDirections_[index].thetaX;
        referenceAngles[2 * index + 1] = referenceDirections_[index].thetaY;
      }
    writeRealDataset(group, "photon_energy_eV",
      std::vector<hsize_t>{static_cast<hsize_t>(photonEnergyEV_.size())},
      photonEnergyEV_, 0, "eV");
    writeRealDataset(group, "omega_rad_per_s",
      std::vector<hsize_t>{static_cast<hsize_t>(omegaDouble.size())},
      omegaDouble, 0, "rad/s");
    writeRealDataset(group, "theta_x_rad",
      std::vector<hsize_t>{static_cast<hsize_t>(thetaX_.size())},
      thetaX_, 0, "rad");
    writeRealDataset(group, "theta_y_rad",
      std::vector<hsize_t>{static_cast<hsize_t>(thetaY_.size())},
      thetaY_, 0, "rad");
    writeRealDataset(group, "reference_angles_rad",
      std::vector<hsize_t>{
        static_cast<hsize_t>(referenceDirections_.size()), 2},
      referenceAngles, 0, "rad; columns=theta_x,theta_y");
    writeRealDataset(group, "window_centers_reduced_observer_time_s",
      std::vector<hsize_t>{static_cast<hsize_t>(windowCenters_.size())},
      windowCenters_, 0, "s");

    writeUnsignedAttribute(group, "format_version", 1);
    writeUnsignedAttribute(group, "samples", sampleCount_);
    writeUnsignedAttribute(group, "shots", config_.shots.size());
    writeUnsignedAttribute(group, "windows_per_shot", windowCenters_.size());
    writeUnsignedAttribute(group, "accumulator_bytes", estimatedBytes_);
    writeDoubleAttribute(group, "interval_start_s",
      config_.timeAverage.startTime);
    writeDoubleAttribute(group, "interval_end_s",
      config_.timeAverage.endTime);
    writeDoubleAttribute(group, "window_duration_s",
      config_.timeAverage.windowDuration);
    writeDoubleAttribute(group, "window_step_s",
      config_.timeAverage.windowStep);
    writeStringAttribute(group, "window", "Hann");
    writeStringAttribute(group, "sample_definition",
      "one Hann-windowed reduced-observer-time spectrum per shot and window");
    writeStringAttribute(group, "reduced_observer_time",
      "u=t_lab-n_dot_r_lab/c; the common propagation delay R/c is omitted");
    writeStringAttribute(group, "spatial_csd_definition",
      "W_ab(ref,target,omega)=mean_samples(conj(E_ref,a)*E_target,b)");
    writeStringAttribute(group, "degree_definition",
      "mu_EM^2=sum_ab(abs(W_ab)^2)/(trace(W_ref_ref)*trace(W_target_target))");
    writeStringAttribute(group, "coherent_mode_weighting",
      "flatten angle and polarization; diagonalize sqrt(dOmega_i)*W_ij*sqrt(dOmega_j)");
    writeStringAttribute(group, "temporal_csd_definition",
      "W_ab(omega1,omega2)=mean_samples(conj(E_a(omega1))*E_b(omega2))");
    writeStringAttribute(group, "stationarity_warning",
      "Interpret time windows as ensemble samples only inside a stationary or quasi-stationary interval");
    H5Gclose(group);
  }
}
