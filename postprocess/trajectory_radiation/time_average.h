#ifndef TRAJECTORY_RADIATION_TIME_AVERAGE_H
#define TRAJECTORY_RADIATION_TIME_AVERAGE_H

#include <cstddef>
#include <vector>

#include "hdf5.h"

#include "radiation_types.h"

namespace radiation
{
  class TimeAverageAccumulator
  {
  public:
    explicit TimeAverageAccumulator(const RadiationConfig& config);

    const std::vector<long double>& omega() const;
    const std::vector<double>& windowCenters() const;
    const std::vector<ObservationDirection>& referenceDirections() const;
    std::size_t estimatedBytes() const;
    unsigned long long sampleCount() const;

    void beginSample(const std::vector<long double>& referenceAmplitude);
    void accumulateAngularBlock(
        const std::vector<long double>& targetAmplitude,
        std::size_t thetaYOffset,
        std::size_t thetaYCount);
    void endSample();
    void write(hid_t file) const;

  private:
    std::size_t spatialIndex(std::size_t reference,
                             std::size_t frequency,
                             std::size_t thetaY,
                             std::size_t thetaX,
                             std::size_t referencePolarization,
                             std::size_t targetPolarization,
                             std::size_t realOrImag) const;
    std::size_t temporalIndex(std::size_t reference,
                              std::size_t firstFrequency,
                              std::size_t secondFrequency,
                              std::size_t firstPolarization,
                              std::size_t secondPolarization,
                              std::size_t realOrImag) const;

    RadiationConfig config_;
    std::vector<double> photonEnergyEV_;
    std::vector<long double> omega_;
    std::vector<double> thetaX_;
    std::vector<double> thetaY_;
    std::vector<double> windowCenters_;
    std::vector<ObservationDirection> referenceDirections_;
    std::vector<long double> currentReferenceField_;
    std::vector<long double> spatialCsd_;
    std::vector<long double> targetIntensity_;
    std::vector<long double> referenceIntensity_;
    std::vector<long double> sampleReferenceIntensity_;
    std::vector<long double> temporalCsd_;
    std::size_t expectedSamples_;
    std::size_t estimatedBytes_;
    unsigned long long sampleCount_;
    bool sampleOpen_;
  };
}

#endif
