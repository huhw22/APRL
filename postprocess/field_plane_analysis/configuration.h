#ifndef FIELD_PLANE_ANALYSIS_CONFIGURATION_H
#define FIELD_PLANE_ANALYSIS_CONFIGURATION_H

#include <cstddef>
#include <string>
#include <vector>

namespace field_analysis
{
  struct ReferenceAngle
  {
    double thetaX;
    double thetaY;
  };

  struct Configuration
  {
    std::string cardPath;
    std::string fieldFile;
    std::string baselineFile;
    bool requireComplete;

    double minimumPhotonEnergyEV;
    double maximumPhotonEnergyEV;
    bool hannTimeWindows;
    double intervalStart;
    double intervalEnd;
    double windowDuration;
    double windowStep;
    std::string transverseWindow;
    std::size_t angularPaddingY;
    std::size_t angularPaddingX;

    std::size_t frequencyBlock;
    std::size_t spatialBatchPoints;
    std::size_t fftThreads;
    std::size_t maximumWorkingMiB;
    double maximumOutputGiB;

    std::vector<double> spatialPhotonEnergyEV;
    std::vector<ReferenceAngle> spatialReferenceAngles;
    std::vector<double> temporalPhotonEnergyEV;
    std::vector<ReferenceAngle> temporalReferenceAngles;

    std::string outputFile;
    unsigned int compression;

    Configuration();
  };

  Configuration loadConfiguration(const std::string& filename);
}

#endif
