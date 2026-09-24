#ifndef FIELD_PLANE_ANALYSIS_H
#define FIELD_PLANE_ANALYSIS_H

#include "configuration.h"

namespace field_analysis
{
  struct AnalysisSummary
  {
    unsigned long long inputSamples;
    unsigned long long ensembleSamples;
    unsigned long long frequencyBins;
    double uniformTimeStep;
    double nyquistPhotonEnergyEV;
    double meanBandEnergy;
    double bandCoherentFraction;
    double estimatedWorkingMiB;
    double estimatedOutputGiB;

    AnalysisSummary();
  };

  AnalysisSummary analyze(const Configuration& config);
}

#endif
