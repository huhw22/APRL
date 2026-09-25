#ifndef FIELD_RECONSTRUCTION_H
#define FIELD_RECONSTRUCTION_H

#include "configuration.h"

namespace reconstruction
{
  struct ReconstructionSummary
  {
    unsigned long long inputRecords;
    unsigned long long acceptedParticles;
    unsigned long long duplicateParticles;
    unsigned long long invalidParticles;
    double outsideChargeFraction;
    double meanGamma;
    double relativeGammaSpread;
    double rmsTransverseBeta;
    double rawSignedEnergy;
    double cleanedSignedEnergy;
    double rawForwardEnergy;
    double cleanedForwardEnergy;
    double relativeForwardEnergyChange;

    ReconstructionSummary();
  };

  ReconstructionSummary reconstructField(
      const ReconstructionConfig& config);
}

#endif
