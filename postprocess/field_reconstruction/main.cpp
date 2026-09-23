#include <exception>
#include <iomanip>
#include <iostream>

#include "configuration.h"
#include "field_reconstruction.h"

int main(int argc, char** argv)
{
  if (argc != 2)
    {
      std::cerr << "Usage: field_reconstruction <configuration.yaml>\n";
      return 2;
    }
  try
    {
      const reconstruction::ReconstructionConfig config =
        reconstruction::loadConfiguration(argv[1]);
      const reconstruction::ReconstructionSummary result =
        reconstruction::reconstructField(config);
      std::cout << std::setprecision(10)
                << "Field reconstruction complete: particles="
                << result.acceptedParticles
                << ", mean_gamma=" << result.meanGamma
                << ", relative_gamma_spread="
                << result.relativeGammaSpread
                << ", rms_transverse_beta="
                << result.rmsTransverseBeta
                << ", outside_charge_fraction="
                << result.outsideChargeFraction << "\n"
                << "Forward energy: raw=" << result.rawForwardEnergy
                << " J, cleaned=" << result.cleanedForwardEnergy
                << " J, relative_change="
                << result.relativeForwardEnergyChange << "\n";
      return 0;
    }
  catch (const std::exception& error)
    {
      std::cerr << "Field reconstruction failure: " << error.what() << "\n";
      return 1;
    }
}
