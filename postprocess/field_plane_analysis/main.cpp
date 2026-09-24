#include <exception>
#include <iomanip>
#include <iostream>

#include "configuration.h"
#include "field_plane_analysis.h"

int main(int argc, char** argv)
{
  if (argc != 2)
    {
      std::cerr << "Usage: field_plane_analysis <analysis.yaml>\n";
      return 2;
    }
  try
    {
      const field_analysis::Configuration config =
        field_analysis::loadConfiguration(argv[1]);
      const field_analysis::AnalysisSummary summary =
        field_analysis::analyze(config);
      std::cout << std::setprecision(10)
                << "Field-plane analysis complete: ensemble_samples="
                << summary.ensembleSamples
                << ", frequency_bins=" << summary.frequencyBins
                << ", dt=" << summary.uniformTimeStep << " s"
                << ", Nyquist=" << summary.nyquistPhotonEnergyEV << " eV"
                << ", mean_band_energy=" << summary.meanBandEnergy << " J"
                << ", band_coherent_fraction="
                << summary.bandCoherentFraction
                << ", working_estimate=" << summary.estimatedWorkingMiB
                << " MiB, output_estimate=" << summary.estimatedOutputGiB
                << " GiB\n";
      if (summary.ensembleSamples == 1)
        std::cout << "Coherence note: one deterministic sample is rank one; "
                     "use Hann windows from a stationary interval to measure "
                     "time-averaged partial coherence.\n";
      return 0;
    }
  catch (const std::exception& error)
    {
      std::cerr << "Field-plane analysis failed: " << error.what() << "\n";
      return 1;
    }
}
