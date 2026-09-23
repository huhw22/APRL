#ifndef FIELD_RECONSTRUCTION_CONFIGURATION_H
#define FIELD_RECONSTRUCTION_CONFIGURATION_H

#include <cstddef>
#include <string>

namespace reconstruction
{
  struct ReconstructionConfig
  {
    std::string cardPath;
    std::string fieldFile;
    std::string particleFile;
    bool requireComplete;
    std::size_t particleReadChunk;

    unsigned int fftThreads;
    std::size_t paddingTime;
    std::size_t paddingY;
    std::size_t paddingX;
    double transverseSmoothing;
    double maximumRelativeGammaSpread;
    double maximumRmsTransverseBeta;
    double maximumOutsideChargeFraction;

    std::string outputFile;
    unsigned int compression;

    ReconstructionConfig();
  };

  ReconstructionConfig loadConfiguration(const std::string& filename);
}

#endif
