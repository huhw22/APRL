#ifndef APRL_RADIATION_SCREEN_OUTPUT_H
#define APRL_RADIATION_SCREEN_OUTPUT_H

#include <mpi.h>
#include <hdf5.h>
#include <string>
#include "boostframe.h"
#include "finite_radiation_screen.h"

namespace aprl
{
  class RadiationScreenOutput
  {
  public:
    RadiationScreenOutput(const RadiationScreenSpec& spec, const std::string& path,
                          bool overwrite, MPI_Comm comm);
    ~RadiationScreenOutput();
    void append(const std::vector<RadiationSegment>& boxSubsteps,
                double boxStop, const BoostFrameTransform& frame);
    void close(bool complete);
  private:
    RadiationScreenOutput(const RadiationScreenOutput&) = delete;
    RadiationScreenOutput& operator=(const RadiationScreenOutput&) = delete;
    void flush();
    FiniteRadiationScreen screen_;
    MPI_Comm communicator_;
    int rank_;
    hid_t file_;
    unsigned long long committed_=0, accelerated_=0, clipped_=0;
    double accumulationSeconds_=0;
    std::vector<RadiationSegment> batch_;
  };
}
#endif
