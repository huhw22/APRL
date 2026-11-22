#ifndef APRL_FINITE_RADIATION_SCREEN_H
#define APRL_FINITE_RADIATION_SCREEN_H

#include <array>
#include <cstddef>
#include <vector>

namespace aprl
{
  // Physical lab endpoints of a committed forced trajectory segment.
  // Position in m, time in s, proper velocity u=p/(m*c), charge in C.
  struct RadiationSegment
  {
    std::array<double,3> x0, x1, u0, u1;
    double t0, t1, charge;
  };

  struct RadiationScreenSpec
  {
    std::size_t nx=0, ny=0;
    double x0=0., y0=0., dx=0., dy=0., z=0.; // x0,y0 are pixel centres
    std::vector<double> frequenciesHz; // optional hints saved for offline processing only
    double phaseOrigin=0.; // absolute lab time reference for offline processing
    std::size_t nt=0;
    double timeStart=0., timeStep=0.; // absolute lab first bin edge and bin width
  };

  class FiniteRadiationScreen
  {
  public:
    explicit FiniteRadiationScreen(const RadiationScreenSpec& spec);
    // Serialize calls on each instance. OpenMP partitions the pixels within a call.
    void accumulate(const std::vector<RadiationSegment>& segments);
    const RadiationScreenSpec& spec() const { return spec_; }
    // [ny,nx,nt,6]: same components, in V/m, mean over each arrival-time bin.
    const std::vector<double>& timeField() const { return timeField_; }
    // Finalization hands the storage to the caller for in-place MPI reduction.
    // Once the output has been taken, accumulation is closed.
    std::vector<double> releaseTimeField();
    std::size_t storageBytes() const;
  private:
    bool released_=false;
    RadiationScreenSpec spec_;
    std::vector<double> timeField_;
  };
}
#endif
