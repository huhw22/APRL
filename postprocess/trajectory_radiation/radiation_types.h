#ifndef TRAJECTORY_RADIATION_TYPES_H
#define TRAJECTORY_RADIATION_TYPES_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace radiation
{
  namespace constants
  {
    const long double pi =
      3.141592653589793238462643383279502884L;
    const long double c = 299792458.0L;
    const long double epsilon0 =
      8.8541878128e-12L;
    const long double mu0 =
      1.25663706212e-6L;
    const long double elementaryCharge =
      1.602176634e-19L;
    const long double hbar =
      1.054571817e-34L;
  }

  struct Vec3
  {
    long double value[3];

    Vec3();
    Vec3(long double x, long double y, long double z);
    long double& operator[](std::size_t index);
    long double operator[](std::size_t index) const;
  };

  Vec3 operator+(const Vec3& left, const Vec3& right);
  Vec3 operator-(const Vec3& left, const Vec3& right);
  Vec3 operator*(long double scalar, const Vec3& vector);
  Vec3 operator/(const Vec3& vector, long double scalar);
  long double dot(const Vec3& left, const Vec3& right);
  Vec3 cross(const Vec3& left, const Vec3& right);
  long double norm(const Vec3& vector);
  Vec3 normalized(const Vec3& vector, const char* description);

  struct AxisDefinition
  {
    double minimum;
    double maximum;
    std::size_t count;
    bool logarithmic;

    AxisDefinition();
    std::vector<double> values() const;
  };

  struct ShotConfig
  {
    std::string name;
    std::vector<std::string> files;
  };

  struct CoherenceConfig
  {
    std::vector<double> spatialPhotonEnergyEV;
    std::vector<Vec3> spatialReferenceAngles;
    std::vector<double> temporalPhotonEnergyEV;
    std::vector<Vec3> temporalReferenceAngles;
  };

  struct TimeAverageConfig
  {
    bool enabled;
    double startTime;
    double endTime;
    double windowDuration;
    double windowStep;
    std::vector<double> photonEnergyEV;
    std::vector<Vec3> referenceAngles;
    std::size_t maximumAccumulatorMiB;

    TimeAverageConfig();
    std::vector<double> windowCenters() const;
  };

  struct RadiationConfig
  {
    std::string cardPath;
    std::vector<ShotConfig> shots;
    bool requireComplete;
    std::size_t readChunkRecords;

    Vec3 observationAxis;
    Vec3 horizontalAxis;
    Vec3 verticalAxis;
    double distanceM;
    AxisDefinition thetaX;
    AxisDefinition thetaY;
    AxisDefinition photonEnergyEV;

    std::size_t frequencyBlock;
    std::size_t thetaYBlock;
    std::size_t minimumRecordsPerParticle;

    std::string outputFile;
    unsigned int compression;
    bool overwrite;
    CoherenceConfig coherence;
    TimeAverageConfig timeAverage;

    RadiationConfig();
  };

  struct TrajectoryRecord
  {
    std::uint64_t particleId;
    std::uint64_t sourceId;
    double time;
    double position[3];
    double properVelocity[3];
    double charge;
    double weight;
    std::uint8_t event;
    std::int8_t boundaryFace;
  };

  struct ParticleTrajectory
  {
    std::uint64_t particleId;
    std::uint64_t sourceId;
    double charge;
    std::vector<TrajectoryRecord> records;
  };

  struct TrajectoryLoadStats
  {
    unsigned long long inputRecords;
    unsigned long long uniqueRecords;
    unsigned long long particles;
    unsigned long long terminalEvents;
    unsigned long long duplicateRecords;

    TrajectoryLoadStats();
  };

  struct ObservationDirection
  {
    Vec3 direction;
    Vec3 horizontal;
    Vec3 vertical;
    double thetaX;
    double thetaY;
  };

  struct ObserverTimeRange
  {
    long double minimum;
    long double maximum;
    unsigned long long internalKnots;

    ObserverTimeRange();
  };

  struct RadiationBlock
  {
    std::size_t frequencyOffset;
    std::size_t frequencyCount;
    std::size_t thetaYOffset;
    std::size_t thetaYCount;
    std::size_t thetaXCount;

    /* [frequency][theta_y][theta_x][polarization][real_or_imag] */
    std::vector<long double> amplitude;

    std::size_t scalarIndex(std::size_t frequency,
                            std::size_t thetaY,
                            std::size_t thetaX,
                            std::size_t polarization,
                            std::size_t realOrImag) const;
  };
}

#endif
