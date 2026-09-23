#ifndef TRAJECTORY_RADIATION_OUTPUT_WRITER_H
#define TRAJECTORY_RADIATION_OUTPUT_WRITER_H

#include <cstddef>
#include <string>
#include <vector>

#include "hdf5.h"

#include "radiation_types.h"

namespace radiation
{
  struct ComplexValue
  {
    double real;
    double imag;
  };

  class FarFieldWriter
  {
  public:
    explicit FarFieldWriter(const RadiationConfig& config);
    ~FarFieldWriter();

    FarFieldWriter(const FarFieldWriter&) = delete;
    FarFieldWriter& operator=(const FarFieldWriter&) = delete;

    void writeShotBlock(std::size_t shot,
                        const RadiationBlock& reducedBlock,
                        const std::vector<long double>& omega);
    void finalizeEnsemble();
    void close();

  private:
    hid_t createComplexType() const;
    hid_t createDataset(const char* name, hid_t type,
                        const std::vector<hsize_t>& dimensions,
                        const std::vector<hsize_t>& chunks) const;
    void writeAxis(const char* name, const std::vector<double>& values,
                   const char* unit) const;
    void writeStringAttribute(hid_t object, const char* name,
                              const std::string& value) const;
    void writeDoubleAttribute(hid_t object, const char* name,
                              double value) const;
    void writeUnsignedAttribute(hid_t object, const char* name,
                                unsigned long long value) const;
    void writeComplexBlock(hid_t dataset, std::size_t shot,
                           const RadiationBlock& block,
                           const std::vector<ComplexValue>& values) const;
    void writeScalarBlock(hid_t dataset, std::size_t shot,
                          const RadiationBlock& block,
                          const std::vector<double>& values,
                          std::size_t trailing = 0) const;
    void readComplexBlock(hid_t dataset, std::size_t shot,
                          const RadiationBlock& block,
                          std::vector<ComplexValue>& values) const;
    void readScalarBlock(hid_t dataset, std::size_t shot,
                         const RadiationBlock& block,
                         std::vector<double>& values) const;
    ComplexValue readFieldPoint(std::size_t shot, std::size_t frequency,
                                std::size_t thetaY, std::size_t thetaX,
                                std::size_t polarization) const;
    void finalizeSpectra();
    void finalizeCoherence();
    std::size_t nearestIndex(const std::vector<double>& values,
                             double requested) const;
    double solidAngleWeight(std::size_t thetaY,
                            std::size_t thetaX) const;

    RadiationConfig config_;
    std::vector<double> photonEnergyEV_;
    std::vector<double> thetaX_;
    std::vector<double> thetaY_;
    std::vector<long double> omega_;
    hid_t file_;
    hid_t group_;
    hid_t complexType_;
    hid_t shotField_;
    hid_t shotIntensity_;
    hid_t shotStokes_;
    hid_t meanField_;
    hid_t meanIntensity_;
    hid_t coherentIntensity_;
    hid_t incoherentIntensity_;
    hid_t meanStokes_;
    hid_t completeDataset_;
    bool open_;
  };
}

#endif
