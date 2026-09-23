#ifndef DIRECT_EB_FIELD_H
#define DIRECT_EB_FIELD_H

#include <cstddef>
#include <memory>
#include <vector>

#include "boostframe.h"
#include "fieldvector.h"

namespace fel
{
  namespace SI
  {
    /* The E/B core is deliberately expressed in SI.  Unit conversion belongs
     * at the input/output boundary, not in Maxwell's update equations. */
    extern const Double c;
    extern const Double epsilon0;
    extern const Double mu0;
    extern const Double elementaryCharge;
    extern const Double electronMass;
  }

  enum class EBMaxwellSolver
  {
    Yee,
    CowanZ
  };

  struct EBCowanCoefficients
  {
    Double alpha[3];
    Double beta[3];
    Double deltaXY;
    Double deltaYZ;
    Double deltaZX;
    Double ratio[3];

    EBCowanCoefficients();
    static EBCowanCoefficients forZDispersion(Double dx, Double dy,
                                               Double dz);
  };

  struct EBGridGeometry
  {
    std::size_t nx;
    std::size_t ny;
    std::size_t nz;
    Double dx;
    Double dy;
    Double dz;
    Double dt;
    EBMaxwellSolver solver;

    EBGridGeometry();
    EBGridGeometry(std::size_t nxValue, std::size_t nyValue,
                   std::size_t nzValue, Double dxValue, Double dyValue,
                   Double dzValue, Double dtValue,
                   EBMaxwellSolver solverValue = EBMaxwellSolver::Yee);
  };

  /* A non-owning pair of one-cell electric-field halos.  Cowan's modified
   * Faraday operator smooths transverse to each derivative, so z-slab MPI
   * needs one E plane beyond each local slab.  Null pointers denote a
   * physical boundary, where the PEC parity extension is used. */
  struct EBElectricHaloView
  {
    const Double* lowerEx;
    const Double* lowerEy;
    const Double* lowerEz;
    const Double* upperEx;
    const Double* upperEy;
    const Double* upperEz;

    EBElectricHaloView();
  };

  /* A compact scalar component on a rectangular Yee lattice.  Components are
   * kept in separate contiguous arrays (SoA) so the hot curl loops stream only
   * the values they need.  x is the unit-stride dimension and a complete z
   * slab remains contiguous for the current MPI decomposition. */
  class YeeComponent
  {
  public:
    YeeComponent();
    YeeComponent(std::size_t nx, std::size_t ny, std::size_t nz);

    void resize(std::size_t nx, std::size_t ny, std::size_t nz);
    void fill(Double value);

    Double& operator()(std::size_t i, std::size_t j, std::size_t k);
    const Double& operator()(std::size_t i, std::size_t j,
                             std::size_t k) const;

    std::size_t nx() const;
    std::size_t ny() const;
    std::size_t nz() const;
    std::size_t size() const;
    std::size_t bytes() const;

    Double* data();
    const Double* data() const;

  private:
    std::size_t index(std::size_t i, std::size_t j, std::size_t k) const;

    std::size_t nx_;
    std::size_t ny_;
    std::size_t nz_;
    std::vector<Double> values_;
  };

  struct RadiationFieldSample
  {
    FieldVector<Double> electric;       /* V/m */
    FieldVector<Double> magnetic;       /* T */
    FieldVector<Double> poynting;       /* W/m^2 */
    Double energyDensity;               /* J/m^3 */

    RadiationFieldSample();
  };

  struct EBMemoryFootprint
  {
    std::size_t electricBytes;
    std::size_t magneticBytes;
    std::size_t currentBytes;
    std::size_t boundaryBytes;

    EBMemoryFootprint();
    std::size_t totalBytes() const;
  };

  class EBFieldGrid;
  class EBConvolutionalPML;

  /* Boundary implementations are deliberately outside the curl kernels.
   * CPML therefore adds only boundary-slab auxiliary arrays without changing
   * the E/B storage or the vacuum update. */
  class EBBoundaryOperator
  {
  public:
    virtual ~EBBoundaryOperator();
    virtual void afterMagneticUpdate(
        EBFieldGrid& fields, const EBElectricHaloView& halo) = 0;
    virtual void afterElectricUpdate(EBFieldGrid& fields) = 0;
    virtual std::size_t memoryBytes() const = 0;
  };

  class PerfectElectricConductorBoundary : public EBBoundaryOperator
  {
  public:
    PerfectElectricConductorBoundary(bool lowerZPhysical = true,
                                     bool upperZPhysical = true);
    virtual void afterMagneticUpdate(
        EBFieldGrid& fields, const EBElectricHaloView& halo);
    virtual void afterElectricUpdate(EBFieldGrid& fields);
    virtual std::size_t memoryBytes() const;

  private:
    bool lowerZPhysical_;
    bool upperZPhysical_;
  };

  class EBFieldGrid
  {
  public:
    explicit EBFieldGrid(const EBGridGeometry& geometry);

    EBFieldGrid(const EBFieldGrid&) = delete;
    EBFieldGrid& operator=(const EBFieldGrid&) = delete;

    const EBGridGeometry& geometry() const;

    YeeComponent& ex();
    YeeComponent& ey();
    YeeComponent& ez();
    YeeComponent& bx();
    YeeComponent& by();
    YeeComponent& bz();
    YeeComponent& jx();
    YeeComponent& jy();
    YeeComponent& jz();

    const YeeComponent& ex() const;
    const YeeComponent& ey() const;
    const YeeComponent& ez() const;
    const YeeComponent& bx() const;
    const YeeComponent& by() const;
    const YeeComponent& bz() const;
    const YeeComponent& jx() const;
    const YeeComponent& jy() const;
    const YeeComponent& jz() const;

    void setBoundary(std::unique_ptr<EBBoundaryOperator> boundary);

    void clearFields();
    void clearCurrent();

    /* Leapfrog Maxwell update in SI:
     *   B^(n+1/2) = B^(n-1/2) - dt curl(E^n)
     *   E^(n+1)   = E^n + dt (curl(B^(n+1/2))/epsilon0 - J/epsilon0)
     * No A/phi value is stored or reconstructed anywhere in this path. */
    void advanceMagnetic();
    void advanceMagnetic(const EBElectricHaloView& halo);
    void advanceElectric();
    void advance();

    Double courantNumber() const;
    Double courantLimit() const;
    const EBCowanCoefficients& cowanCoefficients() const;
    static Double axisPhaseVelocityRatio(Double courantAxis,
                                         Double cellsPerWavelength);
    static Double axisGroupVelocityRatio(Double courantAxis,
                                         Double cellsPerWavelength);
    EBMemoryFootprint memoryFootprint() const;

    /* Return E and B at a cell centre.  The interpolation reconciles the Yee
     * staggering without allocating a second, collocated full-domain copy. */
    RadiationFieldSample radiationSampleCell(
        std::size_t i, std::size_t j, std::size_t k) const;

    /* Transform the sampled field into the lab frame before computing the
     * Poynting vector.  This is the preferred radiation-output boundary. */
    RadiationFieldSample radiationSampleCellLab(
        std::size_t i, std::size_t j, std::size_t k,
        const BoostFrameTransform& frame) const;

    /* Trilinear interpolation at a physical SI position. Component-specific
     * half-cell offsets are handled internally. B is returned at its stored
     * leapfrog half time. */
    void sampleFieldsPosition(
        const FieldVector<Double>& positionSI,
        const FieldVector<Double>& originSI,
        FieldVector<Double>& electric,
        FieldVector<Double>& magnetic) const;

    RadiationFieldSample radiationSamplePositionLab(
        const FieldVector<Double>& positionBoxSI,
        const FieldVector<Double>& originBoxSI,
        const BoostFrameTransform& frame) const;

  private:
    friend class EBConvolutionalPML;

    void validateGeometry() const;
    void advanceMagneticYee();
    void advanceMagneticCowan(const EBElectricHaloView& halo);
    Double electricValue(const YeeComponent& component,
                         unsigned int componentAxis,
                         std::ptrdiff_t i, std::ptrdiff_t j,
                         std::ptrdiff_t k,
                         const EBElectricHaloView& halo) const;
    Double smoothedElectric(const YeeComponent& component,
                            unsigned int componentAxis,
                            unsigned int derivativeAxis,
                            std::ptrdiff_t i, std::ptrdiff_t j,
                            std::ptrdiff_t k,
                            const EBElectricHaloView& halo) const;
    RadiationFieldSample makeRadiationSample(
        const FieldVector<Double>& electric,
        const FieldVector<Double>& magnetic) const;
    Double sampleComponent(const YeeComponent& component,
                           const FieldVector<Double>& positionSI,
                           const FieldVector<Double>& originSI,
                           Double offsetX, Double offsetY,
                           Double offsetZ) const;

    EBGridGeometry geometry_;
    EBCowanCoefficients cowan_;

    YeeComponent ex_;
    YeeComponent ey_;
    YeeComponent ez_;
    YeeComponent bx_;
    YeeComponent by_;
    YeeComponent bz_;

    /* Current density is collocated with the matching electric component. */
    YeeComponent jx_;
    YeeComponent jy_;
    YeeComponent jz_;

    std::unique_ptr<EBBoundaryOperator> boundary_;
  };
}

#endif
