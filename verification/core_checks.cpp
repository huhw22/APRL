#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

#include "boostframe.h"
#include "config.h"
#include "eb_bunch.h"
#include "eb_cpml.h"
#include "eb_deposition.h"
#include "eb_field.h"
#include "eb_initial_field.h"
#include "eb_mpi.h"
#include "eb_particles.h"
#include "longitudinal_mesh.h"
#include "particle_file.h"

namespace
{
  const double pi = 3.141592653589793238462643383279502884;

  void require(bool condition, const std::string& message)
  {
    if (!condition) throw std::runtime_error(message);
  }

  double relativeError(double value, double reference)
  {
    return std::abs(value - reference) /
      std::max(1.0, std::max(std::abs(value), std::abs(reference)));
  }

  void checkRuntimeDefaults(int rank)
  {
    const aprl::RuntimeConfig runtime;
    require(!aprl::InitialSelfFieldConfig().includeTransverseCpmlParticles,
      "Transverse CPML initial sources must remain opt-in");
    require(runtime.shutdownWatchdogEnabled,
      "Shutdown watchdog default unexpectedly disabled");
    require(runtime.shutdownWatchdogTimeoutSeconds == 600,
      "General shutdown watchdog default changed");
    require(runtime.shutdownWatchdogFinalizeTimeoutSeconds == 60,
      "MPI_Finalize watchdog default changed");
    if (rank == 0)
      std::cout << "shutdown_watchdog_timeout_s=600 "
                   "finalize_timeout_s=60" << std::endl;
  }

  void checkInitialFieldTransverseLimit(const std::string& validCard,
                                        const std::string& invalidCard,
                                        MPI_Comm communicator, int rank)
  {
    require(aprl::InitialSelfFieldConfig().maximumMeanTransverseBeta ==
              1.0e-8,
      "Initial self-field transverse limit default changed");
    const aprl::SimulationConfig config =
      aprl::YamlConfigLoader::loadFile(validCard);
    require(config.initialSelfField.maximumMeanTransverseBeta == 1.0e-7,
      "Initial self-field transverse limit was not read from YAML");
    bool invalidRejected = false;
    try
      {
        aprl::YamlConfigLoader::loadFile(invalidCard);
      }
    catch (const std::exception& error)
      {
        invalidRejected = std::string(error.what()).find(
          "maximum_mean_transverse_beta") != std::string::npos;
      }
    require(invalidRejected,
      "Negative initial self-field transverse limit was not rejected");

    const aprl::EBGridGeometry geometry(
      8, 8, 8, 0.1, 0.1, 0.1, 0.05 / aprl::SI::c);
    aprl::EBFieldGrid fields(geometry);
    const aprl::FieldVector<double> origin(0.0);
    const std::size_t guards[3] = {0, 0, 0};
    std::vector<aprl::RelativisticParticleSI> particles(1);
    aprl::RelativisticParticleSI& particle = particles[0];
    particle.position[0] = 0.45;
    particle.position[1] = 0.45;
    particle.position[2] = 0.45;
    particle.properVelocity[0] = 7.5e-8;
    particle.properVelocity[2] = 0.2;
    particle.charge = -aprl::SI::elementaryCharge;
    particle.mass = aprl::SI::electronMass;
    particle.weight = 1.0;

    bool defaultRejected = false;
    try
      {
        aprl::EBGaussFieldInitializer::initialize(
          fields, geometry, 0, origin, particles,
          aprl::InitialSelfFieldModel::RelativisticPoisson,
          1.0e-8, guards, 1.0e-8, 1000, communicator);
      }
    catch (const std::runtime_error& error)
      {
        const std::string message(error.what());
        defaultRejected = message.find(
          "maximum_mean_transverse_beta=1e-08") != std::string::npos;
      }
    require(defaultRejected,
      "Default transverse limit did not reject the tilted beam");

    const aprl::EBGaussInitializationReport report =
      aprl::EBGaussFieldInitializer::initialize(
        fields, geometry, 0, origin, particles,
        aprl::InitialSelfFieldModel::RelativisticPoisson,
        config.initialSelfField.maximumMeanTransverseBeta,
        guards, 1.0e-8, 1000, communicator);
    require(report.meanBetaTransverse > 1.0e-8 &&
            report.meanBetaTransverse < 1.0e-7,
      "Relaxed transverse limit did not admit the expected mean velocity");
    if (rank == 0)
      std::cout << std::setprecision(12)
                << "initial_field_mean_beta_transverse="
                << report.meanBetaTransverse
                << " default_limit=1e-8 relaxed_limit=1e-7"
                << std::endl;
  }

  void checkCowanInitialField(MPI_Comm communicator, int rank, int size)
  {
    require(size == 1,
      "cowan-initial-field verification requires exactly one MPI rank");
    const aprl::EBGridGeometry geometry(
      12, 10, 14, 0.15, 0.12, 0.1, 0.1 / aprl::SI::c,
      aprl::EBMaxwellSolver::CowanZ);
    const aprl::FieldVector<double> origin(0.0);
    const std::size_t guards[3] = {0, 0, 0};
    const double positions[3][3] = {
      {0.66, 0.43, 0.65}, {0.91, 0.71, 0.88}, {0.81, 0.62, 0.44}
    };
    const double velocities[2] = {0.0, 0.23};
    for (double beta : velocities)
      {
        aprl::EBFieldGrid fields(geometry);
        std::vector<aprl::RelativisticParticleSI> particles(3);
        for (std::size_t index = 0; index < particles.size(); ++index)
          {
            for (unsigned int axis = 0; axis < 3; ++axis)
              particles[index].position[axis] = positions[index][axis];
            particles[index].properVelocity[2] =
              beta / std::sqrt(1.0 - beta * beta);
            particles[index].charge = -1.0e-12 * (index + 1);
            particles[index].mass = aprl::SI::electronMass * (index + 1);
            particles[index].weight = index + 1;
            particles[index].id = index + 10;
            particles[index].sourceId = index + 20;
          }
        const std::vector<aprl::RelativisticParticleSI> original(particles);
        const aprl::EBGaussInitializationReport report =
          aprl::EBGaussFieldInitializer::initialize(
            fields, geometry, 0, origin, particles,
            aprl::InitialSelfFieldModel::CowanRelativisticPoisson,
            1.0e-8, guards, 1.0e-11, 2000, communicator);
        require(report.relativeResidual < 4.0e-11,
          "Cowan initial field failed its physical Gauss check");
        require(std::abs(report.meanBetaZ - beta) < 2.0e-15,
          "Cowan initial field did not retain the mean source velocity");
        require(std::isfinite(report.initialElectricEnergy) &&
                report.initialElectricEnergy > 0.0,
          "Cowan initial field has no finite electric energy");
        for (std::size_t index = 0; index < particles.size(); ++index)
          {
            for (unsigned int axis = 0; axis < 3; ++axis)
              require(particles[index].position[axis] ==
                        original[index].position[axis] &&
                      particles[index].properVelocity[axis] ==
                        original[index].properVelocity[axis],
                "Cowan initialization changed a particle position or momentum");
            require(particles[index].charge == original[index].charge &&
                    particles[index].mass == original[index].mass &&
                    particles[index].weight == original[index].weight &&
                    particles[index].id == original[index].id &&
                    particles[index].sourceId == original[index].sourceId &&
                    particles[index].frozen == original[index].frozen &&
                    particles[index].frozenBetaZ == original[index].frozenBetaZ,
              "Cowan initialization changed particle source or held state");
          }

        if (beta == 0.0)
          {
            require(report.initialMagneticEnergy == 0.0,
              "Stationary Cowan source initialized a magnetic field");
            const aprl::YeeComponent* electric[3] = {
              &fields.ex(), &fields.ey(), &fields.ez()
            };
            const aprl::YeeComponent* magnetic[3] = {
              &fields.bx(), &fields.by(), &fields.bz()
            };
            std::vector<double> before[3];
            double electricScale = 0.0;
            for (unsigned int axis = 0; axis < 3; ++axis)
              {
                before[axis].assign(electric[axis]->data(),
                  electric[axis]->data() + electric[axis]->size());
                for (double value : before[axis])
                  electricScale = std::max(electricScale, std::abs(value));
              }
            require(electricScale > 0.0,
              "Stationary Cowan invariant has a vacuous zero-field source");
            // A static charge distribution with J=0 must remain static under
            // the production Maxwell update, including PEC boundary cells.
            // This tests curl compatibility without reproducing its stencil.
            fields.clearCurrent();
            fields.advance();
            double electricChange = 0.0;
            double magneticMaximum = 0.0;
            for (unsigned int axis = 0; axis < 3; ++axis)
              {
                for (std::size_t index = 0;
                     index < electric[axis]->size(); ++index)
                  {
                    require(std::isfinite(electric[axis]->data()[index]),
                      "Stationary Cowan advance produced a non-finite electric field");
                    electricChange = std::max(electricChange,
                      std::abs(electric[axis]->data()[index] -
                               before[axis][index]));
                  }
                for (std::size_t index = 0;
                     index < magnetic[axis]->size(); ++index)
                  {
                    require(std::isfinite(magnetic[axis]->data()[index]),
                      "Stationary Cowan advance produced a non-finite magnetic field");
                    magneticMaximum = std::max(magneticMaximum,
                      std::abs(magnetic[axis]->data()[index]));
                  }
              }
            const double relativeElectricChange = electricChange / electricScale;
            const double relativeMagneticField =
              aprl::SI::c * magneticMaximum / electricScale;
            require(relativeElectricChange < 5.0e-12 &&
                    relativeMagneticField < 5.0e-12,
              "Stationary Cowan initial field changed under zero-current Maxwell advance");
            if (rank == 0)
              std::cout << std::setprecision(12)
                        << "cowan_static_relative_e_change="
                        << relativeElectricChange
                        << " cowan_static_relative_cB="
                        << relativeMagneticField << std::endl;
          }
        else
          require(std::isfinite(report.initialMagneticEnergy) &&
                  report.initialMagneticEnergy > 0.0,
            "Moving Cowan source did not initialize a finite magnetic field");
        if (rank == 0)
          std::cout << std::setprecision(12)
                    << "cowan_initial_beta=" << beta
                    << " physical_gauss_residual=" << report.relativeResidual
                    << " electric_energy=" << report.initialElectricEnergy
                    << " magnetic_energy=" << report.initialMagneticEnergy
                    << std::endl;
      }
  }

  void checkCpmlParticleInitialField(MPI_Comm communicator, int rank)
  {
    const aprl::EBGridGeometry geometry(
      24, 24, 24, 0.1, 0.1, 0.1, 0.05 / aprl::SI::c);
    aprl::EBFieldGrid fields(geometry);
    const aprl::FieldVector<double> origin(0.0);
    const std::size_t guards[3] = {4, 4, 4};
    std::vector<aprl::RelativisticParticleSI> particles(2);
    for (std::size_t index = 0; index < particles.size(); ++index)
      {
        particles[index].position[0] = index == 0 ? 1.2 : 0.15;
        particles[index].position[1] = 1.2;
        particles[index].position[2] = 1.2;
        particles[index].properVelocity[2] = 0.2;
        particles[index].charge = -1.0e-15;
        particles[index].mass = aprl::SI::electronMass;
        particles[index].weight = 1.0;
      }

    std::vector<aprl::RelativisticParticleSI> insideOnly(1, particles[0]);
    const aprl::EBGaussInitializationReport first =
      aprl::EBGaussFieldInitializer::initialize(
        fields, geometry, 0, origin, insideOnly,
        aprl::InitialSelfFieldModel::RelativisticPoisson,
        1.0e-8, guards, 1.0e-8, 2000, communicator, true);
    const double withoutCpml = fields.ex()(5, 12, 12);
    const aprl::EBGaussInitializationReport second =
      aprl::EBGaussFieldInitializer::initialize(
        fields, geometry, 0, origin, particles,
        aprl::InitialSelfFieldModel::RelativisticPoisson,
        1.0e-8, guards, 1.0e-8, 2000, communicator, true);
    require(first.transversePoissonPaddingCells == 2 &&
            second.transversePoissonPaddingCells == 2,
      "CPML-source solve did not pad its temporary Poisson lattice");
    require(std::abs(fields.ex()(5, 12, 12) - withoutCpml) > 1.0e-7,
      "Transverse CPML particle did not affect the physical initial field");
    require(fields.ex()(0, 12, 12) == 0.0 &&
            fields.ey()(0, 12, 12) == 0.0 &&
            fields.bx()(0, 12, 12) == 0.0,
      "CPML evolution fields did not remain zero at initialization");
    require(second.relativeResidual < 4.0e-8,
      "CPML-source initial field failed its physical Gauss check");
    if (rank == 0)
      std::cout << "cpml_particle_initial_delta_ex="
                << fields.ex()(5, 12, 12) - withoutCpml
                << " physical_gauss_residual=" << second.relativeResidual
                << std::endl;
  }

  void checkLongitudinalLayout(int rank)
  {
    const double gamma = 929.13054928585427;
    const double gammaBeta =
      std::sqrt(gamma - 1.0) * std::sqrt(gamma + 1.0);
    const double deltaLab = 6.07998533158;
    const double expectedAdvance = deltaLab / (gamma + gammaBeta);
    const double advance =
      aprl::forwardLightIntersectionAdvance(deltaLab, gamma);
    require(std::abs(advance - expectedAdvance) < 2.0e-18,
      "Stable forward-light intersection advance changed");
    require(aprl::forwardLightIntersectionAdvance(2.5, 1.0) == 2.5,
      "Gamma-one forward-light intersection is not the lab distance");
    const double extremeAdvance =
      aprl::forwardLightIntersectionAdvance(1.0, 1.0e308);
    require(extremeAdvance > 0.0 && extremeAdvance < 1.0e-307,
      "Forward-light intersection overflowed at extreme finite gamma");

    aprl::LongitudinalLayoutInput input;
    input.initialMinimum = -2.0;
    input.initialMaximum = 2.0;
    input.initialMean = 0.0;
    input.initialRms = 0.5;
    input.ballisticMinimum = -3.0;
    input.ballisticMaximum = 5.0;
    input.radiationUpper = 6.0;
    input.dz = 0.25;
    input.cpmlCells = 4;
    input.leftGuardCells = 2;
    input.rightGuardCells = 3;
    input.selfFieldPaddingRms = 8.0;
    input.alignmentCells = 16;
    input.meshCenter = 10.0;
    const aprl::LongitudinalLayout layout =
      aprl::deriveLongitudinalLayout(input);

    require(layout.minimumPhysicalCells == 45,
      "Longitudinal required physical cell count changed");
    require(layout.totalCells == 64 && layout.physicalCells == 56,
      "Longitudinal aligned mesh cell count changed");
    require(layout.alignmentExtraCells == 11,
      "Longitudinal alignment slack changed");
    require(std::abs(layout.extent - 16.0) < 1.0e-14 &&
            std::abs(layout.origin - 2.0) < 1.0e-14,
      "Longitudinal outer mesh placement changed");
    require(std::abs(layout.physicalLower - 3.0) < 1.0e-14 &&
            std::abs(layout.physicalUpper - 17.0) < 1.0e-14,
      "Longitudinal CPML inner boundaries changed");
    require(std::abs(layout.translation - 7.5) < 1.0e-14,
      "Longitudinal bunch translation changed");
    require(std::abs(layout.requiredLower - layout.physicalLower) < 1.0e-14,
      "Longitudinal alignment slack leaked onto the left");
    require(std::abs(layout.requiredUpper - 14.25) < 1.0e-14 &&
            std::abs(layout.rightRoundingPadding - 2.75) < 1.0e-14,
      "Longitudinal right-only alignment padding changed");
    require(std::abs(layout.selfFieldLower - 3.5) < 1.0e-14 &&
            std::abs(layout.selfFieldUpper - 11.5) < 1.0e-14,
      "Longitudinal self-field RMS padding changed");
    require(std::abs(layout.ballisticMinimum - 4.5) < 1.0e-14 &&
            std::abs(layout.ballisticMaximum - 12.5) < 1.0e-14 &&
            std::abs(layout.radiationUpper - 13.5) < 1.0e-14,
      "Translated longitudinal causal boundaries changed");

    aprl::LongitudinalLayoutInput clippedInput(input);
    clippedInput.selfFieldPaddingRms = 0.0;
    const aprl::LongitudinalLayout clippedLayout =
      aprl::deriveLongitudinalLayout(clippedInput);
    require(std::abs(clippedLayout.selfFieldLower -
                     clippedLayout.initialMean) < 1.0e-14 &&
            std::abs(clippedLayout.selfFieldUpper -
                     clippedLayout.initialMean) < 1.0e-14,
      "Zero self-field padding did not collapse the optional RMS envelope");
    require(clippedLayout.minimumPhysicalCells <
              layout.minimumPhysicalCells,
      "Zero self-field padding still enlarges the automatic physical box");
    require(clippedLayout.initialMinimum >=
              clippedLayout.physicalLower &&
            clippedLayout.initialMaximum <=
              clippedLayout.physicalUpper &&
            clippedLayout.ballisticMinimum >=
              clippedLayout.physicalLower &&
            clippedLayout.radiationUpper <=
              clippedLayout.physicalUpper,
      "Box-clipped self-field layout lost sampled or causal content");

    /* Automatic placement is a passive translation of the computational
     * z coordinate.  The box origin must move by the same amount: box time
     * then remains unchanged, box z gains exactly the translation, and the
     * represented laboratory event/worldline is invariant.  Exercise a
     * high-gamma case close to the production FEL cards so a future origin
     * sign error cannot hide behind the generic Lorentz round-trip test. */
    aprl::BoostFrameTransform originalFrame;
    originalFrame.setOriginsFromGamma(
      gamma, aprl::SI::c, 0.0, -0.05, 0.0);
    aprl::BoostFrameTransform translatedFrame;
    translatedFrame.setOriginsFromGamma(
      gamma, aprl::SI::c, 0.0, -0.05, layout.translation);
    const double labTimes[] = {-1.7e-10, 0.0, 2.4e-10};
    const double labPositions[] = {-0.04, 1.25, 3.978};
    double maximumTranslationEventError = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
      {
        double originalTime = 0.0;
        double originalZ = 0.0;
        double translatedTime = 0.0;
        double translatedZ = 0.0;
        originalFrame.labToBox(
          labTimes[i], labPositions[i], originalTime, originalZ);
        translatedFrame.labToBox(
          labTimes[i], labPositions[i], translatedTime, translatedZ);
        maximumTranslationEventError = std::max(
          maximumTranslationEventError,
          std::abs(translatedTime - originalTime));
        maximumTranslationEventError = std::max(
          maximumTranslationEventError,
          std::abs((translatedZ - layout.translation) - originalZ));

        double originalRecoveredTime = 0.0;
        double originalRecoveredZ = 0.0;
        double translatedRecoveredTime = 0.0;
        double translatedRecoveredZ = 0.0;
        originalFrame.boxToLab(
          originalTime, originalZ,
          originalRecoveredTime, originalRecoveredZ);
        translatedFrame.boxToLab(
          translatedTime, translatedZ,
          translatedRecoveredTime, translatedRecoveredZ);
        maximumTranslationEventError = std::max(
          maximumTranslationEventError,
          std::abs(translatedRecoveredTime - originalRecoveredTime));
        maximumTranslationEventError = std::max(
          maximumTranslationEventError,
          std::abs(translatedRecoveredZ - originalRecoveredZ));
      }
    const double sampleBoxTime = 3.0e-11;
    const double sampleBoxZ = -0.0125;
    double originalLabTime = 0.0;
    double originalLabZ = 0.0;
    double translatedLabTime = 0.0;
    double translatedLabZ = 0.0;
    originalFrame.boxToLab(
      sampleBoxTime, sampleBoxZ, originalLabTime, originalLabZ);
    translatedFrame.boxToLab(
      sampleBoxTime, sampleBoxZ + layout.translation,
      translatedLabTime, translatedLabZ);
    maximumTranslationEventError = std::max(
      maximumTranslationEventError,
      std::abs(translatedLabTime - originalLabTime));
    maximumTranslationEventError = std::max(
      maximumTranslationEventError,
      std::abs(translatedLabZ - originalLabZ));
    require(maximumTranslationEventError < 2.0e-12,
      "Automatic longitudinal translation changed a laboratory event");

    aprl::LongitudinalLayoutInput shifted(input);
    shifted.initialMinimum += 100.0;
    shifted.initialMaximum += 100.0;
    shifted.initialMean += 100.0;
    shifted.ballisticMinimum += 100.0;
    shifted.ballisticMaximum += 100.0;
    shifted.radiationUpper += 100.0;
    const aprl::LongitudinalLayout shiftedLayout =
      aprl::deriveLongitudinalLayout(shifted);
    require(shiftedLayout.totalCells == layout.totalCells &&
            std::abs(shiftedLayout.translation -
              (layout.translation - 100.0)) < 1.0e-13 &&
            std::abs(shiftedLayout.requiredUpper -
              layout.requiredUpper) < 1.0e-13,
      "Longitudinal layout is not invariant under input-coordinate shifts");

    bool rejected = false;
    try
      {
        shifted.alignmentCells = 0;
        aprl::deriveLongitudinalLayout(shifted);
      }
    catch (const std::invalid_argument&)
      {
        rejected = true;
      }
    require(rejected, "Longitudinal layout accepted zero-cell alignment");

    if (rank == 0)
      std::cout << std::setprecision(12)
                << "forward_light_advance_m=" << advance
                << " longitudinal_nz=" << layout.totalCells
                << " longitudinal_extent=" << layout.extent
                << " longitudinal_translation=" << layout.translation
                << " translation_event_error="
                << maximumTranslationEventError
                << " longitudinal_right_padding="
                << layout.rightRoundingPadding << std::endl;
  }

  void checkLorentzAndElegant(MPI_Comm communicator, int rank, int size)
  {
    aprl::BoostFrameTransform frame;
    frame.setOriginsFromGamma(37.0, aprl::SI::c,
                              2.5e-12, -0.014, 0.003);
    const double times[] = {-1.0e-11, 2.5e-12, 9.0e-12};
    const double positions[] = {-0.12, -0.014, 0.21};
    double maximumCoordinateError = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
      {
        double boxTime = 0.0;
        double boxZ = 0.0;
        double labTime = 0.0;
        double labZ = 0.0;
        frame.labToBox(times[i], positions[i], boxTime, boxZ);
        frame.boxToLab(boxTime, boxZ, labTime, labZ);
        maximumCoordinateError = std::max(maximumCoordinateError,
          std::max(relativeError(labTime, times[i]),
                   relativeError(labZ, positions[i])));
      }

    const double velocityData[3][3] = {
      {0.2, -0.4, 1173.0},
      {-2.0, 0.5, 40.0},
      {0.01, 0.03, -0.2}
    };
    double maximumMomentumError = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
      {
        aprl::FieldVector<double> lab(0.0);
        for (unsigned int axis = 0; axis < 3; ++axis)
          lab[axis] = velocityData[i][axis];
        aprl::FieldVector<double> box(0.0);
        aprl::FieldVector<double> recovered(0.0);
        frame.properVelocityLabToBox(lab, box);
        frame.properVelocityBoxToLab(box, recovered);
        for (unsigned int axis = 0; axis < 3; ++axis)
          maximumMomentumError = std::max(maximumMomentumError,
            relativeError(recovered[axis], lab[axis]));
      }

    aprl::FieldVector<double> electricLab(0.0);
    aprl::FieldVector<double> magneticLab(0.0);
    electricLab[0] = 1.2e7;
    electricLab[1] = -3.4e6;
    electricLab[2] = 7.0e5;
    magneticLab[0] = 1.0e-3;
    magneticLab[1] = 2.0e-3;
    magneticLab[2] = -4.0e-4;
    aprl::FieldVector<double> electricBox(0.0);
    aprl::FieldVector<double> magneticBox(0.0);
    aprl::FieldVector<double> electricRecovered(0.0);
    aprl::FieldVector<double> magneticRecovered(0.0);
    frame.fieldsLabToBox(electricLab, magneticLab,
                         electricBox, magneticBox);
    frame.fieldsBoxToLab(electricBox, magneticBox,
                         electricRecovered, magneticRecovered);
    double maximumFieldError = 0.0;
    for (unsigned int axis = 0; axis < 3; ++axis)
      {
        maximumFieldError = std::max(maximumFieldError,
          relativeError(electricRecovered[axis], electricLab[axis]));
        maximumFieldError = std::max(maximumFieldError,
          relativeError(magneticRecovered[axis], magneticLab[axis]));
      }

    std::vector<aprl::RelativisticParticleSI> offsetParticles(1);
    aprl::RelativisticParticleSI& offset = offsetParticles[0];
    offset.position[0] = 1.0e-4 * (rank + 1);
    offset.position[1] = -2.0e-4 * (rank + 1);
    offset.position[2] = 2.0e-4 * rank;
    offset.properVelocity[0] = 0.2 + 0.01 * rank;
    offset.properVelocity[1] = -0.1;
    offset.properVelocity[2] = 10.0 + rank;
    offset.weight = 1.0;
    const double originalX = offset.position[0];
    const double distance = 0.02 + offset.position[2];
    const double expectedX = originalX +
      offset.properVelocity[0] / offset.properVelocity[2] * distance;
    const aprl::SILabPlaneProjectionReport offsetReport =
      aprl::SIBunchPreprocessor::projectLabPlaneToSnapshot(
        offsetParticles, 0.0, 0.02, communicator);
    require(offsetReport.particles ==
      static_cast<unsigned long long>(size),
      "Elegant v3 projection lost MPI-distributed records");
    require(relativeError(offsetParticles[0].position[0], expectedX) < 1e-14,
      "Elegant v3 transverse projection is inaccurate");

    std::vector<aprl::RelativisticParticleSI> timedParticles(2);
    for (std::size_t i = 0; i < timedParticles.size(); ++i)
      {
        aprl::RelativisticParticleSI& particle = timedParticles[i];
        particle.position[0] = 2.0e-5 * (rank + 1) * (i + 1);
        particle.position[1] = -1.0e-5 * (i + 1);
        particle.position[2] = (rank * 2.0 + i - 1.0) * 1.0e-13;
        particle.properVelocity[0] = 0.03 * (i + 1);
        particle.properVelocity[1] = -0.02 * (rank + 1);
        particle.properVelocity[2] = 20.0 + rank + i;
        particle.weight = 1.0 + rank + i;
      }
    const double longitudinalOffset = 3.0e-4;
    const aprl::SILabPlaneProjectionReport timedReport =
      aprl::SIBunchPreprocessor::projectLabPlaneEventsToSnapshot(
        timedParticles, 0.0, 0.025, longitudinalOffset, communicator);
    long double localWeight = 0.0L;
    long double localWeightedZ = 0.0L;
    for (std::size_t i = 0; i < timedParticles.size(); ++i)
      {
        localWeight += timedParticles[i].weight;
        localWeightedZ += timedParticles[i].weight *
          timedParticles[i].position[2];
      }
    long double globalWeight = 0.0L;
    long double globalWeightedZ = 0.0L;
    MPI_Allreduce(&localWeight, &globalWeight, 1, MPI_LONG_DOUBLE,
                  MPI_SUM, communicator);
    MPI_Allreduce(&localWeightedZ, &globalWeightedZ, 1, MPI_LONG_DOUBLE,
                  MPI_SUM, communicator);
    const double centroid = static_cast<double>(
      globalWeightedZ / globalWeight);
    require(timedReport.particles ==
      static_cast<unsigned long long>(2 * size),
      "Elegant v4 projection lost MPI-distributed records");
    require(std::abs(centroid - longitudinalOffset) < 2.0e-15,
      "Elegant v4 weighted longitudinal centroid is inaccurate");

    aprl::BoostFrameTransform bunchFrame;
    bunchFrame.setOriginsFromGamma(8.0, aprl::SI::c, 0.0, 0.025, 0.0);
    const aprl::SIBunchBoostReport boostReport =
      aprl::SIBunchPreprocessor::boostLabSnapshotToBoxTimeZero(
        timedParticles, bunchFrame);
    require(boostReport.maximumRelativeMomentumRoundTripError < 2.0e-13,
      "Bunch Lorentz round trip exceeded tolerance");
    require(boostReport.maximumRelativeGammaRoundTripError < 2.0e-13,
      "Bunch gamma round trip exceeded tolerance");
    if (rank == 0)
      std::cout << std::setprecision(12)
                << "lorentz_coordinate_relative_error="
                << maximumCoordinateError
                << " lorentz_momentum_relative_error="
                << maximumMomentumError
                << " lorentz_field_relative_error=" << maximumFieldError
                << " elegant_v4_weighted_z_error_m="
                << std::abs(centroid - longitudinalOffset) << std::endl;
    require(maximumCoordinateError < 2.0e-13,
      "Lorentz coordinate round trip exceeded tolerance");
    require(maximumMomentumError < 2.0e-13,
      "Lorentz momentum round trip exceeded tolerance");
    require(maximumFieldError < 2.0e-12,
      "Lorentz field round trip exceeded tolerance");

  }

  struct BorisResult
  {
    double phaseError;
    double gammaError;
  };

  BorisResult borisOrbit(unsigned int steps)
  {
    aprl::FieldVector<double> properVelocity(0.0);
    properVelocity[0] = 2.0;
    const double initialGamma =
      aprl::BoostFrameTransform::gammaFromProperVelocity(properVelocity);
    aprl::FieldVector<double> electric(0.0);
    aprl::FieldVector<double> magnetic(0.0);
    magnetic[2] = 1.0;
    const double omega = aprl::SI::elementaryCharge /
      (aprl::SI::electronMass * initialGamma);
    const double duration = 2.0 * pi / omega;
    const double dt = duration / static_cast<double>(steps);
    for (unsigned int step = 0; step < steps; ++step)
      aprl::RelativisticBorisPusher::pushMomentum(
        properVelocity, electric, magnetic,
        -aprl::SI::elementaryCharge, aprl::SI::electronMass, dt);
    BorisResult result;
    result.phaseError = std::abs(std::atan2(
      properVelocity[1], properVelocity[0]));
    result.gammaError = std::abs(
      aprl::BoostFrameTransform::gammaFromProperVelocity(properVelocity) -
      initialGamma) / initialGamma;
    return result;
  }

  void checkBorisPhase(int rank)
  {
    const BorisResult coarse = borisOrbit(32);
    const BorisResult fine = borisOrbit(64);
    require(fine.phaseError < 6.0e-3,
      "Fine Boris gyro phase error is too large");
    require(coarse.phaseError / fine.phaseError > 3.7,
      "Boris phase error did not show second-order convergence");
    require(coarse.gammaError < 2.0e-14 && fine.gammaError < 2.0e-14,
      "Magnetic Boris push did not conserve gamma");
    if (rank == 0)
      std::cout << std::setprecision(12)
                << "boris_phase_error_n32=" << coarse.phaseError
                << " boris_phase_error_n64=" << fine.phaseError
                << " convergence_ratio="
                << coarse.phaseError / fine.phaseError
                << " gamma_relative_error=" << fine.gammaError << std::endl;
  }

  void checkChargeContinuity(int rank)
  {
    const aprl::EBGridGeometry geometry(
      9, 8, 10, 0.2, 0.25, 0.15, 0.1 / aprl::SI::c,
      aprl::EBMaxwellSolver::Yee);
    aprl::EBFieldGrid fields(geometry);
    fields.clearCurrent();
    aprl::YeeComponent before(geometry.nx + 1,
                             geometry.ny + 1,
                             geometry.nz + 1);
    aprl::YeeComponent after(geometry.nx + 1,
                            geometry.ny + 1,
                            geometry.nz + 1);
    before.fill(0.0);
    after.fill(0.0);
    const aprl::FieldVector<double> origin(0.0);
    aprl::ChargeConservingCurrentDepositor depositor(fields, origin);
    aprl::FieldVector<double> start(0.0);
    aprl::FieldVector<double> end(0.0);
    start[0] = 0.31; start[1] = 0.42; start[2] = 0.28;
    end[0] = 1.54; end[1] = 1.63; end[2] = 1.21;
    const double charge = 2.3e-9;
    depositor.depositCharge(start, charge, before);
    const aprl::CurrentDepositResult deposit =
      depositor.depositSegment(start, end, charge);
    depositor.depositCharge(end, charge, after);
    const double residual =
      aprl::ChargeConservingCurrentDepositor::maxContinuityResidual(
        fields, before, after);
    double scale = 0.0;
    for (std::size_t k = 0; k < before.nz(); ++k)
      for (std::size_t j = 0; j < before.ny(); ++j)
        for (std::size_t i = 0; i < before.nx(); ++i)
          scale = std::max(scale,
            std::abs(after(i, j, k) - before(i, j, k)) / geometry.dt);
    require(deposit.subsegments >= 8,
      "Continuity trajectory did not exercise multi-cell splitting");
    require(scale > 0.0 && residual / scale < 2.0e-13,
      "Discrete charge continuity residual exceeded round-off tolerance");
    if (rank == 0)
      std::cout << std::setprecision(12)
                << "continuity_subsegments=" << deposit.subsegments
                << " absolute_residual=" << residual
                << " relative_residual=" << residual / scale << std::endl;
  }

  class NoBoundary : public aprl::EBBoundaryOperator
  {
  public:
    virtual void afterMagneticUpdate(
        aprl::EBFieldGrid&, const aprl::EBElectricHaloView&) {}
    virtual void afterElectricUpdate(aprl::EBFieldGrid&) {}
    virtual std::size_t memoryBytes() const { return 0; }
  };

  void checkCowanDispersion(int rank)
  {
    const std::size_t nz = 256;
    const std::size_t ny = 128;
    const double dz = 1.0e-5;
    const aprl::EBGridGeometry geometry(
      4, ny, nz, 4.0 * dz, 5.0 * dz, dz, dz / aprl::SI::c,
      aprl::EBMaxwellSolver::CowanZ);
    aprl::EBFieldGrid fields(geometry);
    fields.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(
      new NoBoundary()));
    const double theta = 2.0 * pi / 16.0;
    for (std::size_t k = 0; k <= nz; ++k)
      for (std::size_t j = 0; j <= geometry.ny; ++j)
        for (std::size_t i = 0; i < geometry.nx; ++i)
          fields.ex()(i, j, k) = std::sin(theta * k);
    for (std::size_t k = 0; k < nz; ++k)
      for (std::size_t j = 0; j <= geometry.ny; ++j)
        for (std::size_t i = 0; i < geometry.nx; ++i)
          fields.by()(i, j, k) =
            std::sin(theta * (k + 1.0)) / aprl::SI::c;

    const unsigned int steps = 48;
    for (unsigned int step = 0; step < steps; ++step)
      fields.advance();
    double maximumError = 0.0;
    for (std::size_t k = steps + 4; k + steps + 4 < nz; ++k)
      maximumError = std::max(maximumError,
        std::abs(fields.ex()(1, ny / 2, k) -
                 std::sin(theta * (static_cast<double>(k) - steps))));

    const aprl::EBCowanCoefficients& coefficient =
      fields.cowanCoefficients();
    const double zNormalization = coefficient.alpha[2] +
      2.0 * coefficient.beta[0] + 2.0 * coefficient.beta[1] +
      4.0 * coefficient.deltaXY;
    const double phase = aprl::EBFieldGrid::axisPhaseVelocityRatio(1.0, 16.0);
    const double group = aprl::EBFieldGrid::axisGroupVelocityRatio(1.0, 16.0);
    if (rank == 0)
      std::cout << std::setprecision(12)
                << "cowan_axial_max_error=" << maximumError
                << " phase_velocity_ratio=" << phase
                << " group_velocity_ratio=" << group
                << " smoothing_normalization=" << zNormalization
                << std::endl;
    require(maximumError < 2.0e-11,
      "Cowan-z axial wave did not translate without phase error");
    require(std::abs(zNormalization - 1.0) < 2.0e-15,
      "Cowan-z smoothing does not preserve an axial plane wave");
    require(std::abs(phase - 1.0) < 2.0e-15 &&
            std::abs(group - 1.0) < 2.0e-15,
      "Cowan-z axial analytical dispersion is not exact");
  }

  void fillFilterMode(aprl::YeeComponent& current,
                      std::size_t localCells, int rank,
                      bool nyquist, double theta)
  {
    for (std::size_t k = 0; k < current.nz(); ++k)
      {
        const std::size_t globalK =
          static_cast<std::size_t>(rank) * localCells + k;
        const double value = nyquist ?
          (globalK % 2 == 0 ? 1.0 : -1.0) :
          std::cos(theta * static_cast<double>(globalK));
        /* Tangential-current planes are shared.  Give the lower-rank copy
         * the complete synthetic value and its neighbour zero so the normal
         * MPI summation reconstructs exactly one global mode. */
        const double localValue = (rank > 0 && k == 0) ? 0.0 : value;
        for (std::size_t j = 0; j < current.ny(); ++j)
          for (std::size_t i = 0; i < current.nx(); ++i)
            current(i, j, k) = localValue;
      }
  }

  void fillTensorFilterMode(aprl::YeeComponent& current,
                            std::size_t localCells, int rank,
                            bool zNodeCentered, double theta)
  {
    for (std::size_t k = 0; k < current.nz(); ++k)
      {
        const std::size_t globalK =
          static_cast<std::size_t>(rank) * localCells + k;
        const bool duplicateLowerNode =
          zNodeCentered && rank > 0 && k == 0;
        for (std::size_t j = 0; j < current.ny(); ++j)
          for (std::size_t i = 0; i < current.nx(); ++i)
            current(i, j, k) = duplicateLowerNode ? 0.0 :
              std::cos(theta * static_cast<double>(i)) *
              std::cos(theta * static_cast<double>(j)) *
              std::cos(theta * static_cast<double>(globalK));
      }
  }

  void checkCurrentFilter(MPI_Comm communicator, int rank, int size)
  {
    require(size == 2,
      "current-filter verification requires exactly two MPI ranks");
    const std::size_t localCells = 8;
    const std::size_t globalCells = localCells *
      static_cast<std::size_t>(size);
    const aprl::EBGridGeometry geometry(
      8, 8, localCells, 4.0e-5, 5.0e-5, 1.0e-5,
      1.0e-5 / aprl::SI::c, aprl::EBMaxwellSolver::CowanZ);
    aprl::EBFieldGrid fields(geometry);
    aprl::EBZSlabHaloExchange halo(communicator);
    const std::size_t i = 4;
    const std::size_t j = 4;

    fillFilterMode(fields.jx(), localCells, rank, true, 0.0);
    halo.synchronizeAndFilterCurrent(
      fields, aprl::EBCurrentFilter::CompensatedBilinear);
    double localNyquistError = 0.0;
    for (std::size_t k = 0; k < fields.jx().nz(); ++k)
      {
        const std::size_t globalK =
          static_cast<std::size_t>(rank) * localCells + k;
        if (globalK >= 2 && globalK + 2 <= globalCells)
          localNyquistError = std::max(localNyquistError,
            std::abs(fields.jx()(i, j, k)));
      }
    double nyquistError = 0.0;
    MPI_Allreduce(&localNyquistError, &nyquistError, 1, MPI_DOUBLE,
                  MPI_MAX, communicator);

    fields.clearCurrent();
    const double theta = 2.0 * pi / 16.0;
    const double response = 1.0 -
      std::pow(std::sin(0.5 * theta), 4.0);
    fillTensorFilterMode(fields.jx(), localCells, rank, true, theta);
    fillTensorFilterMode(fields.jy(), localCells, rank, true, theta);
    fillTensorFilterMode(fields.jz(), localCells, rank, false, theta);
    halo.synchronizeAndFilterCurrent(
      fields, aprl::EBCurrentFilter::CompensatedBilinear);
    double localPassbandError = 0.0;
    const aprl::YeeComponent* component[3] = {
      &fields.jx(), &fields.jy(), &fields.jz()
    };
    const bool zNodeCentered[3] = {true, true, false};
    const std::size_t passbandI = 3;
    const std::size_t passbandJ = 3;
    for (unsigned int field = 0; field < 3; ++field)
      for (std::size_t k = 0; k < component[field]->nz(); ++k)
        {
          const std::size_t globalK =
            static_cast<std::size_t>(rank) * localCells + k;
          const bool hasStencil = globalK >= 2 &&
            (zNodeCentered[field] ? globalK + 2 <= globalCells :
                                    globalK + 2 < globalCells);
          if (hasStencil)
            {
              const double expected = response * response * response *
                std::cos(theta * static_cast<double>(passbandI)) *
                std::cos(theta * static_cast<double>(passbandJ)) *
                std::cos(theta * static_cast<double>(globalK));
              localPassbandError = std::max(localPassbandError,
                std::abs((*component[field])(
                  passbandI, passbandJ, k) - expected));
            }
        }
    double passbandError = 0.0;
    MPI_Allreduce(&localPassbandError, &passbandError, 1, MPI_DOUBLE,
                  MPI_MAX, communicator);

    require(nyquistError < 2.0e-15,
      "Compensated current filter did not cancel the z Nyquist mode");
    require(passbandError < 2.0e-15,
      "Compensated current-filter response or MPI halo is incorrect");
    if (rank == 0)
      std::cout << std::setprecision(12)
                << "current_filter_nyquist_error=" << nyquistError
                << " passband_response=" << response
                << " passband_error=" << passbandError << std::endl;
  }

  struct ReflectionResult
  {
    double incident;
    double reflected;
  };

  ReflectionResult runReflectionCase(bool cpml)
  {
    const std::size_t nz = 200;
    const double dz = 1.0e-5;
    const double dx = 100.0 * dz;
    const double dy = 100.0 * dz;
    const double dt = 0.98 /
      (aprl::SI::c * std::sqrt(1.0 / (dx * dx) +
                              1.0 / (dy * dy) + 1.0 / (dz * dz)));
    const aprl::EBGridGeometry geometry(
      2, 2, nz, dx, dy, dz, dt, aprl::EBMaxwellSolver::Yee);
    aprl::EBFieldGrid fields(geometry);
    if (cpml)
      {
        aprl::EBCPMLParameters parameter;
        parameter.cells[0] = 0;
        parameter.cells[1] = 0;
        parameter.cells[2] = 20;
        parameter.polynomialOrder = 3.0;
        parameter.targetReflection = 1.0e-8;
        parameter.kappaMax = 8.0;
        parameter.alphaFraction = 0.0;
        fields.setBoundary(std::unique_ptr<aprl::EBBoundaryOperator>(
          new aprl::EBConvolutionalPML(
            geometry, nz, 0, parameter, true, true)));
      }

    const double centre = 45.0;
    const double sigma = 10.0;
    const double waveNumber = 2.0 * pi / 20.0;
    const double halfAdvance = 0.5 * aprl::SI::c * dt / dz;
    const auto pulse = [&](double coordinate) {
      const double offset = coordinate - centre;
      return std::exp(-0.5 * offset * offset / (sigma * sigma)) *
             std::cos(waveNumber * offset);
    };
    for (std::size_t k = 0; k <= nz; ++k)
      for (std::size_t j = 0; j <= geometry.ny; ++j)
        for (std::size_t i = 0; i < geometry.nx; ++i)
          fields.ex()(i, j, k) = pulse(static_cast<double>(k));
    for (std::size_t k = 0; k < nz; ++k)
      for (std::size_t j = 0; j <= geometry.ny; ++j)
        for (std::size_t i = 0; i < geometry.nx; ++i)
          fields.by()(i, j, k) = pulse(
            static_cast<double>(k) + 0.5 + halfAdvance) / aprl::SI::c;

    ReflectionResult result = {0.0, 0.0};
    const std::size_t probe = 75;
    for (unsigned int step = 0; step < 360; ++step)
      {
        fields.advance();
        const double value = std::abs(fields.ex()(0, 1, probe));
        if (step < 120) result.incident = std::max(result.incident, value);
        if (step > 190) result.reflected = std::max(result.reflected, value);
      }
    return result;
  }

  void checkCpmlReflection(int rank)
  {
    const ReflectionResult pec = runReflectionCase(false);
    const ReflectionResult cpml = runReflectionCase(true);
    const double pecRatio = pec.reflected / pec.incident;
    const double cpmlRatio = cpml.reflected / cpml.incident;
    require(pec.incident > 0.5 && cpml.incident > 0.5,
      "CPML reflection test did not launch a resolved incident pulse");
    require(pecRatio > 0.5,
      "PEC control did not return a measurable reflection");
    require(cpmlRatio < 0.05,
      "CPML reflected more than five percent of the incident amplitude");
    require(cpml.reflected < 0.1 * pec.reflected,
      "CPML did not suppress reflection by at least one order of magnitude");
    if (rank == 0)
      std::cout << std::setprecision(12)
                << "pec_reflection_ratio=" << pecRatio
                << " cpml_reflection_ratio=" << cpmlRatio
                << " cpml_to_pec_reflected="
                << cpml.reflected / pec.reflected << std::endl;
  }

  struct FileParticleRecord
  {
    double position[3];
    double properVelocity[3];
    std::uint64_t sourceId;
    double macroWeight;
  };

  void h5Require(hid_t handle, const std::string& message)
  {
    if (handle < 0) throw std::runtime_error(message);
  }

  void h5RequireStatus(herr_t status, const std::string& message)
  {
    if (status < 0) throw std::runtime_error(message);
  }

  hid_t particleType(int version, bool fileType)
  {
    const std::size_t f64 = 8;
    const std::size_t bytes = version == 1 ? 7 * f64 : 8 * f64;
    hid_t type = H5Tcreate(H5T_COMPOUND,
      fileType ? bytes : sizeof(FileParticleRecord));
    h5Require(type, "Cannot create particle test datatype");
    const std::size_t positionOffset = fileType ? 0 :
      HOFFSET(FileParticleRecord, position);
    const std::size_t velocityOffset = fileType ? 3 * f64 :
      HOFFSET(FileParticleRecord, properVelocity);
    const std::size_t idOffset = fileType ? 6 * f64 :
      HOFFSET(FileParticleRecord, sourceId);
    const std::size_t weightOffset = fileType ? 7 * f64 :
      HOFFSET(FileParticleRecord, macroWeight);
    hsize_t three[1] = {3};
    hsize_t two[1] = {2};
    hid_t position = H5Tarray_create2(
      fileType ? H5T_IEEE_F64LE : H5T_NATIVE_DOUBLE,
      1, version == 4 ? two : three);
    hid_t velocity = H5Tarray_create2(
      fileType ? H5T_IEEE_F64LE : H5T_NATIVE_DOUBLE, 1, three);
    h5Require(position, "Cannot create particle position type");
    h5Require(velocity, "Cannot create particle velocity type");
    h5RequireStatus(H5Tinsert(type,
      version == 4 ? "plane_position_m" : "position_m",
      positionOffset, position), "Cannot add particle position");
    if (version == 4)
      h5RequireStatus(H5Tinsert(type, "arrival_time_offset_s",
        positionOffset + 2 * sizeof(double),
        fileType ? H5T_IEEE_F64LE : H5T_NATIVE_DOUBLE),
        "Cannot add particle arrival time");
    h5RequireStatus(H5Tinsert(type, "proper_velocity", velocityOffset,
      velocity), "Cannot add particle velocity");
    h5RequireStatus(H5Tinsert(type, "source_id", idOffset,
      fileType ? H5T_STD_U64LE : H5T_NATIVE_UINT64),
      "Cannot add particle source id");
    if (version >= 2)
      h5RequireStatus(H5Tinsert(type, "macro_weight", weightOffset,
        fileType ? H5T_IEEE_F64LE : H5T_NATIVE_DOUBLE),
        "Cannot add particle weight");
    H5Tclose(position);
    H5Tclose(velocity);
    return type;
  }

  void writeVersionFile(const std::string& filename, int version)
  {
    std::vector<FileParticleRecord> records(4);
    for (std::size_t i = 0; i < records.size(); ++i)
      {
        records[i].position[0] = 1.0e-4 * (i + 1);
        records[i].position[1] = -2.0e-4 * (i + 1);
        records[i].position[2] = version == 4 ?
          (static_cast<double>(i) - 1.5) * 1.0e-13 :
          3.0e-4 * (i + 1);
        records[i].properVelocity[0] = 0.01 * i;
        records[i].properVelocity[1] = -0.02 * i;
        records[i].properVelocity[2] = 10.0 + i;
        records[i].sourceId = i + 1;
        records[i].macroWeight = i + 1;
      }
    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,
                           H5P_DEFAULT, H5P_DEFAULT);
    h5Require(file, "Cannot create particle test file");
    hid_t group = H5Gcreate2(file, "/particles", H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    h5Require(group, "Cannot create particle test group");
    hid_t scalar = H5Screate(H5S_SCALAR);
    hid_t attribute = H5Acreate2(group, "format_version", H5T_STD_I32LE,
                                 scalar, H5P_DEFAULT, H5P_DEFAULT);
    h5Require(attribute, "Cannot create particle version attribute");
    h5RequireStatus(H5Awrite(attribute, H5T_NATIVE_INT, &version),
      "Cannot write particle version attribute");
    H5Aclose(attribute);
    H5Sclose(scalar);
    const hsize_t count[1] = {records.size()};
    hid_t space = H5Screate_simple(1, count, NULL);
    hid_t fileRecordType = particleType(version, true);
    hid_t memoryRecordType = particleType(version, false);
    hid_t dataset = H5Dcreate2(group, "records", fileRecordType, space,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    h5Require(dataset, "Cannot create particle test records");
    h5RequireStatus(H5Dwrite(dataset, memoryRecordType, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &records[0]), "Cannot write particle test records");
    H5Dclose(dataset);
    H5Tclose(memoryRecordType);
    H5Tclose(fileRecordType);
    H5Sclose(space);
    H5Gclose(group);
    H5Fclose(file);
  }

  void checkHdf5Versions(const std::string& directory,
                         MPI_Comm communicator, int rank)
  {
    if (rank == 0)
      for (int version = 1; version <= 4; ++version)
        {
          std::ostringstream filename;
          filename << directory << "/particles-v" << version << ".h5";
          writeVersionFile(filename.str(), version);
        }
    MPI_Barrier(communicator);

    for (int version = 1; version <= 4; ++version)
      {
        std::ostringstream filename;
        filename << directory << "/particles-v" << version << ".h5";
        aprl::FieldVector<double> offset(0.0);
        offset[0] = 0.1;
        offset[1] = 0.2;
        offset[2] = 0.3;
        unsigned long long records = 0;
        int detectedVersion = 0;
        const std::vector<aprl::RelativisticParticleSI> particles =
          aprl::ParticleHdf5File::readDistributed(
            filename.str(), 10.0, offset, communicator,
            records, detectedVersion);
        require(records == 4 && detectedVersion == version,
          "Particle HDF5 version or record count changed");
        long double localCharge = 0.0L;
        unsigned long long localRecords = particles.size();
        for (std::size_t i = 0; i < particles.size(); ++i)
          {
            const aprl::RelativisticParticleSI& particle = particles[i];
            localCharge += particle.charge;
            const double expectedWeight = version == 1 ? 1.0 :
              static_cast<double>(particle.sourceId);
            require(std::abs(particle.weight - expectedWeight) < 1.0e-15,
              "Particle HDF5 weight compatibility failed");
            require(particle.sourceId >= 1 && particle.sourceId <= 4,
              "Particle HDF5 source id changed");
            const double source = static_cast<double>(particle.sourceId);
            require(std::abs(particle.position[0] -
              (0.1 + 1.0e-4 * source)) < 1.0e-14,
              "Particle HDF5 transverse offset failed");
            const double expectedZ = version == 4 ?
              (source - 2.5) * 1.0e-13 : 0.3 + 3.0e-4 * source;
            require(std::abs(particle.position[2] - expectedZ) < 1.0e-14,
              "Particle HDF5 longitudinal compatibility failed");
          }
        long double globalCharge = 0.0L;
        unsigned long long globalLocalRecords = 0;
        MPI_Allreduce(&localCharge, &globalCharge, 1, MPI_LONG_DOUBLE,
                      MPI_SUM, communicator);
        MPI_Allreduce(&localRecords, &globalLocalRecords, 1,
                      MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);
        require(globalLocalRecords == records,
          "Distributed HDF5 hyperslabs do not cover every record once");
        const long double expectedCharge =
          -10.0L * static_cast<long double>(aprl::SI::elementaryCharge);
        require(std::abs((globalCharge - expectedCharge) / expectedCharge) <
          2.0e-15L, "Particle HDF5 charge normalization failed");
      }
    if (rank == 0)
      std::cout << "particle_hdf5_versions=1,2,3,4 mpi_distribution=passed"
                << std::endl;
  }

  void checkSaseGeneratorFile(const std::string& filename,
                              MPI_Comm communicator, int rank)
  {
    const double electrons = 1.0e8;
    const unsigned long long expectedRecords = 512;
    aprl::FieldVector<double> offset(0.0);
    unsigned long long records = 0;
    int version = 0;
    const std::vector<aprl::RelativisticParticleSI> particles =
      aprl::ParticleHdf5File::readDistributed(filename, electrons, offset,
        communicator, records, version);
    require(records == expectedRecords && version == 2,
      "SASE generator APRL input contract changed");

    unsigned long long localRecords = particles.size();
    long double localCharge = 0.0L;
    for (std::size_t i = 0; i < particles.size(); ++i)
      {
        const aprl::RelativisticParticleSI& particle = particles[i];
        require(particle.sourceId >= 1 &&
                particle.sourceId <= expectedRecords,
          "SASE generator source ID is outside the record range");
        require(std::abs(particle.weight - 1.0) < 1.0e-15,
          "SASE generator macro weight is not unity");
        require(particle.position[2] >= -4.0e-6 &&
                particle.position[2] < 4.0e-6,
          "SASE generator longitudinal support changed");
        require(particle.properVelocity[2] > 0.0 &&
                std::isfinite(particle.properVelocity[2]),
          "SASE generator produced invalid forward momentum");
        localCharge += particle.charge;
      }
    unsigned long long globalRecords = 0;
    long double globalCharge = 0.0L;
    MPI_Allreduce(&localRecords, &globalRecords, 1,
      MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);
    MPI_Allreduce(&localCharge, &globalCharge, 1,
      MPI_LONG_DOUBLE, MPI_SUM, communicator);
    require(globalRecords == expectedRecords,
      "SASE generator MPI hyperslabs do not cover all records");
    const long double expectedCharge = -static_cast<long double>(electrons) *
      static_cast<long double>(aprl::SI::elementaryCharge);
    require(std::abs((globalCharge - expectedCharge) / expectedCharge) <
      2.0e-15L, "SASE generator represented charge changed");
    if (rank == 0)
      std::cout << "sase_generator_aprl_read=passed records="
                << records << " format_version=" << version << std::endl;
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int localStatus = 0;
  try
    {
      if (argc < 2)
        throw std::invalid_argument(
          "Usage: aprl_core_checks <lorentz-elegant|boris-phase|"
          "longitudinal-layout|charge-continuity|cowan-dispersion|"
          "current-filter|cpml-reflection|runtime-defaults|"
          "initial-field-transverse|cowan-initial-field|"
          "hdf5-versions|sase-generator-file> [path]");
      const std::string check(argv[1]);
      if (check == "lorentz-elegant")
        checkLorentzAndElegant(MPI_COMM_WORLD, rank, size);
      else if (check == "longitudinal-layout")
        checkLongitudinalLayout(rank);
      else if (check == "boris-phase")
        checkBorisPhase(rank);
      else if (check == "charge-continuity")
        checkChargeContinuity(rank);
      else if (check == "cowan-dispersion")
        checkCowanDispersion(rank);
      else if (check == "current-filter")
        checkCurrentFilter(MPI_COMM_WORLD, rank, size);
      else if (check == "cpml-reflection")
        checkCpmlReflection(rank);
      else if (check == "runtime-defaults")
        checkRuntimeDefaults(rank);
      else if (check == "initial-field-transverse")
        {
          if (argc != 4)
            throw std::invalid_argument(
              "initial-field-transverse requires valid and invalid YAML cards");
          checkInitialFieldTransverseLimit(
            argv[2], argv[3], MPI_COMM_WORLD, rank);
        }
      else if (check == "initial-field-cpml-particles")
        checkCpmlParticleInitialField(MPI_COMM_WORLD, rank);
      else if (check == "cowan-initial-field")
        checkCowanInitialField(MPI_COMM_WORLD, rank, size);
      else if (check == "hdf5-versions")
        {
          if (argc != 3)
            throw std::invalid_argument(
              "hdf5-versions requires a writable work directory");
          checkHdf5Versions(argv[2], MPI_COMM_WORLD, rank);
        }
      else if (check == "sase-generator-file")
        {
          if (argc != 3)
            throw std::invalid_argument(
              "sase-generator-file requires a generated HDF5 path");
          checkSaseGeneratorFile(argv[2], MPI_COMM_WORLD, rank);
        }
      else
        throw std::invalid_argument("Unknown core check: " + check);
    }
  catch (const std::exception& error)
    {
      std::cerr << "Core verification failure on rank " << rank << ": "
                << error.what() << std::endl;
      localStatus = 1;
    }
  int globalStatus = 0;
  MPI_Allreduce(&localStatus, &globalStatus, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Finalize();
  return globalStatus;
}
