// One-step support probe for the production Cowan / gather / Boris /
// charge-conserving-current / Ampere operators. This is a numerical dependency
// experiment, not an equilibrated physical beam or a superluminal-signal claim.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

#include "compute_threads.h"
#include "eb_deposition.h"
#include "eb_field.h"
#include "eb_particles.h"
#include "eb_sources.h"

namespace
{
  const std::size_t targetI = 6, targetJ = 6, targetK = 6;

  struct Sample
  {
    bool particle;
    unsigned int substeps;
    unsigned int sourceOffset;
    double fraction;
    double amplitude;
    double byAtNextNode;
    double gatheredEx;
    double gatheredBy;
    double targetJx;
    double targetEx;
    double ampereCurrentContribution;
    double fartherEx;
    double finalUx;
    double finalUz;
    double finalFractionZ;
    double continuityRelative;
    double analyticLinearSlope;
    std::size_t depositSegments;
  };

  Sample experiment(bool withParticle, double fraction,
                    unsigned int substeps, unsigned int sourceOffset,
                    double amplitude)
  {
    const double dz = 1.e-6;
    const aprl::EBGridGeometry geometry(
      12, 12, 16, 4.e-6, 4.e-6, dz, dz / aprl::SI::c,
      aprl::EBMaxwellSolver::CowanZ);
    aprl::EBFieldGrid fields(geometry);
    fields.clearFields();
    fields.clearCurrent();

    // The pulse is uniform transversely near the probe. The two inputs differ
    // only on one Ex z-plane; no initial E/B datum at or next to targetK differs.
    for (std::size_t j = 0; j < fields.ex().ny(); ++j)
      for (std::size_t i = 0; i < fields.ex().nx(); ++i)
        fields.ex()(i, j, targetK + sourceOffset) = amplitude;

    aprl::RelativisticParticleSI particle;
    particle.position[0] = (targetI + .25) * geometry.dx;
    particle.position[1] = (targetJ + .25) * geometry.dy;
    particle.position[2] = (targetK + fraction) * geometry.dz;
    particle.properVelocity[2] = .2;
    particle.charge = -aprl::SI::elementaryCharge;
    particle.mass = aprl::SI::electronMass;
    particle.id = particle.sourceId = 1;
    const aprl::FieldVector<double> origin(0.0);
    const aprl::FieldVector<double> start(particle.position);
    aprl::ChargeConservingCurrentDepositor depositor(fields, origin);
    aprl::YeeComponent rhoBefore(geometry.nx + 1, geometry.ny + 1,
                                 geometry.nz + 1);
    aprl::YeeComponent rhoAfter(geometry.nx + 1, geometry.ny + 1,
                                geometry.nz + 1);
    rhoBefore.fill(0.0);
    rhoAfter.fill(0.0);
    if (withParticle) depositor.depositCharge(start, particle.charge, rhoBefore);

    // Same ordering as Simulation::solve, with no filtering, external fields,
    // boundary crossings or MPI handoff in this interior one-step experiment.
    fields.advanceMagnetic();
    aprl::FieldVector<double> sampledE(0.0), sampledB(0.0);
    fields.sampleFieldsPosition(start, origin, sampledE, sampledB);
    Sample result = {};
    result.particle = withParticle;
    result.substeps = substeps;
    result.sourceOffset = sourceOffset;
    result.fraction = fraction;
    result.amplitude = amplitude;
    result.byAtNextNode = fields.by()(targetI, targetJ, targetK + 1);
    result.gatheredEx = sampledE[0];
    result.gatheredBy = sampledB[1];
    if (withParticle && sourceOffset == 2)
      {
        const double gamma = std::sqrt(1.0 + .2 * .2);
        const double beta = .2 / gamma;
        const double volume = geometry.dx * geometry.dy * geometry.dz;
        const double positionQuadrature =
          (substeps + 1.0) / (2.0 * substeps);
        // Independent first-order estimate: By gathered from the new k+1
        // node is -(f-1/2) amplitude/c. Integrate its transverse Lorentz
        // acceleration using this pusher's end-of-substep position rule,
        // then integrate the linear CIC y/z weights along the short segment.
        result.analyticLinearSlope =
          -particle.charge * particle.charge * geometry.dt * geometry.dt *
          beta * (fraction - .5) /
          (particle.mass * gamma * aprl::SI::epsilon0 * volume) *
          .75 * (1.0 - fraction - .5 * beta) * positionQuadrature;
      }
    if (withParticle)
      {
        const aprl::SIFieldSourceSet emptySources;
        aprl::BoostFrameTransform identityFrame;
        identityFrame.setFromGamma(1.0, aprl::SI::c, 0.0);
        aprl::RelativisticBorisPusher::pushFromGridAndPrescribedLabSubcycled(
          particle, fields, origin, emptySources, identityFrame,
          0.0, geometry.dt, substeps);
        result.depositSegments = depositor.depositSegment(
          start, particle.position, particle.charge).subsegments;
        depositor.depositCharge(particle.position, particle.charge, rhoAfter);
        double scale = 0.0;
        for (std::size_t index = 0; index < rhoBefore.size(); ++index)
          scale = std::max(scale, std::abs(
            rhoAfter.data()[index] - rhoBefore.data()[index]) / geometry.dt);
        const double residual =
          aprl::ChargeConservingCurrentDepositor::maxContinuityResidual(
            fields, rhoBefore, rhoAfter);
        result.continuityRelative = scale == 0.0 ? residual : residual / scale;
      }
    result.targetJx = fields.jx()(targetI, targetJ, targetK);
    result.ampereCurrentContribution =
      -geometry.dt / aprl::SI::epsilon0 * result.targetJx;
    result.finalUx = particle.properVelocity[0];
    result.finalUz = particle.properVelocity[2];
    result.finalFractionZ = particle.position[2] / geometry.dz - targetK;
    fields.advanceElectric();
    result.targetEx = fields.ex()(targetI, targetJ, targetK);
    result.fartherEx = fields.ex()(targetI, targetJ, targetK - 1);
    return result;
  }

  void check(bool condition, const std::string& message,
             std::vector<std::string>& failures)
  {
    if (!condition) failures.push_back(message);
  }

  std::string quote(const std::string& input)
  {
    std::string result = "\"";
    for (char ch : input)
      {
        if (ch == '\\' || ch == '"') result += '\\';
        if (ch == '\n') result += "\\n";
        else result += ch;
      }
    return result + "\"";
  }

  void emit(const Sample& r)
  {
    std::cout << "    {\"particle\":" << (r.particle ? "true" : "false")
      << ",\"substeps\":" << r.substeps
      << ",\"source_offset_cells\":" << r.sourceOffset
      << ",\"initial_z_fraction\":" << r.fraction
      << ",\"source_ex_V_per_m\":" << r.amplitude
      << ",\"by_at_k_plus_1_T\":" << r.byAtNextNode
      << ",\"gathered_ex_V_per_m\":" << r.gatheredEx
      << ",\"gathered_by_T\":" << r.gatheredBy
      << ",\"target_jx_A_per_m2\":" << r.targetJx
      << ",\"target_ex_V_per_m\":" << r.targetEx
      << ",\"ampere_current_contribution_V_per_m\":"
      << r.ampereCurrentContribution
      << ",\"ex_at_k_minus_1_V_per_m\":" << r.fartherEx
      << ",\"final_ux\":" << r.finalUx
      << ",\"final_uz\":" << r.finalUz
      << ",\"final_z_fraction\":" << r.finalFractionZ
      << ",\"continuity_relative_residual\":" << r.continuityRelative
      << ",\"analytic_linear_ex_per_source_ex\":" << r.analyticLinearSlope
      << ",\"deposition_segments\":" << r.depositSegments << '}';
  }
}

int main(int argc, char** argv)
{
  int provided = MPI_THREAD_SINGLE;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  int ranks = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  int status = 1;
  try
    {
      if (ranks != 1 || argc != 1)
        throw std::runtime_error("Use exactly one MPI rank and no CLI arguments");
      aprl::configureComputeThreads(provided, 1);
      std::vector<Sample> samples;
      std::vector<std::string> failures;
      double maximumLinearityError = 0.0;
      double maximumAnalyticError = 0.0;
      for (unsigned int substeps : {1u, 5u})
        for (double fraction : {.6, .7, .8})
          {
            const Sample baseline = experiment(true, fraction, substeps, 2, 0.0);
            samples.push_back(baseline);
            check(baseline.targetEx == 0.0 && baseline.targetJx == 0.0,
              "Zero-perturbation baseline has transverse current or target field", failures);
            double firstResponse = 0.0;
            for (double amplitude : {1.e6, 2.e6, 4.e6})
              {
                const Sample r = experiment(true, fraction, substeps, 2, amplitude);
                samples.push_back(r);
                check(std::isfinite(r.targetEx) && std::abs(r.targetEx) > 1.e-12,
                  "Two-cell PIC dependency was not measurably nonzero", failures);
                check(r.gatheredEx == 0.0 && r.gatheredBy != 0.0,
                  "Expected magnetic gather pathway was not isolated", failures);
                check(r.targetEx == r.ampereCurrentContribution,
                  "Target field is not exactly the local deposited-current contribution", failures);
                check(r.fartherEx == 0.0,
                  "Unexpected propagation to the next farther low-z node", failures);
                check(r.finalUz > 0.0 && r.finalFractionZ < 1.0 &&
                      r.finalFractionZ > fraction && r.depositSegments == 1,
                  "Particle left the intended forward-moving within-cell regime", failures);
                check(r.continuityRelative < 1.e-12,
                  "Charge continuity did not hold to roundoff", failures);
                const double response = r.targetEx / amplitude;
                maximumAnalyticError = std::max(maximumAnalyticError,
                  std::abs(response / r.analyticLinearSlope - 1.0));
                if (amplitude == 1.e6) firstResponse = response;
                else
                  maximumLinearityError = std::max(maximumLinearityError,
                    std::abs(response / firstResponse - 1.0));
                const Sample vacuum = experiment(false, fraction, substeps, 2, amplitude);
                samples.push_back(vacuum);
                check(vacuum.targetEx == 0.0 && vacuum.targetJx == 0.0,
                  "Vacuum control reached a node two cells away in one step", failures);
              }
            const Sample fartherSource = experiment(true, fraction, substeps, 3, 4.e6);
            samples.push_back(fartherSource);
            check(fartherSource.targetEx == 0.0 && fartherSource.targetJx == 0.0 &&
                  fartherSource.gatheredBy == 0.0,
              "Three-cell source negative control affected the local target", failures);
          }
      // q E dz/(m c^2) <= 8e-6; the perturbative Boris response is linear
      // here. This bound allows coordinate-subtraction roundoff in deposition,
      // not physical truncation of a weak nonzero response.
      check(maximumLinearityError < 2.e-5,
        "Amplitude-doubling response is not in the expected linear regime", failures);
      check(maximumAnalyticError < 2.e-5,
        "Measured response differs from the independent first-order slope", failures);
      std::cout << std::setprecision(17)
        << "{\n  \"format\":\"aprl-pic-dependency-probe-v1\",\n"
        << "  \"status\":" << quote(failures.empty() ? "PASS" : "FAIL") << ",\n"
        << "  \"grid_cells\":[12,12,16],\n"
        << "  \"cell_size_m\":[4e-6,4e-6,1e-6],\n"
        << "  \"dt_s\":" << 1.e-6 / aprl::SI::c << ",\n"
        << "  \"target_ijk\":[6,6,6],\n"
        << "  \"maximum_response_linearity_error\":" << maximumLinearityError << ",\n"
        << "  \"maximum_analytic_slope_relative_error\":" << maximumAnalyticError << ",\n"
        << "  \"response_relative_tolerance\":2e-5,\n"
        << "  \"limitations\":["
        << "\"One-step numerical support test, not a propagated physical signal speed.\","
        << "\"The isolated charge is not initialized with a Poisson equilibrium field; charge continuity of its step is checked.\","
        << "\"No MPI halo boundary, CPML crossing, moving window, or production error estimate is tested.\","
        << "\"A nonzero two-cell path disproves a one-cell exact-dependency assumption; it does not establish the maximum dependency radius.\"],\n"
        << "  \"failures\":[";
      for (std::size_t index = 0; index < failures.size(); ++index)
        std::cout << (index ? "," : "") << quote(failures[index]);
      std::cout << "],\n  \"samples\":[\n";
      for (std::size_t index = 0; index < samples.size(); ++index)
        {
          if (index) std::cout << ",\n";
          emit(samples[index]);
        }
      std::cout << "\n  ]\n}\n";
      status = failures.empty() ? 0 : 1;
    }
  catch (const std::exception& error)
    {
      std::cout << "{\"status\":\"FAIL\",\"exception\":"
                << quote(error.what()) << "}\n";
      status = 1;
    }
  MPI_Finalize();
  return status;
}
