#include "moving_window_wave.h"
#include "moving_window_remap.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <sys/stat.h>
#include <vector>

namespace aprl {
namespace verification {
namespace {

constexpr std::size_t NX = 32, NY = 4, FULL_NZ = 512, WINDOW_NZ = 128;
constexpr int FULL_ORIGIN = -192;
constexpr unsigned STEPS = 176, LATE_START = 100;
constexpr double DX = 2.e-6, DY = 2.e-6, DZ = 1.e-6;
constexpr double PULSE_CENTER = 70., PULSE_HALF_SUPPORT = 16., WAVELENGTH = 8.;
constexpr double DETECTOR_INITIAL = 100., WINDOW_BETA = .8;
const double OFFSETS[6][3] = {{.5,0,0},{0,.5,0},{0,0,.5},
                             {0,.5,.5},{.5,0,.5},{.5,.5,0}};
using Clock = std::chrono::steady_clock;

void require(bool ok, const std::string& message) {
  if (!ok) throw std::runtime_error(message);
}

double elapsed(const Clock::time_point& begin) {
  return std::chrono::duration<double>(Clock::now() - begin).count();
}

std::ofstream fresh(const std::string& path) {
  struct stat status;
  require(::stat(path.c_str(), &status) != 0 && errno == ENOENT,
          "Refusing to overwrite wave output: " + path);
  std::ofstream stream(path);
  require(static_cast<bool>(stream), "Cannot create wave output: " + path);
  stream << std::setprecision(17);
  return stream;
}

YeeComponent& component(EBFieldGrid& fields, unsigned c) {
  switch (c) {
    case 0: return fields.ex(); case 1: return fields.ey();
    case 2: return fields.ez(); case 3: return fields.bx();
    case 4: return fields.by(); case 5: return fields.bz();
    case 6: return fields.jx(); case 7: return fields.jy();
    case 8: return fields.jz();
    default: throw std::out_of_range("Wave component index");
  }
}

const YeeComponent& component(const EBFieldGrid& fields, unsigned c) {
  return component(const_cast<EBFieldGrid&>(fields), c);
}

EBGridGeometry geometry(std::size_t nz) {
  return EBGridGeometry(NX, NY, nz, DX, DY, DZ, DZ / SI::c,
                        EBMaxwellSolver::CowanZ);
}

EBCPMLParameters parameters(unsigned zCells) {
  EBCPMLParameters p;
  p.cells[0] = 6;
  p.cells[1] = 0;
  p.cells[2] = zCells;
  return p;
}

void attachCpml(EBFieldGrid& fields, const EBCPMLParameters& p) {
  fields.setBoundary(std::unique_ptr<EBBoundaryOperator>(
      new EBConvolutionalPML(fields.geometry(), fields.geometry().nz,
                             0, p, true, true)));
}

// A_y has units V*s/m. Coordinates and time here are in dz and dt units;
// c*dt=dz. Its compact longitudinal envelope avoids an unmodelled incident
// field at the entering boundary. This is a common manufactured initial
// pulse, not an assertion of an exact one-way Cowan eigenmode.
double potential(double xIndex, double globalZ, double timeStep) {
  if (xIndex <= 0. || xIndex >= NX) return 0.;
  const double xi = globalZ - PULSE_CENTER - timeStep;
  if (std::abs(xi) >= PULSE_HALF_SUPPORT) return 0.;
  const double pi = std::acos(-1.);
  const double x = std::sin(pi * xIndex / NX);
  const double z = std::cos(.5 * pi * xi / PULSE_HALF_SUPPORT);
  return (DZ / SI::c) * (x*x*x*x) * (z*z*z*z) *
         std::cos(2. * pi * xi / WAVELENGTH);
}

void initialize(EBFieldGrid& fields, int origin) {
  fields.clearFields();
  fields.clearCurrent();
  for (unsigned c : {1u, 3u, 5u}) {
    YeeComponent& f = component(fields, c);
    for (std::size_t k = 0; k < f.nz(); ++k)
      for (std::size_t j = 0; j < f.ny(); ++j)
        for (std::size_t i = 0; i < f.nx(); ++i) {
          // Form the world integer before conversion, giving identical
          // floating coordinates on both native lattices at t=0.
          const double z = static_cast<double>(origin + static_cast<int>(k));
          const double x = static_cast<double>(i);
          if (c == 1) {
            // Centered finite-difference c*d_z A at E's integer time.
            f(i,j,k) = SI::c * (potential(x,z+.5,0.) - potential(x,z-.5,0.)) / DZ;
          } else if (c == 3) {
            // Curl A on the native magnetic grid at t=-dt/2.
            f(i,j,k) = -(potential(x,z+1.,-.5) - potential(x,z,-.5)) / DZ;
          } else {
            f(i,j,k) = (potential(x+1.,z,-.5) - potential(x,z,-.5)) / DX;
          }
        }
  }
}

std::uint64_t verifyInitialOverlap(const EBFieldGrid& full, const EBFieldGrid& window,
                                   int fullOrigin, int windowOrigin) {
  std::uint64_t checked = 0;
  for (unsigned c = 0; c < 6; ++c) {
    const YeeComponent& a = component(full,c);
    const YeeComponent& b = component(window,c);
    for (std::size_t k = 0; k < b.nz(); ++k)
      for (std::size_t j = 0; j < b.ny(); ++j)
        for (std::size_t i = 0; i < b.nx(); ++i) {
          const double actual = b(i,j,k);
          const double expected = a(i,j,static_cast<std::size_t>(
              windowOrigin+static_cast<int>(k)-fullOrigin));
          require(std::memcmp(&actual, &expected, sizeof(double)) == 0,
                  "Wave initial overlapping native fields differ");
          ++checked;
        }
  }
  return checked;
}

struct StateAudit {
  std::array<std::uint64_t,12> peakNonzero = {{}};
  std::uint64_t checked = 0;
  bool nonzeroSeen = false;
};

void audit(const EBFieldGrid& fields, StateAudit& result) {
  for (unsigned c = 0; c < 9; ++c) {
    const YeeComponent& f = component(fields,c);
    for (std::size_t i = 0; i < f.size(); ++i) {
      require(std::isfinite(f.data()[i]), "Nonfinite vacuum-wave field/current state");
      ++result.checked;
    }
  }
  const auto& cpml = dynamic_cast<const EBConvolutionalPML&>(fields.boundary());
  require(cpml.checkpointHistoryCount() == result.peakNonzero.size(),
          "Unexpected CPML history contract");
  for (std::size_t h = 0; h < cpml.checkpointHistoryCount(); ++h) {
    std::uint64_t nonzero = 0;
    for (double value : cpml.checkpointHistory(h)) {
      require(std::isfinite(value), "Nonfinite vacuum-wave CPML history");
      nonzero += value != 0.;
      ++result.checked;
    }
    result.peakNonzero[h] = std::max(result.peakNonzero[h], nonzero);
    result.nonzeroSeen = result.nonzeroSeen || nonzero != 0;
  }
}

// Both domains share global interpolation weights. Only the integer address
// changes with their origin; this does not alter or replace any core pusher.
// E(t_n) and B(t_n-dt/2) are sampled at the same z_plane(t_n), explicitly
// retaining the original native-time diagnostic contract.
double sample(const EBFieldGrid& fields, int origin, unsigned c,
              double x, double y, double globalZ) {
  const YeeComponent& f = component(fields,c);
  const double coordinate[3] = {x-OFFSETS[c][0], y-OFFSETS[c][1],
                                globalZ-OFFSETS[c][2]};
  int base[3]; double fraction[3];
  for (unsigned d = 0; d < 3; ++d) {
    base[d] = static_cast<int>(std::floor(coordinate[d]));
    fraction[d] = coordinate[d]-base[d];
  }
  double result = 0.;
  for (int dk = 0; dk < 2; ++dk)
    for (int dj = 0; dj < 2; ++dj)
      for (int di = 0; di < 2; ++di) {
        const double weight = (di ? fraction[0] : 1.-fraction[0]) *
                              (dj ? fraction[1] : 1.-fraction[1]) *
                              (dk ? fraction[2] : 1.-fraction[2]);
        if (weight == 0.) continue;
        const int i = base[0]+di, j = base[1]+dj, k = base[2]+dk-origin;
        require(i >= 0 && j >= 0 && k >= 0 &&
                static_cast<std::size_t>(i) < f.nx() &&
                static_cast<std::size_t>(j) < f.ny() &&
                static_cast<std::size_t>(k) < f.nz(), "Wave sample left native lattice");
        result += weight*f(i,j,k);
      }
  return result;
}

struct Norms {
  long double reference[2] = {0.,0.};
  long double candidate[2] = {0.,0.};
  long double error[2] = {0.,0.};
  double referencePeak = 0.;
  double errorPeak = 0.;
  double peakReferencePlane = 0.;
  double peakCandidatePlane = 0.;
  double peakErrorPlane = 0.;

  long double totalReference() const { return reference[0]+reference[1]; }
  long double totalCandidate() const { return candidate[0]+candidate[1]; }
  long double totalError() const { return error[0]+error[1]; }

  void add(const Norms& other) {
    for (unsigned g = 0; g < 2; ++g) {
      reference[g] += other.reference[g];
      candidate[g] += other.candidate[g];
      error[g] += other.error[g];
    }
    referencePeak = std::max(referencePeak, other.referencePeak);
    errorPeak = std::max(errorPeak, other.errorPeak);
    peakReferencePlane = std::max(peakReferencePlane, other.peakReferencePlane);
    peakCandidatePlane = std::max(peakCandidatePlane, other.peakCandidatePlane);
    peakErrorPlane = std::max(peakErrorPlane, other.peakErrorPlane);
  }
};

double rootRatio(long double numerator, long double denominator) {
  require(denominator > 0., "Vacuum-wave reference norm is zero");
  const double value = static_cast<double>(std::sqrt(numerator/denominator));
  require(std::isfinite(value), "Vacuum-wave statistic overflow");
  return value;
}

Norms plane(const EBFieldGrid& full, const EBFieldGrid& window, int fullOrigin, int origin,
            unsigned step, std::ofstream& csv) {
  Norms result;
  const double z = DETECTOR_INITIAL-WINDOW_BETA*step;
  const std::size_t pixels = (NX-12)*NY;
  for (std::size_t j = 0; j < NY; ++j)
    for (std::size_t i = 6; i < NX-6; ++i)
      for (unsigned c = 0; c < 6; ++c) {
        const double scale = c < 3 ? 1. : SI::c;
        const double a = sample(full,fullOrigin,c,i+.5,j+.5,z)*scale;
        const double b = sample(window,origin,c,i+.5,j+.5,z)*scale;
        const unsigned group = c < 3 ? 0 : 1;
        const long double d = static_cast<long double>(b)-a;
        result.reference[group] += static_cast<long double>(a)*a;
        result.candidate[group] += static_cast<long double>(b)*b;
        result.error[group] += d*d;
        result.referencePeak = std::max(result.referencePeak,std::abs(a));
        result.errorPeak = std::max(result.errorPeak,static_cast<double>(std::abs(d)));
      }
  result.peakReferencePlane = std::sqrt(result.totalReference()/pixels);
  result.peakCandidatePlane = std::sqrt(result.totalCandidate()/pixels);
  result.peakErrorPlane = std::sqrt(result.totalError()/pixels);
  const double refEy = sample(full,fullOrigin,1,NX*.5,NY*.5,z);
  const double winEy = sample(window,origin,1,NX*.5,NY*.5,z);
  const double dt = DZ/SI::c;
  csv << step << ',' << step*dt << ',' << (step-.5)*dt << ',' << origin << ',' << z*DZ;
  for (unsigned g = 0; g < 2; ++g)
    csv << ',' << static_cast<double>(std::sqrt(result.reference[g]/pixels))
        << ',' << static_cast<double>(std::sqrt(result.candidate[g]/pixels))
        << ',' << static_cast<double>(std::sqrt(result.error[g]/pixels));
  csv << ',' << result.peakReferencePlane << ',' << result.peakCandidatePlane
      << ',' << result.peakErrorPlane << ',' << refEy << ',' << winEy
      << ',' << winEy-refEy << '\n';
  return result;
}

void emitArray(std::ofstream& stream, const std::array<std::uint64_t,12>& values) {
  stream << '[';
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) stream << ',';
    stream << values[i];
  }
  stream << ']';
}

} // namespace

void runMovingWindowWave(const std::string& directory, bool staticControl, unsigned zCells,
                         MovingWindowHistoryPolicy historyPolicy, bool extendedObservation) {
  const auto begin = Clock::now();
  require(zCells == 8 || zCells == 20,
          "Wave control supports only original 8 or production-longitudinal 20 CPML cells");
  require(!extendedObservation || !staticControl,
          "Extended observation is available only for a moving candidate");
  require(::mkdir(directory.c_str(),0775) == 0 || errno == EEXIST,
          "Cannot create wave directory");
  const std::size_t candidateNz = (staticControl ? 296 : WINDOW_NZ) + zCells - 8;
  const int initialOrigin = staticControl ? FULL_ORIGIN : 0;
  const unsigned steps = extendedObservation ? 384 : STEPS;
  const std::size_t fullNz = extendedObservation ? 1024 : FULL_NZ;
  const int fullOrigin = extendedObservation ? -512 : FULL_ORIGIN;
  const bool boundaryRelative = historyPolicy == MovingWindowHistoryPolicy::BoundaryRelativeZ;
  require(!staticControl || !boundaryRelative, "Static history-policy comparison is covered by zero-shift state verification");
  const char* mode = staticControl ? (zCells == 20 ? "wave-static20" : "wave-static")
                     : (extendedObservation
                        ? (boundaryRelative ? (zCells == 20 ? "wave-boundary-long20" : "wave-boundary-long")
                                            : (zCells == 20 ? "wave-long20" : "wave-long"))
                        : (boundaryRelative ? (zCells == 20 ? "wave-boundary20" : "wave-boundary")
                                            : (zCells == 20 ? "wave20" : "wave")));
  EBFieldGrid full(geometry(fullNz)), window(geometry(candidateNz));
  const EBCPMLParameters p = parameters(zCells);
  attachCpml(full,p); attachCpml(window,p);
  initialize(full,fullOrigin); initialize(window,initialOrigin);
  const std::uint64_t initialChecks = verifyInitialOverlap(full,window,fullOrigin,initialOrigin);
  std::unique_ptr<MovingWindowRemap> remap;
  if (!staticControl)
    remap.reset(new MovingWindowRemap(window.geometry(),p,std::vector<int>{-1},historyPolicy));
  auto csv = fresh(directory+"/wave-plane.csv");
  csv << "step,time_E_box_s,time_B_box_s,window_origin_global_k,detector_z_box_m,"
      << "E_reference_rms_V_m,E_candidate_rms_V_m,E_difference_rms_V_m,"
      << "cB_reference_rms_V_m,cB_candidate_rms_V_m,cB_difference_rms_V_m,"
      << "combined_reference_rms_V_m,combined_candidate_rms_V_m,combined_difference_rms_V_m,"
      << "Ey_reference_center_V_m,Ey_candidate_center_V_m,Ey_difference_center_V_m\n";
  Norms all, late;
  StateAudit fullAudit, windowAudit;
  double fullAdvance = 0., windowAdvance = 0., remapSeconds = 0.;
  std::uint64_t shifts = 0, moved = 0, incoming = 0, historyMoved = 0, historyZeroed = 0;
  std::uint64_t historyRetained = 0;
  int origin = initialOrigin;
  for (unsigned step = 0; step <= steps; ++step) {
    audit(full,fullAudit); audit(window,windowAudit);
    const Norms observed = plane(full,window,fullOrigin,origin,step,csv);
    all.add(observed);
    if (step >= LATE_START) late.add(observed);
    if (step == steps) break;
    auto start = Clock::now();
    full.advanceMagnetic(); full.advanceElectric();
    fullAdvance += elapsed(start);
    start = Clock::now();
    window.advanceMagnetic(); window.advanceElectric();
    windowAdvance += elapsed(start);
    if (!staticControl) {
      // Exactly floor(-4*n/5), without cumulative floating-origin drift.
      const int nextOrigin = -static_cast<int>((4*(step+1)+4)/5);
      const auto report = remap->shift(window,nextOrigin-origin,nextOrigin*DZ);
      if (nextOrigin != origin) ++shifts;
      moved += report.fieldMovedBytes;
      incoming += report.fieldIncomingBytes;
      historyMoved += report.historyMovedBytes;
      historyZeroed += report.historyZeroedBytes;
      historyRetained += report.historyRetainedBytes;
      remapSeconds += report.totalWallSeconds;
      origin = nextOrigin;
    }
  }
  csv.close();
  require(static_cast<bool>(csv), "Failed writing wave-plane history");
  require(fullAudit.nonzeroSeen && windowAudit.nonzeroSeen,
          "Vacuum pulse failed to exercise CPML histories");
  require(all.referencePeak > 0. && all.peakReferencePlane > 0.,
          "Vacuum pulse did not cross the common detector");
  auto metrics = fresh(directory+"/metrics.json");
  metrics << "{\n\"format\":\"aprl-moving-window-vacuum-wave-v1\",\n"
    << "\"completed\":true,\"mode\":\"" << mode << "\","
    << "\"history_policy\":\"" << (boundaryRelative ? "boundary-relative-z" : "world-position") << "\","
    << "\"physical_status\":\"MEASURED; no physical error threshold or production acceptance imposed\",\n"
    << "\"extended_observation\":" << (extendedObservation ? "true" : "false")
    << ",\"steps\":" << steps << ",\"plane_frames\":" << steps+1
    << ",\"pixels_per_plane\":" << (NX-12)*NY << ",\"nx\":" << NX << ",\"ny\":" << NY
    << ",\"full_nz\":" << fullNz << ",\"window_nz\":" << candidateNz
    << ",\"dx_m\":" << DX << ",\"dy_m\":" << DY << ",\"dz_m\":" << DZ
    << ",\"dt_s\":" << DZ/SI::c << ",\"full_origin_global_k\":" << fullOrigin
    << ",\"window_initial_origin_global_k\":" << initialOrigin
    << ",\"window_final_origin_global_k\":" << origin
    << ",\"window_beta\":" << (staticControl ? 0. : WINDOW_BETA)
    << ",\"detector_beta\":" << WINDOW_BETA
    << ",\"window_initial_cpml_inner_lower_global_k\":" << initialOrigin+static_cast<int>(zCells)
    << ",\"window_initial_cpml_inner_upper_global_k\":" << initialOrigin+static_cast<int>(candidateNz-zCells)
    << ",\"cpml_cells\":[6,0," << zCells << "],\"cpml_polynomial_order\":" << p.polynomialOrder
    << ",\"cpml_target_reflection\":" << p.targetReflection << ",\"cpml_kappa_max\":" << p.kappaMax
    << ",\"cpml_alpha_fraction\":" << p.alphaFraction
    << ",\"cpml_z_matches_production_thickness\":" << (zCells == 20 ? "true" : "false")
    << ",\"complete_production_geometry\":false"
    << ",\"pulse_center_global_k\":" << PULSE_CENTER << ",\"pulse_half_support_cells\":" << PULSE_HALF_SUPPORT
    << ",\"carrier_wavelength_cells\":" << WAVELENGTH
    << ",\"initial_overlap_values_checked\":" << initialChecks << ",\"initial_overlap_bitwise_failures\":0"
    << ",\"finite\":true,\"reference_finite_values_checked\":" << fullAudit.checked
    << ",\"candidate_finite_values_checked\":" << windowAudit.checked
    << ",\"reference_cpml_nonzero_seen\":true,\"candidate_cpml_nonzero_seen\":true,\n"
    << "\"reference_cpml_peak_nonzero_counts\":";
  emitArray(metrics,fullAudit.peakNonzero);
  metrics << ",\"candidate_cpml_peak_nonzero_counts\":";
  emitArray(metrics,windowAudit.peakNonzero);
  metrics << ",\n\"E_relative_l2_error\":" << rootRatio(all.error[0],all.reference[0])
    << ",\"B_relative_l2_error\":" << rootRatio(all.error[1],all.reference[1])
    << ",\"combined_relative_l2_error\":" << rootRatio(all.totalError(),all.totalReference())
    << ",\"reference_peak_component_V_m_equivalent\":" << all.referencePeak
    << ",\"max_component_difference_V_m_equivalent\":" << all.errorPeak
    << ",\"max_component_difference_over_reference_peak\":" << all.errorPeak/all.referencePeak
    << ",\"peak_plane_difference_over_reference_peak_plane\":" << all.peakErrorPlane/all.peakReferencePlane
    << ",\"late_start_step\":" << LATE_START
    << ",\"late_frames\":" << steps-LATE_START+1
    << ",\"late_reference_l2_over_full_reference_l2\":" << rootRatio(late.totalReference(),all.totalReference())
    << ",\"late_candidate_l2_over_full_reference_l2\":" << rootRatio(late.totalCandidate(),all.totalReference())
    << ",\"late_difference_l2_over_full_reference_l2\":" << rootRatio(late.totalError(),all.totalReference())
    << ",\"late_reference_peak_over_reference_peak\":" << late.peakReferencePlane/all.peakReferencePlane
    << ",\"late_candidate_peak_over_reference_peak\":" << late.peakCandidatePlane/all.peakReferencePlane
    << ",\"late_difference_peak_over_reference_peak\":" << late.peakErrorPlane/all.peakReferencePlane
    << ",\"shifts\":" << shifts << ",\"field_moved_bytes\":" << moved
    << ",\"incoming_zero_bytes\":" << incoming << ",\"history_moved_bytes\":" << historyMoved
    << ",\"history_zeroed_bytes\":" << historyZeroed
    << ",\"history_retained_bytes\":" << historyRetained
    << ",\"full_grid_and_cpml_bytes\":" << full.memoryFootprint().totalBytes()
    << ",\"window_grid_and_cpml_bytes\":" << window.memoryFootprint().totalBytes()
    << ",\"remap_helper_bytes\":" << (remap ? remap->memoryBytes() : 0)
    << ",\"reference_advance_wall_s\":" << fullAdvance
    << ",\"candidate_advance_wall_s\":" << windowAdvance
    << ",\"remap_wall_s\":" << remapSeconds << ",\"overall_wall_s\":" << elapsed(begin)
    << ",\n\"field_time_contract\":\"Native box E(t_n), B(t_n-dt/2), both sampled at z_plane(t_n); no simultaneous lab-field transform\","
    << "\"plane_contract\":\"z_box/dz=100-0.8*n, fixed lab plane for boost beta=0.8 gamma=5/3; x cell centers6.5..25.5 and y0.5..3.5\","
    << "\"late_contract\":\"Step" << LATE_START << ".." << steps
    << " fixed before running; late response and difference are not an isolated reflection coefficient\","
    << "\"source_contract\":\"Common compact discrete-A initial fields, zero current; "
    << (staticControl ? "stationary candidate with no remap or layer insertion"
                      : "zero incoming vacuum layers") << "; no full-reference replay\","
    << "\"limitations\":[\"Not an exact one-way analytical mode; both domains share initialization transients and x CPML\","
    << (staticControl
        ? "\"Fixed candidate origin=-192 shares the reference lower CPML and keeps its upper CPML entrance at96; different upper geometry from moving candidate, not a single-variable contrast\","
          "\"The detector samples the upper CPML region during initial steps0..5 before the main pulse arrives\","
        : "\"CPML profiles and memories move only in the candidate; their physical effects remain measured differences\",")
    << "\"20-cell runs match only the existing production longitudinal thickness; x6/y0 remains the common two-dimensional fixture\","
    << "\"Increasing z CPML thickness preserves the upper interface by extending the upper edge; the lower interface moves with thickness\","
    << "\"Single finite pulse and duration do not establish general moving-boundary stability or charged-beam accuracy\","
    << "\"Different-domain timings are workload observations, not a production implementation speedup\"]\n}\n";
  metrics.close();
  require(static_cast<bool>(metrics), "Failed writing vacuum-wave metrics");
}

} // namespace verification
} // namespace aprl
