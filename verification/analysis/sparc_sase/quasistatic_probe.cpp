// Offline only: reconstruct a rigid, common-beta reference at the actual
// detector events. No PIC update, particle push, or production file is changed.
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

#include <mpi.h>
#include <yaml-cpp/yaml.h>

#include "blocked_cowan_reference.h"
#include "compute_threads.h"

using namespace aprl;
using aprl::verification::BlockedCowanReference;

namespace {

struct Event {
  double t, x, y, ux, uy, uz, q, mass, weight;
  std::uint64_t id;
};
struct Request {
  double t_box, t_lab, z_box;
};
static_assert(sizeof(Event) == 80, "Event binary ABI must be 80 bytes");
static_assert(sizeof(Request) == 24, "Request binary ABI must be 24 bytes");

bool exists(const std::string& path) {
  struct stat info;
  if (::stat(path.c_str(), &info) == 0) return true;
  if (errno == ENOENT) return false;
  throw std::runtime_error("Cannot inspect output path: " + path);
}

std::size_t checkedProduct(std::size_t a, std::size_t b) {
  if (b && a > std::numeric_limits<std::size_t>::max() / b)
    throw std::overflow_error("Binary file extent overflows size_t");
  return a * b;
}

template<class T>
std::vector<T> readExact(const std::string& path, std::size_t count) {
  const std::size_t bytes = checkedProduct(count, sizeof(T));
  if (bytes > static_cast<std::size_t>(
        std::numeric_limits<std::streamsize>::max()))
    throw std::overflow_error("Binary file extent exceeds streamsize");
  std::ifstream input(path.c_str(), std::ios::binary | std::ios::ate);
  if (!input || input.tellg() != static_cast<std::streamoff>(bytes))
    throw std::runtime_error("Wrong binary file length: " + path);
  input.seekg(0);
  std::vector<T> values(count);
  if (bytes) input.read(reinterpret_cast<char*>(values.data()), bytes);
  if (!input) throw std::runtime_error("Short binary read: " + path);
  return values;
}

double tolerance(double a, double b = 0., double c = 0.) {
  return 256. * std::numeric_limits<double>::epsilon() *
    std::max(1., std::max(std::abs(a), std::max(std::abs(b), std::abs(c))));
}

void requireFinite(double value, const std::string& name) {
  if (!std::isfinite(value))
    throw std::runtime_error("Nonfinite " + name);
}

struct Configuration {
  YAML::Node global, local;
  explicit Configuration(const std::string& filename, const std::string& name)
    : global(YAML::LoadFile(filename)) {
    if (!global.IsMap()) throw std::runtime_error("Config must be a mapping");
    if (global["cases"]) {
      if (!global["cases"].IsMap() || !global["cases"][name])
        throw std::runtime_error("Missing requested case configuration");
      local = global["cases"][name];
      if (!local.IsMap()) throw std::runtime_error("Case must be a mapping");
    }
  }
  YAML::Node operator()(const std::string& key) const {
    if (local.IsMap() && local[key]) return local[key];
    return global[key];
  }
};

} // namespace

int main(int argc, char** argv) {
  int provided = 0;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  int status = 0;
  try {
    int ranks = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    if (ranks != 1)
      throw std::runtime_error("Offline reference driver requires one MPI rank");
    if (argc != 3)
      throw std::runtime_error("Usage: quasistatic_probe config.yaml CASE");
    const std::string name(argv[2]);
    const Configuration cfg(argv[1], name);
    const int threads = cfg("threads").as<int>(8);
    if (threads <= 0) throw std::runtime_error("Invalid thread count");
    configureComputeThreads(provided, threads);

    const std::string directory = cfg("directory") ?
      cfg("directory").as<std::string>() :
      cfg("root").as<std::string>() + "/" + name;
    const std::string refdir = directory + "/reference";
    const std::string eventFile = cfg("events_file").as<std::string>(
      directory + "/events.bin");
    const std::string requestFile = cfg("requests_file").as<std::string>(
      directory + "/requests.bin");
    for (const char* file : {"fields.bin", "curve.csv", "SOURCE.json",
                            "REFERENCE.json", "COMPLETE.json"})
      if (exists(directory + "/" + file))
        throw std::runtime_error("Refusing existing output: " + directory + "/" + file);
    if (exists(refdir))
      throw std::runtime_error("Refusing existing reference directory");

    const std::size_t nx = cfg("nx").as<std::size_t>();
    const std::size_t ny = cfg("ny").as<std::size_t>();
    const std::size_t nz = cfg("nz").as<std::size_t>();
    const double dx = cfg("dx").as<double>(), dy = cfg("dy").as<double>();
    const double dz = cfg("dz").as<double>(), dt = cfg("dt").as<double>();
    const std::size_t physicalGuard =
      cfg("physical_guard").as<std::size_t>(10);
    if (nx < 3 || ny < 3 || nz < 8 ||
        physicalGuard > (nx - 1) / 2 || physicalGuard > (ny - 1) / 2)
      throw std::runtime_error("Invalid reference/aperture dimensions");
    EBGridGeometry geometry(nx, ny, nz, dx, dy, dz, dt,
                            EBMaxwellSolver::CowanZ);
    const YAML::Node originNode = cfg("origin"), guardsNode = cfg("guards");
    if (!originNode.IsSequence() || originNode.size() != 3 ||
        !guardsNode.IsSequence() || guardsNode.size() != 3)
      throw std::runtime_error("origin and guards must each contain three entries");
    FieldVector<Double> origin(0.);
    std::size_t guards[3];
    for (unsigned a = 0; a < 3; ++a) {
      origin[a] = originNode[a].as<double>();
      requireFinite(origin[a], "grid origin");
      guards[a] = guardsNode[a].as<std::size_t>();
    }
    const double solverBeta = cfg("solver_beta").as<double>();
    const double translationBeta = cfg("translation_beta").as<double>();
    const double referenceTime = cfg("ref_time_box_s").as<double>();
    const std::string rewind = cfg("source_rewind").as<std::string>("common");
    if (!std::isfinite(solverBeta) || std::abs(solverBeta) >= 1. ||
        !std::isfinite(translationBeta) || std::abs(translationBeta) >= 1. ||
        !std::isfinite(referenceTime))
      throw std::runtime_error("Invalid common reference velocity/time");
    if (rewind != "common" && rewind != "individual-full" &&
        rewind != "individual-held")
      throw std::runtime_error("Unknown source_rewind mode");
    BoostFrameTransform frame;
    frame.setOriginsFromGamma(cfg("boost_gamma").as<double>(), SI::c,
      cfg("frame_time_origin_lab_s").as<double>(),
      cfg("frame_z_origin_lab_m").as<double>(),
      cfg("frame_box_origin_m").as<double>());

    const std::size_t particleCount = cfg("n_particles").as<std::size_t>();
    const std::size_t sampleCount = cfg("n_samples").as<std::size_t>();
    if (!particleCount || !sampleCount)
      throw std::runtime_error("Source and sample counts must be nonzero");
    std::vector<Request> requests = readExact<Request>(requestFile, sampleCount);
    const double planeZ = cfg("plane_z_lab_m").as<double>(
      frame.labZFromBoxZT(requests[0].z_box, requests[0].t_box));
    requireFinite(planeZ, "detector laboratory position");
    for (std::size_t s = 0; s < sampleCount; ++s) {
      const Request& r = requests[s];
      requireFinite(r.t_box, "request box time");
      requireFinite(r.t_lab, "request laboratory time");
      requireFinite(r.z_box, "request box position");
      if (s && (!(r.t_box > requests[s-1].t_box) ||
                !(r.t_lab > requests[s-1].t_lab)))
        throw std::runtime_error("Requests must increase in both time coordinates");
      const double expectedZ = frame.boxZFromLabZAndBoxT(planeZ, r.t_box);
      if (std::abs(r.z_box - expectedZ) >
          tolerance(r.z_box, expectedZ, frame.boxZOrigin()))
        throw std::runtime_error("Request is not on the configured fixed lab plane");
      const double expectedTime = frame.labTFromBoxZT(r.z_box, r.t_box);
      const double timeScale = std::max(std::abs(expectedTime),
        std::max(std::abs(r.t_lab), std::abs(frame.gamma() * r.t_box)));
      if (std::abs(r.t_lab - expectedTime) >
          512. * std::numeric_limits<double>::epsilon() *
          std::max(timeScale, dt))
        throw std::runtime_error("Request time coordinates disagree");
    }

    std::vector<Event> events = readExact<Event>(eventFile, particleCount);
    std::vector<RelativisticParticleSI> particles(particleCount);
    long double charge = 0., absCharge = 0., mass = 0., weights = 0.;
    long double fullBetaMass = 0., fullBeta2Mass = 0., heldBetaMass = 0.;
    double zmin = std::numeric_limits<double>::infinity();
    double zmax = -std::numeric_limits<double>::infinity();
    double rewindConsistency = 0.;
    for (std::size_t i = 0; i < particleCount; ++i) {
      const Event& e = events[i];
      const double values[] = {e.t,e.x,e.y,e.ux,e.uy,e.uz,e.q,e.mass,e.weight};
      for (double value : values) requireFinite(value, "particle event");
      if (!e.id || (i && e.id <= events[i-1].id) || e.q == 0. ||
          !(e.mass > 0.) || !(e.weight > 0.))
        throw std::runtime_error("Invalid event ID, charge, mass, or weight");
      RelativisticParticleSI& p = particles[i];
      p.id = p.sourceId = e.id;
      p.charge = e.q;
      p.mass = e.mass;
      p.weight = e.weight;
      FieldVector<Double> labU(0.);
      labU[0] = e.ux; labU[1] = e.uy; labU[2] = e.uz;
      frame.properVelocityLabToBox(labU, p.properVelocity);
      const double fullBeta = p.properVelocity[2] /
        BoostFrameTransform::gammaFromProperVelocity(p.properVelocity);
      const double heldBeta = frame.heldLongitudinalBetaZ(p.properVelocity);
      double tb = 0., zb = 0.;
      frame.labToBox(e.t, planeZ, tb, zb);
      const double sourceBeta = rewind == "common" ? solverBeta :
        (rewind == "individual-full" ? fullBeta : heldBeta);
      p.position[0] = e.x;
      p.position[1] = e.y;
      p.position[2] = std::fma(sourceBeta * SI::c, referenceTime - tb, zb);
      const double recovered = std::fma(solverBeta * SI::c,
                                       tb - referenceTime, p.position[2]);
      const double error = std::abs(recovered - zb);
      rewindConsistency = std::max(rewindConsistency, error);
      if (rewind == "common" && error >
          tolerance(recovered, zb, p.position[2]))
        throw std::runtime_error("Common rewind failed crossing-event consistency");
      zmin = std::min(zmin, static_cast<double>(p.position[2]));
      zmax = std::max(zmax, static_cast<double>(p.position[2]));
      charge += e.q; absCharge += std::abs(e.q);
      mass += e.mass; weights += e.weight;
      fullBetaMass += static_cast<long double>(e.mass) * fullBeta;
      fullBeta2Mass += static_cast<long double>(e.mass) * fullBeta * fullBeta;
      heldBetaMass += static_cast<long double>(e.mass) * heldBeta;
    }
    const double computedBeta = static_cast<double>(fullBetaMass / mass);
    const double computedHeldBeta = static_cast<double>(heldBetaMass / mass);
    const double fullBetaStd = std::sqrt(static_cast<double>(std::max(
      0.L, fullBeta2Mass / mass - (fullBetaMass / mass) * (fullBetaMass / mass))));
    {
      std::ofstream f((directory + "/SOURCE.json").c_str());
      f << std::setprecision(17)
        << "{\n\"particles\":" << particleCount
        << ",\"source_rewind\":\"" << rewind
        << "\",\"source_charge_C\":" << static_cast<double>(charge)
        << ",\"source_absolute_charge_C\":" << static_cast<double>(absCharge)
        << ",\"source_mass_kg\":" << static_cast<double>(mass)
        << ",\"weight_sum\":" << static_cast<double>(weights)
        << ",\"solver_beta\":" << solverBeta
        << ",\"production_translation_beta\":" << translationBeta
        << ",\"full_beta_mass_weighted\":" << computedBeta
        << ",\"full_beta_mass_weighted_std\":" << fullBetaStd
        << ",\"held_beta_mass_weighted\":" << computedHeldBeta
        << ",\"ref_time_box_s\":" << referenceTime
        << ",\"source_z_min_m\":" << zmin
        << ",\"source_z_max_m\":" << zmax
        << ",\"plane_z_lab_m\":" << planeZ
        << ",\"original_event_consistency_zmax_m\":" << rewindConsistency
        << ",\"charge_mass_not_reweighted\":true"
        << ",\"physical_acceptance\":false"
        << ",\"scope\":\"Rigid common-beta finite-domain reference; transverse event coordinates held fixed\"}\n";
      f.close();
      if (!f) throw std::runtime_error("Cannot write source report");
    }
    std::vector<Event>().swap(events);
    std::cerr << std::setprecision(17) << "SOURCE particles=" << particleCount
      << " charge=" << static_cast<double>(charge)
      << " beta=" << solverBeta << " z=" << zmin << "," << zmax
      << " event_error=" << rewindConsistency << std::endl;
    const std::size_t blockLayers = cfg("block_layers").as<std::size_t>(32);
    BlockedCowanReference reference(geometry, origin, particles, solverBeta,
      guards, blockLayers, refdir, MPI_COMM_SELF, false);
    const double chargeError = std::abs(reference.report().depositedCharge -
      static_cast<double>(charge)) / static_cast<double>(absCharge);
    if (!(chargeError < 1e-10))
      throw std::runtime_error("Reference deposition changed source charge");
    std::vector<RelativisticParticleSI>().swap(particles);
    std::cerr << "REFERENCE prepared" << std::endl;

    const std::size_t px = nx - 2 * physicalGuard;
    const std::size_t py = ny - 2 * physicalGuard;
    std::vector<double> cube(checkedProduct(checkedProduct(px, py), 6));
    EBGridGeometry slabGeometry = geometry;
    slabGeometry.nz = 4;
    EBFieldGrid slab(slabGeometry);
    YeeComponent* components[6] =
      {&slab.ex(), &slab.ey(), &slab.ez(),
       &slab.bx(), &slab.by(), &slab.bz()};
    std::ofstream fields((directory + "/fields.bin").c_str(), std::ios::binary);
    std::ofstream curve((directory + "/curve.csv").c_str());
    if (!fields || !curve) throw std::runtime_error("Cannot open waveform outputs");
    curve << "frame,t_box_s,t_lab_s,z_box_m,reference_shift_cells,power_W,"
             "field_energy_per_length_J_per_m\n" << std::setprecision(17);
    for (std::size_t s = 0; s < sampleCount; ++s) {
      const Request& request = requests[s];
      const double relative = (request.z_box - origin[2]) / dz;
      const long off = static_cast<long>(std::floor(relative)) - 1;
      const double shift = solverBeta * SI::c *
        (request.t_box - referenceTime) / dz;
      // The shared reference field() cache is serial. k first keeps all six
      // component queries inside the same bounded native-E/B block.
      for (std::size_t k = 0; k < 5; ++k)
        for (unsigned component = 0; component < 6; ++component) {
          YeeComponent& value = *components[component];
          if (k >= value.nz()) continue;
          const double sourceIndex = static_cast<double>(off) + k - shift;
          const long lo = static_cast<long>(std::floor(sourceIndex));
          const double fraction = sourceIndex - lo;
          const std::size_t nativeNz = nz +
            ((component == 0 || component == 1 || component == 5) ? 1 : 0);
          if (lo < 0 || lo + 1 >= static_cast<long>(nativeNz))
            throw std::runtime_error("Translated reference sample outside source grid");
          for (std::size_t j = 0; j < value.ny(); ++j)
            for (std::size_t i = 0; i < value.nx(); ++i)
              value(i,j,k) =
                (1. - fraction) * reference.field(component, i, j, lo) +
                fraction * reference.field(component, i, j, lo + 1);
        }
      FieldVector<Double> slabOrigin = origin;
      slabOrigin[2] = std::fma(static_cast<double>(off), dz, origin[2]);
      long double power = 0., energyPerLength = 0.;
      int nonfinite = 0;
      #pragma omp parallel for collapse(2) reduction(+:power,energyPerLength) reduction(|:nonfinite) num_threads(threads)
      for (std::size_t j = 0; j < py; ++j)
        for (std::size_t i = 0; i < px; ++i) {
          FieldVector<Double> position(0.), electric(0.), magnetic(0.);
          position[0] = std::fma(i + physicalGuard + .5, dx, origin[0]);
          position[1] = std::fma(j + physicalGuard + .5, dy, origin[1]);
          position[2] = request.z_box;
          slab.sampleFieldsPositionLab(position, slabOrigin, frame,
                                        electric, magnetic);
          const std::size_t pixel = (j * px + i) * 6;
          for (unsigned a = 0; a < 3; ++a) {
            cube[pixel+a] = electric[a];
            cube[pixel+3+a] = magnetic[a];
            nonfinite |= !std::isfinite(electric[a]) || !std::isfinite(magnetic[a]);
          }
          power += (electric[0]*magnetic[1] -
                    electric[1]*magnetic[0]) / SI::mu0;
          energyPerLength += .5 * SI::epsilon0 * electric.norm2() +
                             .5 / SI::mu0 * magnetic.norm2();
        }
      if (nonfinite) throw std::runtime_error("Nonfinite sampled reference fields");
      fields.write(reinterpret_cast<const char*>(cube.data()),
        static_cast<std::streamsize>(cube.size() * sizeof(double)));
      curve << s << "," << request.t_box << "," << request.t_lab << ","
        << request.z_box << "," << shift << ","
        << static_cast<double>(power * dx * dy) << ","
        << static_cast<double>(energyPerLength * dx * dy) << "\n";
      if (!fields || !curve) throw std::runtime_error("Waveform output write failed");
      if ((s+1) % 512 == 0 || s+1 == sampleCount)
        std::cerr << "SAMPLED " << s+1 << "/" << sampleCount << std::endl;
    }
    fields.close(); curve.close();
    if (!fields || !curve) throw std::runtime_error("Waveform output close failed");
    reference.writeReport(directory + "/REFERENCE.json");
    {
      std::ofstream f((directory + "/COMPLETE.json").c_str());
      f << std::setprecision(17)
        << "{\"complete\":true,\"physical_acceptance\":false"
        << ",\"n_particles\":" << particleCount
        << ",\"n_samples\":" << sampleCount
        << ",\"shape\":[" << sampleCount << "," << py << "," << px << ",6]"
        << ",\"field_components\":[\"Ex\",\"Ey\",\"Ez\",\"Bx\",\"By\",\"Bz\"]"
        << ",\"field_dtype\":\"native-endian float64\""
        << ",\"native_field_times\":\"E(0),B(-dt/2) before rigid translation\""
        << ",\"source_rewind\":\"" << rewind << "\""
        << ",\"deposition_relative_charge_error\":" << chargeError
        << ",\"original_event_consistency_zmax_m\":" << rewindConsistency
        << ",\"scope\":\"Offline quasi-static background estimate; not a unique radiation separation\"}\n";
      f.close();
      if (!f) throw std::runtime_error("Cannot write completion report");
    }
    std::cerr << "COMPLETE" << std::endl;
  } catch (const std::exception& error) {
    std::cerr << "ERROR " << error.what() << std::endl;
    status = 1;
  }
  MPI_Finalize();
  return status;
}
