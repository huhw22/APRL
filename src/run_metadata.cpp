#include "run_metadata.h"

#include <algorithm>
#include <chrono>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

#include "runtime_util.h"

#ifndef APRL_SOURCE_REVISION
#define APRL_SOURCE_REVISION "unknown"
#endif

#ifndef APRL_BUILD_DESCRIPTION
#define APRL_BUILD_DESCRIPTION "unknown"
#endif

namespace aprl
{
  namespace
  {
    bool pathExists(const std::string& path)
    {
      struct stat status;
      return ::stat(path.c_str(), &status) == 0;
    }

    std::string parentDirectory(const std::string& path)
    {
      const std::string::size_type separator = path.find_last_of('/');
      if (separator == std::string::npos) return std::string();
      if (separator == 0) return "/";
      return path.substr(0, separator);
    }

    std::string quoteYaml(const std::string& value)
    {
      std::ostringstream result;
      result << '"';
      for (std::size_t i = 0; i < value.size(); ++i)
        {
          const unsigned char character =
            static_cast<unsigned char>(value[i]);
          if (character == '\\' || character == '"')
            result << '\\' << static_cast<char>(character);
          else if (character == '\n') result << "\\n";
          else if (character == '\r') result << "\\r";
          else if (character == '\t') result << "\\t";
          else if (character < 0x20)
            result << "\\x" << std::hex << std::setw(2)
                   << std::setfill('0') << static_cast<unsigned int>(character)
                   << std::dec;
          else result << static_cast<char>(character);
        }
      result << '"';
      return result.str();
    }

    std::string fileIdentity(const std::string& path)
    {
      if (path.empty()) return "generated-gaussian";
      struct stat status;
      if (::stat(path.c_str(), &status) != 0)
        return path + " (stat failed: " + std::strerror(errno) + ")";
      std::ostringstream result;
      result << path << ";size="
             << static_cast<unsigned long long>(status.st_size)
             << ";mtime=" << static_cast<long long>(status.st_mtime);
      return result.str();
    }

    std::string makeRunId()
    {
      const unsigned long long ticks =
        static_cast<unsigned long long>(
          std::chrono::system_clock::now().time_since_epoch().count());
      std::ostringstream result;
      result << std::hex << ticks << '-' << static_cast<unsigned long>(getpid());
      return result.str();
    }

    void broadcastString(std::string& value, int root, MPI_Comm communicator)
    {
      int rank = 0;
      MPI_Comm_rank(communicator, &rank);
      unsigned long long length = rank == root ?
        static_cast<unsigned long long>(value.size()) : 0;
      MPI_Bcast(&length, 1, MPI_UNSIGNED_LONG_LONG, root, communicator);
      if (length > static_cast<unsigned long long>(
            std::numeric_limits<int>::max()))
        throw std::runtime_error("Run metadata string exceeds MPI count range");
      if (rank != root) value.resize(static_cast<std::size_t>(length));
      if (length > 0)
        MPI_Bcast(&value[0], static_cast<int>(length), MPI_CHAR,
                  root, communicator);
    }
  }

  RunMetadata::RunMetadata()
    : runId(), configurationPath(), configurationDigest(), manifestPath(),
      sourceRevision(APRL_SOURCE_REVISION),
      buildDescription(APRL_BUILD_DESCRIPTION), mpiLibrary(),
      inputParticleIdentity(), mpiSize(1)
  {}

  RunMetadata initializeRunMetadata(const SimulationConfig& config,
                                    MPI_Comm communicator,
                                    const LongitudinalLayoutReport*
                                      longitudinalLayout)
  {
    RunMetadata metadata;
    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &metadata.mpiSize);
    metadata.configurationPath = config.configurationPath;
    metadata.configurationDigest = config.configurationDigest;
    metadata.manifestPath = config.output.manifest;
    metadata.inputParticleIdentity = config.beam.type == BeamInputType::Hdf5 ?
      fileIdentity(config.beam.file) : "generated-gaussian";

    char mpiVersion[MPI_MAX_LIBRARY_VERSION_STRING] = {};
    int mpiVersionLength = 0;
    MPI_Get_library_version(mpiVersion, &mpiVersionLength);
    metadata.mpiLibrary.assign(mpiVersion,
      static_cast<std::size_t>(std::max(0, mpiVersionLength)));
    while (!metadata.mpiLibrary.empty() &&
           (metadata.mpiLibrary[metadata.mpiLibrary.size() - 1] == '\0' ||
            metadata.mpiLibrary[metadata.mpiLibrary.size() - 1] == '\n' ||
            metadata.mpiLibrary[metadata.mpiLibrary.size() - 1] == '\r'))
      metadata.mpiLibrary.erase(metadata.mpiLibrary.size() - 1);

    if (rank == 0) metadata.runId = makeRunId();
    broadcastString(metadata.runId, 0, communicator);

    const std::string directory = parentDirectory(metadata.manifestPath);
    if (!directory.empty()) createDirectories(directory, communicator);

    int failed = 0;
    std::string failure;
    if (rank == 0)
      {
        try
          {
            if (!config.output.overwrite && pathExists(metadata.manifestPath))
              throw std::runtime_error(
                "Refusing to overwrite existing run manifest: " +
                metadata.manifestPath +
                "; set output.overwrite: true only for an intentional rerun");
            std::ofstream stream(metadata.manifestPath.c_str(),
              std::ios::out | std::ios::trunc | std::ios::binary);
            if (!stream)
              throw std::runtime_error("Cannot create run manifest: " +
                metadata.manifestPath);
            stream << "format: aprl-run-manifest-v1\n"
                   << "run_id: " << quoteYaml(metadata.runId) << "\n"
                   << "source_revision: "
                   << quoteYaml(metadata.sourceRevision) << "\n"
                   << "build: " << quoteYaml(metadata.buildDescription) << "\n"
                   << "mpi_size: " << metadata.mpiSize << "\n"
                   << "mpi_library: " << quoteYaml(metadata.mpiLibrary) << "\n"
                   << "configuration_path: "
                   << quoteYaml(metadata.configurationPath) << "\n"
                   << "configuration_digest_fnv1a64: "
                   << quoteYaml(metadata.configurationDigest) << "\n"
                   << "particle_input: "
                   << quoteYaml(metadata.inputParticleIdentity) << "\n"
                   << "resolved_runtime:\n"
                   << "  shutdown_watchdog:\n"
                   << "    enabled: "
                   << (config.runtime.shutdownWatchdogEnabled ?
                       "true" : "false") << "\n"
                   << "    timeout_seconds: "
                   << config.runtime.shutdownWatchdogTimeoutSeconds << "\n"
                   << "    finalize_timeout_seconds: "
                   << config.runtime.shutdownWatchdogFinalizeTimeoutSeconds
                   << "\n"
                   << std::setprecision(17)
                   << "resolved_initial_self_field:\n"
                   << "  enabled: "
                   << (config.initialSelfField.enabled ? "true" : "false")
                   << "\n"
                   << "  maximum_mean_transverse_beta: "
                   << config.initialSelfField.maximumMeanTransverseBeta
                   << "\n"
                   << "resolved_mesh:\n"
                   << "  cells: [" << config.mesh.cells[0] << ", "
                   << config.mesh.cells[1] << ", "
                   << config.mesh.cells[2] << "]\n"
                   << "  cell_size_m: [" << config.mesh.cellSize[0]
                   << ", " << config.mesh.cellSize[1] << ", "
                   << config.mesh.cellSize[2] << "]\n"
                   << "  extent_m: [" << config.mesh.extent[0] << ", "
                   << config.mesh.extent[1] << ", "
                   << config.mesh.extent[2] << "]\n"
                   << "  center_m: [" << config.mesh.center[0] << ", "
                   << config.mesh.center[1] << ", "
                   << config.mesh.center[2] << "]\n"
                   << "  ballistic_outside_transverse: "
                   << (config.mesh.ballisticOutsideTransverse ?
                       "true" : "false") << "\n";
            if (longitudinalLayout)
              {
                const LongitudinalLayout& layout =
                  longitudinalLayout->layout;
                stream
                  << "  longitudinal:\n"
                  << "    mode: observation-bounded\n"
                  << "    initial_self_field_enabled: "
                  << (config.initialSelfField.enabled ? "true" : "false")
                  << "\n"
                  << "    self_field_padding_rms: "
                  << config.mesh.longitudinal.selfFieldPaddingRms << "\n"
                  << "    target_lab_z_m: "
                  << longitudinalLayout->observationTargetLab << "\n"
                  << "    initial_head_lab_z_m: "
                  << longitudinalLayout->initialHeadLab << "\n"
                  << "    ballistic_stop_time_s: "
                  << longitudinalLayout->ballisticStopTime << "\n"
                  << "    planned_ballistic_time_s: "
                  << longitudinalLayout->plannedBallisticTime << "\n"
                  << "    radiation_advance_m: "
                  << longitudinalLayout->radiationAdvance << "\n"
                  << "    origin_m: " << layout.origin << "\n"
                  << "    outer_bounds_m: [" << layout.outerLower
                  << ", " << layout.outerUpper << "]\n"
                  << "    physical_bounds_m: [" << layout.physicalLower
                  << ", " << layout.physicalUpper << "]\n"
                  << "    translation_m: " << layout.translation << "\n"
                  << "    initial_bounds_m: [" << layout.initialMinimum
                  << ", " << layout.initialMaximum << "]\n"
                  << "    ballistic_bounds_m: ["
                  << layout.ballisticMinimum << ", "
                  << layout.ballisticMaximum << "]\n"
                  << "    radiation_upper_m: " << layout.radiationUpper
                  << "\n"
                  << "    self_field_bounds_m: [" << layout.selfFieldLower
                  << ", " << layout.selfFieldUpper << "]\n"
                  << "    required_bounds_m: [" << layout.requiredLower
                  << ", " << layout.requiredUpper << "]\n"
                  << "    minimum_physical_cells: "
                  << layout.minimumPhysicalCells << "\n"
                  << "    alignment_extra_cells: "
                  << layout.alignmentExtraCells << "\n"
                  << "    right_rounding_padding_m: "
                  << layout.rightRoundingPadding << "\n";
              }
            stream << "configuration_yaml: |\n";
            std::istringstream card(config.configurationText);
            std::string line;
            while (std::getline(card, line)) stream << "  " << line << "\n";
            if (!stream)
              throw std::runtime_error("Cannot finish run manifest: " +
                metadata.manifestPath);
          }
        catch (const std::exception& error)
          {
            failed = 1;
            failure = error.what();
          }
      }
    MPI_Bcast(&failed, 1, MPI_INT, 0, communicator);
    broadcastString(failure, 0, communicator);
    if (failed) throw std::runtime_error(failure);
    MPI_Barrier(communicator);
    return metadata;
  }
}
