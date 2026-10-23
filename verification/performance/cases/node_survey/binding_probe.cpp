// Compile with the APRL MPI wrapper and OpenMP: mpic++ -std=c++11 -O2 -fopenmp.
// Standalone derivative of cluster_acceptance/binding_probe.cpp.
// Run only after explicit user approval of CPU diagnostics on gnode1 without Slurm.
// The original Slurm probe remains unchanged. This performs no scientific simulation
// or MPI network configuration; its binding and physical-core checks are preserved.
#include <mpi.h>
#include <omp.h>
#include <sched.h>
#include <unistd.h>

#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct ThreadBinding {
  int cpu = -1;
  int package = -1;
  int core = -1;
  bool valid = false;
  std::string affinity;
  std::string reason;
};

bool positive_integer(const char* text, int& value) {
  if (!text || !*text) return false;
  char* end = 0;
  errno = 0;
  const long parsed = std::strtol(text, &end, 10);
  if (errno || !end || *end || parsed < 1 || parsed > INT_MAX) return false;
  value = static_cast<int>(parsed);
  return true;
}

bool topology(int cpu, int& package, int& core) {
  const std::string directory = "/sys/devices/system/cpu/cpu" +
    std::to_string(cpu) + "/topology/";
  std::ifstream package_file(directory + "physical_package_id");
  std::ifstream core_file(directory + "core_id");
  return bool(package_file >> package) && bool(core_file >> core) &&
    package >= 0 && core >= 0;
}

ThreadBinding inspect_thread() {
  ThreadBinding record;
  cpu_set_t mask;
  CPU_ZERO(&mask);
  if (sched_getaffinity(0, sizeof(mask), &mask) != 0) {
    record.reason = std::string("sched_getaffinity failed: ") + std::strerror(errno);
    return record;
  }
  record.cpu = sched_getcpu();
  if (record.cpu < 0 || record.cpu >= CPU_SETSIZE || !CPU_ISSET(record.cpu, &mask)) {
    record.reason = "current CPU is unavailable or outside thread affinity";
    return record;
  }
  if (!topology(record.cpu, record.package, record.core)) {
    record.reason = "cannot read current CPU physical topology";
    return record;
  }
  bool one_core = true;
  std::ostringstream affinity;
  for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    if (!CPU_ISSET(cpu, &mask)) continue;
    if (affinity.tellp() > 0) affinity << ',';
    affinity << cpu;
    int package = -1;
    int core = -1;
    if (!topology(cpu, package, core) || package != record.package || core != record.core)
      one_core = false;
  }
  record.affinity = affinity.str();
  record.valid = one_core && !record.affinity.empty();
  record.reason = record.valid ? "bound to one physical core" :
    "thread affinity spans multiple or unreadable physical cores";
  return record;
}

std::vector<std::string> fields(const std::string& line) {
  std::vector<std::string> result;
  std::istringstream stream(line);
  std::string field;
  while (std::getline(stream, field, '|')) result.push_back(field);
  return result;
}
}  // namespace

int main(int argc, char** argv) {
  const char* authorization = std::getenv("APRL_ALLOW_UNSCHEDULED_CPU_TESTS");
  char approved_host[256] = {};
  if (!authorization || std::strcmp(authorization, "1") != 0 ||
      gethostname(approved_host, sizeof(approved_host) - 1) != 0 ||
      std::strcmp(approved_host, "gnode1") != 0) {
    std::cerr << "ERROR: standalone binding probe requires explicit approval, "
                 "APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1, and hostname exactly gnode1\n";
    return 2;
  }
  int provided = MPI_THREAD_SINGLE;
  if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) != MPI_SUCCESS)
    return 2;
  int rank = 0;
  int ranks = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  int expected_ranks = 0;
  int expected_threads = 0;
  if (argc != 3 || !positive_integer(argv[1], expected_ranks) ||
      !positive_integer(argv[2], expected_threads) || expected_threads > CPU_SETSIZE) {
    if (rank == 0) std::cerr << "Usage: binding_probe EXPECTED_RANKS EXPECTED_THREADS\n";
    MPI_Finalize();
    return 2;
  }
  int failures = 0;
  if (provided < MPI_THREAD_FUNNELED || ranks != expected_ranks) {
    ++failures;
    if (rank == 0) std::cerr << "ERROR: MPI rank count or thread support differs from request\n";
  }
  int environment_threads = 0;
  if (!positive_integer(std::getenv("OMP_NUM_THREADS"), environment_threads) ||
      environment_threads != expected_threads) {
    ++failures;
    std::cerr << "ERROR: rank=" << rank << " OMP_NUM_THREADS differs from request\n";
  }
  char host[MPI_MAX_PROCESSOR_NAME + 1] = {};
  int host_length = 0;
  MPI_Get_processor_name(host, &host_length);
  host[host_length] = '\0';
  const int allowed_processor = std::strcmp(host, "gnode1") == 0 ? 1 : 0;
  int all_processors_allowed = 0;
  MPI_Allreduce(&allowed_processor, &all_processors_allowed, 1, MPI_INT,
                MPI_MIN, MPI_COMM_WORLD);
  if (!all_processors_allowed) {
    if (rank == 0)
      std::cerr << "ERROR: every MPI processor name must be exactly gnode1\n";
    MPI_Finalize();
    return 2;
  }
  omp_set_dynamic(0);
  int actual_threads = 0;
  std::vector<ThreadBinding> records(expected_threads);
#pragma omp parallel num_threads(expected_threads)
  {
#pragma omp single
    actual_threads = omp_get_num_threads();
    records[omp_get_thread_num()] = inspect_thread();
  }
  if (actual_threads != expected_threads) {
    ++failures;
    std::cerr << "ERROR: rank=" << rank << " actual_threads=" << actual_threads
              << " expected_threads=" << expected_threads << '\n';
  }
  std::ostringstream local;
  for (int thread = 0; thread < actual_threads; ++thread) {
    const ThreadBinding& record = records[thread];
    if (!record.valid) ++failures;
    local << host << '|' << record.package << '|' << record.core << '|'
          << rank << '|' << thread << '|' << record.cpu << '|'
          << record.affinity << '|' << (record.valid ? "PASS" : "FAIL") << '|'
          << record.reason << '\n';
  }
  const std::string local_text = local.str();
  const int local_size = static_cast<int>(local_text.size());
  std::vector<int> sizes(rank == 0 ? ranks : 0);
  MPI_Gather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::vector<int> offsets(rank == 0 ? ranks : 0);
  int total_size = 0;
  if (rank == 0) {
    for (int item = 0; item < ranks; ++item) {
      offsets[item] = total_size;
      total_size += sizes[item];
    }
  }
  std::vector<char> gathered(total_size);
  MPI_Gatherv(local_text.data(), local_size, MPI_CHAR, gathered.data(), sizes.data(),
              offsets.data(), MPI_CHAR, 0, MPI_COMM_WORLD);
  if (rank == 0) {
    std::cout << "[binding] ranks=" << ranks << " expected_ranks=" << expected_ranks
              << " threads_per_rank=" << expected_threads
              << " mpi_thread_support=" << provided << '\n';
    std::map<std::string, std::string> owners;
    std::istringstream lines(std::string(gathered.begin(), gathered.end()));
    std::string line;
    while (std::getline(lines, line)) {
      const std::vector<std::string> row = fields(line);
      if (row.size() != 9) {
        ++failures;
        std::cerr << "ERROR: malformed binding record\n";
        continue;
      }
      const std::string owner = "rank=" + row[3] + " thread=" + row[4];
      std::cout << "[binding] " << owner << " host=" << row[0] << " cpu=" << row[5]
                << " package=" << row[1] << " core=" << row[2] << " affinity=" << row[6]
                << " status=" << row[7] << " reason=" << row[8] << '\n';
      if (row[7] != "PASS") continue;
      const std::string key = row[0] + "/package=" + row[1] + "/core=" + row[2];
      const auto insertion = owners.emplace(key, owner);
      if (!insertion.second) {
        ++failures;
        std::cerr << "ERROR: physical core " << key << " is shared by "
                  << insertion.first->second << " and " << owner << '\n';
      }
    }
  }
  int total_failures = 0;
  MPI_Allreduce(&failures, &total_failures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
    std::cout << "[binding] result=" << (total_failures ? "FAIL" : "PASS")
              << " failures=" << total_failures << std::endl;
  MPI_Finalize();
  return total_failures ? 1 : 0;
}
