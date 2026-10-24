// Bounded array-throughput / first-touch diagnostic, not an APRL physics test.
// Reported bytes are two array reads plus one write; they are NOT measured DRAM
// traffic or hardware-counter bandwidth. Compile without fast-math.
#include <omp.h>
#include <sched.h>
#include <unistd.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Clock = std::chrono::steady_clock;

static std::map<int, unsigned long long> numaPages()
{
  std::ifstream file("/proc/self/numa_maps");
  std::map<int, unsigned long long> pages;
  std::string line, word;
  while (std::getline(file, line))
  {
    // Whole-process anonymous pages; includes runtime and allocator overhead.
    if (line.find("anon=") == std::string::npos) continue;
    std::istringstream fields(line);
    while (fields >> word)
      if (word.size() > 3 && word[0] == 'N' && word.find('=') != std::string::npos)
      {
        const auto equal = word.find('=');
        pages[std::stoi(word.substr(1, equal - 1))] += std::stoull(word.substr(equal + 1));
      }
  }
  return pages;
}

int main(int argc, char** argv)
{
  try
  {
    if (argc != 4) throw std::runtime_error("usage: memory_probe TOTAL_KIB THREADS serial|parallel");
    char host[256] = {};
    if (gethostname(host, sizeof(host)) || std::string(host) != "gnode1")
      throw std::runtime_error("Authorized bounded gnode1 CPU diagnostic only");
    const auto kib = std::stoull(argv[1]);
    const int threads = std::stoi(argv[2]);
    const std::string touch(argv[3]);
    if (kib < 64 || kib > 262144 || threads < 1 || threads > 16 ||
        (touch != "serial" && touch != "parallel"))
      throw std::runtime_error("Invalid bounded workload");
    const std::size_t n = kib * 1024 / (3 * sizeof(double));
    double *a = nullptr, *b = nullptr, *out = nullptr;
    if (posix_memalign(reinterpret_cast<void**>(&a), 4096, n * sizeof(double)) ||
        posix_memalign(reinterpret_cast<void**>(&b), 4096, n * sizeof(double)) ||
        posix_memalign(reinterpret_cast<void**>(&out), 4096, n * sizeof(double)))
      throw std::runtime_error("Allocation failed");
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
    std::vector<int> cpus(threads);
    #pragma omp parallel
    cpus[omp_get_thread_num()] = sched_getcpu();
    auto unique = cpus;
    std::sort(unique.begin(), unique.end());
    if (std::adjacent_find(unique.begin(), unique.end()) != unique.end() ||
        unique.front() < 0 || unique.back() >= 16)
      throw std::runtime_error("Threads must bind to distinct verified physical CPUs 0..15");
    if (touch == "parallel")
    {
      #pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < n; ++i) { a[i] = 1.0; b[i] = 2.0; out[i] = 0.0; }
    }
    else
      for (std::size_t i = 0; i < n; ++i) { a[i] = 1.0; b[i] = 2.0; out[i] = 0.0; }
    auto kernel = [&]() {
      #pragma omp parallel for simd schedule(static)
      for (std::size_t i = 0; i < n; ++i) out[i] = a[i] + 3.0 * b[i];
    };
    kernel(); // Unmeasured cache/runtime warm-up; page placement is unchanged.
    const auto pages = numaPages();
    std::vector<double> rates;
    for (int repetition = 0; repetition < 3; ++repetition)
    {
      const auto start = Clock::now();
      double elapsed = 0.0;
      std::size_t iterations = 0;
      do
      {
        kernel();
        ++iterations;
        elapsed = std::chrono::duration<double>(Clock::now() - start).count();
      } while (elapsed < 0.08);
      rates.push_back(3.0 * sizeof(double) * n * iterations / elapsed / 1e9);
    }
    long double checksum = 0;
    for (std::size_t i = 0; i < n; ++i) checksum += out[i];
    if (checksum != 7.0L * n) throw std::runtime_error("Array correctness check failed");
    std::cout << std::setprecision(12) << "{\"array_working_set_bytes\":" << 3 * sizeof(double) * n
      << ",\"threads\":" << threads << ",\"first_touch\":\"" << touch
      << "\",\"cpus\":[";
    for (int i = 0; i < threads; ++i) std::cout << (i ? "," : "") << cpus[i];
    std::cout << "],\"effective_array_GB_per_s\":[";
    for (std::size_t i = 0; i < rates.size(); ++i) std::cout << (i ? "," : "") << rates[i];
    std::cout << "],\"anonymous_numa_pages_whole_process\":{";
    bool comma = false;
    for (const auto& entry : pages)
    {
      std::cout << (comma ? "," : "") << "\"" << entry.first << "\":" << entry.second;
      comma = true;
    }
    std::cout << "},\"page_bytes\":" << sysconf(_SC_PAGESIZE)
      << ",\"correct\":true,\"hardware_counter_bandwidth\":false}\n";
    free(a); free(b); free(out);
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
