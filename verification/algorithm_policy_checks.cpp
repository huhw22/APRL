#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include "config.h"

namespace {
void require(bool ok, const std::string& message) {
  if (!ok) throw std::runtime_error(message);
}
aprl::SimulationConfig load(const std::string& source,
                            const std::string& options) {
  std::string text = source;
  const std::size_t at = text.find("mesh:\n");
  require(at != std::string::npos, "Missing mesh in policy fixture");
  text.insert(at + 6, options);
  std::ofstream out("algorithm-policy.yaml");
  out << text;
  out.close();
  return aprl::YamlConfigLoader::loadFile("algorithm-policy.yaml");
}
void check(const aprl::MeshConfig& mesh, bool quintic) {
  require(mesh.conservativeSourceCoupling, "Coupling fell back to Boris");
  require(!mesh.longitudinalNotch, "Retired seven-point kernel enabled");
  require(mesh.longitudinalQuintic == quintic &&
          mesh.longitudinalBinomial3 == !quintic,
          "Expected exactly one supported kernel");
}
void reject(const std::string& source, const std::string& options,
            const std::string& message) {
  bool rejected = false;
  try { load(source, options); }
  catch (const std::exception& error) {
    rejected = std::string(error.what()).find(message) != std::string::npos;
  }
  require(rejected, "Retired/invalid algorithm did not produce the expected error: " + options);
}
}

int main(int argc, char** argv) {
  try {
    require(argc == 2, "Usage: algorithm_policy_checks <valid-card>");
    std::ifstream in(argv[1]);
    require(bool(in), "Cannot open policy fixture");
    const std::string source((std::istreambuf_iterator<char>(in)),
                             std::istreambuf_iterator<char>());
    check(aprl::MeshConfig(), false);
    check(load(source, "").mesh, false);
    check(load(source, "  particle_coupling: conservative-source\n"
                       "  longitudinal_shape: linear-binomial3\n").mesh, false);
    check(load(source, "  longitudinal_shape: quintic-compensated\n").mesh, true);
    reject(source, "  particle_coupling: boris\n", "no longer supported");
    reject(source, "  particle_coupling: invalid\n", "APRL 2.0 supports only");
    reject(source, "  longitudinal_shape: linear\n", "no longer supported");
    reject(source, "  longitudinal_shape: linear-notch\n", "no longer supported");
    reject(source, "  longitudinal_shape: invalid\n", "APRL 2.0 longitudinal_shape");
    std::string filtered = source;
    const std::size_t at = filtered.find("current_filter: none");
    require(at != std::string::npos, "Missing filter in policy fixture");
    filtered.replace(at, std::string("current_filter: none").size(),
                     "current_filter: compensated-bilinear");
    reject(filtered, "", "APRL 2.0 requires current_filter: none");
    std::cout << "APRL 2.0 algorithm defaults, both supported choices, and retired-option rejection passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
