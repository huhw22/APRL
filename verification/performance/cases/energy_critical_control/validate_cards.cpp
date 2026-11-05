#include "config.h"
#include <exception>
#include <iostream>

// Exercise the real YAML validator without constructing a Simulation or MPI.
int main(int argc, char** argv) {
  if (argc < 2) return 2;
  try {
    for (int i=1; i<argc; ++i) {
      aprl::YamlConfigLoader::loadFile(argv[i]);
      std::cout << "VALID " << argv[i] << '\n';
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
  return 0;
}
