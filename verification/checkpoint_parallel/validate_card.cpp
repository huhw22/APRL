#include "config.h"
#include <exception>
#include <iostream>
int main(int argc, char** argv)
{
  if (argc < 2) return 2;
  try {
    for (int i=1;i<argc;++i) {
      aprl::SimulationConfig config=aprl::YamlConfigLoader::loadFile(argv[i]);
      std::cout << argv[i] << ": valid; cells=" << config.mesh.cells[0]
                << ',' << config.mesh.cells[1] << ',' << config.mesh.cells[2]
                << '\n';
    }
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n'; return 1;
  }
}
