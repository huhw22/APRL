// Parse only: no MPI initialization, particle loading or simulation stepping.
#include <exception>
#include <iostream>
#include "config.h"

int main(int argc, char** argv)
{
  if (argc < 2) return 2;
  try
    {
      for (int i = 1; i < argc; ++i)
        {
          aprl::YamlConfigLoader::loadFile(argv[i]);
          std::cout << "CONFIG_VALID " << argv[i] << '\n';
        }
    }
  catch (const std::exception& error)
    {
      std::cerr << error.what() << '\n';
      return 1;
    }
  return 0;
}
