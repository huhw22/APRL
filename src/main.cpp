#include <exception>
#include <iostream>
#include <string>

#include <mpi.h>

#include "config.h"
#include "runtime_control.h"
#include "runtime_util.h"
#include "simulation.h"

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc != 2)
    {
      if (rank == 0)
        std::cerr << "Usage: simulator <configuration.yaml>" << std::endl;
      MPI_Finalize();
      return 1;
    }

  const double start = MPI_Wtime();
  try
    {
      const fel::SimulationConfig config =
        fel::YamlConfigLoader::loadFile(argv[1]);
      if (config.runtime.interactive())
        fel::RuntimeControl::installSignalHandlers();

      fel::logRoot(MPI_COMM_WORLD,
        "Unnamed FEL simulator (temporary executable name: simulator)");
      fel::Simulation simulation(config, MPI_COMM_WORLD);
      simulation.solve();

      if (rank == 0)
        std::cout << "Elapsed wall time [s]: "
                  << (MPI_Wtime() - start) << std::endl;
    }
  catch (const std::exception& error)
    {
      std::cerr << "Simulation failure on MPI rank " << rank
                << ": " << error.what() << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }

  MPI_Finalize();
  return 0;
}
