#include "mfem.hpp"
#include "logging.hpp"
#include <catch2/catch_session.hpp>
#include <iostream>

const char * DATA_DIR = "../data/";

int
main(int argc, char * argv[])
{

  mfem::OptionsParser args(argc, argv);
  args.AddOption(
      &DATA_DIR, "-dataDir", "--data_directory", "Directory storing input data for tests.");
  args.Parse();
  MPI_Init(&argc, &argv);
  int result = Catch::Session().run(argc, argv);
  MPI_Finalize();

  return result;
}
