#include "mfem.hpp"
#include "logging.hpp"
#include <catch2/catch_session.hpp>
#include <iostream>

const char * DATA_DIR = "../data/";

int
main(int argc, char * argv[])
{

  MPI_Init(&argc, &argv);
  int result = Catch::Session().run(argc, argv);
  MPI_Finalize();

  return result;
}
