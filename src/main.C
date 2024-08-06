#include "PlatypusTestApp.h"
#include "MooseMain.h"
#include "../../NVTX/c/include/nvtx3/nvtx3.hpp"

// Begin the main program.
int
main(int argc, char * argv[])
{
  NVTX3_FUNC_RANGE();
  Moose::main<PlatypusTestApp>(argc, argv);

  return 0;
}