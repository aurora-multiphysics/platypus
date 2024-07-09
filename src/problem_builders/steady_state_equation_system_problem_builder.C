#include "steady_state_equation_system_problem_builder.h"

namespace platypus
{

void
SteadyStateEquationSystemProblemBuilder::InitializeKernels()
{
  ProblemBuilder::InitializeKernels();

  GetEquationSystem()->Init(
      _problem->_gridfunctions, _problem->_fespaces, _problem->_bc_map, _problem->_coefficients);
}

} // namespace platypus
