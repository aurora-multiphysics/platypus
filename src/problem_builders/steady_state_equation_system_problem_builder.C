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

void
SteadyStateEquationSystemProblemBuilder::ConstructOperator()
{
  auto equation_system = std::make_unique<platypus::EquationSystem>();
  auto problem_operator = std::make_shared<platypus::EquationSystemProblemOperator>(
      *_problem, std::move(equation_system));

  _problem_operator = std::move(problem_operator);
}

} // namespace platypus
