#include "steady_state_problem_builder.h"

namespace platypus
{

void
SteadyStateProblemBuilder::SetOperatorGridFunctions()
{
  GetOperator().SetGridFunctions();
}

void
SteadyStateProblemBuilder::ConstructOperator()
{
  _problem_operator.reset();
  _problem_operator = std::make_shared<platypus::ProblemOperator>(*GetProblem());
}

void
SteadyStateProblemBuilder::ConstructState()
{
  GetProblem()->_f =
      std::make_unique<mfem::BlockVector>(GetOperator()._true_offsets); // Vector of dofs
  GetOperator().Init(*(GetProblem()->_f));                              // Set up initial conditions
}
} // namespace platypus
