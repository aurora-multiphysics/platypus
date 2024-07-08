#include "steady_state_problem_builder.h"

namespace platypus
{

void
SteadyStateProblemBuilder::SetOperatorGridFunctions()
{
  GetProblem()->GetOperator()->SetGridFunctions();
}

void
SteadyStateProblemBuilder::ConstructOperator()
{
  GetProblem()->_problem_operator.reset();
  GetProblem()->_problem_operator = std::make_unique<platypus::ProblemOperator>(*GetProblem());
}

void
SteadyStateProblemBuilder::ConstructState()
{
  auto problem_operator = GetProblem()->GetOperator();

  GetProblem()->_f =
      std::make_unique<mfem::BlockVector>(problem_operator->_true_offsets); // Vector of dofs
  problem_operator->Init(*(GetProblem()->_f)); // Set up initial conditions
}
} // namespace platypus
