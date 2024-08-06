#include "steady_state_problem_builder.h"

namespace platypus
{

void
SteadyStateProblemBuilder::SetOperatorGridFunctions()
{
  NVTX3_FUNC_RANGE(); 
  GetProblem()->GetOperator()->SetGridFunctions();
}

void
SteadyStateProblemBuilder::ConstructOperator()
{
  NVTX3_FUNC_RANGE(); 
  GetProblem()->ConstructOperator();
}

void
SteadyStateProblemBuilder::ConstructState()
{
  NVTX3_FUNC_RANGE(); 
  auto problem_operator = GetProblem()->GetOperator();

  GetProblem()->_f =
      std::make_unique<mfem::BlockVector>(problem_operator->_true_offsets); // Vector of dofs
  problem_operator->Init(*(GetProblem()->_f)); // Set up initial conditions
}
} // namespace platypus
