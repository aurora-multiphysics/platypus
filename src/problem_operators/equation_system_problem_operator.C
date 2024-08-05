#include "equation_system_problem_operator.h"

namespace platypus
{
void
EquationSystemProblemOperator::SetGridFunctions()
{
  NVTX3_FUNC_RANGE();
  _trial_var_names = GetEquationSystem()->_trial_var_names;
  ProblemOperator::SetGridFunctions();
}

void
EquationSystemProblemOperator::Init(mfem::Vector & X)
{
  NVTX3_FUNC_RANGE();
  std::cout <<"EquationSystemProblemOperator Init!!" << std::endl;
  ProblemOperator::Init(X);

  GetEquationSystem()->BuildEquationSystem(_problem._bc_map);
}

} // namespace platypus