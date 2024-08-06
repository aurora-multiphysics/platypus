#include "time_domain_problem_operator.h"

namespace platypus
{

std::string
GetTimeDerivativeName(const std::string & name)
{
  NVTX3_FUNC_RANGE();
  return std::string("d") + name + std::string("_dt");
}

std::vector<std::string>
GetTimeDerivativeNames(std::vector<std::string> gridfunction_names)
{
  NVTX3_FUNC_RANGE();
  std::vector<std::string> time_derivative_names;
  for (auto & gridfunction_name : gridfunction_names)
  {
    time_derivative_names.push_back(GetTimeDerivativeName(gridfunction_name));
  }
  return time_derivative_names;
}

void
TimeDomainProblemOperator::SetGridFunctions()
{
  NVTX3_FUNC_RANGE();
  ProblemOperatorInterface::SetGridFunctions();
  width = height = _true_offsets[_trial_variables.size()];
}

} // namespace platypus