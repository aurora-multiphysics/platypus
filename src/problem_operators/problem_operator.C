#include "problem_operator.h"

namespace platypus
{

void
ProblemOperator::SetGridFunctions()
{
  NVTX3_FUNC_RANGE();
  ProblemOperatorInterface::SetGridFunctions();
  width = height = _true_offsets[_trial_variables.size()];
};

} // namespace platypus