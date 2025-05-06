#include "equation_system_problem_operator.h"

namespace platypus
{
void
EquationSystemProblemOperator::SetGridFunctions()
{
  _trial_var_names = GetEquationSystem()->_trial_var_names;
  ProblemOperator::SetGridFunctions();
}

void
EquationSystemProblemOperator::Init(mfem::BlockVector & X)
{
  ProblemOperator::Init(X);

  GetEquationSystem()->BuildEquationSystem(_problem._gridfunctions);
}

void
EquationSystemProblemOperator::Solve(mfem::Vector & X)
{
  GetEquationSystem()->BuildJacobian(_true_x, _true_rhs);
  mfem::Vector zero_vec(_true_rhs.Size());
  zero_vec = 0.0;
  _problem._nonlinear_solver->SetSolver(*_problem._jacobian_solver);
  _problem._nonlinear_solver->SetOperator(*GetEquationSystem());
  _problem._nonlinear_solver->Mult(zero_vec, _true_x);

  GetEquationSystem()->RecoverFEMSolution(_true_x, _problem._gridfunctions);
}

} // namespace platypus
