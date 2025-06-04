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

  GetEquationSystem()->BuildEquationSystem();
}

void
EquationSystemProblemOperator::Solve(mfem::Vector & X)
{
  std::cout << "about to Build Jacobian\n";
  GetEquationSystem()->BuildJacobian(_true_x, _true_rhs);
  std::cout << "Built Jacobian\n";

  _problem._nonlinear_solver->SetSolver(*_problem._jacobian_solver);
  std::cout << "Set Solver\n";
  _problem._nonlinear_solver->SetOperator(*GetEquationSystem());
  std::cout << "Set Operator\n";
  if ( _problem._nonlinear_solver != nullptr ) std::cout << "Non linear solver is not null\n";

  std::cout << "_true_rhs size is " << _true_rhs.Size() << ". _true_x size is " << _true_x.Size() << "\n";

  _problem._nonlinear_solver->Mult(_true_rhs, _true_x);
  std::cout << "Finished MUlt\n";

  GetEquationSystem()->RecoverFEMSolution(_true_x, _problem._gridfunctions);
}

} // namespace platypus
