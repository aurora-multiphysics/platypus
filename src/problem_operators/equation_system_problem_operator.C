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
EquationSystemProblemOperator::AddEstimator(std::shared_ptr<MFEMEstimator> estimator)
{
  _estimator = estimator;
}


void
EquationSystemProblemOperator::SetUpAMR()
{
  _use_amr = true;

  _refiner   = std::make_unique<mfem::ThresholdRefiner>( *_estimator->createEstimator() );
  _refiner->SetTotalErrorFraction(0.7);
}

/*
The refiner will return true if we have met the stopping condition for refinement.
If we don't use AMR, we should just return false.

TODO: should probably call an error if _use_amr is set to false

*/
bool 
EquationSystemProblemOperator::HRefine()
{
  if ( _use_amr ) {
    _refiner->Apply( *_problem._pmesh );
  
    return _refiner->Stop();
  }
  return false;
}


void
EquationSystemProblemOperator::Solve(mfem::Vector & X)
{
  GetEquationSystem()->BuildJacobian(_true_x, _true_rhs);
  
  _problem._nonlinear_solver->SetSolver(*_problem._jacobian_solver);
  _problem._nonlinear_solver->SetOperator(*GetEquationSystem());
  _problem._nonlinear_solver->Mult(_true_rhs, _true_x);

  GetEquationSystem()->RecoverFEMSolution(_true_x, _problem._gridfunctions);
}

} // namespace platypus
