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
EquationSystemProblemOperator::SetUpAMR( std::string estimator_type, std::string estimator_name, InputParameters estimator_params )
{
  _use_amr = true;

  // initialise the MFEMEstimator
  InputParameters base_params = MFEMEstimator::validParams();
  
  // create a MFEMEstimator base object on the stack. That way we have access to the
  // necessary member functions (and the reference to the mfem problem) that means
  // we can constuct the child class we need
  MFEMEstimator base_estimator(base_params);
  
  _estimator = base_estimator.setUp(estimator_type, estimator_name, estimator_params );
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
  // SetUpAMR(); // this should be handled already

  GetEquationSystem()->BuildJacobian(_true_x, _true_rhs);

  _problem._nonlinear_solver->SetSolver(*_problem._jacobian_solver);
  _problem._nonlinear_solver->SetOperator(*GetEquationSystem());
  _problem._nonlinear_solver->Mult(_true_rhs, _true_x);

  GetEquationSystem()->RecoverFEMSolution(_true_x, _problem._gridfunctions);
}

} // namespace platypus
