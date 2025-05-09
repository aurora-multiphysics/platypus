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
EquationSystemProblemOperator::SetUpAMR()
{
  mfem::BilinearFormIntegrator* integ;

  auto eqn_system = _problem._eqn_system;

  // now we need to take a reference to the eqn system and run through each test var name
  // for each test var name, we get a pointer to the domain integrator
  auto test_var_name  = (_problem._eqn_system->_test_var_names)[0];
  auto trial_var_name = (_problem._eqn_system->_trial_var_names)[0];

  // correct blf
  auto blf = eqn_system->_blfs.GetShared(test_var_name);

  // returns an array containing all of our domain integrators
  integ = *(blf->GetDBFI())[0];

  // next make the FE Collections

  // Really bad - just get the FE order from the first FEC in the problem
  auto iter = _problem._fecs.begin();
  int order = iter->second->GetOrder();
  int dim   = _problem._pmesh->Dimension();
  int sdim  = _problem._pmesh->SpaceDimension();

  /*
  Set up error estimator. As per example 6p, we supply a space for the discontinuous
  flux (L2) and a space for the smoothed flux.
  */
  _flux_fec = std::make_unique<mfem::L2_FECollection>(order, dim);
  _flux_fes = std::make_unique<mfem::ParFiniteElementSpace>(_problem._pmesh.get(), _flux_fec.get(), sdim);

  _smooth_flux_fec = std::make_unique<mfem::H1_FECollection>(order, dim);
  _smooth_flux_fes = std::make_unique<mfem::ParFiniteElementSpace>(_problem._pmesh.get(), _smooth_flux_fec.get(), dim);

  
  // only one grid function - it's the concentration!
  auto gridIter = _problem._gridfunctions.begin();
  _error_estimator = std::make_unique<mfem::L2ZienkiewiczZhuEstimator>( *integ, *(gridIter->second), *_flux_fes, *_smooth_flux_fes );

  // next we can init the refiner
  _refiner = std::make_unique<mfem::ThresholdRefiner>(*_error_estimator);
  _refiner->SetTotalErrorFraction(0.7);

}


bool 
EquationSystemProblemOperator::HRefine()
{
  _refiner->Apply( *_problem._pmesh );

  return _refiner->Stop();
}


void
EquationSystemProblemOperator::Solve(mfem::Vector & X)
{
  SetUpAMR();

  GetEquationSystem()->BuildJacobian(_true_x, _true_rhs);

  _problem._nonlinear_solver->SetSolver(*_problem._jacobian_solver);
  _problem._nonlinear_solver->SetOperator(*GetEquationSystem());
  _problem._nonlinear_solver->Mult(_true_rhs, _true_x);

  GetEquationSystem()->RecoverFEMSolution(_true_x, _problem._gridfunctions);
}

} // namespace platypus
