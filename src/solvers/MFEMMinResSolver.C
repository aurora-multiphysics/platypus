#pragma once
#include "MFEMMinResSolver.h"
#include "MFEMProblem.h"

registerMooseObject("PlatypusApp", MFEMMinResSolver);

InputParameters
MFEMMinResSolver::validParams()
{
  InputParameters params = MFEMSolverBase::validParams();
  params.addClassDescription("MFEMMinResSolver solver and preconditioner for the iterative solution "
                             "of MFEM equation systems.");
  params.addParam<double>("rel_tol", 1e-12, "Set the relative tolerance.");
  params.addParam<int>("l_max_its", 2000, "Set the maximum number of iterations.");
  params.addParam<int>("print_level", 2, "Set the solver verbosity.");
  params.addParam<UserObjectName>("preconditioner", "Optional choice of preconditioner to use.");
  return params;
}

MFEMMinResSolver::MFEMMinResSolver(const InputParameters & parameters)
  : MFEMSolverBase(parameters),
    _preconditioner(isParamSetByUser("preconditioner")
      ? getUserObject<MFEMSolverBase>("preconditioner").getSolver()
      : nullptr)
{
  constructSolver(parameters);
}

void
MFEMMinResSolver::constructSolver(const InputParameters & parameters)
{
  _solver = std::make_shared<mfem::MINRESSolver>( getMFEMProblem().mesh().getMFEMParMesh().GetComm() );

  _solver->SetRelTol(getParam<double>("rel_tol"));
  _solver->SetMaxIter(getParam<int>("l_max_its"));

  auto preconditioner = std::dynamic_pointer_cast<mfem::Solver>(_preconditioner);

  if (preconditioner)
    _solver->SetPreconditioner(*preconditioner);
}
