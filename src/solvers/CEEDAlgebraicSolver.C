#pragma once
#include "CEEDAlgebraicSolver.h"
#include "MFEMProblem.h"

registerMooseObject("PlatypusApp", CEEDAlgebraicSolver);

InputParameters
CEEDAlgebraicSolver::validParams()
{
  InputParameters params = MFEMSolverBase::validParams();
  params.addClassDescription("CEED solver for the iterative solution of MFEM equation "
                             "systems using the algebraic multigrid method.");

  params.addParam<double>("l_tol", 1e-5, "Set the relative tolerance.");
  params.addParam<double>("l_abs_tol", 1e-50, "Set the absolute tolerance.");
  params.addParam<int>("l_max_its", 10000, "Set the maximum number of iterations.");
  params.addParam<int>("print_level", 2, "Set the solver verbosity.");
  params.addParam<UserObjectName>("preconditioner", "Optional choice of preconditioner to use.");

  return params;
}

CEEDAlgebraicSolver::CEEDAlgebraicSolver(const InputParameters & parameters)
  : MFEMSolverBase(parameters),
    _preconditioner(isParamSetByUser("preconditioner")
                        ? getUserObject<MFEMSolverBase>("preconditioner").getSolver()
                        : nullptr)
{
  constructSolver(parameters);
}

void
CEEDAlgebraicSolver::constructSolver(const InputParameters & parameters)
{
  auto preconditioner = std::dynamic_pointer_cast<mfem::Solver>(_preconditioner);

  _solver = std::make_shared<mfem::ceed::AlgebraicSolver>(getMFEMProblem().mesh().getMFEMParMesh().GetComm());
  _solver->SetRelTol(getParam<double>("l_tol"));
  _solver->SetAbsTol(getParam<double>("l_abs_tol"));
  _solver->SetMaxIter(getParam<int>("l_max_its"));
  _solver->SetPrintLevel(getParam<int>("print_level"));

  if (preconditioner)
    _solver->SetPreconditioner(*preconditioner);
}
