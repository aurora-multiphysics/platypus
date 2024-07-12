#pragma once
#include "MFEMHypreAMGSolver.h"
#include "MFEMProblem.h"

registerMooseObject("PlatypusApp", MFEMHypreAMGSolver);

InputParameters
MFEMHypreAMGSolver::validParams()
{
  InputParameters params = MFEMSolverBase::validParams();

  params.addParam<double>("tolerance", 1e-16, "Set the tolerance.");
  params.addParam<int>("max_iteration", 1000, "Set the maximum iteration.");

  // NB: - will require option to add a preconditioner somewhere.

  return params;
}

MFEMHypreAMGSolver::MFEMHypreAMGSolver(const InputParameters & parameters)
  : MFEMSolverBase(parameters)
{
  constructSolver(parameters);
}

void
MFEMHypreAMGSolver::constructSolver(const InputParameters & parameters)
{
  _solver = std::make_shared<mfem::HypreBoomerAMG>();

  // Set solver options.
  _solver->SetTol(parameters.get<double>("tolerance"));
  _solver->SetMaxIter(parameters.get<int>("max_iteration"));

  // NB: - set the preconditioner here.
}