#pragma once
#include "MFEMOperatorJacobiSmoother.h"
#include "MFEMProblem.h"

registerMooseObject("PlatypusApp", MFEMOperatorJacobiSmoother);

InputParameters
MFEMOperatorJacobiSmoother::validParams()
{
  InputParameters params = MFEMSolverBase::validParams();
  params.addClassDescription("MFEM solver for performing Jacobi smoothing of the equation system.");
  return params;
}

MFEMOperatorJacobiSmoother::MFEMOperatorJacobiSmoother(const InputParameters & parameters)
  : MFEMSolverBase(parameters)
{
  constructSolver(parameters);
}

void
MFEMOperatorJacobiSmoother::constructSolver(const InputParameters & parameters)
{
  _preconditioner = std::make_shared<mfem::OperatorJacobiSmoother>();
}
