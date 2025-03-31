#include "MFEMVectorFEFunctionBoundaryFluxIntegratedBC.h"

registerMooseObject("PlatypusApp", MFEMVectorFEFunctionBoundaryFluxIntegratedBC);

InputParameters
MFEMVectorFEFunctionBoundaryFluxIntegratedBC::validParams()
{
  InputParameters params = MFEMIntegratedBC::validParams();
  params.addRequiredParam<FunctionName>("function", "The function to be used in the integrated BC");
  return params;
}

// TODO: Currently assumes the vector function coefficient is 3D
MFEMVectorFEFunctionBoundaryFluxIntegratedBC::MFEMVectorFEFunctionBoundaryFluxIntegratedBC(
    const InputParameters & parameters)
  : MFEMIntegratedBC(parameters),
    _coef(getMFEMProblem().getScalarFunctionCoefficient(getParam<FunctionName>("function")))
{
}

// Create a new MFEM integrator to apply to the RHS of the weak form. Ownership managed by the
// caller.
mfem::LinearFormIntegrator *
MFEMVectorFEFunctionBoundaryFluxIntegratedBC::createLFIntegrator()
{
  return new mfem::VectorFEBoundaryFluxLFIntegrator(*_coef);
}

// Create a new MFEM integrator to apply to LHS of the weak form. Ownership managed by the caller.
mfem::BilinearFormIntegrator *
MFEMVectorFEFunctionBoundaryFluxIntegratedBC::createBFIntegrator()
{
  return nullptr;
}
