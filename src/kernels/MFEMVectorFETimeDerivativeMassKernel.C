#include "MFEMVectorFETimeDerivativeMassKernel.h"

registerMooseObject("PlatypusApp", MFEMVectorFETimeDerivativeMassKernel);

InputParameters
MFEMVectorFETimeDerivativeMassKernel::validParams()
{
  InputParameters params = MFEMMassKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the bilinear form "
                             "$(k \\dot{\\vec u}, \\vec v)_\\Omega$ "
                             "arising from the weak form of the operator "
                             "$k \\dot{ \\vec u}$.");
  return params;
}

MFEMVectorFETimeDerivativeMassKernel::MFEMVectorFETimeDerivativeMassKernel(
    const InputParameters & parameters)
  : MFEMVectorFEMassKernel(parameters),
    _var_dot_name(platypus::GetTimeDerivativeName(_test_var_name))
{
}
