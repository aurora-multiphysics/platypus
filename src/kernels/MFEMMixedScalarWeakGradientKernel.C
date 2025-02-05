#include "MFEMMixedScalarWeakGradientKernel.h"
#include "MFEMProblem.h"

registerMooseObject("PlatypusApp", MFEMMixedScalarWeakGradientKernel);

InputParameters
MFEMMixedScalarWeakGradientKernel::validParams()
{
  InputParameters params = MFEMMixedScalarWeakGradientKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the bilinear form "
                             "$-(\\lambda u, \\nabla v)_\\Omega$ "
                             "arising from the weak form of the scalar weak derivative operator "
                             "$\\lambda u$.");
  params.addParam<std::string>("coefficient",
                               "Name of scalar property lambda to multiply the integrator by.");
  return params;
}

MFEMMixedScalarWeakGradientKernel::MFEMMixedScalarWeakGradientKernel(const InputParameters & parameters)
  : MFEMMixedScalarWeakGradientKernel(parameters),
    _coef_name(getParam<std::string>("coefficient")),
    _coef(getMFEMProblem().getProperties().getScalarProperty(_coef_name))
{
}

mfem::MixedScalarWeakDerivativeIntegrator *
MFEMMixedScalarWeakGradientKernel::createIntegrator()
{
  return new mfem::MixedScalarWeakDerivativeIntegrator(_coef);
}
