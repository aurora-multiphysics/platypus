#include "MFEMCurlCurlKernel.h"
#include "MFEMProblem.h"

registerMooseObject("PlatypusApp", MFEMCurlCurlKernel);

InputParameters
MFEMCurlCurlKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription(
      "Adds the domain integrator to an MFEM problem for the bilinear form "
      "$(k\\vec\\nabla \\times \\vec u, \\vec\\nabla \\times \\vec v)_\\Omega$ "
      "arising from the weak form of the curl curl operator "
      "$-k\\vec\\nabla \\times \\vec\\nabla \\times \\vec u$.");
  params.addParam<MFEMScalarCoefficientName>(
      "coefficient", "Name of scalar coefficient k to multiply the integrator by.");
  return params;
}

MFEMCurlCurlKernel::MFEMCurlCurlKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _coef_name(getParam<MFEMScalarCoefficientName>("coefficient")),
    // FIXME: The MFEM bilinear form can also handle vector and matrix
    // coefficients, so ideally we'd handle all three too.
    _coef(getScalarCoefficient(_coef_name))
{
}

mfem::BilinearFormIntegrator *
MFEMCurlCurlKernel::createBFIntegrator()
{
  return new mfem::CurlCurlIntegrator(_coef);
}
