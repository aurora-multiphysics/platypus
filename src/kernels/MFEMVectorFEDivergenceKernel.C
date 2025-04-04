#include "MFEMVectorFEDivergenceKernel.h"

registerMooseObject("PlatypusApp", MFEMVectorFEDivergenceKernel);

InputParameters
MFEMVectorFEDivergenceKernel::validParams()
{
  InputParameters params = MFEMMixedBilinearFormKernel::validParams();
  params.addClassDescription(
      "Adds the domain integrator to an MFEM problem for the mixed bilinear form "
      "$(k\\vec\\nabla u, \\vec v)_\\Omega$ "
      "arising from the weak form of the gradient operator "
      "$k\\vec \\nabla u$.");
  params.addParam<std::string>("coefficient", "Name of property k to use.");
  return params;
}

MFEMVectorFEDivergenceKernel::MFEMVectorFEDivergenceKernel(const InputParameters & parameters)
  : MFEMMixedBilinearFormKernel(parameters),
    _coef_name(getParam<std::string>("coefficient")),
    // FIXME: The MFEM bilinear form can also handle vector and matrix
    // coefficients, so ideally we'd handle all three too.
    _coef(getMFEMProblem().getProperties().getScalarProperty(_coef_name))
{
}

mfem::BilinearFormIntegrator *
MFEMVectorFEDivergenceKernel::createBFIntegrator()
{
  mfem::BilinearFormIntegrator * base_integrator = new mfem::VectorFEDivergenceIntegrator(_coef);
  return createTransposableBFIntegrator(base_integrator);
}
