#include "MFEMVectorFEMassKernel.h"

registerMooseObject("PlatypusApp", MFEMVectorFEMassKernel);

InputParameters
MFEMVectorFEMassKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the bilinear form "
                             "$(k \\vec u, \\vec v)_\\Omega$ "
                             "arising from the weak form of the mass operator "
                             "$k \\vec u$.");
  params.addParam<MFEMScalarCoefficientName>("coefficient",
                                             "Name of property k to multiply the Laplacian by");
  return params;
}

MFEMVectorFEMassKernel::MFEMVectorFEMassKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _coef_name(getParam<MFEMScalarCoefficientName>("coefficient")),
    // FIXME: The MFEM bilinear form can also handle vector and matrix
    // coefficients, so ideally we'd handle all three too.
    _coef(getScalarCoefficient(_coef_name))
{
}

mfem::BilinearFormIntegrator *
MFEMVectorFEMassKernel::createBFIntegrator()
{
  return new mfem::VectorFEMassIntegrator(_coef);
}
