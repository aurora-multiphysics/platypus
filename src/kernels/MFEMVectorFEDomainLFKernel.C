#include "MFEMVectorFEDomainLFKernel.h"

registerMooseObject("PlatypusApp", MFEMVectorFEDomainLFKernel);

InputParameters
MFEMVectorFEDomainLFKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the linear form "
                             "$(\\vec f, \\vec v)_\\Omega$ "
                             "arising from the weak form of the forcing term $\\vec f$.");
  params.addParam<std::string>("vector_coefficient", "Name of MFEM vector coefficient f.");

  return params;
}

MFEMVectorFEDomainLFKernel::MFEMVectorFEDomainLFKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _vec_coef_name(getParam<std::string>("vector_coefficient")),
    _vec_coef(getMFEMProblem().getProblemData()._coefficients._vectors.Get(_vec_coef_name))
{
}

mfem::LinearFormIntegrator *
MFEMVectorFEDomainLFKernel::createIntegrator()
{
  return new mfem::VectorFEDomainLFIntegrator(*_vec_coef);
}
