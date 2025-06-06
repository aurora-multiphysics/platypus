#include "MFEMVectorDomainLFKernel.h"

registerMooseObject("PlatypusApp", MFEMVectorDomainLFKernel);

InputParameters
MFEMVectorDomainLFKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the linear form "
                             "$(\\vec f, \\vec v)_\\Omega$ "
                             "arising from the weak form of the forcing term $\\vec f$.");
  params.addParam<MFEMVectorCoefficientName>("vector_coefficient",
                                             "Name of body force density $\\vec f$.");
  return params;
}

MFEMVectorDomainLFKernel::MFEMVectorDomainLFKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _vec_coef_name(getParam<MFEMVectorCoefficientName>("vector_coefficient")),
    _vec_coef(getVectorCoefficient(_vec_coef_name))
{
}

mfem::LinearFormIntegrator *
MFEMVectorDomainLFKernel::createLFIntegrator()
{
  return new mfem::VectorDomainLFIntegrator(_vec_coef);
}
