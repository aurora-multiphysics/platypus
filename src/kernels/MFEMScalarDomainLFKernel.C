#include "MFEMScalarDomainLFKernel.h"

registerMooseObject("PlatypusApp", MFEMScalarDomainLFKernel);

InputParameters
MFEMScalarDomainLFKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the linear form "
                             "$(\\vec f, \\vec v)_\\Omega$ "
                             "arising from the weak form of the forcing term $\\vec f$.");
  params.addParam<std::string>("coefficient", "Name of the material property $\\vec f$.");
  return params;
}

MFEMScalarDomainLFKernel::MFEMScalarDomainLFKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _coef_name(getParam<std::string>("coefficient")),
    _coef(getMFEMProblem().getProperties().getScalarProperty(_coef_name))
{
}

mfem::LinearFormIntegrator *
MFEMScalarDomainLFKernel::createLFIntegrator()
{
  return new mfem::DomainLFIntegrator(_coef);
}
