#include "MFEMDomainLFKernel.h"

registerMooseObject("PlatypusApp", MFEMDomainLFKernel);

InputParameters
MFEMDomainLFKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the bilinear form "
                             "$(f, v)_\\Omega$ "
                             "arising from the weak form of the mass operator "
                             "$f$.");
  params.addParam<std::string>("coefficient", "Name of function coefficient f.");
  return params;
}

MFEMDomainLFKernel::MFEMDomainLFKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _coef_name(getParam<std::string>("coefficient")),
    _coef(getMFEMProblem().getProperties().getScalarProperty(_coef_name))
{
}

mfem::LinearFormIntegrator *
MFEMDomainLFKernel::createIntegrator()
{
  return new mfem::DomainLFIntegrator(_coef);
}