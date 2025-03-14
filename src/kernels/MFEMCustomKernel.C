#include "MFEMCustomKernel.h"

registerMooseObject("PlatypusApp", MFEMCustomKernel);

InputParameters
MFEMCustomKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the bilinear form "
                             "$(k u, v)_\\Omega$ "
                             "arising from the weak form of the mass operator "
                             "$ku$.");
  params.addParam<std::string>("coefficient", "Name of property for the mass coefficient k.");
  return params;
}

MFEMCustomKernel::MFEMCustomKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _coef_name(getParam<std::string>("coefficient")),
    _coef(getMFEMProblem().getProperties().getScalarProperty(_coef_name))
{
}

mfem::real_t
MFEMCustomKernel::computeQpJacobian()
{
  // the partial derivative of _grad_u is just _grad_phi[_j]
  // return _test[_i][_qp] * (_velocity * _grad_phi[_j][_qp]);

  return _test[_i] * _phi[_j];
}

mfem::BilinearFormIntegrator *
MFEMCustomKernel::createBFIntegrator()
{
  return new platypus::MFEMCustomIntegrator(this);
}