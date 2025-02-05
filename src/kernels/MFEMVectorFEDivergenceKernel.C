#include "MFEMVectorFEDivergenceKernel.h"

registerMooseObject("PlatypusApp", MFEMVectorFEDivergenceKernel);

InputParameters
MFEMVectorFEDivergenceKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the bilinear form "
                             "$-(\\nabla \\cdot \\vec{u}, v)_\\Omega$ "
                             "arising from the weak form of the vector finite element divergence operator.");
  params.addParam<std::string>("coefficient", "Name of property k to multiply the integrator by");
  params.addParam<std::string>(
      "trial_variable",
      "",
      "The trial variable this kernel is acting on and which will be solved for. If empty "
      "(default), it will be the same as the test variable.");
  return params;
}

MFEMVectorFEDivergenceKernel::MFEMVectorFEDivergenceKernel(const InputParameters & parameters)
  : MFEMKernel(parameters), _trial_var_name(getParam<std::string>("trial_variable"))
{
  if (_trial_var_name == "")
  {
    _trial_var_name = _test_var_name;
  }
}

const std::string &
MFEMVectorFEDivergenceKernel::getTrialVariableName() const
{
  return _trial_var_name;
}

MFEMVectorFEDivergenceKernel::MFEMVectorFEDivergenceKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _coef_name(getParam<std::string>("coefficient")),
    _coef(getMFEMProblem().getProperties().getScalarProperty(_coef_name))
{
}

mfem::BilinearFormIntegrator *
MFEMVectorFEDivergenceKernel::createIntegrator()
{
  return new mfem::VectorFEDivergenceIntegrator(_coef);
}
