#include "MFEMLinearElasticityKernel.h"

registerMooseObject("PlatypusApp", MFEMLinearElasticityKernel);

InputParameters
MFEMLinearElasticityKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("blah blah");

  params.addParam<std::string>("lambda",
                               "Name of MFEM Lame costant lambda to multiply the div(u)*I term by");
  params.addParam<std::string>("mu",
                               "Name of MFEM Lame costant mu to multiply the gradients term by");

  return params;
}

MFEMLinearElasticityKernel::MFEMLinearElasticityKernel(const InputParameters & parameters)
  : MFEMKernel(parameters),
    _lambda_name(getParam<std::string>("lambda")),
    _mu_name(getParam<std::string>("mu")),
    _lambda(getMFEMProblem()._coefficients._scalars.Get(_lambda_name)),
    _mu(getMFEMProblem()._coefficients._scalars.Get(_mu_name))
{
}

mfem::BilinearFormIntegrator *
MFEMLinearElasticityKernel::createIntegrator()
{
  return new mfem::ElasticityIntegrator(*_lambda, *_mu);
}
