#include "MFEMVectorFEDomainLFCoupledAuxKernel.h"

registerMooseObject("PlatypusApp", MFEMVectorFEDomainLFCoupledAuxKernel);

InputParameters
MFEMVectorFEDomainLFCoupledAuxKernel::validParams()
{
  InputParameters params = MFEMKernel::validParams();
  params.addClassDescription("Adds the domain integrator to an MFEM problem for the linear form "
                             "$(\\vec f, \\vec v)_\\Omega$ "
                             "arising from the weak form of the forcing term $\\vec f$.");
  params.addRequiredParam<AuxVariableName>("trial_variable", "The name of the coupled auxvariable");
  params.addParam<std::string>("coefficient",
                               "Name of scalar property to multiply the auxvariable by");
  return params;
}

MFEMVectorFEDomainLFCoupledAuxKernel::MFEMVectorFEDomainLFCoupledAuxKernel(
    const InputParameters & parameters)
  : MFEMKernel(parameters),
    _scalar_coef_name(getParam<std::string>("coefficient")),
    _scalar_coef(getMFEMProblem().getProperties().getScalarProperty(_scalar_coef_name)),
    _trial_auxvar(getMFEMProblem().getProblemData()._gridfunctions.Get(
        getParam<AuxVariableName>("trial_variable"))),
    _vec_coef(std::make_shared<mfem::VectorGridFunctionCoefficient>(_trial_auxvar)),
    _source_coef(std::make_shared<mfem::ScalarVectorProductCoefficient>(_scalar_coef, *_vec_coef))
{
}

mfem::LinearFormIntegrator *
MFEMVectorFEDomainLFCoupledAuxKernel::createIntegrator()
{
  return new mfem::VectorFEDomainLFIntegrator(*_source_coef);
}
