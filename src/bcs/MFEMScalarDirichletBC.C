#include "MFEMScalarDirichletBC.h"
#include "MFEMProblem.h"

registerMooseObject("PlatypusApp", MFEMScalarDirichletBC);

InputParameters
MFEMScalarDirichletBC::validParams()
{
  InputParameters params = MFEMEssentialBC::validParams();
  params.addClassDescription("Applies a Dirichlet condition to a scalar variable.");
  params.addRequiredParam<Real>("value", "The scalar value to use in the Dirichlet condition");
  return params;
}

MFEMScalarDirichletBC::MFEMScalarDirichletBC(const InputParameters & parameters)
  : MFEMEssentialBC(parameters),
    _coef(getMFEMProblem().getCoefficients().declareScalar<mfem::ConstantCoefficient>(
        "__ScalarDirichletBC_" + parameters.get<std::string>("_unique_name"),
        getParam<Real>("value")))
{
}

void
MFEMScalarDirichletBC::ApplyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_)
{
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = getBoundaries();
  gridfunc.ProjectBdrCoefficient(_coef, ess_bdrs);
}
