#include "MFEMScalarFunctorDirichletBC.h"

registerMooseObject("PlatypusApp", MFEMScalarFunctorDirichletBC);

InputParameters
MFEMScalarFunctorDirichletBC::validParams()
{
  InputParameters params = MFEMEssentialBC::validParams();
  params.addRequiredParam<MFEMScalarCoefficientName>(
      "coefficient",
      "The coefficient setting the values on the essential boundary. A coefficient can be any of "
      "the "
      "following: a variable, an MFEM material property, a function, or a post-processor.");
  return params;
}

MFEMScalarFunctorDirichletBC::MFEMScalarFunctorDirichletBC(const InputParameters & parameters)
  : MFEMEssentialBC(parameters),
    _coef_name(getParam<MFEMScalarCoefficientName>("coefficient")),
    _coef(getScalarCoefficient(_coef_name))
{
}

void
MFEMScalarFunctorDirichletBC::ApplyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_)
{
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = getBoundaries();
  gridfunc.ProjectBdrCoefficient(_coef, ess_bdrs);
}
