#include "MFEMVectorFunctionNormalDirichletBC.h"

registerMooseObject("PlatypusApp", MFEMVectorFunctionNormalDirichletBC);

InputParameters
MFEMVectorFunctionNormalDirichletBC::validParams()
{
  InputParameters params = MFEMVectorFunctionDirichletBCBase::validParams();
  params.addClassDescription(
      "Applies a Dirichlet condition to the normal components of a vector variable.");
  return params;
}

MFEMVectorFunctionNormalDirichletBC::MFEMVectorFunctionNormalDirichletBC(
    const InputParameters & parameters)
  : MFEMVectorFunctionDirichletBCBase(parameters)
{
}

void
MFEMVectorFunctionNormalDirichletBC::ApplyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_)
{
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = getBoundaries();
  gridfunc.ProjectBdrCoefficientNormal(*_vec_coef, ess_bdrs);
}
