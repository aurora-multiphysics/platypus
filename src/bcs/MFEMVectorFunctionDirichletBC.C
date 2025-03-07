#include "MFEMVectorFunctionDirichletBC.h"

registerMooseObject("PlatypusApp", MFEMVectorFunctionDirichletBC);

MFEMVectorFunctionDirichletBC::MFEMVectorFunctionDirichletBC(const InputParameters & parameters)
  : MFEMVectorFunctionDirichletBCBase(parameters)
{
}

void
MFEMVectorFunctionDirichletBC::ApplyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_)
{
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = getBoundaries();
  gridfunc.ProjectBdrCoefficient(*_vec_coef, ess_bdrs);
}
