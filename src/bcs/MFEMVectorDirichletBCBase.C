#include "MFEMVectorDirichletBCBase.h"

registerMooseObject("PlatypusApp", MFEMVectorDirichletBC);

InputParameters
MFEMVectorDirichletBCBase::validParams()
{
  InputParameters params = MFEMEssentialBC::validParams();
  params.addRequiredParam<UserObjectName>(
      "vector_coefficient", "The vector MFEM coefficient to use in the Dirichlet condition");
  return params;
}

// TODO: Currently assumes the vector function coefficient is 3D
MFEMVectorDirichletBCBase::MFEMVectorDirichletBCBase(const InputParameters & parameters,
        MFEMVectorDirichletBCBase::APPLY_TYPE apply_type)
  : MFEMEssentialBC(parameters),
    _boundary_apply_type(apply_type),
    _vec_coef(const_cast<MFEMVectorCoefficient *>(
        &getUserObject<MFEMVectorCoefficient>("vector_coefficient")))
{
}

void
MFEMVectorDirichletBCBase::ApplyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_)
{
  mfem::Array<int> ess_bdrs(mesh_->bdr_attributes.Max());
  ess_bdrs = GetMarkers(*mesh_);
  if (_vec_coef == nullptr)
  {
    MFEM_ABORT("Boundary condition does not store valid coefficients to specify the "
               "components of the vector at the Dirichlet boundary.");
  }

  switch (_boundary_apply_type)
  {
    case STANDARD:
      gridfunc.ProjectBdrCoefficient(*_vec_coef->getVectorCoefficient(), ess_bdrs);
      break;
    case NORMAL:
      gridfunc.ProjectBdrCoefficientNormal(*_vec_coef->getVectorCoefficient(), ess_bdrs);
      break;
    case TANGENTIAL:
      gridfunc.ProjectBdrCoefficientTangent(*_vec_coef->getVectorCoefficient(), ess_bdrs);
  }
}
