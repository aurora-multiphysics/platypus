#include "MFEMVectorDirichletBCBase.h"

InputParameters
MFEMVectorDirichletBCBase::validParams()
{
  InputParameters params = MFEMBoundaryCondition::validParams();
  params.addRequiredParam<UserObjectName>(
      "vector_coefficient", "The vector MFEM coefficient to use in the Dirichlet condition");
  return params;
}

// TODO: Currently assumes the vector function coefficient is 3D
MFEMVectorDirichletBCBase::MFEMVectorDirichletBCBase(
    const InputParameters & parameters, platypus::VectorDirichletBC::APPLY_TYPE apply_type)
  : MFEMBoundaryCondition(parameters),
    _vec_coef(const_cast<MFEMVectorCoefficient *>(
        &getUserObject<MFEMVectorCoefficient>("vector_coefficient")))
{
  _boundary_condition =
      std::make_shared<platypus::VectorDirichletBC>(getParam<std::string>("variable"),
                                                    bdr_attr,
                                                    _vec_coef->getVectorCoefficient().get(),
                                                    nullptr,
                                                    apply_type);
}
