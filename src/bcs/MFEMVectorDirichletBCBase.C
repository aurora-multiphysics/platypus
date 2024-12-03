#include "MFEMVectorDirichletBCBase.h"

InputParameters
MFEMVectorDirichletBCBase::validParams()
{
  InputParameters params = MFEMEssentialBC::validParams();
  params.addRequiredParam<std::string>(
      "vector_coefficient", "The vector MFEM coefficient to use in the Dirichlet condition");
  return params;
}

// TODO: Currently assumes the vector function coefficient is 3D
MFEMVectorDirichletBCBase::MFEMVectorDirichletBCBase(const InputParameters & parameters)
  : MFEMEssentialBC(parameters),
    _vec_coef_name(getParam<std::string>("vector_coefficient")),
    _vec_coef(getMFEMProblem().getProperties().getVectorProperty(_vec_coef_name))
{
}
