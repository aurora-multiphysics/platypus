#include "MFEMVectorCurlDirichletBC.h"

registerMooseObject("PlatypusApp", MFEMVectorCurlDirichletBC);

// TODO: Currently assumes the vector function coefficient is 3D
MFEMVectorCurlDirichletBC::MFEMVectorCurlDirichletBC(const InputParameters & parameters)
  : MFEMVectorDirichletBCBase(parameters, MFEMVectorDirichletBCBase::APPLY_TYPE::TANGENTIAL)
{
}
