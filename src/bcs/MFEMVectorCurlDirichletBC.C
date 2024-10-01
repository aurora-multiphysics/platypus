#include "MFEMVectorCurlDirichletBC.h"

registerMooseObject("PlatypusApp", MFEMVectorCurlDirichletBC);

// TODO: Currently assumes the vector function coefficient is 3D
MFEMVectorCurlDirichletBC::MFEMVectorCurlDirichletBC(const InputParameters & parameters)
  : MFEMVectorDirichletBCBase(parameters, platypus::VectorDirichletBC::APPLY_TYPE::TANGENTIAL)
{
}
