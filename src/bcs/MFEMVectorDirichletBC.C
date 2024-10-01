#include "MFEMVectorDirichletBC.h"

registerMooseObject("PlatypusApp", MFEMVectorDirichletBC);

// TODO: Currently assumes the vector function coefficient is 3D
MFEMVectorDirichletBC::MFEMVectorDirichletBC(const InputParameters & parameters)
  : MFEMVectorDirichletBCBase(parameters, MFEMVectorDirichletBCBase::APPLY_TYPE::STANDARD)
{
}
