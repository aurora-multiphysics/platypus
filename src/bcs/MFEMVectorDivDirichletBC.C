#include "MFEMVectorDivDirichletBC.h"

registerMooseObject("PlatypusApp", MFEMVectorDivDirichletBC);

// TODO: Currently assumes the vector function coefficient is 3D
MFEMVectorDivDirichletBC::MFEMVectorDivDirichletBC(const InputParameters & parameters)
  : MFEMVectorDirichletBCBase(parameters, MFEMVectorDirichletBCBase::APPLY_TYPE::NORMAL)
{
}
