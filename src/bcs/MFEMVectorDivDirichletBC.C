#include "MFEMVectorDivDirichletBC.h"

registerMooseObject("PlatypusApp", MFEMVectorDivDirichletBC);

// TODO: Currently assumes the vector function coefficient is 3D
MFEMVectorDivDirichletBC::MFEMVectorDivDirichletBC(const InputParameters & parameters)
  : MFEMVectorDirichletBCBase(parameters, platypus::VectorDirichletBC::APPLY_TYPE::NORMAL)
{
}
