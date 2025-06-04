#include "MFEMContactBC.h"

registerMooseObject("PlatypusApp", MFEMContactBC);

InputParameters
MFEMContactBC::validParams()
{
  return MFEMEssentialBC::validParams();
}

MFEMContactBC::MFEMContactBC(const InputParameters & parameters)
  : MFEMBoundaryCondition(parameters)
{
}

// paste in something from the unit test!
void
MFEMContactBC::ApplyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_)
{

}
