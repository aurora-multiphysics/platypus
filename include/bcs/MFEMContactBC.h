#pragma once
#include "MFEMBoundaryCondition.h"

class MFEMContactBC : public MFEMBoundaryCondition
{
public:
  static InputParameters validParams();

  MFEMContactBC(const InputParameters & parameters);
  ~MFEMContactBC() override = default;

  void ApplyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_);

};
