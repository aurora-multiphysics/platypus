#pragma once
#include "MFEMEssentialBC.h"

class MFEMScalarFunctionDirichletBC : public MFEMEssentialBC
{
public:
  static InputParameters validParams();

  MFEMScalarFunctionDirichletBC(const InputParameters & parameters);

  void ApplyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_) override;

protected:
  std::string _coef_name;
  mfem::Coefficient & _coef;
};
