#pragma once

#include "MFEMEssentialBC.h"
#include "MFEMVectorFunctionCoefficient.h"
#include "boundary_conditions.h"

class MFEMVectorDirichletBCBase : public MFEMEssentialBC
{
public:
  static InputParameters validParams();
  void ApplyBC(mfem::GridFunction & gridfunc, mfem::Mesh * mesh_) override;

  ~MFEMVectorDirichletBCBase() override = default;

  enum APPLY_TYPE
  {
    STANDARD,
    TANGENTIAL,
    NORMAL
  };
  APPLY_TYPE _boundary_apply_type;

protected:
  MFEMVectorDirichletBCBase(const InputParameters & parameters,
                            MFEMVectorDirichletBCBase::APPLY_TYPE apply_type);
  MFEMVectorCoefficient * _vec_coef{nullptr};
};
