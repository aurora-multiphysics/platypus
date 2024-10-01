#pragma once

#include "MFEMBoundaryCondition.h"
#include "MFEMVectorFunctionCoefficient.h"
#include "boundary_conditions.h"

class MFEMVectorDirichletBCBase : public MFEMBoundaryCondition
{
public:
  static InputParameters validParams();

  ~MFEMVectorDirichletBCBase() override {}

protected:
  MFEMVectorDirichletBCBase(const InputParameters & parameters,
                            platypus::VectorDirichletBC::APPLY_TYPE apply_type);
  MFEMVectorCoefficient * _vec_coef{nullptr};
};
