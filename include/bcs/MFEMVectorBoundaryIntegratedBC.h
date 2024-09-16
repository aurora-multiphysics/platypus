#pragma once

#include "MFEMBoundaryCondition.h"
#include "MFEMVectorFunctionCoefficient.h"
#include "boundary_conditions.h"

class MFEMVectorBoundaryIntegratedBC : public MFEMBoundaryCondition
{
public:
  static InputParameters validParams();

  MFEMVectorBoundaryIntegratedBC(const InputParameters & parameters);
  ~MFEMVectorBoundaryIntegratedBC() override {}

protected:
  MFEMVectorCoefficient * _vec_coef{nullptr};
};
