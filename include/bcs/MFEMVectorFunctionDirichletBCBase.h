#pragma once

#include "MFEMEssentialBC.h"

class MFEMVectorFunctionDirichletBCBase : public MFEMEssentialBC
{
public:
  static InputParameters validParams();

  ~MFEMVectorFunctionDirichletBCBase() override = default;

protected:
  MFEMVectorFunctionDirichletBCBase(const InputParameters & parameters);
  mfem::VectorCoefficient & _vec_coef;
};
