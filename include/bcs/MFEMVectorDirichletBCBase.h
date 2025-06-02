#pragma once

#include "MFEMEssentialBC.h"

class MFEMVectorDirichletBCBase : public MFEMEssentialBC
{
public:
  static InputParameters validParams();

  ~MFEMVectorDirichletBCBase() override = default;

protected:
  MFEMVectorDirichletBCBase(const InputParameters & parameters);
  std::vector<Real> _vec_value;
  mfem::VectorCoefficient & _vec_coef;
};
