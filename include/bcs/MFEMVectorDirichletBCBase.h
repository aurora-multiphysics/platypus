#pragma once

#include "MFEMEssentialBC.h"
#include "boundary_conditions.h"

class MFEMVectorDirichletBCBase : public MFEMEssentialBC
{
public:
  static InputParameters validParams();

  ~MFEMVectorDirichletBCBase() override = default;

protected:
  MFEMVectorDirichletBCBase(const InputParameters & parameters);
  std::string _vec_coef_name;
  mfem::VectorCoefficient & _vec_coef;
};
