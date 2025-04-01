#pragma once
#include "MFEMKernel.h"

/*
(βu, u')
*/
class MFEMMassKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMMassKernel(const InputParameters & parameters);
  ~MFEMMassKernel() override {}

  virtual mfem::BilinearFormIntegrator * createBFIntegrator() override;

protected:
  const MFEMScalarCoefficientName & _coef_name;
  mfem::Coefficient & _coef;
};
