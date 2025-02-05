#pragma once
#include "MFEMKernel.h"

/*
-(λu, ∇v)
*/
class MixedScalarWeakDerivativeKernel : public MFEMMixedBilinearFormKernel
{
public:
  static InputParameters validParams();

  MixedScalarWeakDerivativeKernel(const InputParameters & parameters);
  ~MixedScalarWeakDerivativeKernel() override {}

  virtual mfem::BilinearFormIntegrator * createIntegrator() override;

protected:
  std::string _coef_name;
  mfem::Coefficient & _coef;
};