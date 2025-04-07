#pragma once
#include "MFEMMixedBilinearFormKernel.h"

/*
(λ∇×u, v)
*/
class MFEMMixedScalarCurlKernel : public MFEMMixedBilinearFormKernel
{
public:
  static InputParameters validParams();

  MFEMMixedScalarCurlKernel(const InputParameters & parameters);
  ~MFEMMixedScalarCurlKernel() override {}

  virtual mfem::BilinearFormIntegrator * createBFIntegrator() override;

protected:
  MFEMScalarCoefficientName _coef_name;
  mfem::Coefficient & _coef;
};
