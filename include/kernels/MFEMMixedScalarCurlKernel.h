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

protected:
  mfem::BilinearFormIntegrator * buildIntegrator() const override;

  std::string _coef_name;
  mfem::Coefficient & _coef;
};
