#pragma once
#include "MFEMMixedBilinearFormKernel.h"

/*
(σ ∇ V, u')
*/
class MFEMMixedVectorGradientKernel : public MFEMMixedBilinearFormKernel
{
public:
  static InputParameters validParams();

  MFEMMixedVectorGradientKernel(const InputParameters & parameters);
  ~MFEMMixedVectorGradientKernel() override = default;

protected:
  mfem::BilinearFormIntegrator * createIntegrator() const override;

  std::string _coef_name;
  mfem::Coefficient & _coef;
};
