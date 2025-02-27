#pragma once
#include "MFEMMixedBilinearFormKernel.h"

/*
(σ ∇ V, u')
*/
class MFEMVectorFEDivergenceKernel : public MFEMMixedBilinearFormKernel
{
public:
  static InputParameters validParams();

  MFEMVectorFEDivergenceKernel(const InputParameters & parameters);
  ~MFEMVectorFEDivergenceKernel() override = default;

protected:
  mfem::BilinearFormIntegrator * buildIntegrator() const override;

  std::string _coef_name;
  mfem::Coefficient & _coef;
};
