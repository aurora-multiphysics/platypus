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

  mfem::BilinearFormIntegrator * createIntegrator() override;

protected:
  std::string _coef_name;
  mfem::Coefficient & _coef;
};
