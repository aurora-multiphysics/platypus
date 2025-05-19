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

  virtual mfem::BilinearFormIntegrator * createBFIntegrator() override;

protected:
  std::string _coef_name;
  mfem::Coefficient & _coef;
};
