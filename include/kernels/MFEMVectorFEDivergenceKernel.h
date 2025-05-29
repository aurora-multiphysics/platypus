#pragma once
#include "MFEMMixedBilinearFormKernel.h"

/*
(k ∇ · u, v)
*/
class MFEMVectorFEDivergenceKernel : public MFEMMixedBilinearFormKernel
{
public:
  static InputParameters validParams();

  MFEMVectorFEDivergenceKernel(const InputParameters & parameters);
  ~MFEMVectorFEDivergenceKernel() override = default;

  virtual mfem::BilinearFormIntegrator * createMBFIntegrator() override;

protected:
  mfem::Coefficient & _coef;
};
