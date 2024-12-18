#pragma once
#include "MFEMKernel.h"

/*
(σ u, ∇ V')
*/
class MFEMVectorFEWeakDivergenceKernel : public MFEMKernel<mfem::BilinearFormIntegrator>
{
public:
  static InputParameters validParams();

  MFEMVectorFEWeakDivergenceKernel(const InputParameters & parameters);
  ~MFEMVectorFEWeakDivergenceKernel() override = default;

  virtual mfem::BilinearFormIntegrator * createIntegrator() override;

protected:
  const platypus::MFEMScalarCoefficientName & _coef_name;
  mfem::Coefficient & _coef;
};
