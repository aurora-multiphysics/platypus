#pragma once
#include "MFEMKernel.h"

/*
(σ u, ∇ V')
*/
class MFEMVectorFEWeakDivergenceKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMVectorFEWeakDivergenceKernel(const InputParameters & parameters);
  ~MFEMVectorFEWeakDivergenceKernel() override = default;

  virtual mfem::BilinearFormIntegrator * createBFIntegrator() override;

protected:
  std::string _coef_name;
  mfem::Coefficient & _coef;
};
