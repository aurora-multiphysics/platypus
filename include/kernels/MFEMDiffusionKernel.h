#pragma once
#include "MFEMKernel.h"

/*
(σ ∇ q, ∇ q')
*/
class MFEMDiffusionKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMDiffusionKernel(const InputParameters & parameters);

  virtual mfem::BilinearFormIntegrator * createJacobianContribution() override;

protected:
  std::string _coef_name;
  mfem::Coefficient & _coef;
};
