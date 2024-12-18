#pragma once
#include "MFEMKernel.h"

/*
(α∇.u, ∇.u')
*/
class MFEMDivDivKernel : public MFEMKernel<mfem::BilinearFormIntegrator>
{
public:
  static InputParameters validParams();

  MFEMDivDivKernel(const InputParameters & parameters);
  ~MFEMDivDivKernel() override {}

  virtual mfem::BilinearFormIntegrator * createIntegrator() override;

protected:
  const platypus::MFEMScalarCoefficientName & _coef_name;
  mfem::Coefficient & _coef;
};
