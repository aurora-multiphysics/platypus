#pragma once
#include "MFEMKernel.h"

/*
(α∇×u, ∇×u')
*/
class MFEMCurlCurlKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMCurlCurlKernel(const InputParameters & parameters);
  ~MFEMCurlCurlKernel() override {}

  virtual mfem::BilinearFormIntegrator * createBFIntegrator() override;

protected:
  const platypus::MFEMScalarCoefficientName & _coef_name;
  mfem::Coefficient & _coef;
};
