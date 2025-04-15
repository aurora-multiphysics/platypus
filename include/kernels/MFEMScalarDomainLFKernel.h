#pragma once
#include "MFEMKernel.h"

/*
(\\vec f, \\vec u')
*/
class MFEMScalarDomainLFKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMScalarDomainLFKernel(const InputParameters & parameters);
  ~MFEMScalarDomainLFKernel() override {}

  virtual mfem::LinearFormIntegrator * createLFIntegrator() override;

protected:
  std::string _coef_name;
  mfem::Coefficient & _coef;
};
