#pragma once
#include "MFEMKernel.h"

/*
(f, u')
*/
class MFEMDomainLFKernel : public MFEMKernel<mfem::LinearFormIntegrator>
{
public:
  static InputParameters validParams();

  MFEMDomainLFKernel(const InputParameters & parameters);
  ~MFEMDomainLFKernel() override {}

  virtual mfem::LinearFormIntegrator * createLFIntegrator() override;

protected:
  std::string _coef_name;
  mfem::Coefficient & _coef;
};
