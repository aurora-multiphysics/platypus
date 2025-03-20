#pragma once
#include "MFEMKernel.h"

/*
(\\vec f, \\vec u')
*/
class MFEMVectorFEDomainLFKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMVectorFEDomainLFKernel(const InputParameters & parameters);
  ~MFEMVectorFEDomainLFKernel() override {}

  virtual mfem::LinearFormIntegrator * createResidualContribution() override;

protected:
  std::shared_ptr<mfem::VectorCoefficient> _vec_coef{nullptr};
};
