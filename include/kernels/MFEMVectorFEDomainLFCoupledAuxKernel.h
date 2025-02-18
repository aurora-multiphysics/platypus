#pragma once
#include "MFEMKernel.h"

/*
(\\vec f, \\vec u')
*/
class MFEMVectorFEDomainLFCoupledAuxKernel : public MFEMKernel<mfem::LinearFormIntegrator>
{
public:
  static InputParameters validParams();

  MFEMVectorFEDomainLFCoupledAuxKernel(const InputParameters & parameters);
  ~MFEMVectorFEDomainLFCoupledAuxKernel() override {}

  virtual mfem::LinearFormIntegrator * createIntegrator() override;

protected:
  std::string _scalar_coef_name;
  mfem::Coefficient & _scalar_coef;
  const mfem::ParGridFunction * _trial_auxvar;
  std::shared_ptr<mfem::VectorCoefficient> _vec_coef{nullptr};
  std::shared_ptr<mfem::VectorCoefficient> _source_coef{nullptr};
};
