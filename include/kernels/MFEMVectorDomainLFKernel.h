#pragma once
#include "MFEMKernel.h"

/*
(\\vec f, \\vec u')
*/
class MFEMVectorDomainLFKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMVectorDomainLFKernel(const InputParameters & parameters);
  ~MFEMVectorDomainLFKernel() override {}

  virtual mfem::LinearFormIntegrator * createResidualContribution() override;

protected:
  std::string _vec_coef_name;
  mfem::VectorCoefficient & _vec_coef;
};
