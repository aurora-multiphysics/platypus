#pragma once
#include "MFEMKernel.h"
#include "MFEMCustomIntegrator.h"
/*
(βu, u')
*/
class MFEMCustomKernel : public MFEMKernel, public platypus::MFEMIntegratorInterface
{
public:
  static InputParameters validParams();

  MFEMCustomKernel(const InputParameters & parameters);
  ~MFEMCustomKernel() override {}

  virtual mfem::real_t computeQpJacobian() override;

  virtual mfem::BilinearFormIntegrator * createBFIntegrator() override;

protected:
  std::string _coef_name;
  mfem::Coefficient & _coef;
};
