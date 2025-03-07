#pragma once
#include "MFEMKernel.h"

/*
(βu, u')
*/
class MFEMVectorFEMassKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMVectorFEMassKernel(const InputParameters & parameters);
  ~MFEMVectorFEMassKernel() override {}

  virtual mfem::BilinearFormIntegrator * createBFIntegrator() override;

protected:
  std::string _coef_name;
  mfem::Coefficient & _coef;
};
