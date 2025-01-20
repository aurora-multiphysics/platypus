#pragma once
#include "MFEMKernel.h"

/*
(βu, u')
*/
class MFEMVectorFEMassKernel : public MFEMKernel<mfem::BilinearFormIntegrator>
{
public:
  static InputParameters validParams();

  MFEMVectorFEMassKernel(const InputParameters & parameters);
  ~MFEMVectorFEMassKernel() override {}

  virtual mfem::BilinearFormIntegrator * createIntegrator() override;

protected:
  const platypus::MFEMScalarCoefficientName & _coef_name;
  mfem::Coefficient & _coef;
};
