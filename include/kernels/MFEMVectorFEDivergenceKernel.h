#pragma once
#include "MFEMKernel.h"

/*
(∇·u, v)
*/
class MFEMVectorFEDivergenceKernel : public MFEMKernel<mfem::BilinearFormIntegrator>
{
public:
  static InputParameters validParams();

  MFEMVectorFEDivergenceKernel(const InputParameters & parameters);
  ~MFEMVectorFEDivergenceKernel() = default;

  // Get name of the trial variable (gridfunction) the kernel acts on.
  // Defaults to the name of the test variable labelling the weak form.
  virtual const std::string & getTrialVariableName() const override;

protected:
  // Name of the trial variable that the kernel is applied to.
  std::string _trial_var_name;
  std::string _coef_name;
  mfem::Coefficient & _coef;
};
