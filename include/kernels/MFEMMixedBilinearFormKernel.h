#pragma once
#include "MFEMKernel.h"

/*
Class to construct an MFEM mixed bilinear form integrator to apply to the equation system.
*/
class MFEMMixedBilinearFormKernel : public MFEMKernel
{
public:
  static InputParameters validParams();

  MFEMMixedBilinearFormKernel(const InputParameters & parameters);
  ~MFEMMixedBilinearFormKernel() = default;

  // Get name of the trial variable (gridfunction) the kernel acts on.
  // Defaults to the name of the test variable labelling the weak form.
  virtual const std::string & getTrialVariableName() const override;

  // Create transposable integrator
  mfem::BilinearFormIntegrator *
  createTransposableIntegrator(mfem::BilinearFormIntegrator * base_integrator);

  // Derived classes must implement this to provide the base integrator
  virtual mfem::BilinearFormIntegrator * createIntegrator() override = 0;

protected:
  // Name of the trial variable that the kernel is applied to.
  std::string _trial_var_name;
  // Bool controlling whether to add the transpose of the integrator to the system
  bool _transpose;
};
