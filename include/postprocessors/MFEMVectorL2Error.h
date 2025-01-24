#pragma once
#include "MFEMPostprocessor.h"
#include "MFEMGeneralUserObject.h"

/*
 * Compute the L2 error for a vector variable.
 */
class MFEMVectorL2Error : public MFEMPostprocessor
{
public:
  static InputParameters validParams();

  MFEMVectorL2Error(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;

  /**
   * Get the L2 Error.
   */
  virtual PostprocessorValue getValue() const override final;

private:
  VariableName _var_name;
  FunctionName _coeff_name;
  std::shared_ptr<mfem::VectorCoefficient> _vec_coeff;
  std::shared_ptr<mfem::GridFunction> _var;
};
