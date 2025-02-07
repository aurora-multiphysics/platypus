#pragma once
#include "mfem/miniapps/common/pfem_extras.hpp"
#include "MFEMAuxKernel.h"

/*
Class to set an H(div) auxvariable to be the curl of an H(curl) vector variable.
*/
class MFEMCurlAux : public MFEMAuxKernel
{
public:
  static InputParameters validParams();

  MFEMCurlAux(const InputParameters & parameters);

  virtual ~MFEMCurlAux() = default;

  // Computes the auxvariable.
  virtual void execute() override;

protected:
  // Name of source MFEMVariable to take the curl of.
  VariableName _source_var_name;
  // Reference to source gridfunction.
  mfem::ParGridFunction & _source_var;
  // Scalar factor to multiply the result by.
  mfem::real_t _scale_factor;
  // FESpaces
  mfem::ParFiniteElementSpace & _hcurl_fespace;
  mfem::ParFiniteElementSpace & _hdiv_fespace;
  // Curl operator
  mfem::common::ParDiscreteCurlOperator _curl;
};
