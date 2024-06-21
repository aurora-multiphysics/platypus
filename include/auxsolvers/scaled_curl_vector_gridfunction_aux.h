#pragma once
#include "scaled_vector_gridfunction_aux.h"

namespace platypus
{

// Scale the curl of a gridfunction in H(Curl) by a scalar Coefficient, and
// store the result. Suitable for solving for H(Div) or H(Curl) conforming
// fields for expressions like E = ρ∇×H
class ScaledCurlVectorGridFunctionAux : public ScaledVectorGridFunctionAux
{
public:
  ScaledCurlVectorGridFunctionAux(
      const std::string & input_gf_name,
      const std::string & scaled_gf_name,
      const std::string & coef_name,
      const double & aConst = 1.0,
      const double & bConst = 1.0,
      const std::string & shift_gf_name = "",
      const platypus::InputParameters & solver_options = platypus::InputParameters());

  ~ScaledCurlVectorGridFunctionAux() override = default;

  void BuildMixedBilinearForm() override;
};
} // namespace platypus
