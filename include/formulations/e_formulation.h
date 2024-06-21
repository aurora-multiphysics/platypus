#pragma once
#include "../common/pfem_extras.hpp"
#include "hcurl_formulation.h"
#include "inputs.h"

namespace platypus
{

class EFormulation : public platypus::HCurlFormulation
{
public:
  EFormulation(const std::string & magnetic_reluctivity_name,
               std::string magnetic_permeability_name,
               const std::string & electric_conductivity_name,
               const std::string & e_field_name);

  ~EFormulation() override = default;

  void RegisterCoefficients() override;

  // Enable auxiliary calculation of J ∈ H(div)
  void RegisterCurrentDensityAux(const std::string & j_field_name,
                                 const std::string & external_j_field_name = "") override;

  // Enable auxiliary calculation of P ∈ L2
  void RegisterJouleHeatingDensityAux(const std::string & p_field_name,
                                      const std::string & e_field_name,
                                      const std::string & j_field_name) override;
  void RegisterJouleHeatingDensityAux(const std::string & p_field_name,
                                      const std::string & e_field_name) override;

protected:
  const std::string _magnetic_permeability_name;
  const std::string & _magnetic_reluctivity_name = platypus::HCurlFormulation::_alpha_coef_name;
  const std::string & _electric_conductivity_name = platypus::HCurlFormulation::_beta_coef_name;
};

} // namespace platypus
