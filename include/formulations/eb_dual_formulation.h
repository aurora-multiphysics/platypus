#pragma once
#include "../common/pfem_extras.hpp"
#include "dual_formulation.h"
#include "inputs.h"

namespace platypus
{

class EBDualFormulation : public platypus::DualFormulation
{
public:
  EBDualFormulation(const std::string & magnetic_reluctivity_name,
                    std::string magnetic_permeability_name,
                    const std::string & electric_conductivity_name,
                    const std::string & e_field_name,
                    const std::string & b_field_name);

  ~EBDualFormulation() override = default;

  // Enable auxiliary calculation of J ∈ H(div)
  void RegisterCurrentDensityAux(const std::string & j_field_name,
                                 const std::string & external_j_field_name = "") override;

  // Enable auxiliary calculation of F ∈ L2
  void RegisterLorentzForceDensityAux(const std::string & f_field_name,
                                      const std::string & b_field_name,
                                      const std::string & j_field_name) override;

  // Enable auxiliary calculation of P ∈ L2
  void RegisterJouleHeatingDensityAux(const std::string & p_field_name,
                                      const std::string & e_field_name,
                                      const std::string & j_field_name) override;
  void RegisterJouleHeatingDensityAux(const std::string & p_field_name,
                                      const std::string & e_field_name) override;

  void RegisterCoefficients() override;

protected:
  const std::string _magnetic_permeability_name;
  const std::string & _magnetic_reluctivity_name = platypus::DualFormulation::_alpha_coef_name;
  const std::string & _electric_conductivity_name = platypus::DualFormulation::_beta_coef_name;
};

} // namespace platypus
