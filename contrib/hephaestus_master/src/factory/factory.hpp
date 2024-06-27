#pragma once
#include "../common/pfem_extras.hpp"
#include "a_formulation.hpp"
#include "av_formulation.hpp"
#include "complex_a_formulation.hpp"
#include "complex_e_formulation.hpp"
#include "e_formulation.hpp"
#include "eb_dual_formulation.hpp"
#include "h_formulation.hpp"
#include "inputs.hpp"
#include "magnetostatic_formulation.hpp"

namespace hephaestus
{

class Factory
{
public:
  static std::shared_ptr<hephaestus::ProblemBuilder>
  CreateProblemBuilder(std::string & formulation_name);

  static std::shared_ptr<hephaestus::TimeDomainEMFormulation>
  CreateTimeDomainEmFormulation(std::string & formulation);

  static std::shared_ptr<hephaestus::FrequencyDomainEMFormulation>
  CreateFrequencyDomainEmFormulation(std::string & formulation);

  static std::shared_ptr<mfem::ParFiniteElementSpace>
  CreateParFESpace(const hephaestus::InputParameters params, mfem::ParMesh & pmesh);
};

} // namespace hephaestus
