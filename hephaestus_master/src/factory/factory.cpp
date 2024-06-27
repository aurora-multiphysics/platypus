#include "factory.hpp"

namespace hephaestus
{

std::shared_ptr<hephaestus::ProblemBuilder>
Factory::CreateProblemBuilder(std::string & formulation_name)
{
  if (formulation_name == "EBForm")
  {
    return std::make_shared<hephaestus::EBDualFormulation>("magnetic_reluctivity",
                                                           "magnetic_permeability",
                                                           "electrical_conductivity",
                                                           "electric_field",
                                                           "magnetic_flux_density");
  }
  else if (formulation_name == "HForm")
  {
    return std::make_shared<hephaestus::HFormulation>("electrical_resistivity",
                                                      "electrical_conductivity",
                                                      "magnetic_permeability",
                                                      "magnetic_field");
  }
  else if (formulation_name == "AForm")
  {
    return std::make_shared<hephaestus::AFormulation>("magnetic_reluctivity",
                                                      "magnetic_permeability",
                                                      "electrical_conductivity",
                                                      "magnetic_vector_potential");
  }
  else if (formulation_name == "EForm")
  {
    return std::make_shared<hephaestus::EFormulation>("magnetic_reluctivity",
                                                      "magnetic_permeability",
                                                      "electrical_conductivity",
                                                      "electric_field");
  }
  else if (formulation_name == "AVForm")
  {
    return std::make_shared<hephaestus::AVFormulation>("magnetic_reluctivity",
                                                       "magnetic_permeability",
                                                       "electrical_conductivity",
                                                       "magnetic_vector_potential",
                                                       "electric_potential");
  }
  else if (formulation_name == "ComplexEForm")
  {
    return std::make_shared<hephaestus::ComplexEFormulation>("magnetic_reluctivity",
                                                             "electrical_conductivity",
                                                             "dielectric_permittivity",
                                                             "frequency",
                                                             "electric_field",
                                                             "electric_field_real",
                                                             "electric_field_imag");
  }
  else if (formulation_name == "ComplexAForm")
  {
    return std::make_shared<hephaestus::ComplexAFormulation>("magnetic_reluctivity",
                                                             "electrical_conductivity",
                                                             "dielectric_permittivity",
                                                             "frequency",
                                                             "magnetic_vector_potential",
                                                             "magnetic_vector_potential_real",
                                                             "magnetic_vector_potential_imag");
  }
  else if (formulation_name == "Custom")
  {
    return std::make_shared<hephaestus::TimeDomainEMFormulation>();
  }
  else
  {
    MFEM_WARNING("Formulation name " << formulation_name << " not recognised.");
    return nullptr;
  }
};

std::shared_ptr<hephaestus::FrequencyDomainEMFormulation>
Factory::CreateFrequencyDomainEmFormulation(std::string & formulation)
{
  if (formulation == "ComplexEForm")
  {
    return std::make_shared<hephaestus::ComplexEFormulation>("magnetic_reluctivity",
                                                             "electrical_conductivity",
                                                             "dielectric_permittivity",
                                                             "frequency",
                                                             "electric_field",
                                                             "electric_field_real",
                                                             "electric_field_imag");
  }
  else if (formulation == "ComplexAForm")
  {
    return std::make_shared<hephaestus::ComplexAFormulation>("magnetic_reluctivity",
                                                             "electrical_conductivity",
                                                             "dielectric_permittivity",
                                                             "frequency",
                                                             "magnetic_vector_potential",
                                                             "magnetic_vector_potential_real",
                                                             "magnetic_vector_potential_imag");
  }
  else
  {
    MFEM_WARNING("Steady formulation name " << formulation << " not recognised.");
  }
  return nullptr;
}

std::shared_ptr<hephaestus::TimeDomainEMFormulation>
Factory::CreateTimeDomainEmFormulation(std::string & formulation)
{
  if (formulation == "EBForm")
  {
    return std::make_shared<hephaestus::EBDualFormulation>("magnetic_reluctivity",
                                                           "magnetic_permeability",
                                                           "electrical_conductivity",
                                                           "electric_field",
                                                           "magnetic_flux_density");
  }
  else if (formulation == "HForm")
  {
    return std::make_shared<hephaestus::HFormulation>("electrical_resistivity",
                                                      "electrical_conductivity",
                                                      "magnetic_permeability",
                                                      "magnetic_field");
  }
  else if (formulation == "AForm")
  {
    return std::make_shared<hephaestus::AFormulation>("magnetic_reluctivity",
                                                      "magnetic_permeability",
                                                      "electrical_conductivity",
                                                      "magnetic_vector_potential");
  }
  else if (formulation == "EForm")
  {
    return std::make_shared<hephaestus::EFormulation>("magnetic_reluctivity",
                                                      "magnetic_permeability",
                                                      "electrical_conductivity",
                                                      "electric_field");
  }
  else if (formulation == "AVForm")
  {
    return std::make_shared<hephaestus::AVFormulation>("magnetic_reluctivity",
                                                       "magnetic_permeability",
                                                       "electrical_conductivity",
                                                       "magnetic_vector_potential",
                                                       "electric_potential");
  }
  else if (formulation == "Custom")
  {
    return std::make_shared<hephaestus::TimeDomainEMFormulation>();
  }
  else
  {
    MFEM_WARNING("Formulation name " << formulation << " not recognised.");
  }
  return nullptr;
}

std::shared_ptr<mfem::ParFiniteElementSpace>
Factory::CreateParFESpace(hephaestus::InputParameters params, mfem::ParMesh & pmesh)
{
  auto fe_type(params.GetParam<std::string>("FESpaceType"));
  int order(params.GetParam<int>("order"));
  int components(params.GetParam<int>("components")); // spatial dimension of mesh. Use
                                                      // FiniteElementCollection::New instead
  if (fe_type == "H1")
  {
    return std::make_shared<mfem::common::H1_ParFESpace>(&pmesh, order, components);
  }
  else if (fe_type == "ND")
  {
    return std::make_shared<mfem::common::ND_ParFESpace>(&pmesh, order, components);
  }
  else if (fe_type == "RT")
  {
    return std::make_shared<mfem::common::RT_ParFESpace>(&pmesh, order, components);
  }
  else if (fe_type == "L2")
  {
    return std::make_shared<mfem::common::L2_ParFESpace>(&pmesh, order, components);
  }
  else
  {
    MFEM_WARNING("FESpaceType " << fe_type << " not recognised.");
  }
  return nullptr;
}

} // namespace hephaestus
