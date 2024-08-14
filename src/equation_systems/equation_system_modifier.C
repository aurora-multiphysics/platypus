#include "equation_system_modifier.h"

namespace platypus
{

bool
EquationSystemModifier::VectorContainsName(const std::vector<std::string> & the_vector,
                                           const std::string & name) const
{

  auto iter = std::find(the_vector.begin(), the_vector.end(), name);

  return (iter != the_vector.end());
}

void
EquationSystemModifier::AddTrialVariableNameIfMissing(const std::string & trial_var_name)
{
  if (!VectorContainsName(getData()->_trial_var_names, trial_var_name))
  {
    getData()->_trial_var_names.push_back(trial_var_name);
  }
}

void
EquationSystemModifier::AddTestVariableNameIfMissing(const std::string & test_var_name)
{
  if (!VectorContainsName(getData()->_test_var_names, test_var_name))
  {
    getData()->_test_var_names.push_back(test_var_name);
  }
}

void
EquationSystemModifier::ApplyBoundaryConditions(platypus::BCMap & bc_map)
{
  getData()->_ess_tdof_lists.resize(getData()->_test_var_names.size());
  for (int i = 0; i < getData()->_test_var_names.size(); i++)
  {
    auto test_var_name = getData()->_test_var_names.at(i);
    // Set default value of gridfunction used in essential BC. Values
    // overwritten in applyEssentialBCs
    *(getData()->_bc_gridfunc.at(i)) = 0.0;
    bc_map.ApplyEssentialBCs(test_var_name,
                             getData()->_ess_tdof_lists.at(i),
                             *(getData()->_bc_gridfunc.at(i)),
                             getData()->_test_pfespaces.at(i)->GetParMesh());
    bc_map.ApplyIntegratedBCs(test_var_name,
                              getData()->_lfs.GetRef(test_var_name),
                              getData()->_test_pfespaces.at(i)->GetParMesh());
  }
}

void
EquationSystemModifier::BuildLinearForms(platypus::BCMap & bc_map)
{
  // Register linear forms
  for (int i = 0; i < getData()->_test_var_names.size(); i++)
  {
    auto test_var_name = getData()->_test_var_names.at(i);
    getData()->_lfs.Register(
        test_var_name, std::make_shared<mfem::ParLinearForm>(getData()->_test_pfespaces.at(i)));
    getData()->_lfs.GetRef(test_var_name) = 0.0;
  }
  // Apply boundary conditions
  ApplyBoundaryConditions(bc_map);

  for (auto & test_var_name : getData()->_test_var_names)
  {
    // Apply kernels
    auto lf = getData()->_lfs.Get(test_var_name);
    if (getData()->_lf_kernels_map.Has(test_var_name))
    {
      auto lf_kernels = getData()->_lf_kernels_map.GetRef(test_var_name);

      for (auto & lf_kernel : lf_kernels)
      {
        lf->AddDomainIntegrator(lf_kernel->createIntegrator());
      }
    }
    lf->Assemble();
  }
}

void
EquationSystemModifier::BuildBilinearForms()
{
  // Register bilinear forms
  for (int i = 0; i < getData()->_test_var_names.size(); i++)
  {
    auto test_var_name = getData()->_test_var_names.at(i);
    getData()->_blfs.Register(
        test_var_name, std::make_shared<mfem::ParBilinearForm>(getData()->_test_pfespaces.at(i)));

    // Apply kernels
    auto blf = getData()->_blfs.Get(test_var_name);
    if (getData()->_blf_kernels_map.Has(test_var_name))
    {
      auto blf_kernels = getData()->_blf_kernels_map.GetRef(test_var_name);

      for (auto & blf_kernel : blf_kernels)
      {
        blf->AddDomainIntegrator(blf_kernel->createIntegrator());
      }
    }
    // Assemble
    blf->Assemble();
  }
}

void
EquationSystemModifier::BuildMixedBilinearForms()
{
  // Register mixed bilinear forms. Note that not all combinations may
  // have a kernel

  // Create mblf for each test/trial pair
  for (int i = 0; i < getData()->_test_var_names.size(); i++)
  {
    auto test_var_name = getData()->_test_var_names.at(i);
    auto test_mblfs = std::make_shared<platypus::NamedFieldsMap<mfem::ParMixedBilinearForm>>();
    for (int j = 0; j < getData()->_test_var_names.size(); j++)
    {
      auto trial_var_name = getData()->_test_var_names.at(j);

      // Register MixedBilinearForm if kernels exist for it, and assemble
      // kernels
      if (getData()->_mblf_kernels_map_map.Has(test_var_name) &&
          getData()->_mblf_kernels_map_map.Get(test_var_name)->Has(trial_var_name))
      {
        auto mblf_kernels =
            getData()->_mblf_kernels_map_map.GetRef(test_var_name).GetRef(trial_var_name);
        auto mblf = std::make_shared<mfem::ParMixedBilinearForm>(getData()->_test_pfespaces.at(j),
                                                                 getData()->_test_pfespaces.at(i));
        // Apply all mixed kernels with this test/trial pair
        for (auto & mblf_kernel : mblf_kernels)
        {
          mblf->AddDomainIntegrator(mblf_kernel->createIntegrator());
        }
        // Assemble mixed bilinear forms
        mblf->Assemble();
        // Register mixed bilinear forms associated with a single trial variable
        // for the current test variable
        test_mblfs->Register(trial_var_name, mblf);
      }
    }
    // Register all mixed bilinear form sets associated with a single test
    // variable
    getData()->_mblfs.Register(test_var_name, test_mblfs);
  }
}

void
EquationSystemModifier::BuildEquationSystem(platypus::BCMap & bc_map)
{
  BuildLinearForms(bc_map);
  BuildBilinearForms();
  BuildMixedBilinearForms();
}

} // namespace platypus