#include "equation_system_assembler.h"

namespace platypus
{

void
EquationSystemAssembler::MakeData()
{
  if (_equation_system_data != nullptr){
    mooseError("_equation_system_data has already been set!")
  }
  else{
    _equation_system_data = std::make_shared<platypus::EquationSystemData>();
  }
}

bool
EquationSystemAssembler::VectorContainsName(const std::vector<std::string> & the_vector,
                                   const std::string & name) const
{

  auto iter = std::find(the_vector.begin(), the_vector.end(), name);

  return (iter != the_vector.end());
}

void
EquationSystemAssembler::AddTrialVariableNameIfMissing(const std::string & trial_var_name)
{
  if (!VectorContainsName(getData()->_trial_var_names, trial_var_name))
  {
    getData()->_trial_var_names.push_back(trial_var_name);
  }
}

void
EquationSystemAssembler::AddTestVariableNameIfMissing(const std::string & test_var_name)
{
  if (!VectorContainsName(getData()->_test_var_names, test_var_name))
  {
    getData()->_test_var_names.push_back(test_var_name);
  }
}

void
EquationSystemAssembler::AddKernel(const std::string & test_var_name,
                          std::shared_ptr<MFEMBilinearFormKernel> blf_kernel)
{
  AddTestVariableNameIfMissing(test_var_name);

  if (!getData()->_blf_kernels_map.Has(test_var_name))
  {
    // 1. Create kernels vector.
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMBilinearFormKernel>>>();

    // 2. Register with map to prevent leaks.
    getData()->_blf_kernels_map.Register(test_var_name, std::move(kernels));
  }

  getData()->_blf_kernels_map.GetRef(test_var_name).push_back(std::move(blf_kernel));
}

void
EquationSystemAssembler::AddKernel(const std::string & test_var_name,
                          std::shared_ptr<MFEMLinearFormKernel> lf_kernel)
{
  AddTestVariableNameIfMissing(test_var_name);

  if (!getData()->_lf_kernels_map.Has(test_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMLinearFormKernel>>>();

    getData()->_lf_kernels_map.Register(test_var_name, std::move(kernels));
  }

  getData()->_lf_kernels_map.GetRef(test_var_name).push_back(std::move(lf_kernel));
}

void
EquationSystemAssembler::AddKernel(const std::string & test_var_name,
                          std::shared_ptr<MFEMNonlinearFormKernel> nlf_kernel)
{
  AddTestVariableNameIfMissing(test_var_name);

  if (!getData()->_nlf_kernels_map.Has(test_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMNonlinearFormKernel>>>();

    getData()->_nlf_kernels_map.Register(test_var_name, std::move(kernels));
  }

  getData()->_nlf_kernels_map.GetRef(test_var_name).push_back(std::move(nlf_kernel));
}

void
EquationSystemAssembler::AddKernel(const std::string & trial_var_name,
                          const std::string & test_var_name,
                          std::shared_ptr<MFEMMixedBilinearFormKernel> mblf_kernel)
{
  if (getData()->_assembly_level != mfem::AssemblyLevel::LEGACY)
    mooseError("Mixed Bilinear Form Kernels are currently only compatible with the LEGACY assembly level.");

  AddTestVariableNameIfMissing(test_var_name);

  // Register new mblf kernels map if not present for this test variable
  if (!getData()->_mblf_kernels_map_map.Has(test_var_name))
  {
    auto kernel_field_map = std::make_shared<
        platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMMixedBilinearFormKernel>>>>();

    getData()->_mblf_kernels_map_map.Register(test_var_name, std::move(kernel_field_map));
  }

  // Register new mblf kernels map if not present for the test/trial variable
  // pair
  if (!getData()->_mblf_kernels_map_map.Get(test_var_name)->Has(trial_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMMixedBilinearFormKernel>>>();

    getData()->_mblf_kernels_map_map.Get(test_var_name)->Register(trial_var_name, std::move(kernels));
  }

  getData()->_mblf_kernels_map_map.GetRef(test_var_name)
      .Get(trial_var_name)
      ->push_back(std::move(mblf_kernel));
}

void
EquationSystemAssembler::ApplyBoundaryConditions(platypus::BCMap & bc_map)
{
  getData()->_ess_tdof_lists.resize(getData()->_test_var_names.size());
  for (int i = 0; i < getData()->_test_var_names.size(); i++)
  {
    auto test_var_name = getData()->_test_var_names.at(i);
    // Set default value of gridfunction used in essential BC. Values
    // overwritten in applyEssentialBCs
    *(getData()->_bc_gridfunc.at(i)) = 0.0;
    bc_map.ApplyEssentialBCs(
        test_var_name, getData()->_ess_tdof_lists.at(i), *(getData()->_bc_gridfunc.at(i)), getData()->_test_pfespaces.at(i)->GetParMesh());
    bc_map.ApplyIntegratedBCs(
        test_var_name, getData()->_lfs.GetRef(test_var_name), getData()->_test_pfespaces.at(i)->GetParMesh());
  }
}


BlockDiagonalEquationSystemAssembler::BlockDiagonalEquationSystemAssembler(mfem::AssemblyLevel assembly_level)
{
    MakeData();
    getData()->_assembly_level = assembly_level;

}

NonBlockDiagonalEquationSystemAssembler::NonBlockDiagonalEquationSystemAssembler()
{
    MakeData();
    getData()->_assembly_level = mfem::AssemblyLevel::LEGACY;

}








/*
bool
EquationSystemAssembler::IsBlockDiagonal()
{
    bool block_diag = true;

  for (int i = 0; i < _test_var_names.size(); ++i)
  {
    for (int j = 0; j < _test_var_names.size(); ++j)
    {
        if (i != j && getData()->_h_blocks(i,j) != nullptr)
        {
            block_diag = false;
            break;
        }
    }
  }

    return block_diag
}
*/

} // namespace platypus
