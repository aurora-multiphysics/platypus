#include "equation_system_assembler.h"

namespace platypus
{

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
    mooseError("Mixed Bilinear Form Kernels are currently only compatible with the LEGACY assembly "
               "level.");

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

    getData()
        ->_mblf_kernels_map_map.Get(test_var_name)
        ->Register(trial_var_name, std::move(kernels));
  }

  getData()
      ->_mblf_kernels_map_map.GetRef(test_var_name)
      .Get(trial_var_name)
      ->push_back(std::move(mblf_kernel));
}

void
EquationSystemAssembler::Init(platypus::GridFunctions & gridfunctions,
                              const platypus::FESpaces & fespaces,
                              platypus::BCMap & bc_map,
                              platypus::Coefficients & coefficients,
                              mfem::AssemblyLevel assembly_level)
{
  for (auto & test_var_name : getData()->_test_var_names)
  {
    getData()->_assembly_level = assembly_level;
    if (!AssemblyIsSupported())
      mooseError("Assembly type chosen is not supported for this system");

    if (!gridfunctions.Has(test_var_name))
    {
      MFEM_ABORT("Test variable " << test_var_name
                                  << " requested by equation system during initialisation was "
                                     "not found in gridfunctions");
    }
    // Store pointers to variable FESpaces
    getData()->_test_pfespaces.push_back(gridfunctions.Get(test_var_name)->ParFESpace());
    // Create auxiliary gridfunctions for applying Dirichlet conditions
    getData()->_bc_gridfunc.emplace_back(
        std::make_unique<mfem::ParGridFunction>(gridfunctions.Get(test_var_name)->ParFESpace()));
  }
}

void
EquationSystemAssembler::BuildJacobian(mfem::BlockVector & trueX, mfem::BlockVector & trueRHS)
{
  getData()->_jacobian_height = trueX.Size();
  getData()->_jacobian_width = trueRHS.Size();
  FormSystem(getData()->_jacobian, trueX, trueRHS);
}

// DiagonalEquationSystemAssembler //////////////////////////

bool
DiagonalEquationSystemAssembler::AssemblyIsSupported()
{
  bool support = true;
  return support;
}

void
DiagonalEquationSystemAssembler::FormSystem(mfem::OperatorHandle & op,
                                            mfem::BlockVector & trueX,
                                            mfem::BlockVector & trueRHS)
{
}

void
DiagonalEquationSystemAssembler::BuildEquationSystem(platypus::BCMap & bc_map)
{
  BuildLinearForms(bc_map);
  BuildBilinearForms();
}

// DenseEquationSystemAssembler //////////////////////////

bool
DenseEquationSystemAssembler::AssemblyIsSupported()
{
  bool support = true;
  if (getData()->_assembly_level != mfem::AssemblyLevel::LEGACY)
    support = false;

  return support;
}

void
DenseEquationSystemAssembler::FormSystem(mfem::OperatorHandle & op,
                                         mfem::BlockVector & trueX,
                                         mfem::BlockVector & trueRHS)
{
  // Allocate block operator
  getData()->_h_blocks.DeleteAll();
  getData()->_h_blocks.SetSize(getData()->_test_var_names.size(),
                               getData()->_test_var_names.size());
  // Form diagonal blocks.
  for (int i = 0; i < getData()->_test_var_names.size(); i++)
  {
    auto & test_var_name = getData()->_test_var_names.at(i);
    auto blf = getData()->_blfs.Get(test_var_name);
    auto lf = getData()->_lfs.Get(test_var_name);
    mfem::Vector aux_x, aux_rhs;
    mfem::HypreParMatrix aux_a;
    blf->FormLinearSystem(getData()->_ess_tdof_lists.at(i),
                          *(getData()->_bc_gridfunc.at(i)),
                          *lf,
                          aux_a,
                          aux_x,
                          aux_rhs);

    getData()->_h_blocks(i, i) = new const mfem::HypreParMatrix(aux_a);
    trueX.GetBlock(i) = aux_x;
    trueRHS.GetBlock(i) = aux_rhs;
  }

  // Form off-diagonal blocks
  for (int i = 0; i < getData()->_test_var_names.size(); i++)
  {
    auto test_var_name = getData()->_test_var_names.at(i);
    for (int j = 0; j < getData()->_test_var_names.size(); j++)
    {
      auto trial_var_name = getData()->_test_var_names.at(j);

      mfem::Vector aux_x, aux_rhs;
      mfem::ParLinearForm aux_lf(getData()->_test_pfespaces.at(i));
      aux_lf = 0.0;
      if (getData()->_mblfs.Has(test_var_name) &&
          getData()->_mblfs.Get(test_var_name)->Has(trial_var_name))
      {
        auto mblf = getData()->_mblfs.Get(test_var_name)->Get(trial_var_name);
        mfem::HypreParMatrix aux_a;
        mblf->FormRectangularLinearSystem(getData()->_ess_tdof_lists.at(j),
                                          getData()->_ess_tdof_lists.at(i),
                                          *(getData()->_bc_gridfunc.at(j)),
                                          aux_lf,
                                          aux_a,
                                          aux_x,
                                          aux_rhs);
        trueRHS.GetBlock(i) += aux_rhs;
        getData()->_h_blocks(i, j) = new const mfem::HypreParMatrix(aux_a);
      }
    }
  }
  // Sync memory
  for (int i = 0; i < getData()->_test_var_names.size(); i++)
  {
    trueX.GetBlock(0).SyncAliasMemory(trueX);
    trueRHS.GetBlock(0).SyncAliasMemory(trueRHS);
  }

  // Create monolithic matrix
  op.Reset(mfem::HypreParMatrixFromBlocks(getData()->_h_blocks));
}

void
DenseEquationSystemAssembler::BuildEquationSystem(platypus::BCMap & bc_map)
{
  BuildLinearForms(bc_map);
  BuildBilinearForms();
  BuildMixedBilinearForms();
}

} // namespace platypus
