#include "equation_system_operator.h"

namespace platypus
{

void
EquationSystemOperator::AddTrialVariableNameIfMissing(const std::string & trial_var_name)
{
  if (!VectorContainsName(DataRead()->_trial_var_names, trial_var_name))
  {
    DataWrite()->_trial_var_names.push_back(trial_var_name);
  }
}

void
EquationSystemOperator::AddTestVariableNameIfMissing(const std::string & test_var_name)
{
  if (!VectorContainsName(DataRead()->_test_var_names, test_var_name))
  {
    DataWrite()->_test_var_names.push_back(test_var_name);
  }
}

void
EquationSystemOperator::ApplyBoundaryConditions(platypus::BCMap & bc_map)
{
  DataWrite()->_ess_tdof_lists.resize(DataRead()->_test_var_names.size());
  for (int i = 0; i < DataRead()->_test_var_names.size(); i++)
  {
    auto test_var_name = DataRead()->_test_var_names.at(i);
    // Set default value of gridfunction used in essential BC. Values
    // overwritten in applyEssentialBCs
    *(DataWrite()->_bc_gridfunc.at(i)) = 0.0;
    bc_map.ApplyEssentialBCs(test_var_name,
                             DataWrite()->_ess_tdof_lists.at(i),
                             *(DataWrite()->_bc_gridfunc.at(i)),
                             DataWrite()->_test_pfespaces.at(i)->GetParMesh());
    bc_map.ApplyIntegratedBCs(test_var_name,
                              DataWrite()->_lfs.GetRef(test_var_name),
                              DataWrite()->_test_pfespaces.at(i)->GetParMesh());
  }
}

void
EquationSystemOperator::BuildLinearForms(platypus::BCMap & bc_map)
{
  // Register linear forms
  for (int i = 0; i < DataRead()->_test_var_names.size(); i++)
  {
    auto test_var_name = DataRead()->_test_var_names.at(i);
    DataWrite()->_lfs.Register(
        test_var_name, std::make_shared<mfem::ParLinearForm>(DataRead()->_test_pfespaces.at(i)));
    DataWrite()->_lfs.GetRef(test_var_name) = 0.0;
  }
  // Apply boundary conditions
  ApplyBoundaryConditions(bc_map);

  for (auto & test_var_name : DataRead()->_test_var_names)
  {
    // Apply kernels
    auto lf = DataRead()->_lfs.Get(test_var_name);
    if (DataRead()->_lf_kernels_map.Has(test_var_name))
    {
      auto lf_kernels = DataRead()->_lf_kernels_map.GetRef(test_var_name);

      for (auto & lf_kernel : lf_kernels)
      {
        lf->AddDomainIntegrator(lf_kernel->createIntegrator());
      }
    }
    lf->Assemble();
  }
}

void
EquationSystemOperator::BuildBilinearForms()
{
  // Register bilinear forms
  for (int i = 0; i < DataRead()->_test_var_names.size(); i++)
  {
    auto test_var_name = DataRead()->_test_var_names.at(i);
    DataWrite()->_blfs.Register(
        test_var_name, std::make_shared<mfem::ParBilinearForm>(DataRead()->_test_pfespaces.at(i)));

    // Apply kernels
    auto blf = DataRead()->_blfs.Get(test_var_name);
    if (DataRead()->_blf_kernels_map.Has(test_var_name))
    {
      auto blf_kernels = DataRead()->_blf_kernels_map.GetRef(test_var_name);

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
EquationSystemOperator::BuildMixedBilinearForms()
{
  // Register mixed bilinear forms. Note that not all combinations may
  // have a kernel

  // Create mblf for each test/trial pair
  for (int i = 0; i < DataRead()->_test_var_names.size(); i++)
  {
    auto test_var_name = DataRead()->_test_var_names.at(i);
    auto test_mblfs = std::make_shared<platypus::NamedFieldsMap<mfem::ParMixedBilinearForm>>();
    for (int j = 0; j < DataRead()->_test_var_names.size(); j++)
    {
      auto trial_var_name = DataRead()->_test_var_names.at(j);

      // Register MixedBilinearForm if kernels exist for it, and assemble
      // kernels
      if (DataRead()->_mblf_kernels_map_map.Has(test_var_name) &&
          DataRead()->_mblf_kernels_map_map.Get(test_var_name)->Has(trial_var_name))
      {
        auto mblf_kernels =
            DataRead()->_mblf_kernels_map_map.GetRef(test_var_name).GetRef(trial_var_name);
        auto mblf = std::make_shared<mfem::ParMixedBilinearForm>(DataRead()->_test_pfespaces.at(j),
                                                                 DataRead()->_test_pfespaces.at(i));
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
    DataWrite()->_mblfs.Register(test_var_name, test_mblfs);
  }
}

void
EquationSystemOperator::AddKernel(const std::string & test_var_name,
                                  std::shared_ptr<MFEMBilinearFormKernel> blf_kernel)
{
  AddTestVariableNameIfMissing(test_var_name);

  if (!DataRead()->_blf_kernels_map.Has(test_var_name))
  {
    // 1. Create kernels vector.
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMBilinearFormKernel>>>();

    // 2. Register with map to prevent leaks.
    DataWrite()->_blf_kernels_map.Register(test_var_name, std::move(kernels));
  }

  DataWrite()->_blf_kernels_map.GetRef(test_var_name).push_back(std::move(blf_kernel));
}

void
EquationSystemOperator::AddKernel(const std::string & test_var_name,
                                  std::shared_ptr<MFEMLinearFormKernel> lf_kernel)
{
  AddTestVariableNameIfMissing(test_var_name);

  if (!DataRead()->_lf_kernels_map.Has(test_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMLinearFormKernel>>>();

    DataWrite()->_lf_kernels_map.Register(test_var_name, std::move(kernels));
  }

  DataWrite()->_lf_kernels_map.GetRef(test_var_name).push_back(std::move(lf_kernel));
}

void
EquationSystemOperator::AddKernel(const std::string & test_var_name,
                                  std::shared_ptr<MFEMNonlinearFormKernel> nlf_kernel)
{
  AddTestVariableNameIfMissing(test_var_name);

  if (!DataRead()->_nlf_kernels_map.Has(test_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMNonlinearFormKernel>>>();

    DataWrite()->_nlf_kernels_map.Register(test_var_name, std::move(kernels));
  }

  DataWrite()->_nlf_kernels_map.GetRef(test_var_name).push_back(std::move(nlf_kernel));
}

void
EquationSystemOperator::AddKernel(const std::string & trial_var_name,
                                  const std::string & test_var_name,
                                  std::shared_ptr<MFEMMixedBilinearFormKernel> mblf_kernel)
{
  if (DataRead()->_assembly_level != mfem::AssemblyLevel::LEGACY)
    mooseError("Mixed Bilinear Form Kernels are currently only compatible with the LEGACY assembly "
               "level.");

  AddTestVariableNameIfMissing(test_var_name);

  // Register new mblf kernels map if not present for this test variable
  if (!DataRead()->_mblf_kernels_map_map.Has(test_var_name))
  {
    auto kernel_field_map = std::make_shared<
        platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMMixedBilinearFormKernel>>>>();

    DataWrite()->_mblf_kernels_map_map.Register(test_var_name, std::move(kernel_field_map));
  }

  // Register new mblf kernels map if not present for the test/trial variable
  // pair
  if (!DataRead()->_mblf_kernels_map_map.Get(test_var_name)->Has(trial_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMMixedBilinearFormKernel>>>();

    DataWrite()
        ->_mblf_kernels_map_map.Get(test_var_name)
        ->Register(trial_var_name, std::move(kernels));
  }

  DataWrite()
      ->_mblf_kernels_map_map.GetRef(test_var_name)
      .Get(trial_var_name)
      ->push_back(std::move(mblf_kernel));
}

void
EquationSystemOperator::Init(platypus::GridFunctions & gridfunctions,
                             const platypus::FESpaces & fespaces,
                             platypus::BCMap & bc_map,
                             platypus::Coefficients & coefficients,
                             mfem::AssemblyLevel assembly_level)
{
  for (auto & test_var_name : DataRead()->_test_var_names)
  {
    DataRead()->_assembly_level = assembly_level;
    if (!AssemblyIsSupported())
      mooseError("Assembly type chosen is not supported for this system");

    if (!gridfunctions.Has(test_var_name))
    {
      MFEM_ABORT("Test variable " << test_var_name
                                  << " requested by equation system during initialisation was "
                                     "not found in gridfunctions");
    }
    // Store pointers to variable FESpaces
    DataWrite()->_test_pfespaces.push_back(gridfunctions.Get(test_var_name)->ParFESpace());
    // Create auxiliary gridfunctions for applying Dirichlet conditions
    DataWrite()->_bc_gridfunc.emplace_back(
        std::make_unique<mfem::ParGridFunction>(gridfunctions.Get(test_var_name)->ParFESpace()));
  }
}

void
EquationSystemOperator::BuildJacobian(mfem::BlockVector & trueX, mfem::BlockVector & trueRHS)
{
  height = trueX.Size();
  width = trueRHS.Size();
  FormSystem(DataWrite()->_jacobian, trueX, trueRHS);
}

void
EquationSystemOperator::Mult(const mfem::Vector & x, mfem::Vector & residual) const
{
  DataRead()->_jacobian->Mult(x, residual);
}

mfem::Operator &
EquationSystemOperator::GetGradient(const mfem::Vector & u) const
{
  return *(DataRead()->_jacobian);
}

void
EquationSystemOperator::RecoverFEMSolution(mfem::BlockVector & trueX,
                                           platypus::GridFunctions & gridfunctions)
{
  for (int i = 0; i < DataRead()->_test_var_names.size(); i++)
  {
    auto & test_var_name = DataRead()->_test_var_names.at(i);
    trueX.GetBlock(i).SyncAliasMemory(trueX);
    gridfunctions.Get(test_var_name)->Distribute(&(trueX.GetBlock(i)));
  }
}

bool
EquationSystemOperator::AssemblyIsSupported()
{
  bool support = true;
  if (DataRead()->_matrix_type == EquationSystemData::MatrixType::DENSE &&
      !(DataRead()->_assembly_level == mfem::AssemblyLevel::LEGACY))
    support = false;

  return support;
}

void
EquationSystemOperator::FormSystem(mfem::OperatorHandle & op,
                                   mfem::BlockVector & trueX,
                                   mfem::BlockVector & trueRHS)
{
  switch (DataRead()->_matrix_type)
  {
    case EquationSystemData::MatrixType::DENSE:
      FormDenseSystem(op, trueX, trueRHS);
      break;

    case EquationSystemData::MatrixType::DIAGONAL:
      FormDiagonalSystem(op, trueX, trueRHS);
      break;
  }
}

void
EquationSystemOperator::FormDenseSystem(mfem::OperatorHandle & op,
                                        mfem::BlockVector & trueX,
                                        mfem::BlockVector & trueRHS)
{
  // Allocate block operator
  DataWrite()->_h_blocks.DeleteAll();
  DataWrite()->_h_blocks.SetSize(DataRead()->_test_var_names.size(),
                                 DataRead()->_test_var_names.size());
  // Form diagonal blocks.
  for (int i = 0; i < DataRead()->_test_var_names.size(); i++)
  {
    auto & test_var_name = DataRead()->_test_var_names.at(i);
    auto blf = DataRead()->_blfs.Get(test_var_name);
    auto lf = DataRead()->_lfs.Get(test_var_name);
    mfem::Vector aux_x, aux_rhs;
    mfem::HypreParMatrix aux_a;
    blf->FormLinearSystem(DataRead()->_ess_tdof_lists.at(i),
                          *(DataRead()->_bc_gridfunc.at(i)),
                          *lf,
                          aux_a,
                          aux_x,
                          aux_rhs);

    DataWrite()->_h_blocks(i, i) = new const mfem::HypreParMatrix(aux_a);
    trueX.GetBlock(i) = aux_x;
    trueRHS.GetBlock(i) = aux_rhs;
  }

  // Form off-diagonal blocks
  for (int i = 0; i < DataRead()->_test_var_names.size(); i++)
  {
    auto test_var_name = DataRead()->_test_var_names.at(i);
    for (int j = 0; j < DataRead()->_test_var_names.size(); j++)
    {
      auto trial_var_name = DataRead()->_test_var_names.at(j);

      mfem::Vector aux_x, aux_rhs;
      mfem::ParLinearForm aux_lf(DataRead()->_test_pfespaces.at(i));
      aux_lf = 0.0;
      if (DataRead()->_mblfs.Has(test_var_name) &&
          DataRead()->_mblfs.Get(test_var_name)->Has(trial_var_name))
      {
        auto mblf = DataRead()->_mblfs.Get(test_var_name)->Get(trial_var_name);
        mfem::HypreParMatrix aux_a;
        mblf->FormRectangularLinearSystem(DataRead()->_ess_tdof_lists.at(j),
                                          DataRead()->_ess_tdof_lists.at(i),
                                          *(DataRead()->_bc_gridfunc.at(j)),
                                          aux_lf,
                                          aux_a,
                                          aux_x,
                                          aux_rhs);
        trueRHS.GetBlock(i) += aux_rhs;
        DataWrite()->_h_blocks(i, j) = new const mfem::HypreParMatrix(aux_a);
      }
    }
  }
  // Sync memory
  for (int i = 0; i < DataRead()->_test_var_names.size(); i++)
  {
    trueX.GetBlock(0).SyncAliasMemory(trueX);
    trueRHS.GetBlock(0).SyncAliasMemory(trueRHS);
  }

  // Create monolithic matrix
  op.Reset(mfem::HypreParMatrixFromBlocks(DataRead()->_h_blocks));
}

void
EquationSystemOperator::FormDiagonalSystem(mfem::OperatorHandle & op,
                                           mfem::BlockVector & trueX,
                                           mfem::BlockVector & trueRHS)
{
  MFEM_ABORT("EquationSystemOperator::FormDiagonalSystem has not yet been implemented!");
}

void
EquationSystemOperator::BuildEquationSystem(platypus::BCMap & bc_map)
{
  BuildLinearForms(bc_map);
  BuildBilinearForms();
  if (DataRead()->_matrix_type == EquationSystemData::MatrixType::DENSE)
    BuildMixedBilinearForms();
}

} // namespace platypus