#include "equation_system_operator_base.h"

namespace platypus
{

// Base methods
bool
EquationSystemOperatorBase::VectorContainsName(const std::vector<std::string> & the_vector,
                                               const std::string & name) const
{

  auto iter = std::find(the_vector.begin(), the_vector.end(), name);

  return (iter != the_vector.end());
}

void
EquationSystemOperatorBase::AddTrialVariableNameIfMissing(const std::string & trial_var_name)
{
  if (!VectorContainsName(GetData()->_trial_var_names, trial_var_name))
  {
    GetData()->_trial_var_names.push_back(trial_var_name);
  }
}

void
EquationSystemOperatorBase::AddTestVariableNameIfMissing(const std::string & test_var_name)
{
  if (!VectorContainsName(GetData()->_test_var_names, test_var_name))
  {
    GetData()->_test_var_names.push_back(test_var_name);
  }
}

void
EquationSystemOperatorBase::ApplyBoundaryConditions(platypus::BCMap & bc_map)
{
  GetData()->_ess_tdof_lists.resize(GetData()->_test_var_names.size());
  for (int i = 0; i < GetData()->_test_var_names.size(); ++i)
  {
    auto test_var_name = GetData()->_test_var_names.at(i);
    // Set default value of gridfunction used in essential BC. Values
    // overwritten in applyEssentialBCs
    *(GetData()->_bc_gridfunc.at(i)) = 0.0;
    bc_map.ApplyEssentialBCs(test_var_name,
                             GetData()->_ess_tdof_lists.at(i),
                             *(GetData()->_bc_gridfunc.at(i)),
                             GetData()->_test_pfespaces.at(i)->GetParMesh());
    bc_map.ApplyIntegratedBCs(test_var_name,
                              GetData()->_lfs.GetRef(test_var_name),
                              GetData()->_test_pfespaces.at(i)->GetParMesh());
  }
}

void
EquationSystemOperatorBase::BuildLinearForms(platypus::BCMap & bc_map)
{
  // Register linear forms
  for (int i = 0; i < GetData()->_test_var_names.size(); ++i)
  {
    auto test_var_name = GetData()->_test_var_names.at(i);
    GetData()->_lfs.Register(
        test_var_name, std::make_shared<mfem::ParLinearForm>(GetData()->_test_pfespaces.at(i)));
    GetData()->_lfs.GetRef(test_var_name) = 0.0;
  }
  // Apply boundary conditions
  ApplyBoundaryConditions(bc_map);

  for (auto & test_var_name : GetData()->_test_var_names)
  {
    // Apply kernels
    auto lf = GetData()->_lfs.Get(test_var_name);
    if (GetData()->_lf_kernels_map.Has(test_var_name))
    {
      auto lf_kernels = GetData()->_lf_kernels_map.GetRef(test_var_name);

      for (auto & lf_kernel : lf_kernels)
      {
        lf->AddDomainIntegrator(lf_kernel->createIntegrator());
      }
    }
    lf->Assemble();
  }
}

void
EquationSystemOperatorBase::BuildBilinearForms()
{
  // Register bilinear forms
  for (int i = 0; i < GetData()->_test_var_names.size(); ++i)
  {
    auto test_var_name = GetData()->_test_var_names.at(i);
    GetData()->_blfs.Register(
        test_var_name, std::make_shared<mfem::ParBilinearForm>(GetData()->_test_pfespaces.at(i)));

    // Apply kernels
    auto blf = GetData()->_blfs.Get(test_var_name);
    if (GetData()->_blf_kernels_map.Has(test_var_name))
    {
      auto blf_kernels = GetData()->_blf_kernels_map.GetRef(test_var_name);

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
EquationSystemOperatorBase::BuildMixedBilinearForms()
{
  // Register mixed bilinear forms. Note that not all combinations may
  // have a kernel

  // Create mblf for each test/trial pair
  for (int i = 0; i < GetData()->_test_var_names.size(); ++i)
  {
    auto test_var_name = GetData()->_test_var_names.at(i);
    auto test_mblfs = std::make_shared<platypus::NamedFieldsMap<mfem::ParMixedBilinearForm>>();
    for (int j = 0; j < GetData()->_test_var_names.size(); ++j)
    {
      auto trial_var_name = GetData()->_test_var_names.at(j);

      // Register MixedBilinearForm if kernels exist for it, and assemble
      // kernels
      if (GetData()->_mblf_kernels_map_map.Has(test_var_name) &&
          GetData()->_mblf_kernels_map_map.Get(test_var_name)->Has(trial_var_name))
      {
        auto mblf_kernels =
            GetData()->_mblf_kernels_map_map.GetRef(test_var_name).GetRef(trial_var_name);
        auto mblf = std::make_shared<mfem::ParMixedBilinearForm>(GetData()->_test_pfespaces.at(j),
                                                                 GetData()->_test_pfespaces.at(i));
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
    GetData()->_mblfs.Register(test_var_name, test_mblfs);
  }
}

void
EquationSystemOperatorBase::AddKernel(const std::string & test_var_name,
                                      std::shared_ptr<MFEMBilinearFormKernel> blf_kernel)
{
  AddTestVariableNameIfMissing(test_var_name);

  if (!GetData()->_blf_kernels_map.Has(test_var_name))
  {
    // 1. Create kernels vector.
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMBilinearFormKernel>>>();

    // 2. Register with map to prevent leaks.
    GetData()->_blf_kernels_map.Register(test_var_name, std::move(kernels));
  }

  GetData()->_blf_kernels_map.GetRef(test_var_name).push_back(std::move(blf_kernel));
}

void
EquationSystemOperatorBase::AddKernel(const std::string & test_var_name,
                                      std::shared_ptr<MFEMLinearFormKernel> lf_kernel)
{
  AddTestVariableNameIfMissing(test_var_name);

  if (!GetData()->_lf_kernels_map.Has(test_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMLinearFormKernel>>>();

    GetData()->_lf_kernels_map.Register(test_var_name, std::move(kernels));
  }

  GetData()->_lf_kernels_map.GetRef(test_var_name).push_back(std::move(lf_kernel));
}

void
EquationSystemOperatorBase::AddKernel(const std::string & test_var_name,
                                      std::shared_ptr<MFEMNonlinearFormKernel> nlf_kernel)
{
  AddTestVariableNameIfMissing(test_var_name);

  if (!GetData()->_nlf_kernels_map.Has(test_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMNonlinearFormKernel>>>();

    GetData()->_nlf_kernels_map.Register(test_var_name, std::move(kernels));
  }

  GetData()->_nlf_kernels_map.GetRef(test_var_name).push_back(std::move(nlf_kernel));
}

void
EquationSystemOperatorBase::AddKernel(const std::string & trial_var_name,
                                      const std::string & test_var_name,
                                      std::shared_ptr<MFEMMixedBilinearFormKernel> mblf_kernel)
{
  if (GetData()->_assembly_level != mfem::AssemblyLevel::LEGACY)
    mooseError("Mixed Bilinear Form Kernels are currently only compatible with the LEGACY assembly "
               "level.");

  // Given we are adding a mixed blf kernel, our block matrix is no longer diagonal
  GetData()->_matrix_type = EquationSystemData::MatrixType::DENSE;

  AddTestVariableNameIfMissing(test_var_name);

  // Register new mblf kernels map if not present for this test variable
  if (!GetData()->_mblf_kernels_map_map.Has(test_var_name))
  {
    auto kernel_field_map = std::make_shared<
        platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMMixedBilinearFormKernel>>>>();

    GetData()->_mblf_kernels_map_map.Register(test_var_name, std::move(kernel_field_map));
  }

  // Register new mblf kernels map if not present for the test/trial variable
  // pair
  if (!GetData()->_mblf_kernels_map_map.Get(test_var_name)->Has(trial_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMMixedBilinearFormKernel>>>();

    GetData()
        ->_mblf_kernels_map_map.Get(test_var_name)
        ->Register(trial_var_name, std::move(kernels));
  }

  GetData()
      ->_mblf_kernels_map_map.GetRef(test_var_name)
      .Get(trial_var_name)
      ->push_back(std::move(mblf_kernel));
}

void
EquationSystemOperatorBase::Init(platypus::GridFunctions & gridfunctions,
                                 const platypus::FESpaces & fespaces,
                                 platypus::BCMap & bc_map,
                                 platypus::Coefficients & coefficients,
                                 mfem::AssemblyLevel assembly_level)
{

  // At first we assume the matrix is diagonal. If we later add a mixed blf kernel, @_matrix_type is
  // set to dense.
  GetData()->_matrix_type = EquationSystemData::MatrixType::DIAGONAL;
  GetData()->_assembly_level = assembly_level;

  for (auto & test_var_name : GetData()->_test_var_names)
  {

    if (!gridfunctions.Has(test_var_name))
    {
      MFEM_ABORT("Test variable " << test_var_name
                                  << " requested by equation system during initialisation was "
                                     "not found in gridfunctions");
    }
    // Store pointers to variable FESpaces
    GetData()->_test_pfespaces.push_back(gridfunctions.Get(test_var_name)->ParFESpace());
    // Create auxiliary gridfunctions for applying Dirichlet conditions
    GetData()->_bc_gridfunc.emplace_back(
        std::make_unique<mfem::ParGridFunction>(gridfunctions.Get(test_var_name)->ParFESpace()));
  }
}

void
EquationSystemOperatorBase::BuildJacobian(mfem::BlockVector & trueX, mfem::BlockVector & trueRHS)
{
  height = trueX.Size();
  width = trueRHS.Size();
  FormSystem(GetData()->_jacobian, trueX, trueRHS);
}

void
EquationSystemOperatorBase::Mult(const mfem::Vector & x, mfem::Vector & residual) const
{
  GetData()->_jacobian->Mult(x, residual);
}

mfem::Operator &
EquationSystemOperatorBase::GetGradient(const mfem::Vector & u) const
{
  return *(GetData()->_jacobian);
}

void
EquationSystemOperatorBase::RecoverFEMSolution(mfem::BlockVector & trueX,
                                               platypus::GridFunctions & gridfunctions)
{
  for (int i = 0; i < GetData()->_test_var_names.size(); ++i)
  {
    auto & test_var_name = GetData()->_test_var_names.at(i);
    trueX.GetBlock(i).SyncAliasMemory(trueX);
    gridfunctions.Get(test_var_name)->Distribute(&(trueX.GetBlock(i)));
  }
}

void
EquationSystemOperatorBase::FormSystem(mfem::OperatorHandle & op,
                                       mfem::BlockVector & trueX,
                                       mfem::BlockVector & trueRHS)
{
  GetData()->_matrix_type = EquationSystemData::MatrixType::DENSE;

  switch (GetData()->_matrix_type)
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
EquationSystemOperatorBase::FormDenseSystem(mfem::OperatorHandle & op,
                                            mfem::BlockVector & trueX,
                                            mfem::BlockVector & trueRHS)
{
  // Allocate block operator
  GetData()->_h_blocks.DeleteAll();
  GetData()->_h_blocks.SetSize(GetData()->_test_var_names.size(),
                               GetData()->_test_var_names.size());
  // Form diagonal blocks.
  for (int i = 0; i < GetData()->_test_var_names.size(); ++i)
  {
    auto & test_var_name = GetData()->_test_var_names.at(i);
    auto blf = GetData()->_blfs.Get(test_var_name);
    auto lf = GetData()->_lfs.Get(test_var_name);
    mfem::Vector aux_x, aux_rhs;
    mfem::HypreParMatrix aux_a;
    blf->FormLinearSystem(GetData()->_ess_tdof_lists.at(i),
                          *(GetData()->_bc_gridfunc.at(i)),
                          *lf,
                          aux_a,
                          aux_x,
                          aux_rhs);

    GetData()->_h_blocks(i, i) = new const mfem::HypreParMatrix(aux_a);
    trueX.GetBlock(i) = aux_x;
    trueRHS.GetBlock(i) = aux_rhs;
  }

  // Form off-diagonal blocks
  for (int i = 0; i < GetData()->_test_var_names.size(); ++i)
  {
    auto test_var_name = GetData()->_test_var_names.at(i);
    for (int j = 0; j < GetData()->_test_var_names.size(); ++j)
    {
      auto trial_var_name = GetData()->_test_var_names.at(j);

      mfem::Vector aux_x, aux_rhs;
      mfem::ParLinearForm aux_lf(GetData()->_test_pfespaces.at(i));
      aux_lf = 0.0;
      if (GetData()->_mblfs.Has(test_var_name) &&
          GetData()->_mblfs.Get(test_var_name)->Has(trial_var_name))
      {
        auto mblf = GetData()->_mblfs.Get(test_var_name)->Get(trial_var_name);
        mfem::HypreParMatrix aux_a;
        mblf->FormRectangularLinearSystem(GetData()->_ess_tdof_lists.at(j),
                                          GetData()->_ess_tdof_lists.at(i),
                                          *(GetData()->_bc_gridfunc.at(j)),
                                          aux_lf,
                                          aux_a,
                                          aux_x,
                                          aux_rhs);
        trueRHS.GetBlock(i) += aux_rhs;
        GetData()->_h_blocks(i, j) = new const mfem::HypreParMatrix(aux_a);
      }
    }
  }
  // Sync memory
  for (int i = 0; i < GetData()->_test_var_names.size(); ++i)
  {
    trueX.GetBlock(0).SyncAliasMemory(trueX);
    trueRHS.GetBlock(0).SyncAliasMemory(trueRHS);
  }

  // Create monolithic matrix
  op.Reset(mfem::HypreParMatrixFromBlocks(GetData()->_h_blocks));
}

void
EquationSystemOperatorBase::FormDiagonalSystem(mfem::OperatorHandle & op,
                                               mfem::BlockVector & trueX,
                                               mfem::BlockVector & trueRHS)
{

  MakeBlockOperator();

  for (int i = 0; i < GetData()->_test_var_names.size(); ++i)
  {
    auto & test_var_name = GetData()->_test_var_names.at(i);

    auto blf = GetData()->_blfs.Get(test_var_name);
    auto lf = GetData()->_lfs.Get(test_var_name);
    blf->SetAssemblyLevel(GetData()->_assembly_level);

    mfem::Vector aux_x, aux_rhs;
    mfem::OperatorPtr * aux_a;
    blf->FormLinearSystem(GetData()->_ess_tdof_lists.at(i),
                          *(GetData()->_bc_gridfunc.at(i)),
                          *lf,
                          *aux_a,
                          aux_x,
                          aux_rhs);

    GetData()->_h_block_op->SetBlock(i, i, aux_a->Ptr());
  }

  op.Reset(GetData()->_h_block_op.get());
}

void
EquationSystemOperatorBase::MakeBlockOperator()
{
  mfem::Array<int> block_offsets(GetData()->_test_var_names.size() + 1);

  for (int i = 0; i < GetData()->_test_var_names.size(); ++i)
    block_offsets[i + 1] =
        GetData()->_blfs.Get(GetData()->_test_var_names.at(i))->FESpace()->GetTrueVSize();

  block_offsets.PartialSum();
  GetData()->_h_block_op = std::make_shared<mfem::BlockOperator>(block_offsets);
}

void
EquationSystemOperatorBase::BuildEquationSystem(platypus::BCMap & bc_map)
{
  BuildLinearForms(bc_map);
  BuildBilinearForms();
  if (GetData()->_matrix_type == EquationSystemData::MatrixType::DENSE)
    BuildMixedBilinearForms();
}

} // namespace platypus