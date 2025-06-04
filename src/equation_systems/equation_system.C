#include "equation_system.h"
#include "/home/kchockalingam/projects/mfem/general/forall.hpp"
namespace platypus
{

EquationSystem::~EquationSystem() { _h_blocks.DeleteAll(); }

bool
EquationSystem::VectorContainsName(const std::vector<std::string> & the_vector,
                                   const std::string & name) const
{

  auto iter = std::find(the_vector.begin(), the_vector.end(), name);

  return (iter != the_vector.end());
}

void
EquationSystem::AddTrialVariableNameIfMissing(const std::string & trial_var_name)
{
  if (!VectorContainsName(_trial_var_names, trial_var_name))
  {
    _trial_var_names.push_back(trial_var_name);
  }
}

void
EquationSystem::AddTestVariableNameIfMissing(const std::string & test_var_name)
{
  if (!VectorContainsName(_test_var_names, test_var_name))
  {
    _test_var_names.push_back(test_var_name);
  }
}

void
EquationSystem::AddKernel(std::shared_ptr<MFEMKernel> kernel)
{
  AddTestVariableNameIfMissing(kernel->getTestVariableName());
  AddTrialVariableNameIfMissing(kernel->getTrialVariableName());
  auto trial_var_name = kernel->getTrialVariableName();
  auto test_var_name = kernel->getTestVariableName();
  if (!_kernels_map.Has(test_var_name))
  {
    auto kernel_field_map =
        std::make_shared<platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>>();
    _kernels_map.Register(test_var_name, std::move(kernel_field_map));
  }
  // Register new kernels map if not present for the test/trial variable
  // pair
  if (!_kernels_map.Get(test_var_name)->Has(trial_var_name))
  {
    auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMKernel>>>();
    _kernels_map.Get(test_var_name)->Register(trial_var_name, std::move(kernels));
  }
  _kernels_map.GetRef(test_var_name).Get(trial_var_name)->push_back(std::move(kernel));
}

void
EquationSystem::AddIntegratedBC(std::shared_ptr<MFEMIntegratedBC> bc)
{
  AddTestVariableNameIfMissing(bc->getTestVariableName());
  AddTrialVariableNameIfMissing(bc->getTrialVariableName());
  auto trial_var_name = bc->getTrialVariableName();
  auto test_var_name = bc->getTestVariableName();
  if (!_integrated_bc_map.Has(test_var_name))
  {
    auto integrated_bc_field_map = std::make_shared<
        platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMIntegratedBC>>>>();
    _integrated_bc_map.Register(test_var_name, std::move(integrated_bc_field_map));
  }
  // Register new integrated bc map if not present for the test/trial variable
  // pair
  if (!_integrated_bc_map.Get(test_var_name)->Has(trial_var_name))
  {
    auto bcs = std::make_shared<std::vector<std::shared_ptr<MFEMIntegratedBC>>>();
    _integrated_bc_map.Get(test_var_name)->Register(trial_var_name, std::move(bcs));
  }
  _integrated_bc_map.GetRef(test_var_name).Get(trial_var_name)->push_back(std::move(bc));
}

void
EquationSystem::AddEssentialBC(std::shared_ptr<MFEMEssentialBC> bc)
{
  AddTestVariableNameIfMissing(bc->getTestVariableName());
  auto test_var_name = bc->getTestVariableName();
  if (!_essential_bc_map.Has(test_var_name))
  {
    auto bcs = std::make_shared<std::vector<std::shared_ptr<MFEMEssentialBC>>>();
    _essential_bc_map.Register(test_var_name, std::move(bcs));
  }
  _essential_bc_map.GetRef(test_var_name).push_back(std::move(bc));
}

void
EquationSystem::ApplyEssentialBCs()
{
  _ess_tdof_lists.resize(_test_var_names.size());
  for (int i = 0; i < _test_var_names.size(); i++)
  {
    auto test_var_name = _test_var_names.at(i);
    // Set default value of gridfunction used in essential BC. Values
    // overwritten in applyEssentialBCs
    mfem::ParGridFunction & trial_gf(*(_xs.at(i)));
    mfem::ParGridFunction & trial_gf_time_derivatives(*(_dxdts.at(i)));
    mfem::ParMesh * pmesh(_test_pfespaces.at(i)->GetParMesh());
    trial_gf = 0.0;
    trial_gf_time_derivatives = 0.0;

    auto bcs = _essential_bc_map.GetRef(test_var_name);
    mfem::Array<int> global_ess_markers(pmesh->bdr_attributes.Max());
    global_ess_markers = 0;
    for (auto & bc : bcs)
    {
      bc->ApplyBC(trial_gf, pmesh);

      mfem::Array<int> ess_bdrs(bc->getBoundaries());
      for (auto it = 0; it != pmesh->bdr_attributes.Max(); ++it)
      {
        global_ess_markers[it] = std::max(global_ess_markers[it], ess_bdrs[it]);
      }
    }
    trial_gf.FESpace()->GetEssentialTrueDofs(global_ess_markers, _ess_tdof_lists.at(i));
   // std::cout << "******** Variable name ***********  " << test_var_name << "  " << i << " **********************" << std::endl;
   // _ess_tdof_lists.at(i).Print(std::cout);

  }
}

void
EquationSystem::FormLinearSystem(mfem::OperatorHandle & op,
                                 mfem::BlockVector & trueX,
                                 mfem::BlockVector & trueRHS) const
{

  switch (_assembly_level)
  {
    case mfem::AssemblyLevel::LEGACY:
      FormLegacySystem(op, trueX, trueRHS);
      break;
    default:
      MFEM_VERIFY(_test_var_names.size() == 1,
                  "Non-legacy assembly is only supported for single-variable systems");
      MFEM_VERIFY(_test_var_names.size() == _trial_var_names.size(),
                  "Non-legacy assembly is only supported for square systems");
      FormSystem(op, trueX, trueRHS);
  }
}

void
EquationSystem::FormSystem(mfem::OperatorHandle & op,
                           mfem::BlockVector & trueX,
                           mfem::BlockVector & trueRHS) const
{
  auto & test_var_name = _test_var_names.at(0);
  auto blf = _blfs.Get(test_var_name);
  auto lf = _lfs.Get(test_var_name);
  mfem::BlockVector aux_x, aux_rhs;
  mfem::OperatorPtr * aux_a = new mfem::OperatorPtr;

  blf->FormLinearSystem(_ess_tdof_lists.at(0), *(_xs.at(0)), *lf, *aux_a, aux_x, aux_rhs);

  trueX.GetBlock(0) = aux_x;
  trueRHS.GetBlock(0) = aux_rhs;
  trueX.SyncFromBlocks();
  trueRHS.SyncFromBlocks();

  _trueRHS.SetSize(width);
  _trueRHS = trueRHS;

  op.Reset(aux_a->Ptr());
}

void
EquationSystem::FormLegacySystem(mfem::OperatorHandle & op,
                                 mfem::BlockVector & trueX,
                                 mfem::BlockVector & trueRHS) const
{

  // Allocate block operator
  _h_blocks.DeleteAll();
  _h_blocks.SetSize(_test_var_names.size(), _test_var_names.size());
  std::cout << "********** I AM HERE ******************************************" << std::endl;
  // Form diagonal blocks.
  for (int i = 0; i < _test_var_names.size(); i++)
  {
    auto & test_var_name = _test_var_names.at(i);
    auto blf = _blfs.Get(test_var_name);
    auto lf = _lfs.Get(test_var_name);
    mfem::Vector aux_x, aux_rhs;
    mfem::HypreParMatrix * aux_a = new mfem::HypreParMatrix;
    // Ownership of aux_a goes to the blf
    blf->FormLinearSystem(_ess_tdof_lists.at(i), *(_xs.at(i)), *lf, *aux_a, aux_x, aux_rhs);
    _h_blocks(i, i) = aux_a;
    trueX.GetBlock(i) = aux_x;
    trueRHS.GetBlock(i) = aux_rhs;
  }

  // Form off-diagonal blocks
  for (int i = 0; i < _test_var_names.size(); i++)
  {
    auto test_var_name = _test_var_names.at(i);
    for (int j = 0; j < _test_var_names.size(); j++)
    {
      auto trial_var_name = _test_var_names.at(j);

      mfem::Vector aux_x, aux_rhs;
      mfem::ParLinearForm aux_lf(_test_pfespaces.at(i));
      aux_lf = 0.0;
      if (_mblfs.Has(test_var_name) && _mblfs.Get(test_var_name)->Has(trial_var_name))
      {
        auto mblf = _mblfs.Get(test_var_name)->Get(trial_var_name);
        mfem::HypreParMatrix * aux_a = new mfem::HypreParMatrix;
        // Ownership of aux_a goes to the blf
        mblf->FormRectangularLinearSystem(_ess_tdof_lists.at(j),
                                          _ess_tdof_lists.at(i),
                                          *(_xs.at(j)),
                                          aux_lf,
                                          *aux_a,
                                          aux_x,
                                          aux_rhs);
        _h_blocks(i, j) = aux_a;
        trueRHS.GetBlock(i) += aux_rhs;
      }
    }
  }
  // Sync memory
  for (int i = 0; i < _test_var_names.size(); i++)
  {
    trueX.GetBlock(0).SyncAliasMemory(trueX);
    trueRHS.GetBlock(0).SyncAliasMemory(trueRHS);
  }
  _trueRHS.SetSize(width);
  _trueRHS = trueRHS;
  // Create monolithic matrix
  op.Reset(mfem::HypreParMatrixFromBlocks(_h_blocks));
}

void
EquationSystem::BuildJacobian(mfem::BlockVector & trueX, mfem::BlockVector & trueRHS)
{
  height = trueX.Size();
  width = trueRHS.Size();
  FormLinearSystem(_jacobian, trueX, trueRHS);
}



void
EquationSystem::UpdateJacobian() const
{
  _h_blocks.DeleteAll();
  _h_blocks.SetSize(_test_var_names.size(), _test_var_names.size());

  for (int i = 0; i < _test_var_names.size(); i++)
    {
      auto & test_var_name = _test_var_names.at(i);
      auto blf = _blfs.Get(test_var_name);
      mfem::OperatorHandle aux_a;
      blf->Update();
      blf->Assemble();
      blf->FormSystemMatrix(_ess_tdof_lists.at(i), aux_a);
      _h_blocks(i, i) = static_cast<mfem::HypreParMatrix*>(aux_a.Ptr());
    }

    // Form off-diagonal blocks
    for (int i = 0; i < _test_var_names.size(); i++)
    {
      auto test_var_name = _test_var_names.at(i);
      for (int j = 0; j < _test_var_names.size(); j++)
      {
        auto trial_var_name = _test_var_names.at(j);
        if (_mblfs.Has(test_var_name) && _mblfs.Get(test_var_name)->Has(trial_var_name))
        {
          auto mblf = _mblfs.Get(test_var_name)->Get(trial_var_name);
          mfem::OperatorHandle aux_a;
          mblf->Update();
          mblf->Assemble();
          mblf->FormRectangularSystemMatrix(empty_tdof,_ess_tdof_lists.at(j), aux_a);
          _h_blocks(i, j) = static_cast<mfem::HypreParMatrix*>(aux_a.Ptr());
        }
      }
    }

    _jacobian.Reset(mfem::HypreParMatrixFromBlocks(_h_blocks));
}

void CopyVec(const mfem::Vector & x, mfem::Vector & y){ y = x;}

void applyDirchValues(const mfem::Vector &k, mfem::Vector &y, mfem::Array<int> dofs)
{
  if(dofs.Size() > 0){ //Only apply if there are constrained DOF's

    const bool use_dev = dofs.UseDevice() || k.UseDevice() || y.UseDevice();
    const int n = dofs.Size();

    // Use read+write access for X - we only modify some of its entries

    auto d_X = y.ReadWrite(use_dev);
    auto d_y = k.Read(use_dev);
    auto d_dofs = dofs.Read(use_dev);

    mfem::forall_switch(use_dev, n, [=] MFEM_HOST_DEVICE (int i)
    {
      const int dof_i = d_dofs[i];

      if (dof_i >= 0)   d_X[dof_i]    =  d_y[dof_i];
      if (!(dof_i >= 0))d_X[-1-dof_i] = -d_y[-1-dof_i];
     });
  }
};
 

void
EquationSystem::Mult(const mfem::Vector & x, mfem::Vector & residual) const
{
  x.HostRead();
  CopyVec(x,_trueBlockX);

  _trueBlockRHS = 0.0;
 
  for (int i = 0; i < _trial_var_names.size(); i++)
    {
      auto & trial_var_name = _trial_var_names.at(i);
      applyDirchValues(*(_xs.at(i)), _trueBlockX.GetBlock(i), _ess_tdof_lists.at(i));
      _gfuncs->Get(trial_var_name)->Distribute(&(_trueBlockX.GetBlock(i)));
    }

  for (int i = 0; i < _test_var_names.size(); i++)
    {
      auto & test_var_name = _test_var_names.at(i);
      auto lf = _lfs.GetShared(test_var_name);
      lf->Assemble();
     // lf->ParallelAssemble(_trueBlockRHS.GetBlock(i));
     // _trueBlockRHS.GetBlock(i).SetSubVector(_ess_tdof_lists.at(i) , 0.00);
    }

  //UpdateJacobian();
   FormLinearSystem(_jacobian, _trueBlockX, _trueBlockRHS);
  _jacobian->Mult(_trueBlockX, residual);
  residual.HostRead();
  residual -= _trueBlockRHS;

  //residual -= _trueRHS;
}

void
TimeDependentEquationSystem::update_old_state()
{
    // Update solution values on Dirichlet values to be in terms of du/dt instead of u
  for (int i = 0; i < _test_var_names.size(); i++)
    {
     auto & test_var_name = _test_var_names.at(i);
     //Below is a horrible hack - this has to be before the mult is called
     CopyVec(*_trial_variables.Get(test_var_name), _trueBlockX_Old.GetBlock(i));
    }
}

void
TimeDependentEquationSystem::Mult(const mfem::Vector & truedXdt, mfem::Vector & residual) const
{
  truedXdt.HostRead();
  CopyVec(truedXdt, _trueBlockdXdt);

  _trueBlockRHS = 0.0;

  for (int i = 0; i < _trial_var_names.size(); i++)
    {
      auto & trial_var_name = _trial_var_names.at(i);
      std::cout << " trial_var_name  " << trial_var_name << std::endl;
      //_trueBlockdXdt.GetBlock(i).SetSubVector(_ess_tdof_lists.at(i) , 0.00);
      _gfuncs->Get(trial_var_name)->Distribute(&(_trueBlockdXdt.GetBlock(i)));
    }
    
    for (int i = 0; i < _test_var_names.size(); i++)
    {
      auto & test_var_name = _test_var_names.at(i);
      applyDirchValues(*(_xs.at(i)), *(_gfuncs->Get(test_var_name)),_ess_tdof_lists.at(i));
    }

  for (int i = 0; i < _test_var_names.size(); i++)
    {
      auto & test_var_name = _test_var_names.at(i);
      auto lf = _lfs.GetShared(test_var_name);
      lf->Assemble();
      lf->ParallelAssemble(_trueBlockRHS.GetBlock(i));
    }

  //UpdateJacobian();
   FormLinearSystem(_jacobian, _trueBlockdXdt, _trueBlockRHS);
   
  _jacobian->Mult(_trueBlockdXdt, residual);
  residual.HostRead();
  residual -= _trueBlockRHS;
  
  //residual -= _trueRHS;
}

mfem::Operator &
EquationSystem::GetGradient(const mfem::Vector & u) const
{
  return *_jacobian;
}

void
EquationSystem::RecoverFEMSolution(mfem::BlockVector & trueX,
                                   platypus::GridFunctions & gridfunctions)
{
  for (int i = 0; i < _trial_var_names.size(); i++)
  {
    auto & trial_var_name = _trial_var_names.at(i);
    trueX.GetBlock(i).SyncAliasMemory(trueX);
    gridfunctions.Get(trial_var_name)->Distribute(&(trueX.GetBlock(i)));
  }
}

void
EquationSystem::Init(platypus::GridFunctions & gridfunctions,
                     const platypus::FESpaces & fespaces,
                     mfem::AssemblyLevel assembly_level)
{
  _assembly_level = assembly_level;

  for (auto & test_var_name : _test_var_names)
  {
    if (!gridfunctions.Has(test_var_name))
    {
      MFEM_ABORT("Test variable " << test_var_name
                                  << " requested by equation system during initialisation was "
                                     "not found in gridfunctions");
    }
    // Store pointers to variable FESpaces
    _test_pfespaces.push_back(gridfunctions.Get(test_var_name)->ParFESpace());
    // Create auxiliary gridfunctions for applying Dirichlet conditions
    _xs.emplace_back(
        std::make_unique<mfem::ParGridFunction>(gridfunctions.Get(test_var_name)->ParFESpace()));
    _dxdts.emplace_back(
        std::make_unique<mfem::ParGridFunction>(gridfunctions.Get(test_var_name)->ParFESpace()));
    _trial_variables.Register(test_var_name, gridfunctions.GetShared(test_var_name));
  }
}

void
EquationSystem::BuildLinearForms()
{
  // Register linear forms
  for (int i = 0; i < _test_var_names.size(); i++)
  {
    auto test_var_name = _test_var_names.at(i);
    _lfs.Register(test_var_name, std::make_shared<mfem::ParLinearForm>(_test_pfespaces.at(i)));
    _lfs.GetRef(test_var_name) = 0.0;
  }
  // Apply boundary conditions
  ApplyEssentialBCs();

  for (auto & test_var_name : _test_var_names)
  {
    // Apply kernels
    auto lf = _lfs.GetShared(test_var_name);
    ApplyDomainLFIntegrators(test_var_name, lf, _kernels_map);
    ApplyBoundaryLFIntegrators(test_var_name, lf, _integrated_bc_map);
    lf->Assemble(); //<-------  Don't need it????
  }
}

void
EquationSystem::BuildBilinearForms()
{
  // Register bilinear forms
  for (int i = 0; i < _test_var_names.size(); i++)
  {
    auto test_var_name = _test_var_names.at(i);
    _blfs.Register(test_var_name, std::make_shared<mfem::ParBilinearForm>(_test_pfespaces.at(i)));

    // Apply kernels
    auto blf = _blfs.GetShared(test_var_name);
    blf->SetAssemblyLevel(_assembly_level);
    ApplyBoundaryBLFIntegrators<mfem::ParBilinearForm>(
        test_var_name, test_var_name, blf, _integrated_bc_map);
    ApplyDomainBLFIntegrators<mfem::ParBilinearForm>(
        test_var_name, test_var_name, blf, _kernels_map);
    // Assemble
    blf->Assemble(); //<-------  Don't need it????
  }
}

void
EquationSystem::BuildMixedBilinearForms()
{
  // Register mixed bilinear forms. Note that not all combinations may
  // have a kernel

  // Create mblf for each test/trial pair
  for (int i = 0; i < _test_var_names.size(); i++)
  {
    auto test_var_name = _test_var_names.at(i);
    auto test_mblfs = std::make_shared<platypus::NamedFieldsMap<mfem::ParMixedBilinearForm>>();
    for (int j = 0; j < _test_var_names.size(); j++)
    {
      auto trial_var_name = _test_var_names.at(j);
      auto mblf = std::make_shared<mfem::ParMixedBilinearForm>(_test_pfespaces.at(j),
                                                               _test_pfespaces.at(i));
      // Register MixedBilinearForm if kernels exist for it, and assemble
      // kernels
      if (_kernels_map.Has(test_var_name) && _kernels_map.Get(test_var_name)->Has(trial_var_name) &&
          test_var_name != trial_var_name)
      {
        mblf->SetAssemblyLevel(_assembly_level);
        // Apply all mixed kernels with this test/trial pair
        ApplyDomainBLFIntegrators<mfem::ParMixedBilinearForm>(
            trial_var_name, test_var_name, mblf, _kernels_map);
        // Assemble mixed bilinear forms
        mblf->Assemble(); //<-------  Don't need it????
        // Register mixed bilinear forms associated with a single trial variable
        // for the current test variable
        test_mblfs->Register(trial_var_name, mblf);
      }
    }
    // Register all mixed bilinear form sets associated with a single test
    // variable
    _mblfs.Register(test_var_name, test_mblfs);
  }
}

void
EquationSystem::BuildEquationSystem(platypus::GridFunctions & gridfunctions, mfem::Array<int> & btoffsets)
{
  _gfuncs = &gridfunctions;
  _block_true_offsets = &btoffsets;
  _trueBlockX.Update(*_block_true_offsets);
  _trueBlockRHS.Update(*_block_true_offsets);
  _trueBlockdXdt.Update(*_block_true_offsets);
  _trueBlockX_Old.Update(*_block_true_offsets);
  _h_blocks.DeleteAll();
  _h_blocks.SetSize(_test_var_names.size(), _test_var_names.size());
  BuildBilinearForms();
  BuildMixedBilinearForms();
  BuildLinearForms();
}

TimeDependentEquationSystem::TimeDependentEquationSystem() : _dt_coef(1.0) {}

void
TimeDependentEquationSystem::AddTrialVariableNameIfMissing(const std::string & var_name)
{
  // The TimeDependentEquationSystem operator expects to act on a vector of variable time
  // derivatives
  std::string var_time_derivative_name = GetTimeDerivativeName(var_name);

  if (!VectorContainsName(_trial_var_names, var_time_derivative_name))
  {
    _trial_var_names.push_back(var_time_derivative_name);
    _trial_var_time_derivative_names.push_back(var_time_derivative_name);
  }
}

void
TimeDependentEquationSystem::SetTimeStep(double dt)
{
  if (fabs(dt - _dt_coef.constant) > 1.0e-12 * dt)
  {
    _dt_coef.constant = dt;
    for (auto test_var_name : _test_var_names)
    {
      auto blf = _blfs.Get(test_var_name);
      blf->Update();
      blf->Assemble();
    }
  }
}

void
TimeDependentEquationSystem::AddKernel(std::shared_ptr<MFEMKernel> kernel)
{
  if (kernel->getTrialVariableName() == GetTimeDerivativeName(kernel->getTestVariableName()))
  {
    auto trial_var_name = kernel->getTrialVariableName();
    auto test_var_name = kernel->getTestVariableName();
    AddTestVariableNameIfMissing(test_var_name);
    AddTrialVariableNameIfMissing(test_var_name);
    if (!_td_kernels_map.Has(test_var_name))
    {
      auto kernel_field_map =
          std::make_shared<platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>>();
      _td_kernels_map.Register(test_var_name, std::move(kernel_field_map));
    }
    // Register new kernels map if not present for the test variable
    if (!_td_kernels_map.Get(test_var_name)->Has(test_var_name))
    {
      auto kernels = std::make_shared<std::vector<std::shared_ptr<MFEMKernel>>>();
      _td_kernels_map.Get(test_var_name)->Register(test_var_name, std::move(kernels));
    }
    _td_kernels_map.GetRef(test_var_name).Get(test_var_name)->push_back(std::move(kernel));
  }
  else
  {
    EquationSystem::AddKernel(kernel);
  }
}

void
TimeDependentEquationSystem::BuildBilinearForms()
{
  EquationSystem::BuildBilinearForms();

  // Build and assemble bilinear forms acting on time derivatives
  for (int i = 0; i < _test_var_names.size(); i++)
  {
    auto test_var_name = _test_var_names.at(i);

    _td_blfs.Register(test_var_name,
                      std::make_shared<mfem::ParBilinearForm>(_test_pfespaces.at(i)));

    // Apply kernels to td_blf
    auto td_blf = _td_blfs.GetShared(test_var_name);
    td_blf->SetAssemblyLevel(_assembly_level);
    ApplyBoundaryBLFIntegrators<mfem::ParBilinearForm>(
        test_var_name, test_var_name, td_blf, _integrated_bc_map);
    ApplyDomainBLFIntegrators<mfem::ParBilinearForm>(
        test_var_name, test_var_name, td_blf, _td_kernels_map);

    // Recover and scale integrators from blf. This is to apply the dt*du/dt contributions from the
    // operator on the trial variable in the implicit integration scheme
    auto blf = _blfs.Get(test_var_name);
    auto integs = blf->GetDBFI();
    auto b_integs = blf->GetBBFI();
    auto markers = blf->GetBBFI_Marker();

    mfem::SumIntegrator * sum = new mfem::SumIntegrator;
    ScaleIntegrator * scaled_sum = new ScaleIntegrator(sum, _dt_coef.constant, false);

    for (int i = 0; i < integs->Size(); ++i)
    {
      sum->AddIntegrator(*integs[i]);
    }

    for (int i = 0; i < b_integs->Size(); ++i)
    {
      td_blf->AddBoundaryIntegrator(new ScaleIntegrator(*b_integs[i], _dt_coef.constant, false),
                                    *(*markers[i]));
    }

    // scaled_sum is owned by td_blf
    td_blf->AddDomainIntegrator(scaled_sum);

    // Assemble form
    td_blf->Assemble();
  }
}

void
TimeDependentEquationSystem::FormLegacySystem(mfem::OperatorHandle & op,
                                              mfem::BlockVector & truedXdt,
                                              mfem::BlockVector & trueRHS) const
{

  // Allocate block operator
  _h_blocks.DeleteAll();
  _h_blocks.SetSize(_test_var_names.size(), _test_var_names.size());
  // Form diagonal blocks.
    std::cout << "********** TIME: I AM HERE ******************************************" << std::endl;
  for (int i = 0; i < _test_var_names.size(); i++)
  {
    auto & test_var_name = _test_var_names.at(i);
    auto td_blf = _td_blfs.Get(test_var_name);
    auto blf = _blfs.Get(test_var_name);
    auto lf = _lfs.Get(test_var_name);
    // if implicit, add contribution to linear form from terms involving state
    // variable at previous timestep: {
    // blf->AddMult(*_trial_variables.Get(test_var_name), *lf, -1.0);
     blf->AddMult(_trueBlockX_Old.GetBlock(i), *lf, -1.0);
    // }
    mfem::Vector aux_x, aux_rhs;
    // Update solution values on Dirichlet values to be in terms of du/dt instead of u
    mfem::Vector bc_x = *(_xs.at(i).get());
    //bc_x -= *_trial_variables.Get(test_var_name);
    bc_x -= _trueBlockX_Old.GetBlock(i);
    bc_x /= _dt_coef.constant;

    // Form linear system for operator acting on vector of du/dt
    mfem::HypreParMatrix * aux_a = new mfem::HypreParMatrix;
    // Ownership of aux_a goes to the blf
    td_blf->FormLinearSystem(_ess_tdof_lists.at(i), bc_x, *lf, *aux_a, aux_x, aux_rhs);
    _h_blocks(i, i) = aux_a;
    truedXdt.GetBlock(i) = aux_x;
    trueRHS.GetBlock(i) = aux_rhs;
  }

  truedXdt.SyncFromBlocks();
  trueRHS.SyncFromBlocks();

  _trueRHS.SetSize(width);
  _trueRHS = trueRHS;

  // Create monolithic matrix
  op.Reset(mfem::HypreParMatrixFromBlocks(_h_blocks));
}

void
TimeDependentEquationSystem::FormSystem(mfem::OperatorHandle & op,
                                        mfem::BlockVector & truedXdt,
                                        mfem::BlockVector & trueRHS) const
{
  std::cout << "********** TIME FormSystem: I AM HERE ******************************************" << std::endl;
  auto & test_var_name = _test_var_names.at(0);
  auto td_blf = _td_blfs.Get(test_var_name);
  auto blf = _blfs.Get(test_var_name);
  auto lf = _lfs.Get(test_var_name);
  // if implicit, add contribution to linear form from terms involving state
  // variable at previous timestep: {

  // The AddMult method in mfem::BilinearForm is not defined for non-legacy assembly
  mfem::Vector lf_prev(lf->Size());
  //blf->Mult(*_trial_variables.Get(test_var_name), lf_prev);
  blf->Mult(_trueBlockX_Old.GetBlock(0), lf_prev);
  *lf -= lf_prev;
  // }
  mfem::Vector aux_x, aux_rhs;
  // Update solution values on Dirichlet values to be in terms of du/dt instead of u
  mfem::Vector bc_x = *(_xs.at(0).get());
  //bc_x -= *_trial_variables.Get(test_var_name);
  bc_x -= _trueBlockX_Old.GetBlock(0);
  bc_x /= _dt_coef.constant;

  // Form linear system for operator acting on vector of du/dt
  mfem::OperatorPtr * aux_a = new mfem::OperatorPtr;
  // Ownership of aux_a goes to the blf
  td_blf->FormLinearSystem(_ess_tdof_lists.at(0), bc_x, *lf, *aux_a, aux_x, aux_rhs);

  truedXdt.GetBlock(0) = aux_x;
  trueRHS.GetBlock(0) = aux_rhs;
  truedXdt.SyncFromBlocks();
  trueRHS.SyncFromBlocks();

  _trueRHS.SetSize(width);
  _trueRHS = trueRHS;

  // Create monolithic matrix
  op.Reset(aux_a->Ptr());
}

void
TimeDependentEquationSystem::UpdateEquationSystem()
{
  BuildBilinearForms();
  BuildMixedBilinearForms();
  BuildLinearForms();
}

} // namespace platypus
