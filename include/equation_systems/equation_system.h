#pragma once
#include "mfem/miniapps/common/pfem_extras.hpp"
#include "MFEMIntegratedBC.h"
#include "MFEMEssentialBC.h"
#include "MFEMContainers.h"
#include "MFEMKernel.h"
#include "MFEMMixedBilinearFormKernel.h"
#include "ScaleIntegrator.h"

namespace platypus
{

/*
Class to store weak form components (bilinear and linear forms, and optionally
mixed and nonlinear forms) and build methods
*/
class EquationSystem : public mfem::Operator
{
public:
  EquationSystem() = default;
  ~EquationSystem() override;

  // Test variables are associated with LinearForms,
  // whereas trial variables are associated with gridfunctions.

  // Names of all variables corresponding to gridfunctions. This may differ
  // from test_var_names when time derivatives are present.
  std::vector<std::string> _trial_var_names;
  // Pointers to trial variables.
  platypus::GridFunctions _trial_variables;
  // Names of all test variables corresponding to linear forms in this equation
  // system
  std::vector<std::string> _test_var_names;
  std::vector<mfem::ParFiniteElementSpace *> _test_pfespaces;

  // Components of weak form. // Named according to test variable
  platypus::NamedFieldsMap<mfem::ParBilinearForm> _blfs;
  platypus::NamedFieldsMap<mfem::ParLinearForm> _lfs;
  platypus::NamedFieldsMap<mfem::ParNonlinearForm> _nlfs;
  platypus::NamedFieldsMap<platypus::NamedFieldsMap<mfem::ParMixedBilinearForm>>
      _mblfs; // named according to trial variable

  // add test variable to EquationSystem;
  virtual void AddTestVariableNameIfMissing(const std::string & test_var_name);
  virtual void AddTrialVariableNameIfMissing(const std::string & trial_var_name);

  // Add kernels.
  virtual void AddKernel(std::shared_ptr<MFEMKernel> kernel);
  virtual void AddIntegratedBC(std::shared_ptr<MFEMIntegratedBC> kernel);
  virtual void AddEssentialBC(std::shared_ptr<MFEMEssentialBC> bc);
  virtual void ApplyEssentialBCs();

  // Build forms
  virtual void Init(platypus::GridFunctions & gridfunctions,
                    const platypus::FESpaces & fespaces,
                    mfem::AssemblyLevel assembly_level);
  virtual void BuildLinearForms();
  virtual void BuildBilinearForms();
  virtual void BuildMixedBilinearForms();
  virtual void BuildEquationSystem(platypus::GridFunctions & gridfunctions, mfem::Array<int> & btoffsets);

  // Form linear system, with essential boundary conditions accounted for
  virtual void FormLinearSystem(mfem::OperatorHandle & op,
                                mfem::BlockVector & trueX,
                                mfem::BlockVector & trueRHS);

  virtual void
  FormSystem(mfem::OperatorHandle & op, mfem::BlockVector & trueX, mfem::BlockVector & trueRHS);
  virtual void FormLegacySystem(mfem::OperatorHandle & op,
                                mfem::BlockVector & trueX,
                                mfem::BlockVector & trueRHS);

  // Build linear system, with essential boundary conditions accounted for
  virtual void BuildJacobian(mfem::BlockVector & trueX, mfem::BlockVector & trueRHS);

  /// Compute residual y = Mu
  void Mult(const mfem::Vector & u, mfem::Vector & residual) const override;

  /// Compute J = M + grad_H(u)
  mfem::Operator & GetGradient(const mfem::Vector & u) const override;

  // Update variable from solution vector after solve
  virtual void RecoverFEMSolution(mfem::BlockVector & trueX,
                                  platypus::GridFunctions & gridfunctions);

  std::vector<mfem::Array<int>> _ess_tdof_lists;

protected:
  bool VectorContainsName(const std::vector<std::string> & the_vector,
                          const std::string & name) const;

  /**
   * Template method for applying BilinearFormIntegrators on domains from kernels to a BilinearForm,
   * or MixedBilinearForm
   */
  template <class FormType>
  void ApplyDomainBLFIntegrators(
      const std::string & trial_var_name,
      const std::string & test_var_name,
      std::shared_ptr<FormType> form,
      platypus::NamedFieldsMap<platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>> &
          kernels_map)
  {
    if (kernels_map.Has(test_var_name) && kernels_map.Get(test_var_name)->Has(trial_var_name))
    {
      auto kernels = kernels_map.GetRef(test_var_name).GetRef(trial_var_name);
      for (auto & kernel : kernels)
      {
        mfem::BilinearFormIntegrator * integ = kernel->createBFIntegrator();
        if (integ != nullptr)
        {
          kernel->isSubdomainRestricted()
              ? form->AddDomainIntegrator(std::move(integ), kernel->getSubdomains())
              : form->AddDomainIntegrator(std::move(integ));
        }
      }
    }
  }

  void ApplyDomainLFIntegrators(
      const std::string & test_var_name,
      std::shared_ptr<mfem::ParLinearForm> form,
      platypus::NamedFieldsMap<platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>> &
          kernels_map)
  {
    if (kernels_map.Has(test_var_name))
    {
      auto kernels = kernels_map.GetRef(test_var_name).GetRef(test_var_name);
      for (auto & kernel : kernels)
      {
        mfem::LinearFormIntegrator * integ = kernel->createLFIntegrator();
        if (integ != nullptr)
        {
          kernel->isSubdomainRestricted()
              ? form->AddDomainIntegrator(std::move(integ), kernel->getSubdomains())
              : form->AddDomainIntegrator(std::move(integ));
        }
      }
    }
  }

  template <class FormType>
  void ApplyBoundaryBLFIntegrators(
      const std::string & trial_var_name,
      const std::string & test_var_name,
      std::shared_ptr<FormType> form,
      platypus::NamedFieldsMap<
          platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMIntegratedBC>>>> &
          integrated_bc_map)
  {
    if (integrated_bc_map.Has(test_var_name) &&
        integrated_bc_map.Get(test_var_name)->Has(trial_var_name))
    {
      auto bcs = integrated_bc_map.GetRef(test_var_name).GetRef(trial_var_name);
      for (auto & bc : bcs)
      {
        mfem::BilinearFormIntegrator * integ = bc->createBFIntegrator();
        if (integ != nullptr)
        {
          bc->isBoundaryRestricted()
              ? form->AddBoundaryIntegrator(std::move(integ), bc->getBoundaries())
              : form->AddBoundaryIntegrator(std::move(integ));
        }
      }
    }
  }

  void ApplyBoundaryLFIntegrators(
      const std::string & test_var_name,
      std::shared_ptr<mfem::ParLinearForm> form,
      platypus::NamedFieldsMap<
          platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMIntegratedBC>>>> &
          integrated_bc_map)
  {
    if (integrated_bc_map.Has(test_var_name))
    {
      auto bcs = integrated_bc_map.GetRef(test_var_name).GetRef(test_var_name);
      for (auto & bc : bcs)
      {
        mfem::LinearFormIntegrator * integ = bc->createLFIntegrator();
        if (integ != nullptr)
        {
          bc->isBoundaryRestricted()
              ? form->AddBoundaryIntegrator(std::move(integ), bc->getBoundaries())
              : form->AddBoundaryIntegrator(std::move(integ));
        }
      }
    }
  }

  // gridfunctions for setting Dirichlet BCs
  std::vector<std::unique_ptr<mfem::ParGridFunction>> _xs;
  std::vector<std::unique_ptr<mfem::ParGridFunction>> _dxdts;

  mfem::Array2D<const mfem::HypreParMatrix *> _h_blocks;

  // Arrays to store kernels to act on each component of weak form. Named
  // according to test variable
  platypus::NamedFieldsMap<platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>>
      _kernels_map;
  platypus::NamedFieldsMap<platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMIntegratedBC>>>>
      _integrated_bc_map;
  platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMEssentialBC>>> _essential_bc_map;

  mutable mfem::OperatorHandle _jacobian;
  mutable mfem::Vector _trueRHS;
  mutable mfem::BlockVector _trueBlockX, _trueBlockRHS, _trueBlockdXdt, _trueBlockX_Old;
  mfem::AssemblyLevel _assembly_level;
  platypus::GridFunctions * _gfuncs;
  mfem::Array<int> * _block_true_offsets;
  mfem::Array<int> empty_tdof;
};

/*
Class to store weak form components for time dependent PDEs
*/
class TimeDependentEquationSystem : public EquationSystem
{
public:
  TimeDependentEquationSystem();
  ~TimeDependentEquationSystem() override = default;

  void AddTrialVariableNameIfMissing(const std::string & trial_var_name) override;

  virtual void SetTimeStep(double dt);
  virtual void UpdateEquationSystem();
  mfem::ConstantCoefficient _dt_coef; // Coefficient for timestep scaling
  std::vector<std::string> _trial_var_time_derivative_names;

  platypus::NamedFieldsMap<platypus::NamedFieldsMap<std::vector<std::shared_ptr<MFEMKernel>>>>
      _td_kernels_map;
  // Container to store contributions to weak form of the form (F du/dt, v)
  platypus::NamedFieldsMap<mfem::ParBilinearForm> _td_blfs;

  virtual void AddKernel(std::shared_ptr<MFEMKernel> kernel) override;
  virtual void BuildBilinearForms() override;
  virtual void FormLegacySystem(mfem::OperatorHandle & op,
                                mfem::BlockVector & truedXdt,
                                mfem::BlockVector & trueRHS) override;
  virtual void FormSystem(mfem::OperatorHandle & op,
                          mfem::BlockVector & truedXdt,
                          mfem::BlockVector & trueRHS) override;
};

} // namespace platypus
