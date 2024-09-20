#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.h"
#include "named_fields_map.h"
#include "MFEMKernel.h"
#include "equation_system_data.h"

namespace platypus
{

/*
Base class to store weak form components (bilinear and linear forms, and optionally
mixed and nonlinear forms) and build methods for static problems
*/
class EquationSystemOperatorBase : public mfem::Operator
{

public:
  using MFEMBilinearFormKernel = MFEMKernel<mfem::BilinearFormIntegrator>;
  using MFEMLinearFormKernel = MFEMKernel<mfem::LinearFormIntegrator>;
  using MFEMNonlinearFormKernel = MFEMKernel<mfem::NonlinearFormIntegrator>;
  using MFEMMixedBilinearFormKernel = MFEMKernel<mfem::BilinearFormIntegrator>;

  EquationSystemOperatorBase() = default;
  ~EquationSystemOperatorBase() = default;

  // add test variable to EquationSystem;
  virtual void AddTestVariableNameIfMissing(const std::string & test_var_name);
  virtual void AddTrialVariableNameIfMissing(const std::string & trial_var_name);

  // Add kernels.
  virtual void AddKernel(const std::string & test_var_name,
                         std::shared_ptr<MFEMBilinearFormKernel> blf_kernel);

  void AddKernel(const std::string & test_var_name,
                 std::shared_ptr<MFEMLinearFormKernel> lf_kernel);

  void AddKernel(const std::string & test_var_name,
                 std::shared_ptr<MFEMNonlinearFormKernel> nlf_kernel);

  void AddKernel(const std::string & trial_var_name,
                 const std::string & test_var_name,
                 std::shared_ptr<MFEMMixedBilinearFormKernel> mblf_kernel);

  virtual void ApplyBoundaryConditions(platypus::BCMap & bc_map);

  // Build forms
  virtual void Init(platypus::GridFunctions & gridfunctions,
                    const platypus::FESpaces & fespaces,
                    platypus::BCMap & bc_map,
                    mfem::AssemblyLevel assembly_level);
  virtual void BuildLinearForms(platypus::BCMap & bc_map);
  virtual void BuildBilinearForms();
  virtual void BuildMixedBilinearForms();
  virtual void BuildEquationSystem(platypus::BCMap & bc_map);

  // Forms the BlockOperatorSolver that holds the equation system and respective solvers
  virtual void MakeBlockOperatorSolver();

  // Form linear system, with essential boundary conditions accounted for
  virtual void FormSystem(mfem::OperatorHandle & op,
                                mfem::BlockVector & trueX,
                                mfem::BlockVector & trueRHS);
  virtual void FormDenseSystem(mfem::OperatorHandle & op,
                                mfem::BlockVector & trueX,
                                mfem::BlockVector & trueRHS);
  virtual void FormDiagonalSystem(mfem::OperatorHandle & op,
                                mfem::BlockVector & trueX,
                                mfem::BlockVector & trueRHS);



  // Build linear system, with essential boundary conditions accounted for
  virtual void BuildJacobian(mfem::BlockVector & trueX, mfem::BlockVector & trueRHS);

  /// Compute residual y = Mu
  void Mult(const mfem::Vector & u, mfem::Vector & residual) const override;

  /// Compute J = M + grad_H(u)
  virtual mfem::Operator & GetGradient(const mfem::Vector & u) const override;

  // Update variable from solution vector after solve
  virtual void RecoverFEMSolution(mfem::BlockVector & trueX,
                                  platypus::GridFunctions & gridfunctions);

  /**
   * Template method for storing kernels.
   */
  template <class T>
  void addKernelToMap(std::shared_ptr<T> kernel,
                      platypus::NamedFieldsMap<std::vector<std::shared_ptr<T>>> & kernels_map)
  {
    auto test_var_name = kernel->getTestVariableName();
    if (!kernels_map.Has(test_var_name))
    {
      // 1. Create kernels vector.
      auto kernels = std::make_shared<std::vector<std::shared_ptr<T>>>();
      // 2. Register with map to prevent leaks.
      kernels_map.Register(test_var_name, std::move(kernels));
    }
    kernels_map.GetRef(test_var_name).push_back(std::move(kernel));
  }

  virtual EquationSystemData * GetData() const = 0;

protected:
  bool VectorContainsName(const std::vector<std::string> & the_vector,
                          const std::string & name) const;

private:
  std::shared_ptr<EquationSystemData> _equation_system_data{nullptr};
};

/*
Class to store weak form components (bilinear and linear forms, and optionally
mixed and nonlinear forms) and build methods for static problems
*/
class EquationSystemOperator : public EquationSystemOperatorBase
{

public:
  EquationSystemOperator();
  ~EquationSystemOperator();

  // Data retrieval for reading/writing to @_equation_system_data
  EquationSystemData * GetData() const override
  {
    if (!_equation_system_data)
    {
      MFEM_ABORT("platypus::EquationSystemData instance is NULL.");
    }
    return _equation_system_data.get();
  }

private:
  std::shared_ptr<EquationSystemData> _equation_system_data{nullptr};
};

/*
Class to store weak form components (bilinear and linear forms, and optionally
mixed and nonlinear forms) and build methods for time dependent problems
*/
class TimeDependentEquationSystemOperator : public EquationSystemOperatorBase
{

public:
  TimeDependentEquationSystemOperator();
  ~TimeDependentEquationSystemOperator();

  void AddTrialVariableNameIfMissing(const std::string & trial_var_name) override;

  virtual void SetTimeStep(double dt);
  virtual void UpdateEquationSystem(platypus::BCMap & bc_map);

  virtual void AddKernel(const std::string & test_var_name,
                         std::shared_ptr<MFEMBilinearFormKernel> blf_kernel) override;
  virtual void BuildBilinearForms() override;
  virtual void FormDenseSystem(mfem::OperatorHandle & op,
                                mfem::BlockVector & truedXdt,
                                mfem::BlockVector & trueRHS) override;

  // Data retrieval for reading/writing to @_equation_system_data
  TimeDependentEquationSystemData * GetData() const override
  {
    if (!_equation_system_data)
    {
      MFEM_ABORT("platypus::TimeDependentEquationSystemData instance is NULL.");
    }
    return _equation_system_data.get();
  }

private:
  std::shared_ptr<TimeDependentEquationSystemData> _equation_system_data{nullptr};
};

} // namespace platypus