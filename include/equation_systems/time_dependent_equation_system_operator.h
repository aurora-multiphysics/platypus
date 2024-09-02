#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.h"
#include "named_fields_map.h"
#include "MFEMKernel.h"
#include "equation_system_data.h"
#include "equation_system_operator_base.h"

namespace platypus
{

// Time dependent operator

class TimeDependentEquationSystemOperator : public EquationSystemOperatorBase
{

public:

  // Constructor
  TimeDependentEquationSystemOperator(std::shared_ptr<TimeDependentEquationSystemData> data)  : _equation_system_data{data}
  {
    mfem::ConstantCoefficient dt(1.0);
    DataWrite()->_dt_coef = dt;
  }

  void AddTrialVariableNameIfMissing(const std::string & trial_var_name) override;
  void AddTestVariableNameIfMissing(const std::string & test_var_name) override;
  void ApplyBoundaryConditions(platypus::BCMap & bc_map) override;
  bool AssemblyIsSupported() override;

  // Forms methods
  void BuildLinearForms(platypus::BCMap & bc_map) override;
  void BuildBilinearForms() override;
  void BuildMixedBilinearForms() override;
  void BuildEquationSystem(platypus::BCMap & bc_map) override;

  // Add Kernel methods
  void AddKernel(const std::string & test_var_name,
                 std::shared_ptr<MFEMBilinearFormKernel> blf_kernel) override;
  void AddKernel(const std::string & test_var_name,
                 std::shared_ptr<MFEMLinearFormKernel> lf_kernel) override;
  void AddKernel(const std::string & test_var_name,
                 std::shared_ptr<MFEMNonlinearFormKernel> nlf_kernel) override;
  void AddKernel(const std::string & trial_var_name,
                 const std::string & test_var_name,
                 std::shared_ptr<MFEMMixedBilinearFormKernel> mblf_kernel) override;

  void Init(platypus::GridFunctions & gridfunctions,
            const platypus::FESpaces & fespaces,
            platypus::BCMap & bc_map,
            platypus::Coefficients & coefficients,
            mfem::AssemblyLevel assembly_level) override;

  // Operator methods
  void Mult(const mfem::Vector & u, mfem::Vector & residual) const override;
  mfem::Operator & GetGradient(const mfem::Vector & u) const override;
  void RecoverFEMSolution(mfem::BlockVector & trueX,
                                  platypus::GridFunctions & gridfunctions) override;

  // Form system methods
  void
  FormSystem(mfem::OperatorHandle & op, mfem::BlockVector & trueX, mfem::BlockVector & trueRHS) override;
  void
  FormDenseSystem(mfem::OperatorHandle & op,
                  mfem::BlockVector & trueX,
                  mfem::BlockVector & trueRHS) override;
  void
  FormDiagonalSystem(mfem::OperatorHandle & op,
                     mfem::BlockVector & trueX,
                     mfem::BlockVector & trueRHS) override;

  // Build linear system, with essential boundary conditions accounted for
  void BuildJacobian(mfem::BlockVector & trueX, mfem::BlockVector & trueRHS) override;



  void SetTimeStep(double dt);
  void UpdateEquationSystem(platypus::BCMap & bc_map);

  static std::string GetTimeDerivativeName(std::string name)
  {
    return std::string("d") + name + std::string("_dt");
  }

  // Data retrieval for writing to @_equation_system_data
  std::shared_ptr<TimeDependentEquationSystemData> DataWrite()
  {
    if (!_equation_system_data)
    {
      MFEM_ABORT("platypus::TimeDependentEquationSystemData instance is NULL.");
    }
    return _equation_system_data;
  }

  // Data retrieval for reading @_equation_system_data
  std::shared_ptr<TimeDependentEquationSystemData> DataRead() const
  {
    if (!_equation_system_data)
    {
      MFEM_ABORT("platypus::EquationSystemData instance is NULL.");
    }
    return _equation_system_data;
  }


private:

  std::shared_ptr<TimeDependentEquationSystemData> _equation_system_data{nullptr};

};


} // namespace platypus
