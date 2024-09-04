#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.h"
#include "named_fields_map.h"
#include "MFEMKernel.h"
#include "equation_system_data.h"

namespace platypus
{

/*
Base class containing the common methods to access and modify an EquationSystemData
*/

class EquationSystemOperatorBase : public mfem::Operator
{

public:
  EquationSystemOperatorBase() = default;
  ~EquationSystemOperatorBase() = default;

  bool VectorContainsName(const std::vector<std::string> & the_vector,
                          const std::string & name) const;

  // Methods to be defined in derived classes

  virtual void AddTrialVariableNameIfMissing(const std::string & trial_var_name);
  virtual void AddTestVariableNameIfMissing(const std::string & test_var_name);
  virtual void ApplyBoundaryConditions(platypus::BCMap & bc_map);

  // Forms methods
  virtual void BuildLinearForms(platypus::BCMap & bc_map);
  virtual void BuildBilinearForms();
  virtual void BuildMixedBilinearForms();
  virtual void BuildEquationSystem(platypus::BCMap & bc_map);

  // Add Kernel methods
  virtual void AddKernel(const std::string & test_var_name,
                         std::shared_ptr<MFEMBilinearFormKernel> blf_kernel);
  virtual void AddKernel(const std::string & test_var_name,
                         std::shared_ptr<MFEMLinearFormKernel> lf_kernel);
  virtual void AddKernel(const std::string & test_var_name,
                         std::shared_ptr<MFEMNonlinearFormKernel> nlf_kernel);
  virtual void AddKernel(const std::string & trial_var_name,
                         const std::string & test_var_name,
                         std::shared_ptr<MFEMMixedBilinearFormKernel> mblf_kernel);

  virtual void Init(platypus::GridFunctions & gridfunctions,
                    const platypus::FESpaces & fespaces,
                    platypus::BCMap & bc_map,
                    platypus::Coefficients & coefficients,
                    mfem::AssemblyLevel assembly_level);

  // Operator methods
  virtual void Mult(const mfem::Vector & u, mfem::Vector & residual) const;
  virtual mfem::Operator & GetGradient(const mfem::Vector & u) const;
  virtual void RecoverFEMSolution(mfem::BlockVector & trueX,
                                  platypus::GridFunctions & gridfunctions);

  // Form system methods
  virtual void
  FormSystem(mfem::OperatorHandle & op, mfem::BlockVector & trueX, mfem::BlockVector & trueRHS);
  virtual void FormDenseSystem(mfem::OperatorHandle & op,
                               mfem::BlockVector & trueX,
                               mfem::BlockVector & trueRHS);
  virtual void FormDiagonalSystem(mfem::OperatorHandle & op,
                                  mfem::BlockVector & trueX,
                                  mfem::BlockVector & trueRHS);

  virtual void MakeBlockOperator();

  // Build linear system, with essential boundary conditions accounted for
  virtual void BuildJacobian(mfem::BlockVector & trueX, mfem::BlockVector & trueRHS);

  // Data retrieval for writing to @_equation_system_data
  virtual EquationSystemData * GetData() const = 0;

private:
  std::shared_ptr<EquationSystemData> _equation_system_data{nullptr};
};

} // namespace platypus
