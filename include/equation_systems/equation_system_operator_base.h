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

  virtual void AddTrialVariableNameIfMissing(const std::string & trial_var_name) = 0;
  virtual void AddTestVariableNameIfMissing(const std::string & test_var_name) = 0;
  virtual void ApplyBoundaryConditions(platypus::BCMap & bc_map) = 0;
  virtual bool AssemblyIsSupported() = 0;

  // Forms methods
  virtual void BuildLinearForms(platypus::BCMap & bc_map) = 0;
  virtual void BuildBilinearForms() = 0;
  virtual void BuildMixedBilinearForms() = 0;
  virtual void BuildEquationSystem(platypus::BCMap & bc_map) = 0;

  // Add Kernel methods
  virtual void AddKernel(const std::string & test_var_name,
                 std::shared_ptr<MFEMBilinearFormKernel> blf_kernel) = 0;
  virtual void AddKernel(const std::string & test_var_name,
                 std::shared_ptr<MFEMLinearFormKernel> lf_kernel) = 0;
  virtual void AddKernel(const std::string & test_var_name,
                 std::shared_ptr<MFEMNonlinearFormKernel> nlf_kernel) = 0;
  virtual void AddKernel(const std::string & trial_var_name,
                 const std::string & test_var_name,
                 std::shared_ptr<MFEMMixedBilinearFormKernel> mblf_kernel) = 0;

  virtual void Init(platypus::GridFunctions & gridfunctions,
            const platypus::FESpaces & fespaces,
            platypus::BCMap & bc_map,
            platypus::Coefficients & coefficients,
            mfem::AssemblyLevel assembly_level) = 0;

  // Operator methods
  virtual void Mult(const mfem::Vector & u, mfem::Vector & residual) const = 0;
  virtual mfem::Operator & GetGradient(const mfem::Vector & u) const = 0;
  virtual void RecoverFEMSolution(mfem::BlockVector & trueX,
                                  platypus::GridFunctions & gridfunctions) = 0;

  // Form system methods
  virtual void
  FormSystem(mfem::OperatorHandle & op, mfem::BlockVector & trueX, mfem::BlockVector & trueRHS) = 0;
  virtual void
  FormDenseSystem(mfem::OperatorHandle & op,
                                            mfem::BlockVector & trueX,
                                            mfem::BlockVector & trueRHS) = 0;
  virtual void
  FormDiagonalSystem(mfem::OperatorHandle & op,
                                            mfem::BlockVector & trueX,
                                            mfem::BlockVector & trueRHS) = 0;

  // Build linear system, with essential boundary conditions accounted for
  virtual void BuildJacobian(mfem::BlockVector & trueX, mfem::BlockVector & trueRHS) = 0;

};


} // namespace platypus
