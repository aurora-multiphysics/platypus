#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.h"
#include "named_fields_map.h"
#include "MFEMKernel.h"
#include "equation_system_data.h"
#include "equation_system_modifier.h"

namespace platypus
{

/*
Class to store weak form components (bilinear and linear forms, and optionally
mixed and nonlinear forms) and build methods
*/
class EquationSystemAssembler : public EquationSystemModifier
{
public:
  using MFEMBilinearFormKernel = MFEMKernel<mfem::BilinearFormIntegrator>;
  using MFEMLinearFormKernel = MFEMKernel<mfem::LinearFormIntegrator>;
  using MFEMNonlinearFormKernel = MFEMKernel<mfem::NonlinearFormIntegrator>;
  using MFEMMixedBilinearFormKernel = MFEMKernel<mfem::BilinearFormIntegrator>;

  EquationSystemAssembler(std::shared_ptr<EquationSystemData> data);
  ~EquationSystemAssembler() = default;

  // Add Kernels
  void AddKernel(const std::string & test_var_name,
                 std::shared_ptr<MFEMBilinearFormKernel> blf_kernel);
  void AddKernel(const std::string & test_var_name,
                 std::shared_ptr<MFEMLinearFormKernel> lf_kernel);
  void AddKernel(const std::string & test_var_name,
                 std::shared_ptr<MFEMNonlinearFormKernel> nlf_kernel);
  void AddKernel(const std::string & trial_var_name,
                 const std::string & test_var_name,
                 std::shared_ptr<MFEMMixedBilinearFormKernel> mblf_kernel);

  void ApplyBoundaryConditions(platypus::BCMap & bc_map);

  void Init(platypus::GridFunctions & gridfunctions,
            const platypus::FESpaces & fespaces,
            platypus::BCMap & bc_map,
            platypus::Coefficients & coefficients,
            mfem::AssemblyLevel assembly_level);

  virtual bool AssemblyIsSupported() = 0;
  virtual void
  FormSystem(mfem::OperatorHandle & op, mfem::BlockVector & trueX, mfem::BlockVector & trueRHS) = 0;

  // Build linear system, with essential boundary conditions accounted for
  virtual void BuildJacobian(mfem::BlockVector & trueX, mfem::BlockVector & trueRHS);
};

class DiagonalEquationSystemAssembler : public EquationSystemAssembler
{
public:
  DiagonalEquationSystemAssembler() = default;
  ~DiagonalEquationSystemAssembler() = default;

  bool AssemblyIsSupported() override;
  void FormSystem(mfem::OperatorHandle & op,
                  mfem::BlockVector & trueX,
                  mfem::BlockVector & trueRHS) override;
  void BuildEquationSystem(platypus::BCMap & bc_map) override;
};

class DenseEquationSystemAssembler : public EquationSystemAssembler
{
public:
  DenseEquationSystemAssembler() = default;
  ~DenseEquationSystemAssembler() = default;

  bool AssemblyIsSupported() override;
  void FormSystem(mfem::OperatorHandle & op,
                  mfem::BlockVector & trueX,
                  mfem::BlockVector & trueRHS) override;
  void BuildEquationSystem(platypus::BCMap & bc_map) override;
};

} // namespace platypus
