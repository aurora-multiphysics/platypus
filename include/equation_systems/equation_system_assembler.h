#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.h"
#include "named_fields_map.h"
#include "MFEMKernel.h"
#include "equation_system_data.h"

namespace platypus
{

/*
Class to store weak form components (bilinear and linear forms, and optionally
mixed and nonlinear forms) and build methods
*/
class EquationSystemAssembler
{
public:
  using MFEMBilinearFormKernel = MFEMKernel<mfem::BilinearFormIntegrator>;
  using MFEMLinearFormKernel = MFEMKernel<mfem::LinearFormIntegrator>;
  using MFEMNonlinearFormKernel = MFEMKernel<mfem::NonlinearFormIntegrator>;
  using MFEMMixedBilinearFormKernel = MFEMKernel<mfem::BilinearFormIntegrator>;

  EquationSystemAssembler(std::shared_ptr<EquationSystemData> data);
  ~EquationSystemAssembler() = default;

  void MakeData();

  void AddTrialVariableNameIfMissing(const std::string & trial_var_name);
  void AddTestVariableNameIfMissing(const std::string & test_var_name);

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
  virtual void BuildLinearForms(platypus::BCMap & bc_map);
  virtual void BuildBilinearForms();
  virtual void BuildMixedBilinearForms();

  virtual bool AssemblyIsSupported() = 0;
  virtual void
  FormSystem(mfem::OperatorHandle & op, mfem::BlockVector & trueX, mfem::BlockVector & trueRHS) = 0;
  virtual void BuildEquationSystem(platypus::BCMap & bc_map) = 0;

  [[nodiscard]] EquationSystemData * getData()
  {
    if (!_equation_system_data)
    {
      MFEM_ABORT("platypus::EquationSystemData instance is NULL.");
    }

    return _equation_system_data.get();
  };

private:
  std::shared_ptr<EquationSystemData> _equation_system_data{nullptr};

protected:
  bool VectorContainsName(const std::vector<std::string> & the_vector,
                          const std::string & name) const;
};

class DiagonalEquationSystemAssembler : EquationSystemAssembler
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

class DenseEquationSystemAssembler : EquationSystemAssembler
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
