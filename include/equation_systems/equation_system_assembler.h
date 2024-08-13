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

  EquationSystemAssembler(mfem::AssemblyLevel assembly_level);
  ~EquationSystemAssembler() = default;

  void MakeData();

  bool VectorContainsName(const std::vector<std::string> & the_vector, const std::string & name) const;

  void AddTrialVariableNameIfMissing(const std::string & trial_var_name);
  void AddTestVariableNameIfMissing(const std::string & test_var_name);

  void AddKernel(const std::string & test_var_name, std::shared_ptr<MFEMBilinearFormKernel> blf_kernel);
  void AddKernel(const std::string & test_var_name, std::shared_ptr<MFEMLinearFormKernel> lf_kernel);
  void AddKernel(const std::string & test_var_name, std::shared_ptr<MFEMNonlinearFormKernel> nlf_kernel);
  void AddKernel(const std::string & trial_var_name, const std::string & test_var_name, std::shared_ptr<MFEMMixedBilinearFormKernel> mblf_kernel);

  void ApplyBoundaryConditions(platypus::BCMap & bc_map);


  [[nodiscard]] EquationSystemData * getData()
  {
    if (!_equation_system_data)
    {
      MFEM_ABORT("platypus::EquationSystemData instance is NULL.");
    }

    return _equation_system_data.get();

  };

private:

  std::shared_ptr<platypus::EquationSystemData> _equation_system_data{nullptr};

};

class BlockDiagonalEquationSystemAssembler : EquationSystemAssembler
{
public:

  BlockDiagonalEquationSystemAssembler(mfem::AssemblyLevel assembly_level) override;
  ~BlockDiagonalEquationSystemAssembler() override = default;

};

class NonBlockDiagonalEquationSystemAssembler : EquationSystemAssembler
{
  NonBlockDiagonalEquationSystemAssembler() override;
  ~NonBlockDiagonalEquationSystemAssembler() override = default;
};


// Can be derived from Block-Diagonal or Non-Block-Diagonal assemblers
template<class T>
class TimeDependentEquationSystemAssembler : public T
{
public:
  TimeDependentEquationSystemAssembler();
  ~TimeDependentEquationSystemAssembler() override = default;

  static std::string GetTimeDerivativeName(std::string name)
  {
    return std::string("d") + name + std::string("_dt");
  }

  void AddTrialVariableNameIfMissing(const std::string & trial_var_name) override;

  virtual void SetTimeStep(double dt);
  virtual void UpdateEquationSystem(platypus::BCMap & bc_map);
};

} // namespace platypus
