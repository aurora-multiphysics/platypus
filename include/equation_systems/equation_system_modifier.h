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

class EquationSystemModifier
{

public:
  EquationSystemModifier(std::shared_ptr<EquationSystemData> data) : _equation_system_data{data} {};
  ~EquationSystemModifier() = default;

  bool VectorContainsName(const std::vector<std::string> & the_vector,
                                   const std::string & name) const;
  void AddTrialVariableNameIfMissing(const std::string & trial_var_name);
  void AddTestVariableNameIfMissing(const std::string & test_var_name);
  void ApplyBoundaryConditions(platypus::BCMap & bc_map);

  virtual void BuildLinearForms(platypus::BCMap & bc_map);
  virtual void BuildBilinearForms();
  virtual void BuildMixedBilinearForms();
  virtual void BuildEquationSystem(platypus::BCMap & bc_map);

  [[nodiscard]] EquationSystemData * getData()
  {
    if (!_equation_system_data)
    {
      MFEM_ABORT("platypus::EquationSystemData instance is NULL.");
    }

    return _equation_system_data.get();
  }

  protected:

  std::shared_ptr<EquationSystemData> _equation_system_data{nullptr};

};

} // namespace platypus

