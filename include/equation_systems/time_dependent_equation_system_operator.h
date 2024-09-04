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
  TimeDependentEquationSystemOperator(std::shared_ptr<TimeDependentEquationSystemData> data)
    : _equation_system_data{data}
  {
    mfem::ConstantCoefficient dt(1.0);
    GetData()->_dt_coef = dt;
  }

  void SetTimeStep(double dt);
  void UpdateEquationSystem(platypus::BCMap & bc_map);
  void AddTrialVariableNameIfMissing(const std::string & var_name) override;

  static std::string GetTimeDerivativeName(std::string name)
  {
    return std::string("d") + name + std::string("_dt");
  }

  // Data retrieval for writing to @_equation_system_data
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
