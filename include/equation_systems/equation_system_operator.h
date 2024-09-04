#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.h"
#include "named_fields_map.h"
#include "MFEMKernel.h"
#include "equation_system_data.h"
#include "equation_system_operator_base.h"

namespace platypus
{

// Static operator

class EquationSystemOperator : public EquationSystemOperatorBase
{

public:
  // Constructor
  EquationSystemOperator(std::shared_ptr<EquationSystemData> data) : _equation_system_data{data} {}

  // Data retrieval for writing to @_equation_system_data
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

} // namespace platypus
