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
Class to use in time-dependent systems. It updates the data in EquationSystemData for the next
timestep.
*/

class EquationSystemUpdater : public EquationSystemModifier
{
public:
  EquationSystemUpdater(std::shared_ptr<EquationSystemData> data);
  ~EquationSystemUpdater() = default;

  void SetTimeStep(double dt);
  void UpdateEquationSystem(platypus::BCMap & bc_map);

  static std::string GetTimeDerivativeName(std::string name)
  {
    return std::string("d") + name + std::string("_dt");
  }
};

} // namespace platypus
