#include "equation_system_updater.h"

namespace platypus
{

EquationSystemUpdater::EquationSystemUpdater(std::shared_ptr<EquationSystemData> data) : EquationSystemModifier(data)
{
    mfem::ConstantCoefficient dt(1.0);
    getData()->_dt_coef = dt;
}

void
EquationSystemUpdater::SetTimeStep(double dt)
{
  if (fabs(dt - getData()->_dt_coef.constant) > 1.0e-12 * dt)
  {
    getData()->_dt_coef.constant = dt;
    for (auto test_var_name : getData()->_test_var_names)
    {
      auto blf = getData()->_blfs.Get(test_var_name);
      blf->Update();
      blf->Assemble();
    }
  }
}

void
EquationSystemUpdater::UpdateEquationSystem(platypus::BCMap & bc_map)
{
    BuildEquationSystem(bc_map);
}


} // namespace platypus
