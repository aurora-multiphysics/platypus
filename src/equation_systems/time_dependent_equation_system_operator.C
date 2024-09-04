#include "time_dependent_equation_system_operator.h"

namespace platypus
{

void
TimeDependentEquationSystemOperator::AddTrialVariableNameIfMissing(const std::string & var_name)
{

  if (!VectorContainsName(GetData()->_trial_var_names, var_name))
  {
    GetData()->_trial_var_names.push_back(var_name);
  }

  std::string var_time_derivative_name = GetTimeDerivativeName(var_name);
  if (std::find(GetData()->_trial_var_time_derivative_names.begin(),
                GetData()->_trial_var_time_derivative_names.end(),
                var_time_derivative_name) == GetData()->_trial_var_time_derivative_names.end())
  {
    GetData()->_trial_var_time_derivative_names.push_back(var_time_derivative_name);
  }
}

void
TimeDependentEquationSystemOperator::SetTimeStep(double dt)
{
  if (fabs(dt - GetData()->_dt_coef.constant) > 1.0e-12 * dt)
  {
    GetData()->_dt_coef.constant = dt;
    for (auto test_var_name : GetData()->_test_var_names)
    {
      auto blf = GetData()->_blfs.Get(test_var_name);
      blf->Update();
      blf->Assemble();
    }
  }
}

void
TimeDependentEquationSystemOperator::UpdateEquationSystem(platypus::BCMap & bc_map)
{
  BuildEquationSystem(bc_map);
}

} // namespace platypus