#include "time_domain_problem_builder.h"

namespace platypus
{

std::vector<mfem::ParGridFunction *>
TimeDomainProblemBuilder::RegisterTimeDerivatives(std::vector<std::string> gridfunction_names,
                                                  platypus::GridFunctions & gridfunctions)
{
  std::vector<mfem::ParGridFunction *> time_derivatives;

  for (auto & gridfunction_name : gridfunction_names)
  {
    gridfunctions.Register(GetTimeDerivativeName(gridfunction_name),
                           std::make_shared<mfem::ParGridFunction>(
                               gridfunctions.Get(gridfunction_name)->ParFESpace()));

    time_derivatives.push_back(gridfunctions.Get(GetTimeDerivativeName(gridfunction_name)));
  }

  return time_derivatives;
}

void
TimeDomainProblemBuilder::RegisterGridFunctions()
{
  std::vector<std::string> gridfunction_names;
  for (auto const & [name, gf] : _problem->_gridfunctions)
  {
    gridfunction_names.push_back(name);
  }
  RegisterTimeDerivatives(gridfunction_names, _problem->_gridfunctions);
}

void
TimeDomainProblemBuilder::SetOperatorGridFunctions()
{
  GetOperator().SetGridFunctions();
}

void
TimeDomainProblemBuilder::ConstructOperator()
{

  _problem_operator.reset();
  _problem_operator = std::make_shared<platypus::TimeDomainProblemOperator>(*_problem);
}

void
TimeDomainProblemBuilder::ConstructState()
{
  // Vector of dofs.
  _problem->_f = std::make_unique<mfem::BlockVector>(GetOperator()._true_offsets);
  *(_problem->_f) = 0.0;               // give initial value
  GetOperator().Init(*(_problem->_f)); // Set up initial conditions
  GetOperator().SetTime(0.0);
}

void
TimeDomainProblemBuilder::ConstructTimestepper()
{
  _problem->_ode_solver = std::make_unique<mfem::BackwardEulerSolver>();
  _problem->_ode_solver->Init(GetOperator());
}

} // namespace platypus
