#include "problem_builder.h"

namespace platypus
{

void
ProblemBuilder::ConstructNonlinearSolver()
{
  auto nl_solver = std::make_shared<mfem::NewtonSolver>(_problem->_pmesh->GetComm());

  // Defaults to one iteration, without further nonlinear iterations
  nl_solver->SetRelTol(0.0);
  nl_solver->SetAbsTol(0.0);
  nl_solver->SetMaxIter(1);

  _problem->_nonlinear_solver = nl_solver;
}

void
ProblemBuilder::InitializeKernels()
{
}

void
ProblemBuilder::InitializeOutputs()
{
  _problem->_outputs.Init(_problem->_gridfunctions);
}

void
ProblemBuilder::FinalizeProblem()
{
  RegisterGridFunctions();

  InitializeKernels();
  SetOperatorGridFunctions();

  ConstructNonlinearSolver();

  ConstructState();
  ConstructTimestepper();
  InitializeOutputs();
}

} // namespace platypus
