#include "problem_builder.h"

namespace platypus
{

Problem::~Problem()
{
  // Ensure that all owned memory is properly freed!
  _f.reset();
  _ode_solver.reset();
}

void
ProblemBuilder::SetMesh(std::shared_ptr<mfem::ParMesh> pmesh)
{
  // logger.info("Setting Mesh");
  GetProblem()->_pmesh = pmesh;
  GetProblem()->_comm = pmesh->GetComm();
  MPI_Comm_size(pmesh->GetComm(), &(GetProblem()->_num_procs));
  MPI_Comm_rank(pmesh->GetComm(), &(GetProblem()->_myid));
}

void
ProblemBuilder::SetFESpaces(platypus::FESpaces & fespaces)
{
  // logger.info("Setting FE Spaces");
  GetProblem()->_fespaces = fespaces;
}

void
ProblemBuilder::SetGridFunctions(platypus::GridFunctions & gridfunctions)
{
  // logger.info("Setting GridFunctions");
  GetProblem()->_gridfunctions = gridfunctions;
}

void
ProblemBuilder::AddFESpace(std::string fespace_name, std::string fec_name, int vdim, int ordering)
{
  // logger.info("Adding {} FE Space to problem", fespace_name);
  if (GetProblem()->_fespaces.Has(fespace_name))
  {
    const std::string error_message = "A fespace with the name " + fespace_name +
                                      " has already been added to the problem fespaces.";
    mfem::mfem_error(error_message.c_str());
  }
  if (!GetProblem()->_fecs.Has(fec_name))
  {
    auto fec = std::shared_ptr<mfem::FiniteElementCollection>(
        mfem::FiniteElementCollection::New(fec_name.c_str()));
    GetProblem()->_fecs.Register(fec_name, fec);
  }

  if (!GetProblem()->_fespaces.Has(fespace_name))
  {
    mfem::ParMesh * pmesh = GetProblem()->_pmesh.get();
    if (pmesh == nullptr)
    {
      MFEM_ABORT("ParMesh not found when trying to add " << fespace_name << " to fespaces.");
    }
    auto pfes = std::make_shared<mfem::ParFiniteElementSpace>(
        GetProblem()->_pmesh.get(), GetProblem()->_fecs.Get(fec_name), vdim, ordering);

    GetProblem()->_fespaces.Register(fespace_name, std::move(pfes));
  }
}

void
ProblemBuilder::AddGridFunction(std::string gridfunction_name, std::string fespace_name)
{
  if (GetProblem()->_gridfunctions.Has(gridfunction_name))
  {
    const std::string error_message = "A gridfunction with the name " + gridfunction_name +
                                      " has already been added to the problem gridfunctions.";
    mfem::mfem_error(error_message.c_str());
  }

  if (!GetProblem()->_fespaces.Has(fespace_name))
  {
    MFEM_ABORT("FESpace " << fespace_name << " not found in fespaces when trying to add "
                          << gridfunction_name
                          << " associated with it into gridfunctions. Please add " << fespace_name
                          << " to fespaces before adding this gridfunction.");
  }

  auto fespace = GetProblem()->_fespaces.Get(fespace_name);

  auto gridfunc = std::make_shared<mfem::ParGridFunction>(fespace);
  *gridfunc = 0.0;

  GetProblem()->_gridfunctions.Register(gridfunction_name, std::move(gridfunc));
}

void
ProblemBuilder::ConstructNonlinearSolver()
{
  auto nl_solver = std::make_shared<mfem::NewtonSolver>(GetProblem()->_comm);

  // Defaults to one iteration, without further nonlinear iterations
  nl_solver->SetRelTol(0.0);
  nl_solver->SetAbsTol(0.0);
  nl_solver->SetMaxIter(1);

  GetProblem()->_nonlinear_solver = nl_solver;
}

void
ProblemBuilder::InitializeKernels()
{
}

void
ProblemBuilder::InitializeOutputs()
{
  GetProblem()->_outputs.Init(GetProblem()->_gridfunctions);
}

void
ProblemBuilder::FinalizeProblem(bool build_operator)
{
  RegisterFESpaces();
  RegisterGridFunctions();
  RegisterCoefficients();

  if (build_operator)
  {
    ConstructOperator();
  }

  InitializeKernels();
  SetOperatorGridFunctions();

  ConstructNonlinearSolver();

  ConstructState();
  ConstructTimestepper();
  InitializeOutputs();
}

} // namespace platypus
