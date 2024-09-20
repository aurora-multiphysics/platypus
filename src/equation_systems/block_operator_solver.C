#include "block_operator_solver.h"

namespace platypus
{

void BlockOperatorSolver::AddPreconditioner(std::shared_ptr<mfem::Solver> prec)
{
  _preconditioners.push_back(prec);
}

void BlockOperatorSolver::AddSolver(std::shared_ptr<mfem::Solver> solver)
{
  _solvers.push_back(solver);
}

void BlockOperatorSolver::SolveBlocks(const mfem::BlockVector & x, mfem::BlockVector & residual)
{
  MFEM_VERIFY(_solvers.size() == NumRowBlocks(), "Mismatched number of solvers and blocks!");

  for (int i=0; i<NumRowBlocks(); ++i)
  {
    _solvers.at(i)->SetOperator(GetBlock(i,i));
    _solvers.at(i)->Mult(x.GetBlock(i), residual.GetBlock(i));
  }
}


} //namespace platypus