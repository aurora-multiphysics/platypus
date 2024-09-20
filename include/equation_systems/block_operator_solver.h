#pragma once
#include "../common/pfem_extras.hpp"
#include "inputs.h"
#include "named_fields_map.h"
#include "MFEMKernel.h"



namespace platypus
{

// This is a wrapper class for MFEM's block operator which solves independent blocks separately
class BlockOperatorSolver : public mfem::BlockOperator
{

public:

  using BlockOperator::BlockOperator;
  ~BlockOperatorSolver() = default;

  void SolveBlocks(const mfem::BlockVector & x, mfem::BlockVector & residual);
  void AddPreconditioner(std::shared_ptr<mfem::Solver> prec);
  void AddSolver(std::shared_ptr<mfem::Solver> solver);

private:

  std::vector<std::shared_ptr<mfem::Solver>> _preconditioners;
  std::vector<std::shared_ptr<mfem::Solver>> _solvers;

};

} // namespace platypus