#pragma once
#include "MFEMSolverBase.h"
#include "mfem.hpp"
#include <memory>

/**
 * Wrapper for mfem::HypreBoomerAMG solver.
 */
class MFEMHypreAMGSolver : public MFEMSolverBase
{
public:
  static InputParameters validParams();

  MFEMHypreAMGSolver(const InputParameters & parameters);

  std::shared_ptr<mfem::Solver> getSolver() const override { return _solver; }

protected:
  void constructSolver(const InputParameters & parameters) override;

private:
  std::shared_ptr<mfem::HypreBoomerAMG> _solver{nullptr};
};
