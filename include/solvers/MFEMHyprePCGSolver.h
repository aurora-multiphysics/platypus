#pragma once
#include "MFEMSolverBase.h"
#include "mfem.hpp"
#include <memory>

/**
 * Wrapper for mfem::HyprePCG solver.
 */
class MFEMHyprePCGSolver : public MFEMSolverBase
{
public:
  static InputParameters validParams();

  MFEMHyprePCGSolver(const InputParameters & parameters);

  std::shared_ptr<mfem::Solver> getSolver() const override { return _solver; }

protected:
  void constructSolver(const InputParameters & parameters) override;

private:
  std::shared_ptr<mfem::HyprePCG> _solver{nullptr};
};
