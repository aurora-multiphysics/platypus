#pragma once
#include "MFEMSolverBase.h"
#include "mfem.hpp"
#include <memory>

/**
 * Wrapper for mfem::CGSolver.
 */
class MFEMCGSolver : public MFEMSolverBase
{
public:
  static InputParameters validParams();

  MFEMCGSolver(const InputParameters & parameters);

  /// Returns a shared pointer to the instance of the Solver derived-class.
  std::shared_ptr<mfem::Solver> getSolver() const override { return _solver; }

protected:
  void constructSolver(const InputParameters & parameters) override;

private:
  std::shared_ptr<mfem::Solver> _preconditioner{nullptr};
  std::shared_ptr<mfem::CGSolver> _solver{nullptr};
};
