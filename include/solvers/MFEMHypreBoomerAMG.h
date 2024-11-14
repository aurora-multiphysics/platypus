#pragma once
#include "MFEMSolverBase.h"
#include "mfem.hpp"
#include <memory>

/**
 * Wrapper for mfem::HypreBoomerAMG solver.
 */
class MFEMHypreBoomerAMG : public MFEMSolverBase
{
public:
  static InputParameters validParams();

  MFEMHypreBoomerAMG(const InputParameters &);

  /// Returns a shared pointer to the instance of the Solver derived-class.
  std::shared_ptr<mfem::Solver> getSolver() const override { return _solver; }

protected:
  void constructSolver(const InputParameters & parameters) override;

private:
  std::shared_ptr<mfem::ParFiniteElementSpace> _mfem_fespace{nullptr};
  std::shared_ptr<mfem::HypreBoomerAMG> _solver{nullptr};
};
