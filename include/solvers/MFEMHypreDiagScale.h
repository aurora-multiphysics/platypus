#pragma once
#include "MFEMSolverBase.h"
#include "mfem.hpp"
#include <memory>

/**
 * Wrapper for mfem::HypreBoomerAMG solver.
 */
 class MFEMHypreDiagScale : public MFEMSolverBase
 {
 public:
   static InputParameters validParams();
 
   MFEMHypreDiagScale(const InputParameters &);
 
   /// Returns a shared pointer to the instance of the Solver derived-class.
   std::shared_ptr<mfem::Solver> getSolver() const override { return _solver; }
 
 protected:
   void constructSolver(const InputParameters & parameters) override;
 
 private:
   std::shared_ptr<mfem::ParFiniteElementSpace> _mfem_fespace{nullptr};
   mfem::real_t _strength_threshold;
   std::shared_ptr<mfem::HypreDiagScale> _solver{nullptr};
 };
 

