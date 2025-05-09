#pragma once
#include "problem_operator.h"
#include "problem_operator_interface.h"
#include "equation_system_interface.h"

namespace platypus
{
/// Steady-state problem operator with an equation system.
class EquationSystemProblemOperator : public ProblemOperator, public EquationSystemInterface
{
public:
  EquationSystemProblemOperator(MFEMProblemData & problem)
    : ProblemOperator(problem), _equation_system(problem._eqn_system)
  {
  }

  void SetGridFunctions() override;
  void Init(mfem::BlockVector & X) override;
  virtual void Solve(mfem::Vector & X) override;
  void SetUpAMR();
  bool HRefine();

  ~EquationSystemProblemOperator() override = default;

  [[nodiscard]] platypus::EquationSystem * GetEquationSystem() const override
  {
    if (!_equation_system)
    {
      MFEM_ABORT("No equation system has been added to ProblemOperator.");
    }

    return _equation_system.get();
  }

private:
  std::shared_ptr<platypus::EquationSystem>     _equation_system{nullptr};
  std::shared_ptr<mfem::ErrorEstimator>         _error_estimator;
  std::unique_ptr<mfem::ThresholdRefiner>       _refiner;
  std::unique_ptr<mfem::H1_FECollection>        _smooth_flux_fec;
  std::unique_ptr<mfem::L2_FECollection>        _flux_fec;
  std::unique_ptr<mfem::ParFiniteElementSpace>  _smooth_flux_fes;
  std::unique_ptr<mfem::ParFiniteElementSpace>  _flux_fes;
};

} // namespace platypus