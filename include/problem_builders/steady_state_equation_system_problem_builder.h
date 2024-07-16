#pragma once
#include "equation_system_problem_operator.h"
#include "problem_builder_base.h"
#include "steady_state_problem_builder.h"
#include "equation_system_interface.h"

namespace platypus
{

/// Problem-builder for SteadyStateEquationSystemProblem.
class SteadyStateEquationSystemProblemBuilder : public SteadyStateProblemBuilder,
                                                public EquationSystemProblemBuilderInterface
{
public:
  SteadyStateEquationSystemProblemBuilder() = default;

  ~SteadyStateEquationSystemProblemBuilder() override = default;

  /// NB: use of final! This calls ProblemBuilder::InitializeKernels and also ensures that the
  /// equation system is initialized.
  void InitializeKernels() final;

  void ConstructOperator() override;

protected:
  [[nodiscard]] platypus::EquationSystemProblemOperator & GetOperator() const
  {
    return static_cast<platypus::EquationSystemProblemOperator &>(*_problem_operator);
  }

  [[nodiscard]] platypus::EquationSystem * GetEquationSystem() const override
  {
    return GetOperator().GetEquationSystem();
  }
};

} // namespace platypus
