#pragma once
#include "problem_builder_base.h"
#include "time_domain_problem_builder.h"
#include "time_domain_equation_system_problem_operator.h"
#include "equation_system_interface.h"

namespace platypus
{

// Problem-builder for TimeDomainEquationSystemProblem.
class TimeDomainEquationSystemProblemBuilder : public TimeDomainProblemBuilder,
                                               public EquationSystemProblemBuilderInterface
{
public:
  TimeDomainEquationSystemProblemBuilder() = default;

  ~TimeDomainEquationSystemProblemBuilder() override = default;

  /// NB: - note use of final. Ensure that the equation system is initialized.
  void InitializeKernels() final;

  void ConstructOperator() override;

protected:
  [[nodiscard]] TimeDomainEquationSystemProblemOperator & GetOperator() const
  {
    return static_cast<TimeDomainEquationSystemProblemOperator &>(*_problem_operator);
  }

  [[nodiscard]] platypus::TimeDependentEquationSystem * GetEquationSystem() const override
  {
    return GetOperator().GetEquationSystem();
  }
};

} // namespace platypus
