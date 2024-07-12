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
  /// NB: set "_problem" member variable in parent class.
  TimeDomainEquationSystemProblemBuilder() : TimeDomainProblemBuilder() {}

  ~TimeDomainEquationSystemProblemBuilder() override = default;

  /// NB: - note use of final. Ensure that the equation system is initialized.
  void InitializeKernels() final;

  void ConstructOperator() override
  {
    auto equation_system = std::make_unique<platypus::TimeDependentEquationSystem>();
    auto problem_operator = std::make_shared<platypus::TimeDomainEquationSystemProblemOperator>(
        *_problem, std::move(equation_system));

    _problem_operator = std::move(problem_operator);
  }

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
