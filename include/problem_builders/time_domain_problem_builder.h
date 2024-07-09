#pragma once
#include "problem_builder_base.h"
#include "time_domain_problem_operator.h"

namespace platypus
{
/// Problem-builder for TimeDomainProblem.
class TimeDomainProblemBuilder : public ProblemBuilder
{
public:
  TimeDomainProblemBuilder() : ProblemBuilder() {}

  ~TimeDomainProblemBuilder() override = default;

  static std::vector<mfem::ParGridFunction *>
  RegisterTimeDerivatives(std::vector<std::string> gridfunction_names,
                          platypus::GridFunctions & gridfunctions);

  void RegisterGridFunctions() override;

  void SetOperatorGridFunctions() override;

  void ConstructOperator() override;

  void ConstructState() override;

  void ConstructTimestepper() override;

protected:
  [[nodiscard]] platypus::TimeDomainProblemOperator & GetOperator() const
  {
    return static_cast<TimeDomainProblemOperator &>(*_problem_operator);
  }
};

} // namespace platypus
