#pragma once
#include "problem_operator.h"
#include "problem_builder_base.h"
namespace platypus
{

class SteadyStateProblemBuilder : public ProblemBuilder
{
public:
  SteadyStateProblemBuilder() : ProblemBuilder(new platypus::MFEMProblemData) {}

  ~SteadyStateProblemBuilder() override = default;

  void RegisterGridFunctions() override {}

  void SetOperatorGridFunctions() override;

  void ConstructOperator() override;

  void ConstructState() override;

  void ConstructTimestepper() override {}

protected:
  // NB: constructor for derived classes.
  SteadyStateProblemBuilder(platypus::MFEMProblemData * problem) : ProblemBuilder(problem) {}

  [[nodiscard]] platypus::ProblemOperator & GetOperator() const
  {
    return static_cast<platypus::ProblemOperator &>(*_problem_operator);
  }
};

} // namespace platypus
