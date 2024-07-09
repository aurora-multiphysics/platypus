#pragma once
#include "problem_operator.h"
#include "problem_builder_base.h"
namespace platypus
{

/// Class for steady-state problems with no equation system.
class SteadyStateProblem : public Problem
{
};

class SteadyStateProblemBuilder : public ProblemBuilder
{
public:
  SteadyStateProblemBuilder() : ProblemBuilder(new platypus::SteadyStateProblem) {}

  ~SteadyStateProblemBuilder() override = default;

  void RegisterGridFunctions() override {}

  void SetOperatorGridFunctions() override;

  void ConstructOperator() override;

  void ConstructState() override;

  void ConstructTimestepper() override {}

protected:
  // NB: constructor for derived classes.
  SteadyStateProblemBuilder(platypus::SteadyStateProblem * problem) : ProblemBuilder(problem) {}

  [[nodiscard]] platypus::ProblemOperator & GetOperator() const
  {
    return static_cast<platypus::ProblemOperator &>(*_problem_operator);
  }

  [[nodiscard]] platypus::SteadyStateProblem * GetProblem() const override
  {
    return ProblemBuilder::GetProblem<platypus::SteadyStateProblem>();
  }
};

} // namespace platypus
