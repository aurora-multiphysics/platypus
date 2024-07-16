#pragma once
#include "executioner_base.h"
#include "steady_state_problem_builder.h"
#include "problem_operator.h"

namespace platypus
{

class SteadyExecutioner : public Executioner
{
public:
  SteadyExecutioner() = default;
  explicit SteadyExecutioner(const platypus::InputParameters & params);

  void Solve() const override;

  void Execute() const override;

private:
  platypus::MFEMProblemData * _problem{nullptr};
  platypus::ProblemOperator * _problem_operator{nullptr};
};

} // namespace platypus
