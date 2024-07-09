#pragma once
#include "../common/pfem_extras.hpp"
#include "problem_builder_base.h"
#include "problem_operator_interface.h"
#include "problem_operator_base.h"

namespace platypus
{

std::string GetTimeDerivativeName(const std::string & name);

std::vector<std::string> GetTimeDerivativeNames(std::vector<std::string> gridfunction_names);

/// Problem operator for time-dependent problems with no equation system. The user will need to subclass this since the solve is not
/// implemented.
class TimeDomainProblemOperator : public mfem::TimeDependentOperator,
                                  public ProblemOperatorInterface,
                                  public ProblemOperatorBase
{
public:
  TimeDomainProblemOperator(platypus::Problem & problem) : ProblemOperatorInterface(problem) {}
  ~TimeDomainProblemOperator() override = default;

  void SetGridFunctions() override;

  void ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt) override {}
};

} // namespace platypus