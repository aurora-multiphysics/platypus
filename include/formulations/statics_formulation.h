#pragma once
#include "../common/pfem_extras.hpp"
#include "formulation.h"
#include "inputs.h"
#include "sources.h"

namespace platypus
{

class StaticsFormulation : public SteadyStateEMFormulation
{
public:
  StaticsFormulation(std::string alpha_coef_name, std::string h_curl_var_name);

  ~StaticsFormulation() override = default;

  void ConstructJacobianPreconditioner() override;

  void ConstructJacobianSolver() override;

  void ConstructOperator() override;

  void RegisterGridFunctions() override;

  void RegisterCoefficients() override;

protected:
  const std::string _alpha_coef_name;
  const std::string _h_curl_var_name;
};

class StaticsOperator : public ProblemOperator
{
public:
  StaticsOperator(platypus::Problem & problem,
                  std::string h_curl_var_name,
                  std::string stiffness_coef_name);

  ~StaticsOperator() override = default;

  void SetGridFunctions() override;
  void Init(mfem::Vector & X) override;
  void Solve(mfem::Vector & X) override;

private:
  std::string _h_curl_var_name, _stiffness_coef_name;

  mfem::Coefficient * _stiff_coef{nullptr}; // Stiffness Material Coefficient
};

} // namespace platypus
