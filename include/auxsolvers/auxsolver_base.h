#pragma once
#include "coefficients.h"
#include "gridfunctions.h"
#include "hephaestus_solvers.h"
#include "inputs.h"
#include "mfem.hpp"

// Specify classes that perform auxiliary calculations on GridFunctions or
// Coefficients.
namespace platypus
{

class AuxSolver;

class AuxSolver
{
public:
  AuxSolver() = default;

  // NB: require virtual destructor to avoid leaks.
  virtual ~AuxSolver() = default;

  virtual void Init(const platypus::GridFunctions & gridfunctions,
                    platypus::Coefficients & coefficients) = 0;

  virtual void Solve(double t = 0.0) = 0;

  // Set priority. Lower values are evaluated first.
  void SetPriority(const int priority) { _priority = priority; };

  [[nodiscard]] inline int Priority() const { return _priority; }

  static bool PriorityComparator(
      const std::pair<std::shared_ptr<platypus::AuxSolver>, std::string> & first_comp,
      const std::pair<std::shared_ptr<platypus::AuxSolver>, std::string> & second_comp)
  {
    return (first_comp.first->Priority() < second_comp.first->Priority());
  }

private:
  int _priority{0};
};

} // namespace platypus
