#pragma once
#include "equation_system.h"
#include "problem_operator_base.h"
#include "gridfunctions.h"
#include "inputs.h"
#include <fstream>
#include <iostream>
#include <memory>

namespace platypus
{

/// Base problem class.
class Problem
{
public:
  Problem() = default;
  virtual ~Problem();

  std::shared_ptr<mfem::ParMesh> _pmesh{nullptr};
  platypus::BCMap _bc_map;
  platypus::Coefficients _coefficients;
  platypus::Outputs _outputs;

  std::unique_ptr<mfem::ODESolver> _ode_solver{nullptr};
  std::unique_ptr<mfem::BlockVector> _f{nullptr};

  std::shared_ptr<mfem::Solver> _jacobian_preconditioner{nullptr};
  std::shared_ptr<mfem::Solver> _jacobian_solver{nullptr};
  std::shared_ptr<mfem::NewtonSolver> _nonlinear_solver{nullptr};

  platypus::FECollections _fecs;
  platypus::FESpaces _fespaces;
  platypus::GridFunctions _gridfunctions;

  /// Returns a pointer to the operator. See derived classes.
  [[nodiscard]] virtual mfem::Operator * GetOperator() const = 0;
};

/// ProblemBuilder base class.
class ProblemBuilder
{
public:
  /// NB: delete empty constructor to allow only derived classes to be constructed.
  ProblemBuilder() = delete;

  // Virtual destructor required to prevent leaks.
  virtual ~ProblemBuilder() = default;

  virtual void RegisterGridFunctions() = 0;

  virtual void SetOperatorGridFunctions() = 0;
  virtual void ConstructNonlinearSolver();
  virtual void ConstructOperator() = 0;
  virtual void ConstructState() = 0;
  virtual void ConstructTimestepper() = 0;

  virtual void InitializeKernels();

  void InitializeOutputs();

  /// @brief Call @a FinalizeProblem to setup a problem. The operator should
  /// have already been constructed by calling ConstructOperator().
  void FinalizeProblem();

  /// Returns a shared pointer to the problem.
  [[nodiscard]] std::shared_ptr<platypus::Problem> ReturnProblem() const { return _problem; }

  /// Returns a shared pointer to the problem operator.
  [[nodiscard]] std::shared_ptr<platypus::ProblemOperatorBase> ReturnOperator() const
  {
    return _problem_operator;
  }

protected:
  /// Protected constructor. Derived classes must call this constructor.
  ProblemBuilder(platypus::Problem * problem)
    : _problem(std::shared_ptr<platypus::Problem>(problem))
  {
  }

  /// Overridden in derived classes.
  [[nodiscard]] virtual platypus::Problem * GetProblem() const = 0;

  /// Helper template getter with safety check.
  template <class TDerivedProblem>
  [[nodiscard]] TDerivedProblem * GetProblem() const
  {
    if (!_problem)
    {
      MFEM_ABORT("platypus::Problem instance is NULL.");
    }

    return static_cast<TDerivedProblem *>(_problem.get());
  }

  /// Coefficient used in some derived classes.
  mfem::ConstantCoefficient _one_coef{1.0};

private:
  std::shared_ptr<platypus::Problem> _problem{nullptr};
  std::shared_ptr<platypus::ProblemOperatorBase> _problem_operator{nullptr};
};

/// Interface for problem builders that are constructing problems with an equation system.
class EquationSystemProblemBuilderInterface
{
public:
  EquationSystemProblemBuilderInterface() = default;
  virtual ~EquationSystemProblemBuilderInterface() = default;

  /// Add a kernel to the problem's equation system.
  template <class T>
  void AddKernel(std::string var_name, std::shared_ptr<platypus::Kernel<T>> kernel)
  {
    GetEquationSystem()->AddTrialVariableNameIfMissing(var_name);
    GetEquationSystem()->AddKernel(var_name, std::move(kernel));
  }

protected:
  /// Implemented in derived classes. Returns a pointer to the problem operator's equation system.
  [[nodiscard]] virtual platypus::EquationSystem * GetEquationSystem() const = 0;
};

} // namespace platypus