#pragma once
#include "MFEMExecutioner.h"
#include "equation_system_problem_operator.h"

class MFEMContact : public MFEMExecutioner
{
public:
  static InputParameters validParams();

  MFEMContact() = default;
  explicit MFEMContact(const InputParameters & params);
  ~MFEMContact() override = default;

  void constructProblemOperator() override;
  virtual void init() override;
  virtual void execute() override;

  void initTribol();

protected:
  // Time variables used for consistency with MOOSE, needed for outputs.
  // Important for future synchronisation of solves in MultiApps
  Real   _system_time;
  int &  _time_step;
  Real & _time;

  /// Iteration number obtained from the main application
  unsigned int _output_iteration_number;

private:
  bool _last_solve_converged;
  std::unique_ptr<platypus::ProblemOperator> _problem_operator{nullptr};

  // Stuff that we need for the problem setup that we copy from the test case
  std::set<int>  mortar_attrs    = {4};
  std::set<int>  nonmortar_attrs = {5};
};