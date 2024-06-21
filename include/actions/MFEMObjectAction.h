#pragma once

#include "MFEMProblem.h"

// MOOSE includes
#include "MooseObjectAction.h"

/*
 * This class adds a getMFEMProblem method.
 */
class MFEMObjectAction : public MooseObjectAction
{
public:
  static InputParameters validParams();

  MFEMObjectAction(const InputParameters & parameters);

  /// Returns a reference to the MFEMProblem instance.
  MFEMProblem & getMFEMProblem() const { return _mfem_problem; }

private:
  MFEMProblem & _mfem_problem;
};