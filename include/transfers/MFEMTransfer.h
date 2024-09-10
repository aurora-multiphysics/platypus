#pragma once

// MOOSE includes
#include "Transfer.h"

// Forwards declaration.
class MFEMProblem;

/*
 * This class adds a getMFEMProblem method.
 */
class MFEMTransfer : public Transfer
{
public:
  static InputParameters validParams();

  MFEMTransfer(const InputParameters & parameters);

  /// Returns a reference to the MFEMProblem instance.
  MFEMProblem & getMFEMProblem() const { return _mfem_problem; }

private:
  MFEMProblem & _mfem_problem;
};