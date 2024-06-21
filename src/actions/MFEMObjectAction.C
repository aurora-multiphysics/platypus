#include "MFEMObjectAction.h"

InputParameters
MFEMObjectAction::validParams()
{
  return MooseObjectAction::validParams();
}

MFEMObjectAction::MFEMObjectAction(const InputParameters & parameters)
  : MooseObjectAction(parameters), _mfem_problem(static_cast<MFEMProblem &>(*_problem))
{
}
