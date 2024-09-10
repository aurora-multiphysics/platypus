#include "MFEMTransfer.h"
#include "MFEMProblem.h"

InputParameters
MFEMTransfer::validParams()
{
  return Transfer::validParams();
}

MFEMTransfer::MFEMTransfer(const InputParameters & parameters)
  : Transfer(parameters), _mfem_problem(static_cast<MFEMProblem &>(_fe_problem))
{
}
