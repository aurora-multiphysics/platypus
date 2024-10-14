#include "MFEMFunctionCoefficient.h"

registerMooseObject("PlatypusApp", MFEMFunctionCoefficient);

InputParameters
MFEMFunctionCoefficient::validParams()
{
  InputParameters params = MFEMCoefficient::validParams();
  params.addClassDescription(
      "Class for defining an mfem::FunctionCoefficient object to add to an MFEMProblem.");
  params.addParam<FunctionName>("function", 0, "The function to associated with the Dirichlet BC");
  return params;
}

MFEMFunctionCoefficient::MFEMFunctionCoefficient(const InputParameters & parameters)
  : MFEMCoefficient(parameters),
    _func(getFunction("function")),
    _coefficient(std::make_shared<mfem::FunctionCoefficient>(
        [&](const mfem::Vector & p, double t) { return _func.value(t, PointFromMFEMVector(p)); }))
{
}

MFEMFunctionCoefficient::~MFEMFunctionCoefficient() {}
