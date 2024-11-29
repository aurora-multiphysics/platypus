#include "MFEMScalarIC.h"
#include "MFEMProblem.h"
#include <libmesh/libmesh_common.h>
#include <mfem.hpp>

registerMooseObject("PlatypusApp", MFEMScalarIC);

InputParameters
MFEMScalarIC::validParams()
{
  auto params = MFEMGeneralUserObject::validParams();
  params.addRequiredParam<std::string>("variable",
                                       "The variable to apply the initial condition for");
  params.addRequiredParam<std::string>("coefficient", "The scalar coefficient");
  params.registerBase("InitialCondition");
  // We cannot generally execute this at construction time since the coefficient may depend on, for
  // instance, a MOOSE function which is not itself setup until its initialSetup is called.
  // UserObject initial execution occurs after function initialSetup
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL};
  params.suppressParameter<ExecFlagEnum>("execute_on");
  return params;
}

MFEMScalarIC::MFEMScalarIC(const InputParameters & params) : MFEMGeneralUserObject(params) {}

void
MFEMScalarIC::execute()
{
  auto coeff = getMFEMProblem().getCoefficient(getParam<std::string>("coefficient"));
  auto grid_function = getMFEMProblem().getGridFunction(getParam<std::string>("variable"));
  grid_function->ProjectCoefficient(*coeff);
}
