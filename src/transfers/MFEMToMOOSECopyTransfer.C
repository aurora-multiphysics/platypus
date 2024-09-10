#include "MFEMToMOOSECopyTransfer.h"
#include "MFEMProblem.h"

registerMooseObject("PlatypusApp", MFEMToMOOSECopyTransfer);

InputParameters
MFEMToMOOSECopyTransfer::validParams()
{
  InputParameters params = MFEMTransfer::validParams();
  params.addParam<MultiAppName>("from_multi_app", "The name of the MultiApp to receive data from");
  params.addParam<MultiAppName>("to_multi_app", "The name of the MultiApp to transfer the data to");

  // MultiAppTransfers by default will execute with their associated MultiApp. These flags will be
  // added by FEProblemBase when the transfer is added.
  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum.addAvailableFlags(EXEC_SAME_AS_MULTIAPP);
  exec_enum = EXEC_SAME_AS_MULTIAPP;
  params.setDocString("execute_on", exec_enum.getDocString());

  params.addParam<bool>(
      "check_multiapp_execute_on",
      true,
      "When false the check between the multiapp and transfer execute on flags is not performed.");
  params.addParam<bool>("displaced_source_mesh",
                        false,
                        "Whether or not to use the displaced mesh for the source mesh.");
  params.addParam<bool>("displaced_target_mesh",
                        false,
                        "Whether or not to use the displaced mesh for the target mesh.");

  params.addRequiredParam<std::vector<AuxVariableName>>(
      "variable", "The auxiliary variable to store the transferred values in.");
  params.addRequiredParam<std::vector<VariableName>>("source_variable",
                                                     "The variable to transfer from.");

  params.addClassDescription(
      "Copies variables (nonlinear and auxiliary) between multiapps that have identical meshes.");
  return params;
}

MFEMToMOOSECopyTransfer::MFEMToMOOSECopyTransfer(const InputParameters & parameters)
  : MFEMTransfer(parameters)
{
  // Get the multiapps from their names
  if (!isParamValid("multi_app"))
  {
    if (isParamValid("from_multi_app"))
    {
      _from_multi_app = getMFEMProblem().getMultiApp(getParam<MultiAppName>("from_multi_app"));
    }
    if (isParamValid("to_multi_app"))
    {
      _to_multi_app = getMFEMProblem().getMultiApp(getParam<MultiAppName>("to_multi_app"));
    }
    if (!isParamValid("direction") && !isParamValid("from_multi_app") &&
        !isParamValid("to_multi_app"))
      mooseError("from_multi_app and/or to_multi_app must be specified");
  }

  // if (getParam<bool>("check_multiapp_execute_on"))
  //   checkMultiAppExecuteOn();

  // Fill direction attributes, for backward compatibility but also convenience
  if (!isParamValid("direction"))
  {
    if (_from_multi_app && (!_to_multi_app || _from_multi_app == _to_multi_app))
      _directions.setAdditionalValue("from_multiapp");
    if (_to_multi_app && (!_from_multi_app || _from_multi_app == _to_multi_app))
      _directions.setAdditionalValue("to_multiapp");
    if (_from_multi_app && _to_multi_app && _from_multi_app != _to_multi_app)
      _directions.setAdditionalValue("between_multiapp");

    // So it's available in the next constructors
    _direction = _directions[0];
    _current_direction = _directions[0];
  }

  // Handle deprecated parameters
  if (parameters.isParamSetByUser("direction"))
  {
    if (!isParamValid("multi_app"))
      paramError("direction",
                 "The deprecated 'direction' parameter is meant to be used in conjunction with the "
                 "'multi_app' parameter");
    if (isParamValid("to_multi_app") || isParamValid("from_multi_app"))
      paramError("direction",
                 "The deprecated 'direction' parameter is not meant to be used in conjunction with "
                 "the 'from_multi_app' or 'to_multi_app' parameters");
  }
}

void
MFEMToMOOSECopyTransfer::initialSetup()
{
  // // Check for siblings transfer support
  // if (_to_multi_app && _from_multi_app)
  //   checkSiblingsTransferSupported();

  // getAppInfo();

  // if (_from_multi_app)
  //   _from_multi_app->addAssociatedTransfer(*this);
  // if (_to_multi_app)
  //   _to_multi_app->addAssociatedTransfer(*this);

  const FEProblemBase * from_problem;
  const FEProblemBase * to_problem;

  if (_current_direction == FROM_MULTIAPP)
  {
    // Subdomain and variable type information is shared on all subapps
    from_problem = &getFromMultiApp()->appProblemBase(getFromMultiApp()->firstLocalApp());
    to_problem = &getFromMultiApp()->problemBase();
  }
  else if (_current_direction == TO_MULTIAPP)
  {
    from_problem = &getToMultiApp()->problemBase();
    to_problem = &getToMultiApp()->appProblemBase(getToMultiApp()->firstLocalApp());
  }
  else
  {
    from_problem = &getFromMultiApp()->appProblemBase(getFromMultiApp()->firstLocalApp());
    to_problem = &getToMultiApp()->appProblemBase(getToMultiApp()->firstLocalApp());
  }
}

void
MFEMToMOOSECopyTransfer::execute()
{
  // TIME_SECTION("MFEMToMOOSECopyTransfer::execute()", 5, "Copies variables");

  if (_current_direction == TO_MULTIAPP)
  {
    FEProblemBase & from_problem = getToMultiApp()->problemBase();
    for (unsigned int i = 0; i < getToMultiApp()->numGlobalApps(); i++)
      if (getToMultiApp()->hasLocalApp(i))
        transfer(getToMultiApp()->appProblemBase(i), from_problem);
  }

  else if (_current_direction == FROM_MULTIAPP)
  {
    FEProblemBase & to_problem = getFromMultiApp()->problemBase();
    for (unsigned int i = 0; i < getFromMultiApp()->numGlobalApps(); i++)
      if (getFromMultiApp()->hasLocalApp(i))
        transfer(to_problem, getFromMultiApp()->appProblemBase(i));
  }
  // else if (_current_direction == BETWEEN_MULTIAPP)
  // {
  //   int transfers_done = 0;
  //   for (unsigned int i = 0; i < getFromMultiApp()->numGlobalApps(); i++)
  //   {
  //     if (getFromMultiApp()->hasLocalApp(i))
  //     {
  //       if (getToMultiApp()->hasLocalApp(i))
  //       {
  //         transfer(getToMultiApp()->appProblemBase(i), getFromMultiApp()->appProblemBase(i));
  //         transfers_done++;
  //       }
  //     }
  //   }
  //   if (!transfers_done)
  //     mooseError("BETWEEN_MULTIAPP transfer not supported if there is not at least one subapp "
  //                "per multiapp involved on each rank");
  // }
}

void
MFEMToMOOSECopyTransfer::transfer(FEProblemBase & to_problem, FEProblemBase & from_problem)
{
}