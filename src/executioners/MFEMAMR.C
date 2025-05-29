#include "MFEMAMR.h"
#include "MFEMProblem.h"
#include "MFEMEstimator.h"


registerMooseObject("PlatypusApp", MFEMAMR);

/**
  TODO: Support multiple aux kernels in a block

  Note the way that estimators are set up - by fetching a single grid function
  via a single variable/test_var_name. We may need to modify this class to support
  estimators on multiple variables.
*/

InputParameters
MFEMAMR::validParams()
{
  InputParameters params = MFEMExecutioner::validParams();
  params.addClassDescription("Executioner for steady state MFEM problems with AMR.");
  params.addParam<Real>("time", 0.0, "System time");

  params.addRequiredParam<std::string>("variable", "Variable to perform amr with");
  params.addRequiredParam<std::string>("kernel", "Kernel to perform amr with");
  params.addRequiredParam<std::string>("fe_space", "FESpace to set order from");

  return params;
}


MFEMAMR::MFEMAMR(const InputParameters & params)
  : MFEMExecutioner(params),
    _system_time(getParam<Real>("time")),
    _time_step(_mfem_problem.timeStep()),
    _time(_mfem_problem.time()),
    _output_iteration_number(0),
    _test_var_name(getParam<std::string>("variable")),
    _kernel_name(getParam<std::string>("kernel")),
    _fe_space_name(getParam<std::string>("fe_space"))
{
  _time = _system_time;
}


void
MFEMAMR::constructProblemOperator()
{
  _problem_data._eqn_system = std::make_shared<platypus::EquationSystem>();
  auto problem_operator = std::make_unique<platypus::EquationSystemProblemOperator>(_problem_data);

  _problem_operator.reset();
  _problem_operator = std::move(problem_operator);

}


void
MFEMAMR::init()
{
  _mfem_problem.execute(EXEC_PRE_MULTIAPP_SETUP);
  _mfem_problem.initialSetup();
  
  // Set up initial conditions
  _problem_data._eqn_system->Init(
  _problem_data._gridfunctions,
  _problem_data._fespaces,
  getParam<MooseEnum>("assembly_level").getEnum<mfem::AssemblyLevel>());
  
  _problem_operator->SetGridFunctions();
  _problem_operator->Init(_problem_data._f);

  // remove the hardcoding here
  InputParameters estimator_params = _factory.getValidParams("MFEMZienkiewiczZhuEstimator");
  
  // bad practice - this prevent hardcoding but it is not nice
  estimator_params.addParam<std::string>("variable", _test_var_name, "Variable to perform amr with");
  estimator_params.addParam<std::string>("kernel", _kernel_name, "Kernel to perform amr with");
  estimator_params.addParam<std::string>("fe_space", _fe_space_name, "FESpace to set order from");
  
  _mfem_problem.addUserObject( "MFEMZienkiewiczZhuEstimator", "MFEMZienkiewiczZhuEstimator", estimator_params );

  const UserObject * est_uo = &(_mfem_problem.getUserObjectBase("MFEMZienkiewiczZhuEstimator"));
  if ( dynamic_cast<const MFEMEstimator *>(est_uo) != nullptr )
  {
    std::shared_ptr<MooseObject> object_ptr  = _mfem_problem.getUserObject<MFEMEstimator>("MFEMZienkiewiczZhuEstimator").getSharedPtr();
    std::shared_ptr<MFEMEstimator> estimator = std::dynamic_pointer_cast<MFEMEstimator>(object_ptr);

    // construct the estimator
    _problem_operator->AddEstimator(estimator);
  }
  
  else
  {
    mooseError("Cannot add estimator :()");
  }
}


void
MFEMAMR::execute()
{
  _problem_operator->SetUpAMR();
  if (_app.isRecovering())
  {
    _console << "\nCannot recover steady solves!\nExiting...\n" << std::endl;
    return;
  }

  _time_step = 0;
  _time = _time_step;
  _mfem_problem.outputStep(EXEC_INITIAL);
  _time = _system_time;

  preExecute();

  _mfem_problem.advanceState();

  // first step in any steady state solve is always 1 (preserving backwards compatibility)
  _time_step = 1;
  _mfem_problem.timestepSetup();

  // Solve equation system.
  if (_mfem_problem.shouldSolve())
  {
    _problem_operator->Solve(_problem_data._f);

    _problem_operator->HRefine();
    _problem_operator->UpdateAfterRefinement();

    _problem_operator->Solve(_problem_data._f);
  }

  // Displace mesh, if required
  _mfem_problem.displaceMesh();

  _mfem_problem.computeIndicators();
  _mfem_problem.computeMarkers();

  // need to keep _time in sync with _time_step to get correct output
  _time = _time_step;
  // Execute user objects at timestep end
  _mfem_problem.execute(EXEC_TIMESTEP_END);
  _mfem_problem.outputStep(EXEC_TIMESTEP_END);
  _time = _system_time;

  {
    TIME_SECTION("final", 1, "Executing Final Objects")
    _mfem_problem.execMultiApps(EXEC_FINAL);
    _mfem_problem.finalizeMultiApps();
    _mfem_problem.postExecute();
    _mfem_problem.execute(EXEC_FINAL);
    _time = _time_step;
    _mfem_problem.outputStep(EXEC_FINAL);
    _time = _system_time;
  }

  postExecute();

}
