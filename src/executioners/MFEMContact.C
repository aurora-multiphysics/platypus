#include "MFEMContact.h"
#include "MFEMProblem.h"

#include "axom/slic.hpp"
#include "tribol/interface/tribol.hpp"
#include "tribol/interface/mfem_tribol.hpp"

/*
  We should be able to inherit all of this
*/

registerMooseObject("PlatypusApp", MFEMContact);

InputParameters
MFEMContact::validParams()
{
  InputParameters params = MFEMExecutioner::validParams();
  params.addClassDescription("Executioner for steady state MFEM problems.");
  params.addParam<Real>("time", 0.0, "System time");
  return params;
}

MFEMContact::MFEMContact(const InputParameters & params)
  : MFEMExecutioner(params),
    _system_time(getParam<Real>("time")),
    _time_step(_mfem_problem.timeStep()),
    _time(_mfem_problem.time()),
    _output_iteration_number(0)
{
  std::cout << "Initialising the Executioner\n";
  _time = _system_time;
}


void
MFEMContact::constructProblemOperator()
{
  _problem_data._eqn_system = std::make_shared<platypus::EquationSystem>();
  auto problem_operator = std::make_unique<platypus::EquationSystemProblemOperator>(_problem_data);
  
  _problem_operator.reset();
  _problem_operator = std::move(problem_operator);
}


/*
  1. Perform all tribol initialisaiton first (even if it means doing some redundant stuff)
  2. That means we can register the pressure  grid function with the problem data early
  3. Before calling eqn system Init, (but after(!) registering the grid function), we can add pressure to test/trial var names
      - Otherwise the sizes will be all messed up

*/
void
MFEMContact::init()
{
  _mfem_problem.execute(EXEC_PRE_MULTIAPP_SETUP);

  _mfem_problem.initialSetup();

  initTribol();

  // we can add params to mark which are the ess_t_dofs etc

  // the coords nodal grid function is set in mfemmesh:66

  // need to make sure coupling scheme is set before doing this!
  
  // pressure needs to be registered as a grid function so we can get it into the block vector later
  mfem::ParGridFunction& pressure = tribol::getMfemPressure(0);
  
  // need a shared pointer here
  std::shared_ptr<mfem::ParGridFunction> pressure_ptr = std::make_shared<mfem::ParGridFunction>(pressure);
  if (!pressure_ptr) MFEM_ABORT("Creating pressure shared pointer failed");
  _problem_data._gridfunctions.Register("pressure", pressure_ptr);
  
  _problem_data._eqn_system->AddTestVariableNameIfMissing("pressure");
  _problem_data._eqn_system->AddTrialVariableNameIfMissing("pressure");

  // Set up initial conditions
  // makes and registers some pargrid functions
  _problem_data._eqn_system->Init(
    _problem_data._gridfunctions,
    _problem_data._fespaces,
    getParam<MooseEnum>("assembly_level").getEnum<mfem::AssemblyLevel>()
  );
    
  _problem_operator->SetGridFunctions();

  // this should set the boundary conditions!
  _problem_operator->Init(_problem_data._f);


}

void
MFEMContact::initTribol()
{
  // silence warning
  axom::slic::initialize();

  // Create a Tribol coupling scheme: defines contact surfaces and enforcement
  int coupling_scheme_id = 0;
  
  const int dimensions = _mfem_problem.mesh().dimension();
  tribol::initialize(dimensions, _problem_data._comm);

  /* WARNING - HARCODING BELOW */

  // While there is a single mfem ParMesh for this problem, Tribol
  // defines a mortar and a nonmortar contact mesh, each with a unique mesh ID.
  // The Tribol mesh IDs for each contact surface are defined here.
  int mesh1_id = 0;
  int mesh2_id = 1;

  mfem::ParGridFunction* coords = dynamic_cast<mfem::ParGridFunction*>( _mfem_problem.getProblemData()._pmesh->GetNodes() );
  if ( !coords ) MFEM_ABORT("Failed to cast from GridFunction to ParGridFunction");

  // take a reference to the pmesh
  mfem::ParMesh& pmesh = *(_mfem_problem.getProblemData()._pmesh);

  tribol::registerMfemCouplingScheme(
    coupling_scheme_id, mesh1_id, mesh2_id,
    pmesh, *coords, mortar_attrs, nonmortar_attrs,
    tribol::SURFACE_TO_SURFACE,
    tribol::NO_CASE,
    tribol::SINGLE_MORTAR,
    tribol::FRICTIONLESS,
    tribol::LAGRANGE_MULTIPLIER,
    tribol::BINNING_GRID
  );

  // Set Tribol options for Lagrange multiplier enforcement
  tribol::setLagrangeMultiplierOptions(
    coupling_scheme_id,
    tribol::ImplicitEvalMode::MORTAR_RESIDUAL_JACOBIAN
  );
 
  // #4: Update contact mesh decomposition so the on-rank Tribol meshes
  // coincide with the current configuration of the mesh. This must be called
  // before tribol::update().
  tribol::updateMfemParallelDecomposition();
  
  // #5: Update contact gaps, forces, and tangent stiffness contributions
  int cycle = 1;   // pseudo cycle
  mfem::real_t t = 1.0;  // pseudo time
  mfem::real_t dt = 1.0; // pseudo dt
  tribol::update(cycle, t, dt);

}

void
MFEMContact::execute()
{
  _time_step = 0;
  _time = _time_step;
  _mfem_problem.outputStep(EXEC_INITIAL);
  _time = _system_time;

  preExecute();

  _mfem_problem.advanceState();

  // first step in any steady state solve is always 1 (preserving backwards compatibility)
  _time_step = 1;
  _mfem_problem.timestepSetup();
  // // Solve equation system.
  // if (_mfem_problem.shouldSolve())
  // {
    // }
    
  _problem_operator->Solve(_problem_data._f);

  // all this is shutdown stuff

  // Update displacement and coords grid functions
  // mfem::ParFiniteElementSpace* fespace = _problem_data._fespaces.Get(fespaceName);
  
  // that only works because we should have placed the contact stuff in the 0,0 block
  // auto& displacement_true = _problem_operator->_true_x.GetBlock(0);
  // fespace->GetProlongationMatrix()->Mult(displacement_true, displacement);


}

/*
PROBLEM OPERATOR IS AN "ALIAS" FOR EQN SYSTEM

MOST IMPORTANT STUFF HAPPENS IN _problem_operator->Solve(_problem_data._f);
THIS CALLS FormLegacySystem

_problem_operator points to platypus::ProblemOperator
  seems like an alias to EquationSystemProblemOperator

we can 

*/



