#include "MFEMObjectUnitTest.h"
#include "MFEMScalarDirichletBC.h"
#include "MFEMDiffusionKernel.h"
#include "MFEMHypreBoomerAMG.h"
#include "equation_system_problem_operator.h"
#include "MFEMParaViewDataCollection.h"
#include "MFEMSteady.h"

class MFEMMeshRefinementTest : public MFEMObjectUnitTest
{
public:
  MFEMMeshRefinementTest() : MFEMObjectUnitTest("PlatypusApp") {}
};

/// Test based on MFEM example 6
TEST_F(MFEMMeshRefinementTest, DiffusionRefinement)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // _mfem_mesh_ptr->buildMesh();
  mfem::ParMesh & pmesh = _mfem_mesh_ptr->getMFEMParMesh();
  int dim = pmesh.Dimension();
  int sdim = pmesh.SpaceDimension();

  int order = 1;

  // FE space
  InputParameters fespace_params = _factory.getValidParams("MFEMFESpace");
  fespace_params.set<MooseEnum>("fec_type") = "H1";
  fespace_params.set<MooseEnum>("fec_order") = "FIRST";
  _mfem_problem->addFESpace("MFEMFESpace", "H1FESpace", fespace_params);

  // bottom BC
  InputParameters bc_bottom_params = _factory.getValidParams("MFEMScalarDirichletBC");
  bc_bottom_params.set<std::string>("variable") = "diffused";
  bc_bottom_params.set<Real>("value") = 1.;
  bc_bottom_params.set<std::vector<BoundaryName>>("boundary") = {"1"};
  auto & bc_bottom =
      addObject<MFEMScalarDirichletBC>("MFEMScalarDirichletBC", "bottom", bc_bottom_params);

  // low terminal BC
  InputParameters bc_low_terminal_params = _factory.getValidParams("MFEMScalarDirichletBC");
  bc_low_terminal_params.set<std::string>("variable") = "diffused";
  bc_low_terminal_params.set<Real>("value") = 0.;
  bc_low_terminal_params.set<std::vector<BoundaryName>>("boundary") = {"2"};
  auto & bc_low_terminal = addObject<MFEMScalarDirichletBC>(
      "MFEMScalarDirichletBC", "low_terminal", bc_low_terminal_params);

  // diffusion coefficient material
  InputParameters coef_params = _factory.getValidParams("MFEMGenericConstantMaterial");
  coef_params.set<std::vector<std::string>>("prop_names") = {"diffusivity"};
  coef_params.set<std::vector<double>>("prop_values") = {1.0};
  _mfem_problem->addMaterial("MFEMGenericConstantMaterial", "material1", coef_params);

  // diffusion kernel
  InputParameters kernel_params = _factory.getValidParams("MFEMDiffusionKernel");
  kernel_params.set<std::string>("variable") = "diffused";
  kernel_params.set<std::string>("coefficient") = "diffusivity";
  auto & kernel = addObject<MFEMDiffusionKernel>("MFEMDiffusionKernel", "kernel1", kernel_params);

  // construct solver
  InputParameters solver_params = _factory.getValidParams("MFEMHypreBoomerAMG");
  solver_params.set<double>("l_tol") = 1e-7; // HypreBoomerAMG cannot set absolute tolerance
  auto & solver = addObject<MFEMHypreBoomerAMG>("MFEMHypreBoomerAMG", "solver1", solver_params);

  // construct problem operator
  auto & problem_data = _mfem_problem->getProblemData();
  problem_data._eqn_system = std::make_shared<platypus::EquationSystem>();
  auto problem_operator = std::make_unique<platypus::EquationSystemProblemOperator>(problem_data);
  problem_operator.reset();

  // construct data collection
  InputParameters data_collection_params = _factory.getValidParams("MFEMParaViewDataCollection");
  data_collection_params.set<std::string>("file_base") = "/dev/null";
  _mfem_problem->addOutput(
      "MFEMParaViewDataCollection", "ParaviewDataCollection", data_collection_params);

  // construct MFEMExecutioner
  InputParameters executioner_params = _factory.getValidParams("MFEMSteady");
  executioner_params.set<std::string>("device") = "cpu";
  std::shared_ptr<Executioner> executioner =
      _factory.create<Executioner>("MFEMSteady", "MFEMSteady", executioner_params);

  // setup MOOSE quadrature rules
  _mfem_problem->createQRules(QuadratureType::QGAUSS, Order::INVALID_ORDER);
  // call prepare to setup things
  // done in SetupMeshCompleteAction normally
  _mfem_mesh_ptr->prepare(nullptr);

  // init equation system
  _mfem_problem->execute(EXEC_PRE_MULTIAPP_SETUP);
  _mfem_problem->initialSetup();

  // set up initial conditions
  problem_data._eqn_system->Init(problem_data._gridfunctions,
                                 problem_data._fespaces,
                                 problem_data._bc_map,
                                 mfem::AssemblyLevel::FULL);

  problem_operator->SetGridFunctions();
  problem_operator->Init(problem_data._f);

  // execute
  _mfem_problem->advanceState();

  _mfem_problem->timestepSetup();

  // Solve equation system.
  problem_operator->Solve(problem_data._f);

  // Displace mesh, if required
  _mfem_problem->displaceMesh();

  _mfem_problem->computeIndicators();
  _mfem_problem->computeMarkers();

  // TODO: rest of test
  EXPECT_EQ(1 + 1, 2);
  ////////////////////////////////
  ////////////////////////////////
  // mfem::Array<int> ess_tdof_list, ess_bdr;
  // if (pmesh.bdr_attributes.Size())
  // {
  //   ess_bdr.SetSize(pmesh.bdr_attributes.Max());
  //   ess_bdr = 1;
  //   fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
  // }
  // mfem::FunctionCoefficient f(fexact);
  // mfem::ParLinearForm b(&fespace);
  // b.AddDomainIntegrator(new mfem::DomainLFIntegrator(f));
  // b.Assemble();
  //
  // mfem::ParBilinearForm a(&fespace);
  // mfem::ConstantCoefficient one(1.0);
  // a.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
  // a.Assemble();
  //
  // mfem::ParGridFunction x(&fespace);
  // mfem::FunctionCoefficient uex(uexact);
  // x = 0.0;
  // x.ProjectBdrCoefficient(uex, ess_bdr);
  //
  // mfem::OperatorPtr A;
  // mfem::Vector B, X;
  // a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
  //
  // solver.SetOperator(*A);
  // solver.Mult(B, X);
  //
  // mfem::Vector Y(X.Size());
  // A->Mult(X, Y);
  // Y -= B;
  // ASSERT_LE(Y.Norml2(), tol);
}
