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
  int dim               = pmesh.Dimension();
  int sdim              = pmesh.SpaceDimension();

  int order = 1;

  // construct problem operator
  auto & problem_data      = _mfem_problem->getProblemData();
  problem_data._eqn_system = std::make_shared<platypus::EquationSystem>();
  auto eqn_system          = problem_data._eqn_system;

  // FE space
  InputParameters fespace_params             = _factory.getValidParams("MFEMScalarFESpace");
  fespace_params.set<MooseEnum>("fec_type")  = "H1";
  fespace_params.set<MooseEnum>("fec_order") = "FIRST";
  _mfem_problem->addFESpace("MFEMScalarFESpace", "H1FESpace", fespace_params);

  // bottom BC
  InputParameters bc_bottom_params              = _factory.getValidParams("MFEMScalarDirichletBC");
  bc_bottom_params.set<std::string>("variable") = "diffused";
  bc_bottom_params.set<Real>("value")           = 1.;
  bc_bottom_params.set<std::vector<BoundaryName>>("boundary") = {"1"};
  _mfem_problem->addBoundaryCondition("MFEMScalarDirichletBC", "bottom", bc_bottom_params);
  
  // low terminal BC
  InputParameters bc_low_terminal_params              = _factory.getValidParams("MFEMScalarDirichletBC");
  bc_low_terminal_params.set<std::string>("variable") = "diffused";
  bc_low_terminal_params.set<Real>("value")           = 0.;
  bc_low_terminal_params.set<std::vector<BoundaryName>>("boundary") = {"2"};
  _mfem_problem->addBoundaryCondition("MFEMScalarDirichletBC", "low_terminal", bc_low_terminal_params);
  

  // diffusion coefficient material
  InputParameters coef_params = _factory.getValidParams("MFEMGenericConstantMaterial");
  coef_params.set<std::vector<std::string>>("prop_names") = {"diffusivity"};
  coef_params.set<std::vector<double>>("prop_values")     = {1.0};
  _mfem_problem->addMaterial("MFEMGenericConstantMaterial", "material1", coef_params);

  // diffusion kernel
  InputParameters kernel_params                   = _factory.getValidParams("MFEMDiffusionKernel");
  kernel_params.set<std::string>("variable")      = "diffused";
  kernel_params.set<std::string>("coefficient")   = "diffusivity";
  // auto & kernel                                   = addObject<MFEMDiffusionKernel>("MFEMDiffusionKernel", "kernel1", kernel_params);
  _mfem_problem->addKernel("MFEMDiffusionKernel", "kernel1", kernel_params);

  // construct precon
  InputParameters precon_params      = _factory.getValidParams("MFEMHypreBoomerAMG");
  precon_params.set<double>("l_tol") = 1e-7; // HypreBoomerAMG cannot set absolute tolerance
  _mfem_problem->addMFEMPreconditioner("MFEMHypreBoomerAMG", "precon1", precon_params);

  // construct solver
  InputParameters solver_params          = _factory.getValidParams("MFEMHypreGMRES");
  solver_params.set<double>("l_tol")     = 1e-7;
  solver_params.set<double>("l_abs_tol") = 1e-5;
  _mfem_problem->addMFEMSolver("MFEMHypreGMRES", "solver1", solver_params);
  

  // next, the variables
  InputParameters var_params                = _factory.getValidParams("MFEMVariable");
  var_params.set<UserObjectName>("fespace") = "H1FESpace";
  _mfem_problem->addVariable("MFEMVariable", "diffused", var_params);

  _mfem_problem->addMFEMNonlinearSolver();

  // TODO: rest of test
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
  // a.FormLinearSystem(ess_tdof_list, x, b, A, X, B); <- call this in equation system
  //
  mfem::BlockVector    X;

  eqn_system->Init(
    problem_data._gridfunctions,
    problem_data._fespaces,
    mfem::AssemblyLevel::LEGACY
  );

  // problem operator
  auto problem_operator = std::make_unique<platypus::EquationSystemProblemOperator>(problem_data);
  problem_operator->SetGridFunctions();
  problem_operator->Init( X );
  
  // solve!
  problem_operator->Solve( X );

  ASSERT_EQ(1+1,2);
}


