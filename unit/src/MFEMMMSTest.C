// largely copying the mfemmeshtest here!

/*

  Could also try to inherit from the MFEMSolverTest here, but make
  sure to bypass its constructor
*/

#include "MFEMObjectUnitTest.h" // should have all the headers i need
#include "MFEMMesh.h"

#include "MFEMDiffusionKernel.h"
#include "MFEMScalarDirichletBC.h"
#include "MFEMHypreBoomerAMG.h"

/**
 * Performs a MMS test on deformed(?) MFEM mesh.
 * Inspiration: https://mooseframework.inl.gov/python/mms.html
 * 
 * We assume an exact solution of sin(2*pi*x) * sin(2*pi*y)
 * and an exact forcing function of 8pi^2 u
 */
class MFEMMMSTest : public MFEMObjectUnitTest
{
public:
  MFEMMMSTest(): MFEMObjectUnitTest("PlatypusApp", "data/mug.e") {}

  void buildMFEMMesh(MeshFileName filename, int serial_ref = 0, int parallel_ref = 0);
  void testDiffusionSolve(mfem::Solver & solver, mfem::real_t tol);

  // void SetUp() override;

  /**
   * Computes f based on given x vector
   */
  static double fexact(const mfem::Vector& x)
  {
    return 8 * M_PI * M_PI * sin( 2* M_PI * x(0) ) * sin( 2* M_PI * x(1) ) * sin( 2* M_PI * x(2) );
  }

  /**
   * Computes u based on given x vector
   */
  static double uexact(const mfem::Vector& x)
  {
    return sin( 2* M_PI * x(0) ) * sin( 2* M_PI * x(1) ) * sin( 2* M_PI * x(2) );
  }

  mfem::GridFunction _scalar_gridfunc;

protected:
  std::string _mesh_type;

};

// DIRECTLY COPIED FROM MFEMSolverTest - consider simply inheriting here
void
MFEMMMSTest::testDiffusionSolve(mfem::Solver & solver, mfem::real_t tol)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int ne = 4;
  int dim = 3;
  mfem::Mesh mesh;
  mesh = mfem::Mesh::MakeCartesian3D(ne, ne, ne, mfem::Element::HEXAHEDRON, 1.0, 1.0, 1.0);

  // create parmesh from regular mesh
  mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);
  mesh.Clear();
  int order = 3;
  mfem::H1_FECollection fec(order, dim);
  mfem::ParFiniteElementSpace fespace(&pmesh, &fec);
  mfem::Array<int> ess_tdof_list, ess_bdr;
  if (pmesh.bdr_attributes.Size())
  {
    ess_bdr.SetSize(pmesh.bdr_attributes.Max());
    ess_bdr = 1;
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
  }
  mfem::FunctionCoefficient f(fexact);
  mfem::ParLinearForm b(&fespace);
  b.AddDomainIntegrator(new mfem::DomainLFIntegrator(f));
  b.Assemble();

  mfem::ParBilinearForm a(&fespace);
  mfem::ConstantCoefficient one(1.0);
  a.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
  a.Assemble();

  mfem::ParGridFunction x(&fespace);
  mfem::FunctionCoefficient uex(uexact);
  x = 0.0;
  x.ProjectBdrCoefficient(uex, ess_bdr);

  mfem::OperatorPtr A;
  mfem::Vector B, X;
  a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

  solver.SetOperator(*A);
  solver.Mult(B, X);

  mfem::Vector Y(X.Size());
  A->Mult(X, Y);
  Y -= B;
  ASSERT_LE(Y.Norml2(), tol);
}

TEST_F(MFEMMMSTest, SpatialConvergenceTest)
{
  // buildMFEMMesh("data/mug.e");

  // set the boundary conditions
  // top:
  InputParameters bc_params_top                        = _factory.getValidParams("MFEMScalarDirichletBC");
  bc_params_top.set<std::string>("variable")               = "diffused";
  bc_params_top.set<std::vector<BoundaryName>>("boundary") = {"2"};
  bc_params_top.set<Real>("value")                         = 0.;
  MFEMScalarDirichletBC & dirichlet_bc_top             = addObject<MFEMScalarDirichletBC>(
    "MFEMScalarDirichletBC", "top", bc_params_top
  );


  // bottom:
  InputParameters bc_params_bottom                     = _factory.getValidParams("MFEMScalarDirichletBC");
  bc_params_bottom.set<std::string>("variable")               = "diffused";
  bc_params_bottom.set<std::vector<BoundaryName>>("boundary") = {"1"};
  bc_params_bottom.set<Real>("value")           = {1.};
  MFEMScalarDirichletBC & dirichlet_bc_btm          = addObject<MFEMScalarDirichletBC>(
    "MFEMScalarDirichletBC", "bottom", bc_params_bottom
  );
  
  // material...
  InputParameters material_params                         = _factory.getValidParams("MFEMGenericConstantMaterial");
  material_params.set<std::vector<std::string>>("prop_names") = {"diffusivity"};
  material_params.set<std::vector<double>>("prop_values")     = {1.0};
  _mfem_problem->addMaterial("MFEMGenericConstantMaterial", "material1", material_params);


  // kernel params!
  InputParameters kernel_params                 =  _factory.getValidParams("MFEMDiffusionKernel");
  kernel_params.set<std::string>("variable")    = "diffused";
  kernel_params.set<std::string>("coefficient") = "diffusivity";
  MFEMDiffusionKernel & kernel                  = addObject<MFEMDiffusionKernel>(
    "MFEMDiffusionKernel", "diffKernel", kernel_params
  );


  // also need to set the FESpace params in the solver
  InputParameters fespace_params = _factory.getValidParams("MFEMFESpace");
  fespace_params.set<MooseEnum>("fec_type")  = "H1";
  fespace_params.set<MooseEnum>("fec_order") = "FIRST";
  MFEMFESpace & fespace = addObject<MFEMFESpace>("MFEMFESpace", "H1FESpace", fespace_params);
  
  // solver params!
  InputParameters solver_params          = _factory.getValidParams("MFEMHypreBoomerAMG");
  solver_params.set<double>("l_tol")     = 1e-16; // HypreBoomerAMG cannot set absolute tolerance
  solver_params.set<int>("l_max_its")    = 1000;
  MFEMHypreBoomerAMG & solver =
      addObject<MFEMHypreBoomerAMG>("MFEMHypreBoomerAMG", "solver1", solver_params);
  
  // Test MFEMSolver returns an solver of the expected type
  auto solver_downcast = std::dynamic_pointer_cast<mfem::HypreBoomerAMG>(solver.getSolver());
  solver_downcast->SetErrorMode(mfem::HypreSolver::ErrorMode::IGNORE_HYPRE_ERRORS);
  ASSERT_NE(solver_downcast.get(), nullptr);
  
  // here we go...
  testDiffusionSolve(*solver_downcast.get(), 1e-5);
}
