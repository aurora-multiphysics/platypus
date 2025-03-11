#include "MFEMObjectUnitTest.h"

#include "axom/slic.hpp"

#include "tribol/interface/tribol.hpp"
#include "tribol/interface/mfem_tribol.hpp"

class MFEMContactTest : public MFEMObjectUnitTest
{
public:
  MFEMContactTest() : MFEMObjectUnitTest("PlatypusApp") {}

  bool Run();

};

bool
MFEMContactTest::Run()
{
  // silence warning
  axom::slic::initialize();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  std::cout << "Inside unit test on rank " << rank << "\n";

  int ref_levels = 2;             // number of times to uniformly refine the serial mesh
  mfem::real_t lambda = 50.0;     // Lame parameter lambda
  mfem::real_t mu     = 50.0;     // Lame parameter mu (shear modulus)
  
  std::string mesh_file = "data/two-hex.mesh";
  constexpr int dim     = 3;
  constexpr int order   = 1;

  std::set<int>  mortar_attrs({4});
  std::set<int>  nonmortar_attrs({5});

  std::vector<std::set<int>> fixed_attrs(dim);
  fixed_attrs[0] = {1}; // x=0 plane of both blocks
  fixed_attrs[1] = {2}; // y=0 plane of both blocks
  fixed_attrs[2] = {3, 6};
  mfem::Mesh serial_mesh(mesh_file);
  for (int i = 0; i < ref_levels; ++i)
  {
     serial_mesh.UniformRefinement();
  }

  mfem::ParMesh mesh(MPI_COMM_WORLD, serial_mesh);

  serial_mesh.Clear();


  // Create an H1 finite element space on the mesh for displacements/forces
  mfem::H1_FECollection fec(order, dim);
  mfem::ParFiniteElementSpace fespace(&mesh, &fec, dim);
  auto n_displacement_dofs = fespace.GlobalTrueVSize();

  std::cout << "Number of displacement unknowns: " << n_displacement_dofs << "\n";

  // Create coordinate and displacement grid functions
  mfem::ParGridFunction coords(&fespace);
  mesh.SetNodalGridFunction(&coords);
  mfem::ParGridFunction displacement(&fespace);
  displacement = 0.0;

  // Find true dofs with homogeneous Dirichlet BCs
  mfem::Array<int> ess_tdof_list;
  {
    mfem::Array<int> ess_vdof_marker(fespace.GetVSize());
    ess_vdof_marker = 0;
    for (int i = 0; i < dim; ++i)
    {
        mfem::Array<int> ess_bdr(mesh.bdr_attributes.Max());
        ess_bdr = 0;
        for (auto xfixed_attr : fixed_attrs[i])
        {
          ess_bdr[xfixed_attr-1] = 1;
        }
        mfem::Array<int> new_ess_vdof_marker;
        fespace.GetEssentialVDofs(ess_bdr, new_ess_vdof_marker, i);
        for (int j = 0; j < new_ess_vdof_marker.Size(); ++j)
        {
          ess_vdof_marker[j] = ess_vdof_marker[j] || new_ess_vdof_marker[j];
        }
    }
    mfem::Array<int> ess_tdof_marker;
    fespace.GetRestrictionMatrix()->BooleanMult(ess_vdof_marker, ess_tdof_marker);
    mfem::FiniteElementSpace::MarkerToList(ess_tdof_marker, ess_tdof_list);
  }

  // Create small deformation linear elastic bilinear form
  mfem::ParBilinearForm a(&fespace);
  mfem::ConstantCoefficient lambda_coeff(lambda);
  mfem::ConstantCoefficient mu_coeff(mu);
  a.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda_coeff, mu_coeff));

  // Compute elasticity contribution to tangent stiffness matrix
  a.Assemble();
  std::unique_ptr<mfem::HypreParMatrix> A(new mfem::HypreParMatrix);
  a.FormSystemMatrix(ess_tdof_list, *A);


  tribol::initialize(dim, MPI_COMM_WORLD);

  // #2: Create a Tribol coupling scheme: defines contact surfaces and enforcement
  int coupling_scheme_id = 0;
  // NOTE: While there is a single mfem ParMesh for this problem, Tribol
  // defines a mortar and a nonmortar contact mesh, each with a unique mesh ID.
  // The Tribol mesh IDs for each contact surface are defined here.
  int mesh1_id = 0;
  int mesh2_id = 1;
  tribol::registerMfemCouplingScheme(
    coupling_scheme_id, mesh1_id, mesh2_id,
    mesh, coords, mortar_attrs, nonmortar_attrs,
    tribol::SURFACE_TO_SURFACE,
    tribol::NO_CASE,
    tribol::SINGLE_MORTAR,
    tribol::FRICTIONLESS,
    tribol::LAGRANGE_MULTIPLIER,
    tribol::BINNING_GRID
  );

  // #3: Set additional options/access pressure grid function on contact surfaces
  // Access Tribol's pressure grid function (on the contact surface). The
  // pressure ParGridFunction is created upon calling
  // registerMfemCouplingScheme(). It's lifetime coincides with the lifetime of
  // the coupling scheme, so the host code can reference and update it as
  // needed.
  auto& pressure = tribol::getMfemPressure(coupling_scheme_id); 
  
  std::cout << "Number of pressure unknowns: " <<
    pressure.ParFESpace()->GlobalTrueVSize() << std::endl;
  
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


  // #6a: Return contact contribution to the tangent stiffness matrix as a
  // block operator. See documentation for getMfemBlockJacobian() for block
  // definitions.
  auto A_blk = tribol::getMfemBlockJacobian(coupling_scheme_id);
  // Add elasticity contribution to the block operator returned by
  // getMfemBlockJacobian(). Note the (0,0) block Tribol returns is empty, so
  // we can simply set the (0,0) block to A.
  A_blk->SetBlock(0, 0, A.release());

  // Convert block operator to a single HypreParMatrix
  mfem::Array2D<const mfem::HypreParMatrix*> hypre_blocks(2, 2);
  for (int i{0}; i < 2; ++i)
  {
    for (int j{0}; j < 2; ++j)
    {
        if (A_blk->GetBlock(i, j).Height() != 0 && A_blk->GetBlock(i, j).Width() != 0)
        {
          hypre_blocks(i, j) =
              dynamic_cast<const mfem::HypreParMatrix*>(&A_blk->GetBlock(i, j));
        }
        else
        {
          hypre_blocks(i, j) = nullptr;
        }
    }
  }
  auto A_hyprePar = std::unique_ptr<mfem::HypreParMatrix>(
                      mfem::HypreParMatrixFromBlocks(hypre_blocks)
                    );
  
  // Create block RHS vector holding forces and gaps at tdofs
  mfem::BlockVector B_blk(A_blk->RowOffsets());
  B_blk = 0.0;

  // Fill with initial nodal gaps.
  // Note forces from contact are currently zero since pressure is zero prior
  // to first solve.
  mfem::Vector gap;
  // #6b: Return computed gap constraints on the contact surfaces
  tribol::getMfemGap(coupling_scheme_id, gap); // gap on ldofs
  auto& P_submesh = *pressure.ParFESpace()->GetProlongationMatrix();
  auto& gap_true = B_blk.GetBlock(1); // gap tdof vectorParFESpace()
  // gap is a dual vector, so (gap tdof vector) = P^T * (gap ldof vector)
  P_submesh.MultTranspose(gap, gap_true);

  // Create block solution vector holding displacements and pressures at tdofs
  mfem::BlockVector X_blk(A_blk->ColOffsets());
  X_blk = 0.0;

  mfem::MINRESSolver solver(MPI_COMM_WORLD);
  solver.SetRelTol(1.0e-12);
  solver.SetMaxIter(2000);
  solver.SetPrintLevel(3);
  solver.SetOperator(*A_hyprePar);
  mfem::HypreDiagScale prec(*A_hyprePar);
  prec.SetErrorMode(mfem::HypreSolver::IGNORE_HYPRE_ERRORS);
  solver.SetPreconditioner(prec);

  solver.Mult(B_blk, X_blk);

  // Update displacement and coords grid functions
  auto& displacement_true = X_blk.GetBlock(0);
  fespace.GetProlongationMatrix()->Mult(displacement_true, displacement);
  displacement.Neg();
  coords += displacement;

  // Update the pressure grid function
  auto& pressure_true = X_blk.GetBlock(1);
  P_submesh.Mult(pressure_true, pressure);

  // Verify the forces are in equilibrium, i.e. f_int = A*u = -f_contact = -B^T*p
  // This should be true if the solver converges.
  mfem::Vector f_int_true(fespace.GetTrueVSize());
  f_int_true = 0.0;
  mfem::Vector f_contact_true(f_int_true);
  A_blk->GetBlock(0, 0).Mult(displacement_true, f_int_true);
  A_blk->GetBlock(0, 1).Mult(pressure_true, f_contact_true);
  mfem::Vector resid_true(f_int_true);
  resid_true += f_contact_true;
  for (int i{0}; i < ess_tdof_list.Size(); ++i)
  {
     resid_true[ess_tdof_list[i]] = 0.0;
  }
  auto resid_linf = resid_true.Normlinf();

  // before we would test mfem::Mpi::Root()
  if (rank==0)
  {
    MPI_Reduce(MPI_IN_PLACE, &resid_linf, 1, MPI_FLOAT, MPI_MAX, 0,
              MPI_COMM_WORLD);
    std::cout << "|| force residual ||_(infty) = " << resid_linf << std::endl;
  }
  else
  {
    MPI_Reduce(&resid_linf, &resid_linf, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
  }

  // Verify the gap is closed by the displacements, i.e. B*u = gap
  // This should be true if the solver converges.
  mfem::Vector gap_resid_true(gap_true.Size());
  gap_resid_true = 0.0;
  A_blk->GetBlock(1, 0).Mult(displacement_true, gap_resid_true);
  gap_resid_true -= gap_true;
  auto gap_resid_linf = gap_resid_true.Normlinf();

  if (rank==0)
  {
     MPI_Reduce(MPI_IN_PLACE, &gap_resid_linf, 1, MPI_FLOAT, MPI_MAX, 0,
                MPI_COMM_WORLD);
     std::cout << "|| gap residual ||_(infty) = " << gap_resid_linf << std::endl;
  }
  else
  {
     MPI_Reduce(&gap_resid_linf, &gap_resid_linf, 1, MPI_FLOAT, MPI_MAX, 0,
                MPI_COMM_WORLD);
  }

  // Update the Tribol mesh based on deformed configuration
  // NOTE: This is done to update the contact submesh based on the current
  // deformed shape of the mesh.
  tribol::updateMfemParallelDecomposition();

  // Save data in VisIt format
  mfem::VisItDataCollection visit_vol_dc("contact-patch-test-volume", &mesh);
  visit_vol_dc.RegisterField("coordinates", &coords);
  visit_vol_dc.RegisterField("displacement", &displacement);
  visit_vol_dc.Save();
  mfem::VisItDataCollection visit_surf_dc("contact-patch-test-surface",
                                          pressure.ParFESpace()->GetMesh());
  visit_surf_dc.RegisterField("pressure", &pressure);
  visit_surf_dc.Save();

  // #7: Tribol cleanup: deletes coupling schemes and clears associated memory
  tribol::finalize();

  return true;
}


TEST_F(MFEMContactTest, Run)
{
  // Build required kernel inputs
  ASSERT_TRUE( Run() );
}


