#include "MFEMObjectUnitTest.h"
#include "MFEMMHDTest.h"

class MFEMMHDTest : public MFEMObjectUnitTest
{
public:
  MFEMMHDTest() : MFEMObjectUnitTest("PlatypusApp") {}
};

void
LiMHDPreconditioner::Assemble()
{
}

void
LiMHDPreconditioner::ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt)
{
  Assemble();
}

void
u_exact(const mfem::Vector & x, mfem::Vector & f)
{
  double u_max(1.0);
  double y_max(1.0);
  double z_max(1.0);
  double y(x(1));
  double z(x(2));
  f(0) = (9.0 / 4.0) * u_max * (1 - (y * y) / (y_max * y_max)) * (1 - (z * z) / (z_max * z_max));
  f(1) = 0.0;
  f(2) = 0.0;
}

/**
 * Test MFEMMHDSolver
 */
TEST_F(MFEMMHDTest, MFEMMHDTestSolve)
{
  MFEMProblemData problem_data;
  LiMHDPreconditioner preconditioner(problem_data);

  mfem::Mesh mesh("data/square_Ha_20.e", 1, 1);
  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);

  mesh.EnsureNodes();
  int dim = mesh.Dimension();

  auto fec_h1 = std::make_shared<mfem::H1_FECollection>(1, pmesh->Dimension());
  auto fec_h1_o2 = std::make_shared<mfem::H1_FECollection>(2, pmesh->Dimension());
  auto fec_nd = std::make_shared<mfem::ND_FECollection>(1, pmesh->Dimension());
  auto fec_rt = std::make_shared<mfem::RT_FECollection>(1, pmesh->Dimension());
  auto fec_l2 = std::make_shared<mfem::L2_FECollection>(1, pmesh->Dimension());

  mfem::ParFiniteElementSpace d_fespace(pmesh.get(), fec_rt.get());
  mfem::ParFiniteElementSpace phi_fespace(pmesh.get(), fec_l2.get());
  mfem::ParFiniteElementSpace v_fespace(pmesh.get(), fec_h1_o2.get(), 3);
  mfem::ParFiniteElementSpace q_fespace(pmesh.get(), fec_h1.get());

  mfem::ParGridFunction j_n(&d_fespace);
  mfem::ParGridFunction phi_n(&phi_fespace);
  mfem::ParGridFunction u_n(&v_fespace);
  mfem::ParGridFunction u_n_1(&v_fespace);   // u_{n-1}
  mfem::ParGridFunction u_n_2(&v_fespace);   // u_{n-2}
  mfem::ParGridFunction ubar_n(&v_fespace);  // 0.5(u_n + u_{n-1})
  mfem::ParGridFunction ustar_n(&v_fespace); // 0.5(3u_{n-1} - u_{n-2})
  mfem::ParGridFunction p_n(&q_fespace);
  mfem::ParGridFunction xi(&q_fespace);
  mfem::ParGridFunction eta(&q_fespace);

  mfem::ConstantCoefficient zero(0.0);
  mfem::ConstantCoefficient one(1.0);
  mfem::ConstantCoefficient negone(1.0);
  mfem::ConstantCoefficient half(0.5);
  mfem::ConstantCoefficient reciprocal_Re(1 / 100.0);
  mfem::ConstantCoefficient alpha(0.5);
  mfem::VectorGridFunctionCoefficient ustar_coef(&ustar_n);
  mfem::VectorFunctionCoefficient ucoef(3, u_exact);
  mfem::Vector zero_vec(3);
  zero_vec = 0.0;
  mfem::VectorConstantCoefficient zero_vec_coef(zero_vec);

  // Initial conditions
  u_n.ProjectCoefficient(ucoef);
  u_n_1.ProjectCoefficient(ucoef);
  u_n_2.ProjectCoefficient(ucoef);
  ubar_n.ProjectCoefficient(ucoef);
  ustar_n.ProjectCoefficient(ucoef);
  j_n.ProjectCoefficient(zero_vec_coef);
  p_n.ProjectCoefficient(zero);
  phi_n.ProjectCoefficient(zero);
  xi.ProjectCoefficient(zero);
  eta.ProjectCoefficient(zero);

  mfem::ParLinearForm r_j(&d_fespace);
  r_j.Assemble();

  mfem::ParLinearForm r_phi(&phi_fespace);
  r_phi.Assemble();

  mfem::ParLinearForm r_u(&v_fespace);
  r_u.Assemble();

  mfem::ParLinearForm r_p(&q_fespace);
  r_p.Assemble();

  // $(q, q')$ where $q, q'$ is in $H1$
  mfem::ParBilinearForm blf_M_p(&q_fespace);
  blf_M_p.AddDomainIntegrator(new mfem::MassIntegrator());
  blf_M_p.Assemble();

  // (Q \nabla q, \nabla q')$ where $q, q'$ is in $H1$
  mfem::ParBilinearForm blf_S_p(&q_fespace);
  blf_S_p.AddDomainIntegrator(new mfem::DiffusionIntegrator());
  blf_S_p.Assemble();

  // $(\phi, \phi')$ where $\phi, \phi'$ is in $L2$
  mfem::ParBilinearForm blf_M_phi(&phi_fespace);
  blf_M_phi.AddDomainIntegrator(new mfem::MassIntegrator());
  blf_M_phi.Assemble();

  mfem::ParBilinearForm blf_F_kappa(&v_fespace);
  // $(2/tau)*(v, v')$
  blf_F_kappa.AddDomainIntegrator(new mfem::VectorMassIntegrator());
  // A_AL = $(1/Re)*(\nabla w, \nabla v) + \alpha (\nabla \cdot w, \nabla \cdot v) $
  blf_F_kappa.AddDomainIntegrator(new mfem::VectorDiffusionIntegrator(reciprocal_Re));
  // nlf_F_kappa.AddDomainIntegrator(new mfem::ElasticityIntegrator());
  // O = $.5*(u \cdot \nabla v, w ) - .5*(u \cdot \nabla w, v )$
  // nlf_F_kappa.AddDomainIntegrator(new mfem::SkewSymmetricVectorConvectionNLFIntegrator());
  blf_F_kappa.AddDomainIntegrator(new mfem::ConvectionIntegrator(ustar_coef, 0.5));
  blf_F_kappa.AddDomainIntegrator(new mfem::ConservativeConvectionIntegrator(ustar_coef, 0.5));
  // $\kappa*(Kv, Kv')$
  // <custom integrator here?>

  // $(Q \nabla \cdot v, q)$ where $v=(v_1,\cdots,v_n)$ in H1^n, $v$ in H1.
  mfem::ParMixedBilinearForm blf_B(&v_fespace, &q_fespace);
  blf_B.AddDomainIntegrator(new mfem::VectorDivergenceIntegrator());

  // $(d, d')+(\nabla \cdot d, \nabla \cdot d')$ where $d, d'$ is in $H(div)$
  mfem::ParBilinearForm blf_D_j(&d_fespace);
  blf_D_j.AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
  blf_D_j.AddDomainIntegrator(new mfem::DivDivIntegrator());
  blf_D_j.Assemble();

  // -(Q \phi, \nabla \cdot d)$ where $\phi$ is in $L_2$, and $d$ is in $H(div)$
  mfem::ParMixedBilinearForm blf_G(&phi_fespace, &d_fespace);
  blf_G.AddDomainIntegrator(new mfem::MixedScalarWeakGradientIntegrator());
  blf_G.Assemble();

  // (d, B_n \times v')$ where $v$ is in $H1^n$, and $d$ is in $H(div)$
  mfem::ParMixedBilinearForm blf_K(&d_fespace, &v_fespace);
  // <custom integrator here?>

  // preconditioner solve Py = r

  // Form the linear systems for both
  //       M_p xi = r_p,
  //       S_p eta = r_p.
  {
    mfem::Array<int> pressure_ess_dofs;
    q_fespace.GetBoundaryTrueDofs(pressure_ess_dofs);

    // Solve the M_p xi = r_p system using PCG with Jacobi preconditioner.
    {
      mfem::HypreParMatrix Mp;
      mfem::Vector Bmp, Xi;
      blf_M_p.FormLinearSystem(pressure_ess_dofs, xi, r_p, Mp, Xi, Bmp);
      mfem::CGSolver M_solver(MPI_COMM_WORLD);
      mfem::HypreSmoother M_prec;

      M_solver.iterative_mode = false;
      M_solver.SetRelTol(1e-8);
      M_solver.SetAbsTol(0.0);
      M_solver.SetMaxIter(10);
      M_solver.SetPrintLevel(0);
      M_prec.SetType(mfem::HypreSmoother::Jacobi); // Works for now, but check if diagonal...
      M_solver.SetPreconditioner(M_prec);
      M_solver.SetOperator(Mp);

      M_solver.Mult(Bmp, Xi);
      blf_M_p.RecoverFEMSolution(Xi, r_p, xi);
    }

    // Solve the S_p eta = r_p system using two iterations of HypreBoomerAMG.
    {
      mfem::HypreParMatrix Sp;
      mfem::Vector Bsp, Eta;
      blf_S_p.FormLinearSystem(pressure_ess_dofs, eta, r_p, Sp, Eta, Bsp);
      mfem::HypreBoomerAMG amg(Sp);
      mfem::HyprePCG pcg(Sp);
      pcg.SetTol(1e-12);
      pcg.SetMaxIter(2);
      pcg.SetPrintLevel(2);
      pcg.SetPreconditioner(amg);
      pcg.Mult(Bsp, Eta);
      blf_S_p.RecoverFEMSolution(Eta, r_p, eta);
    }
    // Do the sum y_p = -alpha1 xi - eta.  Note: above alpha1 must be defined negative.
    p_n.Add(alpha.constant, xi);
    p_n.Add(-1.0, eta);
  }

  // Solve the M_phi phi_n = -r_phi system using PCG with Jacobi preconditioner.
  mfem::Array<int> phi_ess_dofs;
  phi_fespace.GetBoundaryTrueDofs(phi_ess_dofs);
  {
    r_phi *= -1.0;
    mfem::HypreParMatrix Mphi;
    mfem::Vector Bmpphi, Phi_n;
    blf_M_phi.FormLinearSystem(phi_ess_dofs, phi_n, r_phi, Mphi, Phi_n, Bmpphi);

    mfem::CGSolver M_phi_solver(MPI_COMM_WORLD);
    mfem::HypreSmoother M_phi_prec;

    M_phi_solver.iterative_mode = false;
    M_phi_solver.SetRelTol(1e-8);
    M_phi_solver.SetAbsTol(0.0);
    M_phi_solver.SetMaxIter(10);
    M_phi_solver.SetPrintLevel(0);
    M_phi_prec.SetType(mfem::HypreSmoother::Jacobi); // Works for now, but check if diagonal...
    M_phi_solver.SetPreconditioner(M_phi_prec);
    M_phi_solver.SetOperator(Mphi);

    M_phi_solver.Mult(Bmpphi, Phi_n);
    blf_M_phi.RecoverFEMSolution(Phi_n, r_phi, phi_n);
  }

  // Solve the F_k u_n = r_u - B^T p_n system with an additive Schwarz preconditioner
  mfem::Array<int> u_ess_dofs;
  v_fespace.GetBoundaryTrueDofs(u_ess_dofs);
  // u_n.ProjectBdrCoefficient(*_vec_coef->getVectorCoefficient(), ess_bdrs);
  // mfem::common::AttrToMarker(mesh.bdr_attributes.Max(), _bdr_attributes, ess_bdrs);

  {
    blf_B.AddMultTranspose(p_n, r_u, -1.0); // currently failing
    mfem::HypreParMatrix F_kappa;
    mfem::Vector Bfk, Un;
    blf_F_kappa.FormLinearSystem(u_ess_dofs, u_n, r_u, F_kappa, Un, Bfk);
    mfem::HypreBoomerAMG additive_schwartz(F_kappa);
    mfem::GMRESSolver gmres(MPI_COMM_WORLD);
    gmres.SetRelTol(1e-3);
    gmres.SetMaxIter(100);
    gmres.SetPrintLevel(2);
    gmres.SetPreconditioner(additive_schwartz);
    gmres.SetOperator(F_kappa);
    gmres.Mult(Bfk, Un);
    blf_F_kappa.RecoverFEMSolution(Un, r_u, u_n);
  }

  // Solve the D_j J_n = r_j - 2G^T phi_n - 2 K^T u_n system using 5 iterations of HypreADS.
  mfem::Array<int> j_ess_dofs;
  d_fespace.GetBoundaryTrueDofs(j_ess_dofs);
  {
    blf_G.AddMult(phi_n, r_j, -2.0);
    //  blf_K->AddMultTranspose(u_n, r_j, -2.0);
    mfem::HypreParMatrix Dj;
    mfem::Vector Bdj, Jn;
    blf_D_j.FormLinearSystem(j_ess_dofs, j_n, r_j, Dj, Jn, Bdj);
    mfem::HypreADS ads(Dj, &d_fespace);
    mfem::HyprePCG pcg(Dj);
    pcg.SetTol(1e-12);
    pcg.SetMaxIter(5);
    pcg.SetPrintLevel(2);
    pcg.SetPreconditioner(ads);
    pcg.Mult(Bdj, Jn);
    blf_D_j.RecoverFEMSolution(Jn, r_j, j_n);
  }
}
