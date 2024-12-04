#include "MFEMObjectUnitTest.h"
#include "MFEMMHDTest.h"

class MFEMMHDTest : public MFEMObjectUnitTest
{
public:
  MFEMMHDTest() : MFEMObjectUnitTest("PlatypusApp") {}
};

void
LiMHDPreconditioner::ImplicitSolve(const double dt, const mfem::Vector & X, mfem::Vector & dX_dt)
{
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
  mfem::ParFiniteElementSpace v_fespace(pmesh.get(), fec_h1.get(), 3);
  mfem::ParFiniteElementSpace q_fespace(pmesh.get(), fec_h1_o2.get());

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

  mfem::ParNonlinearForm nlf_F_kappa(&v_fespace);
  // $(2/tau)*(v, v')$
  nlf_F_kappa.AddDomainIntegrator(new mfem::MassIntegrator());
  // A_AL = $(1/Re)*(\nabla w, \nabla v) + \alpha (\nabla \cdot w, \nabla \cdot v) $
  // nlf_F_kappa.AddDomainIntegrator(new mfem::ElasticityIntegrator());
  // O = $.5*(u \cdot \nabla v, w ) - .5*(u \cdot \nabla w, v )$
  nlf_F_kappa.AddDomainIntegrator(new mfem::SkewSymmetricVectorConvectionNLFIntegrator());
  // $\kappa*(Kv, Kv')$
  // <custom integrator here?>

  // $(Q \nabla \cdot v, q)$ where $v=(v_1,\cdots,v_n)$ in H1^n, $v$ in H1.
  mfem::ParMixedBilinearForm blf_B(&v_fespace, &q_fespace);
  blf_B.AddDomainIntegrator(new mfem::VectorDivergenceIntegrator());

  // $(d, d')+(\nabla \cdot d, \nabla \cdot d')$ where $d, d'$ is in $H(div)$
  mfem::ParBilinearForm blf_D_j(&d_fespace);
  blf_D_j.AddDomainIntegrator(new mfem::MassIntegrator());
  blf_D_j.AddDomainIntegrator(new mfem::DivDivIntegrator());
  blf_D_j.Assemble();

  // -(Q \phi, \nabla \cdot d)$ where $\phi$ is in $L_2$, and $d$ is in $H(div)$
  mfem::ParMixedBilinearForm blf_G(&phi_fespace, &d_fespace);
  blf_G.AddDomainIntegrator(new mfem::MixedScalarWeakGradientIntegrator());
  blf_G.Assemble();

  // (d, B_n \times v')$ where $v$ is in $H1^n$, and $d$ is in $H(div)$
  mfem::ParMixedBilinearForm blf_K(&d_fespace, &v_fespace);
  // <custom integrator here?>
}

/**
 * Test MFEMMHDSolver
 */
TEST_F(MFEMMHDTest, MFEMMHDTestSolve)
{
  MFEMProblemData problem_data;
  LiMHDPreconditioner preconditioner(problem_data);
}
