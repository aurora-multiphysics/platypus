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
  d_fespace;
  phi_fespace;
  v_fespace;
  q_fespace;

  // $(q, q')$ where $q, q'$ is in $H1$
  mfem::BilinearForm blf_M_p(q_fespace);
  blf_M_p->AddDomainIntegrator(new mfem::MassIntegrator());
  blf_M_p.Assemble();

  // (Q \nabla q, \nabla q')$ where $q, q'$ is in $H1$
  mfem::BilinearForm blf_S_p(q_fespace);
  blf_S_p->AddDomainIntegrator(new mfem::DiffusionIntegrator());
  blf_S_p.Assemble();

  // $(\phi, \phi')$ where $\phi, \phi'$ is in $L2$
  mfem::BilinearForm blf_M_phi(phi_fespace);
  blf_M_phi->AddDomainIntegrator(new mfem::MassIntegrator());
  blf_M_phi.Assemble();

  mfem::NonlinearForm nlf_F_kappa(v_fespace);
  // $(2/tau)*(v, v')$
  nlf_F_kappa->AddDomainIntegrator(new mfem::MassIntegrator());
  // A_AL = $(1/Re)*(\nabla w, \nabla v) + \alpha (\nabla \cdot w, \nabla \cdot v) $
  nlf_F_kappa->AddDomainIntegrator(new mfem::ElasticityIntegrator());
  // O = $.5*(u \cdot \nabla v, w ) - .5*(u \cdot \nabla w, v )$
  nlf_F_kappa->AddDomainIntegrator(new mfem::SkewSymmetricVectorConvectionNLFIntegrator());
  // $\kappa*(Kv, Kv')$
  // <custom integrator here?>

  // $(Q \nabla \cdot v, q)$ where $v=(v_1,\cdots,v_n)$ in H1^n, $v$ in H1.
  mfem::MixedBilinearForm blf_B(v_fespace, q_fespace);
  blf_B->AddDomainIntegrator(new mfem::VectorDivergenceIntegrator());

  // $(d, d')+(\nabla \cdot d, \nabla \cdot d')$ where $d, d'$ is in $H(div)$
  mfem::BilinearForm blf_D_j(d_fespace);
  blf_D_j->AddDomainIntegrator(new mfem::MassIntegrator());
  blf_D_j->AddDomainIntegrator(new mfem::DivDivIntegrator());
  blf_D_j.Assemble();

  // -(Q \phi, \nabla \cdot d)$ where $\phi$ is in $L_2$, and $d$ is in $H(div)$
  mfem::MixedBilinearForm blf_G(phi_fespace, d_fespace);
  blf_G->AddDomainIntegrator(new mfem::MixedScalarWeakGradientIntegrator());
  blf_G.Assemble();

  // (d, B_n \times v')$ where $v$ is in $H1^n$, and $d$ is in $H(div)$
  mfem::MixedBilinearForm blf_K(d_fespace, v_fespace);
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
