#include "gtest/gtest.h"
#include "mfem.hpp"

double
f(double u)
{
  return u;
}
double
df(double u)
{
  return 1.0;
}

// Define a coefficient that, given a grid function u, function func, returns func(u)
class NonlinearGridFunctionCoefficient : public mfem::Coefficient
{
  mfem::GridFunction & gf;
  std::function<double(double)> func;

public:
  NonlinearGridFunctionCoefficient(mfem::GridFunction & gf_, std::function<double(double)> func_)
    : gf(gf_), func(func_)
  {
  }
  double Eval(mfem::ElementTransformation & T, const mfem::IntegrationPoint & ip)
  {
    return func(gf.GetValue(T, ip));
  }
};

// Define a nonlinear integrator that computes (f(u), v) and its linearized
// operator, (u df(u), v).
//
// Note that the action (f(u), v) can be computed using DomainLFIntegrator
// and the Jacobian matrix linearized operator can be computed using
// MassIntegrator with the appropriate coefficients.
// This integrator is supported for H1 and L2 fespaces
class NonlinearMassIntegrator : public mfem::NonlinearFormIntegrator
{
  mfem::FiniteElementSpace & fes;
  mfem::GridFunction gf;
  mfem::Array<int> dofs;
  std::function<double(double)> func;
  std::function<double(double)> dfunc;

public:
  NonlinearMassIntegrator(mfem::FiniteElementSpace & fes_,
                          std::function<double(double)> func_,
                          std::function<double(double)> dfunc_)
    : fes(fes_), gf(&fes), func(func_), dfunc(dfunc_)
  {
  }

  virtual void AssembleElementVector(const mfem::FiniteElement & el,
                                     mfem::ElementTransformation & Tr,
                                     const mfem::Vector & elfun,
                                     mfem::Vector & elvect)
  {
    fes.GetElementDofs(Tr.ElementNo, dofs);
    gf.SetSubVector(dofs, elfun);
    NonlinearGridFunctionCoefficient coeff(gf, func);
    mfem::DomainLFIntegrator integ(coeff);
    integ.AssembleRHSElementVect(el, Tr, elvect);
  }

  virtual void AssembleElementGrad(const mfem::FiniteElement & el,
                                   mfem::ElementTransformation & Tr,
                                   const mfem::Vector & elfun,
                                   mfem::DenseMatrix & elmat)
  {
    fes.GetElementDofs(Tr.ElementNo, dofs);
    gf.SetSubVector(dofs, elfun);
    NonlinearGridFunctionCoefficient coeff(gf, dfunc);
    mfem::MassIntegrator integ(coeff);
    integ.AssembleElementMatrix(el, Tr, elmat);
  }
};

class NLOperator : public mfem::Operator
{
private:
  // The FE-spaces and operators
  int nBlocks;
  mfem::ParFiniteElementSpace * feSpace = NULL; // FE-Spaces
  mfem::Operator * NLForm = NULL;               // Operators used for E(U)
  mfem::BilinearForm * a = NULL;
  mfem::OperatorPtr A;

  // Forms used for Forcing vector (b)
  mfem::ParLinearForm * LForm = NULL; // Linear forms

  // Array of constrainedNodes
  mfem::Array<int> bdr_tDofs;

  // vectors
  mutable mfem::Vector b_vec, x_Dirch; // BC vectors
  mutable mfem::Vector z_vec, w_vec;   // Temporary vectors

public:
  // Creates my block non-linear form
  NLOperator(mfem::ParFiniteElementSpace * feSpace_,
             mfem::ParNonlinearForm * NLForm_,
             mfem::ParLinearForm * LForms_,
             const mfem::MemoryType deviceMT);

  // Destroys my block linear form
  ~NLOperator();

  // The residual operator
  virtual void Mult(const mfem::Vector & x, mfem::Vector & y) const override;

  // The Jacobian operator (Can be used for preconditioning)
  mfem::Operator & GetGradient(const mfem::Vector & x) const override;

  // A null and void function that inherits from
  // the Operator class
  virtual void SetOperator(const mfem::Operator & op) {};
};

NLOperator::NLOperator(mfem::ParFiniteElementSpace * feSpace_,
                             mfem::ParNonlinearForm * NLForm_,
                             mfem::ParLinearForm * LForms_,
                             const mfem::MemoryType deviceMT)
  : Operator(feSpace_->TrueVSize()), NLForm(NLForm_), LForm(LForms_)
{
  feSpace = new mfem::ParFiniteElementSpace(*feSpace_);
  b_vec = mfem::Vector(feSpace_->TrueVSize(), deviceMT);
  b_vec = 0.0;
  x_Dirch = mfem::Vector(feSpace_->TrueVSize(), deviceMT);
  x_Dirch = 0.0;
  z_vec = mfem::Vector(feSpace_->TrueVSize(), deviceMT);
  z_vec = 0.0;
  w_vec = mfem::Vector(feSpace_->TrueVSize(), deviceMT);
  w_vec = 0.0;
}

NLOperator::~NLOperator() {};

void
NLOperator::Mult(const mfem::Vector & x, mfem::Vector & y) const
{
  NLForm->Mult(x, y);
};

mfem::Operator &
NLOperator::GetGradient(const mfem::Vector & x) const
{
  return (NLForm->GetGradient(x));
};

TEST(CheckData, TestNonLinearDiffusionIntegratorInhomogenous)
{

  // 1. Parse command line options

  int order = 1;
  bool nonzero_rhs = true;

  // 2. Read the mesh from the given mesh file, and refine once uniformly.
  mfem::Mesh mesh =
      mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::Type::QUADRILATERAL, true, 2.0, 2.0);
  mesh.UniformRefinement();
  mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);

  // 3. Define a finite element space on the mesh. Here we use H1 continuous
  //    high-order Lagrange finite elements of the given order.
  mfem::H1_FECollection fec(order, mesh.Dimension());
  mfem::ParFiniteElementSpace fespace(&pmesh, &fec);

  // 4. Extract the list of all the boundaries. These will be marked as
  //    Dirichlet in order to enforce zero boundary conditions.
  mfem::Array<int> ess_bdr(pmesh.bdr_attributes.Max());
  ess_bdr = 1;
  mfem::Array<int> ess_tdof_list;
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  // 5. Define the solution x as a finite element grid function in fespace. Set
  //    the initial guess to zero, which also sets the boundary conditions.
  mfem::ParGridFunction u1(&fespace), u2(&fespace);
  u1 = 0.0;
  u2 = 0.0;

  // Solve as a non-linear problem.

  // 6. Set up the nonlinear form n(u,v) = (grad u, grad v) + (f(u), v)
  mfem::ParNonlinearForm n(&fespace);
  n.AddDomainIntegrator(new NonlinearMassIntegrator(fespace, f, df));
  n.AddDomainIntegrator(new mfem::DiffusionIntegrator);

  // 7. Set up the the right-hand side. For simplicitly, we just use a zero
  //    vector. Because of the form of the nonlinear function f, it is still
  //    nontrivial to solve n(u,v) = 0.
  mfem::ParLinearForm b(&fespace);
  b = 0.0;
  if (nonzero_rhs)
  {
    mfem::ConstantCoefficient five(5.0);
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(five));
    b.Assemble();
  }

  // 8. Get true dof vectors and set essential BCs on rhs.
  mfem::Vector X(fespace.GetTrueVSize()), B(fespace.GetTrueVSize());
  u1.GetTrueDofs(X);
  b.ParallelAssemble(B);
  n.SetEssentialBC(ess_bdr, &B);

  // 9. Set up the Newton solver. Each Newton iteration requires a linear
  //    solve. Here we use UMFPack as a direct solver for these systems.
  mfem::CGSolver solver(MPI_COMM_WORLD);
  mfem::NewtonSolver newton(MPI_COMM_WORLD);
  newton.SetOperator(n);
  newton.SetSolver(solver);
  newton.SetPrintLevel(1);
  newton.SetRelTol(1e-10);
  newton.SetMaxIter(20);

  // 10. Solve the nonlinear system.
  newton.Mult(B, X);
  u1.Distribute(X);

  // Solve as a linear problem.
  {
    mfem::BilinearForm a(&fespace);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator);
    a.AddDomainIntegrator(new mfem::MassIntegrator);
    a.Assemble();

    mfem::OperatorPtr A;
    mfem::Vector C, Y;
    a.FormLinearSystem(ess_tdof_list, u2, b, A, Y, C);
    mfem::GSSmoother M((mfem::SparseMatrix &)(*A));
    mfem::PCG(*A, M, C, Y, 1, 500, 1e-12, 0.0);
    a.RecoverFEMSolution(Y, b, u2);
  }

  u1 -= u2;

  EXPECT_NEAR(u1.Norml2(), 0, 1e-5);
}
