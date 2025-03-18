#include "gtest/gtest.h"
#include "mfem.hpp"

// # Steepness of profile in the Dirichlet boundary condition
const mfem::real_t a(6.0);
// # Offset of the centre of the profile in the Dirichlet boundary condition from the centre of hte
// mesh
const mfem::real_t offset(0.0);
// # Angle of the magnetic field below the x-axis
const mfem::real_t theta(2.0);
// # Ratio of the perpendicular diffusivity to the parallel diffusivity
const mfem::real_t eps(1e-3);
// # Square root of the total number of elements in the mesh
const int meshres = 64;
// # How much the resolution should be distributed in the (roughly) perpendicular direction isntead
// of the parallel direction
const int mesh_anisotropy = 4;

mfem::real_t alpha = tan(theta * M_PI / 180.0);
const mfem::real_t k_par(10.0);
mfem::real_t k_per = eps * k_par;

mfem::real_t
T_exact_func(const mfem::Vector & x, mfem::real_t t)
{
  mfem::real_t T_func = 0.5 + 0.5 * tanh(a * (x[1] - offset - 77. + alpha * x[0])) *
                                  tanh(a * (23. - x[1] + offset - alpha * x[0]));
  return T_func;
}

mfem::real_t
forcing_func(const mfem::Vector & x, mfem::real_t t)
{
  mfem::real_t arg_py = a * (x[1] - offset + alpha * x[0] - 77);
  mfem::real_t arg_my = a * (-x[1] + offset - alpha * x[0] + 23);
  mfem::real_t f =
      a * a *
      (tanh(arg_my) * tanh(arg_py) * (pow(1.0 / cosh(arg_my), 2) + pow(1.0 / cosh(arg_py), 2)) +
       pow(1.0 / (cosh(arg_my) * cosh(arg_py)), 2)) *
      (1 + alpha * alpha) * k_per;
  return f;
}

// Anisotropic diffusion coefficient
void
kc_func_2D(const mfem::Vector & x, mfem::DenseMatrix & kc)
{
  mfem::Vector bhat({cos(theta * M_PI / 180.0), -sin(theta * M_PI / 180.0)});
  mfem::real_t k_dif(k_par - k_per);
  kc.SetSize(2);
  kc(0, 0) = k_dif * bhat[0] * bhat[0] + k_per;
  kc(0, 1) = k_dif * bhat[0] * bhat[1];
  kc(1, 0) = kc(0, 1);
  kc(1, 1) = k_dif * bhat[1] * bhat[1] + k_per;
}

TEST(CheckData, NCDiffusionTest)
{
  // 1. Build mesh
  mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(meshres / mesh_anisotropy + 1,
                                                meshres * mesh_anisotropy + 1,
                                                mfem::Element::Type::QUADRILATERAL,
                                                true,
                                                100.0,
                                                100.0);
  for (int i = 0; i < mesh.GetNBE(); i++)
  {
    mfem::Element * be = mesh.GetBdrElement(i);
    mfem::Array<int> vertices;
    be->GetVertices(vertices);
    be->SetAttribute(1);
  }
  mesh.SetAttributes();
  mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);

  // 2. Define a finite element space on the mesh. Here we use H1 continuous
  //    Lagrange finite elements of the given order.
  int order = 4;
  mfem::H1_FECollection fec(order, mesh.Dimension());
  mfem::ParFiniteElementSpace fespace(&pmesh, &fec);

  // 4. Extract the list of all the boundaries. These will be marked as
  //    Dirichlet in order to enforce zero boundary conditions.
  mfem::Array<int> ess_bdr(pmesh.bdr_attributes.Max());
  ess_bdr = 1;
  mfem::Array<int> ess_tdof_list;
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  // 5. Define the temperature T as a finite element grid function in fespace. Set
  //    the initial guess to the exact function.
  // mfem::FunctionCoefficient bcL_coef(boundary_forcing_func);
  mfem::FunctionCoefficient T_exact_coef(T_exact_func);
  mfem::ParGridFunction T(&fespace), T_exact(&fespace);
  T = 0.0;
  T.ProjectBdrCoefficient(T_exact_coef, ess_bdr);
  T_exact.ProjectCoefficient(T_exact_coef);

  // 6. Set up the the right-hand side.
  mfem::ParLinearForm b(&fespace);
  b = 0.0;
  mfem::FunctionCoefficient forcing_coef(forcing_func);
  b.AddDomainIntegrator(new mfem::DomainLFIntegrator(forcing_coef));
  b.Assemble();

  mfem::MatrixFunctionCoefficient kc_2D_coef(2, kc_func_2D);
  // Solve as a linear problem.
  {
    mfem::BilinearForm a(&fespace);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(kc_2D_coef));
    a.Assemble();

    mfem::OperatorPtr A;
    mfem::Vector C, Y;
    a.FormLinearSystem(ess_tdof_list, T, b, A, Y, C);

    mfem::GSSmoother M((mfem::SparseMatrix &)(*A));
    mfem::PCG(*A, M, C, Y, 1, 5000, 1e-12, 0.0);

    a.RecoverFEMSolution(Y, b, T);
  }

  bool debug_output = true;
  if (debug_output)
  {
    mfem::ParaViewDataCollection paraview_dc("NCDiffusion", &pmesh);
    paraview_dc.SetPrefixPath("ParaView");
    paraview_dc.SetLevelsOfDetail(order);
    paraview_dc.SetCycle(0);
    paraview_dc.SetDataFormat(mfem::VTKFormat::BINARY);
    paraview_dc.SetHighOrderOutput(true);
    paraview_dc.SetTime(0.0); // set the time
    paraview_dc.RegisterField("T", &T);
    paraview_dc.RegisterField("T_exact", &T_exact);
    paraview_dc.Save();
  }

  T -= T_exact;

  EXPECT_NEAR(T.Norml2(), 0, 1e-5);
}
