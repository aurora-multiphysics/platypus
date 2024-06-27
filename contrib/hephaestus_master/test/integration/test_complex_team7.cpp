#include "auxsolvers.hpp"
#include "steady_executioner.hpp"

#include "hephaestus.hpp"

#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

class TestComplexTeam7
{
protected:
  static void SourceCurrent(const mfem::Vector & xv, mfem::Vector & J)
  {
    double x0(194e-3);  // Coil centre x coordinate
    double y0(100e-3);  // Coil centre y coordinate
    double a(50e-3);    // Coil thickness
    double i0(2742);    // Coil current in Ampere-turns
    double s(2.5e-3);   // Coil cross sectional area
    double freq(200.0); // Frequency in Hz

    double x = xv(0);
    double y = xv(1);

    // Current density magnitude
    double jmag = (i0 / s);

    // Calculate x component of current density unit vector
    if (abs(x - x0) < a)
    {
      J(0) = -(y - y0) / abs(y - y0);
    }
    else if (abs(y - y0) < a)
    {
      J(0) = 0.0;
    }
    else
    {
      J(0) =
          -(y - (y0 + a * ((y - y0) / abs(y - y0)))) /
          hypot(x - (x0 + a * ((x - x0) / abs(x - x0))), y - (y0 + a * ((y - y0) / abs(y - y0))));
    }

    // Calculate y component of current density unit vector
    if (abs(y - y0) < a)
    {
      J(1) = (x - x0) / abs(x - x0);
    }
    else if (abs(x - x0) < a)
    {
      J(1) = 0.0;
    }
    else
    {
      J(1) =
          (x - (x0 + a * ((x - x0) / abs(x - x0)))) /
          hypot(x - (x0 + a * ((x - x0) / abs(x - x0))), y - (y0 + a * ((y - y0) / abs(y - y0))));
    }

    // Calculate z component of current density unit vector
    J(2) = 0.0;

    // Scale by current density magnitude
    J *= jmag;
  }

  hephaestus::InputParameters TestParams()
  {
    hephaestus::Subdomain air("air", 1);
    air._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1.0));
    hephaestus::Subdomain plate("plate", 2);
    plate._scalar_coefficients.Register("electrical_conductivity",
                                        std::make_shared<mfem::ConstantCoefficient>(3.526e7));
    hephaestus::Subdomain coil1("coil1", 3);
    coil1._scalar_coefficients.Register("electrical_conductivity",
                                        std::make_shared<mfem::ConstantCoefficient>(1.0));
    hephaestus::Subdomain coil2("coil2", 4);
    coil2._scalar_coefficients.Register("electrical_conductivity",
                                        std::make_shared<mfem::ConstantCoefficient>(1.0));
    hephaestus::Subdomain coil3("coil3", 5);
    coil3._scalar_coefficients.Register("electrical_conductivity",
                                        std::make_shared<mfem::ConstantCoefficient>(1.0));
    hephaestus::Subdomain coil4("coil4", 6);
    coil4._scalar_coefficients.Register("electrical_conductivity",
                                        std::make_shared<mfem::ConstantCoefficient>(1.0));

    hephaestus::Coefficients coefficients(
        std::vector<hephaestus::Subdomain>({air, plate, coil1, coil2, coil3, coil4}));

    coefficients._scalars.Register("frequency", std::make_shared<mfem::ConstantCoefficient>(200.0));
    coefficients._scalars.Register("magnetic_permeability",
                                   std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));
    coefficients._scalars.Register("dielectric_permittivity",
                                   std::make_shared<mfem::ConstantCoefficient>(0.0));

    hephaestus::BCMap bc_map;

    mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./team7.g")).c_str(), 1, 1);

    hephaestus::Outputs outputs;
    outputs.Register("VisItDataCollection",
                     std::make_shared<mfem::VisItDataCollection>("ComplexMaxwellTeam7VisIt"));
    outputs.Register("ParaViewDataCollection",
                     std::make_shared<mfem::ParaViewDataCollection>("ComplexMaxwellTeam7ParaView"));

    hephaestus::Sources sources;

    // NB: needs to live to end of program so register to keep non-zero reference count.
    auto j_src_coef = std::make_shared<mfem::VectorFunctionCoefficient>(3, SourceCurrent);
    coefficients._vectors.Register("source_coefficient", j_src_coef);

    mfem::Array<mfem::VectorCoefficient *> sourcecoefs(4);
    sourcecoefs[0] = j_src_coef.get();
    sourcecoefs[1] = j_src_coef.get();
    sourcecoefs[2] = j_src_coef.get();
    sourcecoefs[3] = j_src_coef.get();

    mfem::Array<int> coilsegments(4);
    coilsegments[0] = 3;
    coilsegments[1] = 4;
    coilsegments[2] = 5;
    coilsegments[3] = 6;

    auto j_src_restricted =
        std::make_shared<mfem::PWVectorCoefficient>(3, coilsegments, sourcecoefs);
    coefficients._vectors.Register("source", j_src_restricted);

    hephaestus::InputParameters current_solver_options;
    current_solver_options.SetParam("Tolerance", float(1.0e-12));
    current_solver_options.SetParam("MaxIter", (unsigned int)200);

    sources.Register(
        "source",
        std::make_shared<hephaestus::DivFreeSource>(
            "source", "source", "HCurl", "H1", "electric_potential", current_solver_options));

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-16));
    solver_options.SetParam("MaxIter", (unsigned int)1000);

    hephaestus::InputParameters params;
    params.SetParam("UseGLVis", true);

    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);
    params.SetParam("SolverOptions", solver_options);
    hephaestus::logger.info("Created params ");
    return params;
  }
};

TEST_CASE_METHOD(TestComplexTeam7, "TestComplexTeam7", "[CheckRun]")
{
  hephaestus::InputParameters params(TestParams());
  auto pmesh = std::make_shared<mfem::ParMesh>(params.GetParam<mfem::ParMesh>("Mesh"));
  mfem::H1_FECollection fecm(1, 3);
  mfem::ParFiniteElementSpace pfespace(pmesh.get(), &fecm, 3);
  // Necessary, in case the nodal FE space is not set on the pmesh because it is lowest order.
  pmesh->SetNodalFESpace(&pfespace);

  auto bc_map(params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
  auto coefficients(params.GetParam<hephaestus::Coefficients>("Coefficients"));
  auto sources(params.GetParam<hephaestus::Sources>("Sources"));
  auto outputs(params.GetParam<hephaestus::Outputs>("Outputs"));
  auto solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
      "SolverOptions", hephaestus::InputParameters()));

  auto problem_builder =
      std::make_unique<hephaestus::ComplexAFormulation>("magnetic_reluctivity",
                                                        "electrical_conductivity",
                                                        "dielectric_permittivity",
                                                        "frequency",
                                                        "magnetic_vector_potential",
                                                        "magnetic_vector_potential_real",
                                                        "magnetic_vector_potential_imag");
  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace("HCurl", "ND_3D_P1");
  problem_builder->AddFESpace("HDiv", "RT_3D_P0");
  problem_builder->AddFESpace("H1", "H1_3D_P1");
  problem_builder->AddGridFunction("magnetic_vector_potential_real", "HCurl");
  problem_builder->AddGridFunction("magnetic_vector_potential_imag", "HCurl");
  problem_builder->AddGridFunction("magnetic_flux_density_real", "HDiv");
  problem_builder->AddGridFunction("magnetic_flux_density_imag", "HDiv");
  problem_builder->RegisterMagneticFluxDensityAux("magnetic_flux_density_real",
                                                  "magnetic_flux_density_imag");

  problem_builder->SetBoundaryConditions(bc_map);
  problem_builder->SetCoefficients(coefficients);
  problem_builder->SetSources(sources);
  problem_builder->SetOutputs(outputs);
  problem_builder->SetSolverOptions(solver_options);

  // Call LineSampler to save values
  std::string gridfunction_name("magnetic_flux_density_real");
  std::string csv_name("SimulatedA1B1Transect.csv");
  const int num_pts = 100;
  // Mesh bounding box (for the full serial mesh).
  mfem::Vector pos_min, pos_max;
  pmesh->GetBoundingBox(pos_min, pos_max, 1);
  pos_min(1) = 0.072;
  pos_max(1) = 0.072;
  pos_min(2) = 0.034;
  pos_max(2) = 0.034;
  std::shared_ptr<hephaestus::LineSamplerAux> linesamplerwriter =
      std::make_shared<hephaestus::LineSamplerAux>(
          gridfunction_name,
          pos_min,
          pos_max,
          num_pts,
          csv_name,
          "t (s), x (m), y (m), z (m), B_x (T), B_y (T), B_z (T)");
  linesamplerwriter->SetPriority(5);
  problem_builder->AddPostprocessor("LineSamplerWriter", linesamplerwriter);

  problem_builder->FinalizeProblem();

  auto problem = problem_builder->ReturnProblem();

  hephaestus::InputParameters exec_params;
  exec_params.SetParam("Problem", static_cast<hephaestus::SteadyStateProblem *>(problem.get()));

  auto executioner = std::make_unique<hephaestus::SteadyExecutioner>(exec_params);

  hephaestus::logger.info("Created exec ");

  executioner->Execute();
}
