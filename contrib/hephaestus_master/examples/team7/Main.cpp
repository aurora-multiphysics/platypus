// Implements the TEAM Problem 7 benchmark in the time domain
// available from https://www.compumag.org/wp/team/.
// Reference results available from
// Fujiwara, K. and Nakata, T. (1990), Results for Benchmark Problem 7
// (asymmetrical conductor with a hole), COMPEL, Vol. 9 No. 3, pp. 137-154.
// https://doi.org/10.1108/eb010071

#include "hephaestus.hpp"

const char * DATA_DIR = "../../data/";

static void
source_current(const mfem::Vector & xv, double t, mfem::Vector & J)
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
  double jmag = (i0 / s) * sin(2 * M_PI * freq * t);

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
    J(0) = -(y - (y0 + a * ((y - y0) / abs(y - y0)))) /
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
    J(1) = (x - (x0 + a * ((x - x0) / abs(x - x0)))) /
           hypot(x - (x0 + a * ((x - x0) / abs(x - x0))), y - (y0 + a * ((y - y0) / abs(y - y0))));
  }

  // Calculate z component of current density unit vector
  J(2) = 0.0;

  // Scale by current density magnitude
  J *= jmag;
}

hephaestus::Coefficients
defineCoefficients()
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
  coefficients._scalars.Register("magnetic_permeability",
                                 std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  auto j_src_coef = std::make_shared<mfem::VectorFunctionCoefficient>(3, source_current);

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

  // Register to prevent j_src_coef being destroyed when it goes out of scope.
  coefficients._vectors.Register("source_coefficient", std::move(j_src_coef));

  auto j_src_restricted = std::make_shared<mfem::PWVectorCoefficient>(3, coilsegments, sourcecoefs);
  coefficients._vectors.Register("source", j_src_restricted);

  return coefficients;
}

hephaestus::Sources
defineSources()
{
  hephaestus::InputParameters current_solver_options;
  current_solver_options.SetParam("Tolerance", float(1.0e-12));
  current_solver_options.SetParam("MaxIter", (unsigned int)200);

  hephaestus::Sources sources;
  sources.Register(
      "source",
      std::make_shared<hephaestus::DivFreeSource>(
          "source", "source", "HCurl", "H1", "_source_potential", current_solver_options, false));

  return sources;
}

hephaestus::Outputs
defineOutputs()
{
  hephaestus::Outputs outputs;
  outputs.Register("ParaViewDataCollection",
                   std::make_shared<mfem::ParaViewDataCollection>("Team7ParaView"));
  return outputs;
}

int
main(int argc, char * argv[])
{
  mfem::OptionsParser args(argc, argv);
  args.AddOption(
      &DATA_DIR, "-dataDir", "--data_directory", "Directory storing input data for tests.");
  args.Parse();
  MPI_Init(&argc, &argv);

  // Create Formulation
  auto problem_builder = std::make_unique<hephaestus::AFormulation>("magnetic_reluctivity",
                                                                    "magnetic_permeability",
                                                                    "electrical_conductivity",
                                                                    "magnetic_vector_potential");
  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./team7.g")).c_str(), 1, 1);
  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);
  mfem::H1_FECollection fecm(1, 3);
  mfem::ParFiniteElementSpace pfespace(pmesh.get(), &fecm, 3);
  // Necessary, in case the nodal FE space is not set on the pmesh because it is lowest order.
  pmesh->SetNodalFESpace(&pfespace);

  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace("H1", "H1_3D_P1");
  problem_builder->AddFESpace("HCurl", "ND_3D_P1");
  problem_builder->AddFESpace("HDiv", "RT_3D_P0");
  problem_builder->AddFESpace("Scalar_L2", "L2_3D_P0");
  problem_builder->AddFESpace("Vector_L2", "L2_3D_P0", 3);
  problem_builder->AddGridFunction("magnetic_vector_potential", "HCurl");

  problem_builder->AddGridFunction("magnetic_flux_density", "HDiv");
  problem_builder->RegisterMagneticFluxDensityAux("magnetic_flux_density");

  problem_builder->AddGridFunction("current_density", "HDiv");
  problem_builder->RegisterCurrentDensityAux("current_density");

  problem_builder->AddGridFunction("electric_field", "HCurl");
  problem_builder->RegisterElectricFieldAux("electric_field");

  problem_builder->AddGridFunction("lorentz_force_density", "Vector_L2");
  problem_builder->RegisterLorentzForceDensityAux(
      "lorentz_force_density", "magnetic_flux_density", "current_density");

  problem_builder->AddGridFunction("joule_heating_density", "Scalar_L2");
  problem_builder->RegisterJouleHeatingDensityAux(
      "joule_heating_density", "electric_field", "current_density");

  hephaestus::Coefficients coefficients = defineCoefficients();
  problem_builder->SetCoefficients(coefficients);

  hephaestus::Sources sources = defineSources();
  problem_builder->SetSources(sources);

  hephaestus::Outputs outputs = defineOutputs();
  problem_builder->SetOutputs(outputs);
  {
    // Call LineSampler to save values
    std::string gridfunction_name("magnetic_flux_density");
    std::string csv_name("SimulatedA1B1Transect.csv");
    const int num_pts = 100;
    // Mesh bounding box (for the full serial mesh).
    mfem::Vector pos_min, pos_max;
    mesh.GetBoundingBox(pos_min, pos_max, 1);
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

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-16));
    solver_options.SetParam("MaxIter", (unsigned int)1000);
    problem_builder->SetSolverOptions(solver_options);

    problem_builder->FinalizeProblem();

    auto problem = problem_builder->ReturnProblem();
    hephaestus::InputParameters exec_params;
    exec_params.SetParam("TimeStep", float(0.001));
    exec_params.SetParam("StartTime", float(0.00));
    exec_params.SetParam("EndTime", float(0.002));
    exec_params.SetParam("VisualisationSteps", int(1));
    exec_params.SetParam("Problem", static_cast<hephaestus::TimeDomainProblem *>(problem.get()));

    auto executioner = std::make_unique<hephaestus::TransientExecutioner>(exec_params);

    hephaestus::logger.info("Created executioner");
    executioner->Execute();
  }
  MPI_Finalize();
}
