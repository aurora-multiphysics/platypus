// Implements the TEAM Problem 4 benchmark (FELIX brick)
// available from https://www.compumag.org/wp/team/.
// Reference results available from
// A. Kameari, Results for benchmark calculations of problem 4
// (the felix brick), COMPEL, Vol 7. Num 1 (1988).

#include "hephaestus.hpp"
const char * DATA_DIR = "../../data/";

static void
external_db_dt(const mfem::Vector & xv, double t, mfem::Vector & db_dt)
{
  // External magnetic flux density B = B0*exp(-t/tau)
  double b0(0.1);     // Initial external magnetic field (T)
  double tau(0.0119); // Time constant (s)

  // Calculate rate of change of external flux density dB/dt
  db_dt(0) = 0.0;
  db_dt(1) = 0.0;
  db_dt(2) = -(b0 / tau) * exp(-t / tau);
}

static double
external_dpsi_dt(const mfem::Vector & xv, double t)
{
  // Calculate magnetic potential for external source field
  // External magnetic flux density B = B0*exp(-t/tau)
  mfem::Vector dbext_dt(3);
  external_db_dt(xv, t, dbext_dt);

  return xv(2) * dbext_dt(2);
}

static void
boundary_dh_dt(const mfem::Vector & x, double t, mfem::Vector & dh_dt)
{
  // Assumes boundary is sufficiently far from sources
  dh_dt(0) = 0.0;
  dh_dt(1) = 0.0;
  dh_dt(2) = 0.0;
}

hephaestus::Coefficients
defineCoefficients()
{
  hephaestus::Subdomain brick("brick", 1);
  brick._scalar_coefficients.Register("electric_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(2.538e7));
  hephaestus::Subdomain vacuum("vacuum", 2);
  vacuum._scalar_coefficients.Register("electric_conductivity",
                                       std::make_shared<mfem::ConstantCoefficient>(1.0));

  hephaestus::Coefficients coefficients(std::vector<hephaestus::Subdomain>({brick, vacuum}));
  coefficients._scalars.Register("magnetic_permeability",
                                 std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  auto dpsi_dt_coef = std::make_shared<mfem::FunctionCoefficient>(external_dpsi_dt);

  // Register to prevent dpsi_dt_coef being destroyed when it goes out of scope.
  coefficients._scalars.Register("magnetic_potential_time_derivative", std::move(dpsi_dt_coef));

  auto dh_dt_vec_coef = std::make_shared<mfem::VectorFunctionCoefficient>(3, boundary_dh_dt);
  coefficients._vectors.Register("surface_tangential_dHdt", dh_dt_vec_coef);

  return coefficients;
}

hephaestus::Sources
defineSources()
{
  hephaestus::InputParameters source_solver_options;
  source_solver_options.SetParam("Tolerance", float(1.0e-20));
  source_solver_options.SetParam("MaxIter", (unsigned int)2000);

  hephaestus::Sources sources;
  sources.Register("source",
                   std::make_shared<hephaestus::ScalarPotentialSource>("dhext_dt",
                                                                       "dmagnetic_potential_dt",
                                                                       "HCurl",
                                                                       "H1",
                                                                       "_one",
                                                                       -1.0,
                                                                       source_solver_options));

  return sources;
}

hephaestus::Outputs
defineOutputs()
{
  hephaestus::Outputs outputs;
  outputs.Register("ParaViewDataCollection",
                   std::make_shared<mfem::ParaViewDataCollection>("Team4ParaView"));
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
  auto problem_builder = std::make_unique<hephaestus::HFormulation>(
      "electric_resistivity", "electric_conductivity", "magnetic_permeability", "magnetic_field");
  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./team4_symmetrized.g")).c_str(), 1, 1);
  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);

  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace("H1", "H1_3D_P1");
  problem_builder->AddFESpace("HCurl", "ND_3D_P1");
  problem_builder->AddFESpace("HDiv", "RT_3D_P0");
  problem_builder->AddFESpace("Scalar_L2", "L2_3D_P0");
  problem_builder->AddFESpace("Vector_L2", "L2_3D_P0", 3);
  problem_builder->AddGridFunction("magnetic_field", "HCurl");

  problem_builder->AddGridFunction("dmagnetic_potential_dt", "H1");
  problem_builder->AddGridFunction("magnetic_flux_density", "HDiv");
  problem_builder->RegisterMagneticFluxDensityAux("magnetic_flux_density");

  problem_builder->AddGridFunction("current_density", "HDiv");
  problem_builder->RegisterCurrentDensityAux("current_density");

  problem_builder->AddGridFunction("electric_field", "HCurl");
  problem_builder->RegisterElectricFieldAux("electric_field");

  hephaestus::Coefficients coefficients = defineCoefficients();
  problem_builder->SetCoefficients(coefficients);

  hephaestus::Sources sources = defineSources();
  problem_builder->SetSources(sources);

  hephaestus::Outputs outputs = defineOutputs();
  problem_builder->SetOutputs(outputs);

  problem_builder->AddBoundaryCondition("tangential_dhdt_bc",
                                        std::make_shared<hephaestus::VectorDirichletBC>(
                                            "dmagnetic_field_dt",
                                            mfem::Array<int>({1, 2, 5, 6}),
                                            coefficients._vectors.Get("surface_tangential_dHdt")));
  problem_builder->AddBoundaryCondition(
      "magnetic_potential_bc",
      std::make_shared<hephaestus::ScalarDirichletBC>(
          "dmagnetic_potential_dt",
          mfem::Array<int>({1, 2}),
          coefficients._scalars.Get("magnetic_potential_time_derivative")));

  auto fluxmonitor = std::make_shared<hephaestus::FluxMonitorAux>("current_density", 3);
  fluxmonitor->SetPriority(2);
  problem_builder->AddPostprocessor("FluxMonitor", fluxmonitor);

  hephaestus::InputParameters solver_options;
  solver_options.SetParam("AbsTolerance", float(1.0e-20));
  solver_options.SetParam("Tolerance", float(1.0e-20));
  solver_options.SetParam("MaxIter", (unsigned int)500);
  problem_builder->SetSolverOptions(solver_options);

  problem_builder->FinalizeProblem();

  auto problem = problem_builder->ReturnProblem();
  hephaestus::InputParameters exec_params;
  exec_params.SetParam("TimeStep", float(0.001));
  exec_params.SetParam("StartTime", float(0.00));
  exec_params.SetParam("EndTime", float(0.02));
  exec_params.SetParam("VisualisationSteps", int(1));
  exec_params.SetParam("Problem", static_cast<hephaestus::TimeDomainProblem *>(problem.get()));

  auto executioner = std::make_unique<hephaestus::TransientExecutioner>(exec_params);

  executioner->Execute();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double current;
  double t;
  for (std::size_t i = 0; i < fluxmonitor->_times.Size(); ++i)
  {
    if (rank == 0)
    {
      current = -2 * fluxmonitor->_fluxes[i];
      t = fluxmonitor->_times[i];
      hephaestus::logger.info("t = {} s, I = {} A", t, current);
    }
  }

  MPI_Finalize();
}
