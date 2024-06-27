#include "hephaestus.hpp"

const char * DATA_DIR = "../../data/";

hephaestus::Coefficients
defineCoefficients()
{
  hephaestus::Subdomain air("air", 1);
  air._scalar_coefficients.Register("electrical_conductivity",
                                    std::make_shared<mfem::ConstantCoefficient>(1.0));
  air._scalar_coefficients.Register("magnetic_permeability",
                                    std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  hephaestus::Subdomain plate("plate", 2);
  plate._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(3.526e7));
  plate._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  hephaestus::Subdomain coil1("coil1", 3);
  coil1._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1.0));
  coil1._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  hephaestus::Subdomain coil2("coil2", 4);
  coil2._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1.0));
  coil2._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  hephaestus::Subdomain coil3("coil3", 5);
  coil3._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1.0));
  coil3._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  hephaestus::Subdomain coil4("coil4", 6);
  coil4._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1.0));
  coil4._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  hephaestus::Coefficients coefficients(
      std::vector<hephaestus::Subdomain>({air, plate, coil1, coil2, coil3, coil4}));
  coefficients._scalars.Register("frequency", std::make_shared<mfem::ConstantCoefficient>(200.0));
  coefficients._scalars.Register("dielectric_permittivity",
                                 std::make_shared<mfem::ConstantCoefficient>(8.854e-12));

  coefficients._scalars.Register("I", std::make_shared<mfem::ConstantCoefficient>(2742));

  return coefficients;
}

hephaestus::Sources
defineSources()
{
  hephaestus::InputParameters div_free_source_params;

  // This vector of subdomains will form the coil that we pass to
  // ClosedCoilSolver
  int order = 1;
  int electrode_attr = 7;
  std::string coil_attr = "3 4 5 6";
  mfem::Array<int> coil_domains;
  std::stringstream ss(coil_attr);
  int att;

  while (ss >> att)
    coil_domains.Append(att);

  hephaestus::Sources sources;
  sources.Register("source",
                   std::make_shared<hephaestus::ClosedCoilSolver>("source_grad_phi",
                                                                  "HCurl",
                                                                  "H1",
                                                                  "I",
                                                                  "electrical_conductivity",
                                                                  coil_domains,
                                                                  electrode_attr,
                                                                  true));
  return sources;
}

hephaestus::Outputs
defineOutputs()
{
  hephaestus::Outputs outputs;
  outputs.Register("ParaViewDataCollection",
                   std::make_shared<mfem::ParaViewDataCollection>("ClosedCoilParaView"));
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

  hephaestus::logger.set_level(spdlog::level::info);

  // Create Formulation
  auto problem_builder = std::make_unique<hephaestus::MagnetostaticFormulation>(
      "magnetic_reluctivity", "magnetic_permeability", "magnetic_vector_potential");

  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./team7.g")).c_str(), 1, 1);
  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);

  int par_ref_lvl = 1;
  for (int l = 0; l < par_ref_lvl; ++l)
    pmesh->UniformRefinement();

  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P1"));
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P1"));
  problem_builder->AddFESpace(std::string("HDiv"), std::string("RT_3D_P0"));
  problem_builder->AddGridFunction(std::string("magnetic_vector_potential"), std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("source_grad_phi"), std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("magnetic_flux_density"), std::string("HDiv"));
  problem_builder->RegisterMagneticFluxDensityAux("magnetic_flux_density");

  hephaestus::Coefficients coefficients = defineCoefficients();
  problem_builder->SetCoefficients(coefficients);

  hephaestus::Sources sources = defineSources();
  problem_builder->SetSources(sources);

  hephaestus::Outputs outputs = defineOutputs();
  problem_builder->SetOutputs(outputs);

  hephaestus::InputParameters solver_options;
  solver_options.SetParam("Tolerance", float(1.0e-13));
  solver_options.SetParam("AbsTolerance", float(1.0e-16));
  solver_options.SetParam("MaxIter", (unsigned int)500);
  problem_builder->SetSolverOptions(solver_options);

  problem_builder->FinalizeProblem();

  auto problem = problem_builder->ReturnProblem();
  hephaestus::InputParameters exec_params;
  exec_params.SetParam("VisualisationSteps", int(1));
  exec_params.SetParam("UseGLVis", true);
  exec_params.SetParam("Problem", static_cast<hephaestus::SteadyStateProblem *>(problem.get()));

  auto executioner = std::make_unique<hephaestus::SteadyExecutioner>(exec_params);

  hephaestus::logger.info("Created executioner");
  executioner->Execute();

  MPI_Finalize();
}