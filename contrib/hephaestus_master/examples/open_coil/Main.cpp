#include "hephaestus.hpp"

const char * DATA_DIR = "../../data/";

static void
constVec(const mfem::Vector & x, mfem::Vector & V)
{
  V = 1.0;
}

hephaestus::Coefficients
defineCoefficients(double Itotal)
{

  hephaestus::Coefficients coefficients;
  coefficients._scalars.Register("magnetic_permeability",
                                 std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));

  // Electrical conductivity
  coefficients._scalars.Register("electrical_conductivity",
                                 std::make_shared<mfem::ConstantCoefficient>(1.0));

  // Time-dependent current
  coefficients._scalars.Register("I", std::make_shared<mfem::ConstantCoefficient>(Itotal));

  return coefficients;
}

hephaestus::Sources
defineSources(std::pair<int, int> elec, mfem::Array<int> coil_domains)
{
  hephaestus::Sources sources;
  sources.Register("source",
                   std::make_shared<hephaestus::OpenCoilSolver>(
                       "grad_phi", "phi", "I", "electrical_conductivity", coil_domains, elec));
  return sources;
}

hephaestus::Outputs
defineOutputs()
{
  hephaestus::Outputs outputs;
  outputs.Register("ParaViewDataCollection",
                   std::make_shared<mfem::ParaViewDataCollection>("OpenCoilParaView"));
  return outputs;
}

int
main(int argc, char * argv[])
{

  // Refinement and order
  int par_ref_lvl = -1;
  int order = 1;

  // Total electrical current going around the coil. Must be nonzero, can be
  // changed later.
  double itotal = 10;

  // Attribute that defines the internal faces over which we apply the potential
  // difference
  std::pair<int, int> elec_attrs;

  // Mesh file
  std::string mesh_filename = "coil.gen";

  // Domain attributes of the coil to be solved and boundary attributes of the
  // electrodes
  std::string coil_attr = "1";
  std::string elec_bdr_attr = "1 2";

  mfem::OptionsParser args(argc, argv);
  args.AddOption(
      &DATA_DIR, "-dataDir", "--data_directory", "Directory storing input data for tests.");
  args.AddOption(&par_ref_lvl, "-ref", "--parallel-refinement", "Parallel refinement level.");
  args.AddOption(&order, "-o", "--order", "Base functions order");
  args.AddOption(&itotal, "-I", "--Itotal", "Total electrical current.");
  args.AddOption(&mesh_filename, "-f", "--mesh-filename", "Mesh file name");
  args.AddOption(&coil_attr,
                 "-cd",
                 "--coil-domains",
                 "List of coil domain attributes separated by spaces, e.g. \'1 3 4\'");
  args.AddOption(&elec_bdr_attr,
                 "-e",
                 "--electrode-attrs",
                 "List of electrode attributes separated by spaces, e.g. \'1 "
                 "2\'. Must be two values.");

  args.Parse();

  MPI_Init(&argc, &argv);

  // Set Mesh
  mfem::Mesh mesh((std::string(DATA_DIR) + mesh_filename).c_str(), 1, 1);
  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);

  for (int l = 0; l < par_ref_lvl; ++l)
    pmesh->UniformRefinement();

  // This vector of subdomains will form the coil that we pass to
  // OpenCoilSolver
  mfem::Array<int> coil_domains;

  // This vector of attributes will form the electrodes that we pass to
  // OpenCoilSolver
  mfem::Array<int> elec_bdr_array;

  // Parsing the string of attributes
  std::stringstream ss(coil_attr);
  int att;
  while (ss >> att)
    coil_domains.Append(att);

  std::stringstream ss_bdr(elec_bdr_attr);
  while (ss_bdr >> att)
    elec_bdr_array.Append(att);

  if (elec_bdr_array.Size() != 2)
    mfem::mfem_error("Electrode boundary attribute list must contain two attributes.");

  elec_attrs.first = elec_bdr_array[0];
  elec_attrs.second = elec_bdr_array[1];

  // Create Formulation
  auto problem_builder = std::make_unique<hephaestus::MagnetostaticFormulation>(
      "magnetic_reluctivity", "magnetic_permeability", "magnetic_vector_potential");

  // Set Mesh
  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P1"));
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P1"));
  problem_builder->AddFESpace(std::string("HDiv"), std::string("RT_3D_P0"));
  problem_builder->AddGridFunction(std::string("magnetic_vector_potential"), std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("phi"), std::string("H1"));
  problem_builder->AddGridFunction(std::string("grad_phi"), std::string("HCurl"));
  problem_builder->AddGridFunction(std::string("magnetic_flux_density"), std::string("HDiv"));
  problem_builder->RegisterMagneticFluxDensityAux("magnetic_flux_density");
  hephaestus::Coefficients coefficients = defineCoefficients(itotal);

  mfem::Array<int> a_dbc_bdr(3);
  a_dbc_bdr[0] = 1;
  a_dbc_bdr[1] = 2;
  a_dbc_bdr[2] = 4;

  auto magnetic_vector_func_coefficient =
      std::make_shared<mfem::VectorFunctionCoefficient>(3, constVec);
  coefficients._vectors.Register("magnetic_vector_func", magnetic_vector_func_coefficient);

  auto a_dbc = std::make_shared<hephaestus::VectorDirichletBC>(
      "magnetic_vector_potential", a_dbc_bdr, magnetic_vector_func_coefficient.get());

  problem_builder->AddBoundaryCondition("a_dbc", a_dbc);

  problem_builder->SetCoefficients(coefficients);

  hephaestus::Sources sources = defineSources(elec_attrs, coil_domains);
  problem_builder->SetSources(sources);

  hephaestus::Outputs outputs = defineOutputs();
  problem_builder->SetOutputs(outputs);

  hephaestus::InputParameters solver_options;
  solver_options.SetParam("Tolerance", float(1.0e-10));
  solver_options.SetParam("AbsTolerance", float(1.0e-10));
  solver_options.SetParam("MaxIter", (unsigned int)1000);
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