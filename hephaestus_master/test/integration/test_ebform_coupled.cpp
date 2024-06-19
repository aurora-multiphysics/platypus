#include "hephaestus.hpp"
#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

class CopperConductivityCoefficient : public hephaestus::CoupledCoefficient
{
public:
  CopperConductivityCoefficient(const hephaestus::InputParameters & params)
    : hephaestus::CoupledCoefficient(params){};
  double Eval(mfem::ElementTransformation & T, const mfem::IntegrationPoint & ip) override
  {
    return 2.0 * M_PI * 10 + 0.001 * (_gf->GetValue(T, ip));
  }
};

class TestEBFormCoupled
{
protected:
  static double PotentialHigh(const mfem::Vector & x, double t)
  {
    double wj(2.0 * M_PI / 60.0);
    return 2 * cos(wj * t);
  }
  static double PotentialGround(const mfem::Vector & x, double t) { return 0.0; }
  static void EdotBc(const mfem::Vector & x, mfem::Vector & E) { E = 0.0; }
  static void JSrc(const mfem::Vector & x, double t, mfem::Vector & j)
  {
    double wj(2.0 * M_PI / 60.0);
    j[0] = 0.0;
    j[1] = 0.0;
    j[2] = 2 * sin(wj * t);
  }

  hephaestus::InputParameters TestParams()
  {
    double sigma = 2.0 * M_PI * 10;

    double sigma_air;

    sigma_air = 1.0e-6 * sigma;

    hephaestus::FESpaces fespaces;
    hephaestus::GridFunctions gridfunctions;

    // materialCopper instances and get property coefs? init can be for
    // all...
    // CoupledCoefficients must also be added to AuxSolvers
    hephaestus::InputParameters copper_conductivity_params;
    copper_conductivity_params.SetParam("CoupledVariableName", std::string("temperature"));

    std::shared_ptr<hephaestus::CoupledCoefficient> wire_conductivity =
        std::make_shared<CopperConductivityCoefficient>(copper_conductivity_params);

    hephaestus::Subdomain wire("wire", 1);
    wire._scalar_coefficients.Register("electrical_conductivity", wire_conductivity);

    hephaestus::Subdomain air("air", 2);
    air._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(sigma_air));

    hephaestus::Coefficients coefficients(std::vector<hephaestus::Subdomain>({wire, air}));

    // coefficients.scalars.Register(
    //     "electrical_conductivity",
    //     new mfem::PWCoefficient(coefficients.getGlobalScalarProperty(
    //         std::string("electrical_conductivity"))),
    //     true);

    hephaestus::AuxSolvers preprocessors;
    preprocessors.Register("CoupledCoefficient", wire_conductivity);

    hephaestus::BCMap bc_map;

    auto edot_vec_coef = std::make_shared<mfem::VectorFunctionCoefficient>(3, EdotBc);
    coefficients._vectors.Register("surface_tangential_dEdt", edot_vec_coef);

    bc_map.Register("tangential_dEdt",
                    std::make_shared<hephaestus::VectorDirichletBC>(std::string("electric_field"),
                                                                    mfem::Array<int>({1, 2, 3}),
                                                                    edot_vec_coef.get()));
    coefficients._scalars.Register("magnetic_permeability",
                                   std::make_shared<mfem::ConstantCoefficient>(1.0));

    mfem::Array<int> high_terminal(1);
    high_terminal[0] = 1;

    auto potential_src = std::make_shared<mfem::FunctionCoefficient>(PotentialHigh);
    coefficients._scalars.Register("source_potential", potential_src);

    bc_map.Register("high_potential",
                    std::make_shared<hephaestus::ScalarDirichletBC>(
                        std::string("electric_potential"), high_terminal, potential_src.get()));

    auto potential_ground = std::make_shared<mfem::FunctionCoefficient>(PotentialGround);
    coefficients._scalars.Register("ground_potential", potential_ground);

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 2;
    bc_map.Register(
        "ground_potential",
        std::make_shared<hephaestus::ScalarDirichletBC>(
            std::string("electric_potential"), ground_terminal, potential_ground.get()));

    mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./cylinder-hex-q2.gen")).c_str(), 1, 1);

    hephaestus::Outputs outputs;
    outputs.Register("VisItDataCollection",
                     std::make_shared<mfem::VisItDataCollection>("EBFormVisIt"));
    outputs.Register("ParaViewDataCollection",
                     std::make_shared<mfem::ParaViewDataCollection>("EBFormParaView"));

    hephaestus::AuxSolvers postprocessors;
    hephaestus::Sources sources;

    hephaestus::InputParameters current_solver_options;
    current_solver_options.SetParam("Tolerance", float(1.0e-9));
    current_solver_options.SetParam("MaxIter", (unsigned int)1000);
    sources.Register("source",
                     std::make_shared<hephaestus::ScalarPotentialSource>("source",
                                                                         "electric_potential",
                                                                         "HCurlSource",
                                                                         "H1Source",
                                                                         "electrical_conductivity",
                                                                         -1,
                                                                         current_solver_options));

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-9));
    solver_options.SetParam("MaxIter", (unsigned int)1000);

    hephaestus::InputParameters params;
    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("FESpaces", fespaces);
    params.SetParam("GridFunctions", gridfunctions);
    params.SetParam("PreProcessors", preprocessors);
    params.SetParam("PostProcessors", postprocessors);
    params.SetParam("Sources", sources);
    params.SetParam("Outputs", outputs);
    params.SetParam("SolverOptions", solver_options);

    return params;
  }
};

TEST_CASE_METHOD(TestEBFormCoupled, "TestEBFormCoupled", "[CheckRun]")
{
  hephaestus::InputParameters params(TestParams());

  auto problem_builder = std::make_unique<hephaestus::EBDualFormulation>("magnetic_reluctivity",
                                                                         "magnetic_permeability",
                                                                         "electrical_conductivity",
                                                                         "electric_field",
                                                                         "magnetic_flux_density");

  auto bc_map(params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
  auto coefficients(params.GetParam<hephaestus::Coefficients>("Coefficients"));
  auto preprocessors(params.GetParam<hephaestus::AuxSolvers>("PreProcessors"));
  auto postprocessors(params.GetParam<hephaestus::AuxSolvers>("PostProcessors"));
  auto sources(params.GetParam<hephaestus::Sources>("Sources"));
  auto outputs(params.GetParam<hephaestus::Outputs>("Outputs"));
  auto solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
      "SolverOptions", hephaestus::InputParameters()));

  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(params.GetParam<mfem::ParMesh>("Mesh"));
  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("L2"), std::string("L2_3D_P1"));
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P1"));
  problem_builder->AddFESpace(std::string("HDiv"), std::string("RT_3D_P0"));
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P1"));
  problem_builder->AddFESpace(std::string("H1Source"), std::string("H1_3D_P2"));
  problem_builder->AddFESpace(std::string("HCurlSource"), std::string("ND_3D_P2"));
  problem_builder->AddGridFunction(std::string("current_density"), std::string("L2"));
  problem_builder->AddGridFunction(std::string("joule_heating_load"), std::string("L2"));
  problem_builder->AddGridFunction(std::string("temperature"), std::string("H1"));
  problem_builder->RegisterCurrentDensityAux("current_density");
  problem_builder->RegisterJouleHeatingDensityAux(
      "joule_heating_load", "electric_field", "current_density");
  problem_builder->SetBoundaryConditions(bc_map);
  problem_builder->SetAuxSolvers(preprocessors);
  problem_builder->SetCoefficients(coefficients);
  problem_builder->SetPostprocessors(postprocessors);
  problem_builder->SetSources(sources);
  problem_builder->SetOutputs(outputs);
  problem_builder->SetSolverOptions(solver_options);

  problem_builder->FinalizeProblem();

  auto problem = problem_builder->ReturnProblem();
  hephaestus::InputParameters exec_params;
  exec_params.SetParam("TimeStep", float(0.5));
  exec_params.SetParam("StartTime", float(0.00));
  exec_params.SetParam("EndTime", float(2.5));
  exec_params.SetParam("VisualisationSteps", int(1));
  exec_params.SetParam("Problem", static_cast<hephaestus::TimeDomainProblem *>(problem.get()));

  auto executioner = std::make_unique<hephaestus::TransientExecutioner>(exec_params);

  executioner->Execute();
}
