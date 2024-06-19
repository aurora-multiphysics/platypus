#include "hephaestus.hpp"
#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

class TestComplexAFormRod
{
protected:
  static double PotentialHigh(const mfem::Vector & x, double t)
  {
    // double wj_(2.0 * M_PI / 60.0);
    // return 2 * cos(wj_ * t);
    return 2.0;
  }
  static double PotentialGround(const mfem::Vector & x, double t) { return 0.0; }
  static void ABcR(const mfem::Vector & x, mfem::Vector & A)
  {
    A.SetSize(3);
    A = 0.0;
  }
  static void ABcI(const mfem::Vector & x, mfem::Vector & A)
  {
    A.SetSize(3);
    A = 0.0;
  }

  hephaestus::InputParameters TestParams()
  {
    double sigma = 2.0 * M_PI * 10;

    double sigma_air;

    sigma_air = 1.0e-6 * sigma;

    hephaestus::Subdomain wire("wire", 1);
    wire._scalar_coefficients.Register("electrical_conductivity",
                                       std::make_shared<mfem::ConstantCoefficient>(sigma));

    hephaestus::Subdomain air("air", 2);
    air._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(sigma_air));

    hephaestus::Coefficients coefficients(std::vector<hephaestus::Subdomain>({wire, air}));

    coefficients._scalars.Register("frequency",
                                   std::make_shared<mfem::ConstantCoefficient>(1.0 / 60.0));
    coefficients._scalars.Register("dielectric_permittivity",
                                   std::make_shared<mfem::ConstantCoefficient>(0.0));
    coefficients._scalars.Register("magnetic_permeability",
                                   std::make_shared<mfem::ConstantCoefficient>(1.0));

    hephaestus::BCMap bc_map;

    coefficients._vectors.Register("ABcR",
                                   std::make_shared<mfem::VectorFunctionCoefficient>(3, ABcR));
    coefficients._vectors.Register("ABcI",
                                   std::make_shared<mfem::VectorFunctionCoefficient>(3, ABcI));

    bc_map.Register(
        "tangential_A",
        std::make_shared<hephaestus::VectorDirichletBC>(std::string("magnetic_vector_potential"),
                                                        mfem::Array<int>({1, 2, 3}),
                                                        coefficients._vectors.Get("ABcR"),
                                                        coefficients._vectors.Get("ABcI")));

    mfem::Array<int> high_terminal(1);
    high_terminal[0] = 1;

    auto potential_src = std::make_shared<mfem::FunctionCoefficient>(PotentialHigh);
    coefficients._scalars.Register("source_potential", potential_src);

    bc_map.Register("high_potential",
                    std::make_shared<hephaestus::ScalarDirichletBC>(
                        std::string("electric_potential"), high_terminal, potential_src.get()));

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 2;

    auto potential_ground = std::make_shared<mfem::FunctionCoefficient>(PotentialGround);
    coefficients._scalars.Register("ground_potential", potential_ground);

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

    hephaestus::GridFunctions gridfunctions;
    hephaestus::AuxSolvers preprocessors;
    hephaestus::AuxSolvers postprocessors;
    hephaestus::Sources sources;
    hephaestus::InputParameters current_solver_options;
    current_solver_options.SetParam("Tolerance", float(1.0e-9));
    current_solver_options.SetParam("MaxIter", (unsigned int)1000);
    sources.Register("source",
                     std::make_shared<hephaestus::ScalarPotentialSource>("source",
                                                                         "electric_potential",
                                                                         "HCurl",
                                                                         "H1",
                                                                         "electrical_conductivity",
                                                                         -1,
                                                                         current_solver_options));

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-9));
    solver_options.SetParam("MaxIter", (unsigned int)1000);

    hephaestus::InputParameters params;
    params.SetParam("UseGLVis", true);
    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("GridFunctions", gridfunctions);
    params.SetParam("PreProcessors", preprocessors);
    params.SetParam("PostProcessors", postprocessors);
    params.SetParam("Sources", sources);
    params.SetParam("Outputs", outputs);
    params.SetParam("SolverOptions", solver_options);

    return params;
  }
};

TEST_CASE_METHOD(TestComplexAFormRod, "TestComplexAFormRod", "[CheckRun]")
{
  hephaestus::InputParameters params(TestParams());
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(params.GetParam<mfem::ParMesh>("Mesh"));

  auto problem_builder =
      std::make_unique<hephaestus::ComplexAFormulation>("magnetic_reluctivity",
                                                        "electrical_conductivity",
                                                        "dielectric_permittivity",
                                                        "frequency",
                                                        "magnetic_vector_potential",
                                                        "magnetic_vector_potential_real",
                                                        "magnetic_vector_potential_imag");

  auto bc_map(params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
  auto coefficients(params.GetParam<hephaestus::Coefficients>("Coefficients"));
  auto preprocessors(params.GetParam<hephaestus::AuxSolvers>("PreProcessors"));
  auto postprocessors(params.GetParam<hephaestus::AuxSolvers>("PostProcessors"));
  auto sources(params.GetParam<hephaestus::Sources>("Sources"));
  auto outputs(params.GetParam<hephaestus::Outputs>("Outputs"));
  auto solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
      "SolverOptions", hephaestus::InputParameters()));

  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P1"));
  problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P1"));
  problem_builder->AddGridFunction("magnetic_vector_potential_real", "HCurl");
  problem_builder->AddGridFunction("magnetic_vector_potential_imag", "HCurl");
  problem_builder->SetBoundaryConditions(bc_map);
  problem_builder->SetAuxSolvers(preprocessors);
  problem_builder->SetCoefficients(coefficients);
  problem_builder->SetPostprocessors(postprocessors);

  problem_builder->AddFESpace("HDiv", "RT_3D_P0");
  problem_builder->AddFESpace("Scalar_L2", "L2_3D_P0");

  problem_builder->AddGridFunction("electric_field_real", "HCurl");
  problem_builder->AddGridFunction("electric_field_imag", "HCurl");
  problem_builder->RegisterElectricFieldAux("electric_field_real", "electric_field_imag");

  problem_builder->AddGridFunction("magnetic_flux_density_real", "HDiv");
  problem_builder->AddGridFunction("magnetic_flux_density_imag", "HDiv");
  problem_builder->RegisterMagneticFluxDensityAux("magnetic_flux_density_real",
                                                  "magnetic_flux_density_imag");

  problem_builder->AddGridFunction("current_density_real", "HDiv");
  problem_builder->AddGridFunction("current_density_imag", "HDiv");
  problem_builder->RegisterCurrentDensityAux("current_density_real", "current_density_imag");

  problem_builder->AddGridFunction("joule_heating_density", "Scalar_L2");
  problem_builder->RegisterJouleHeatingDensityAux("joule_heating_density",
                                                  "electric_field_real",
                                                  "electric_field_imag",
                                                  "current_density_real",
                                                  "current_density_imag");

  problem_builder->SetSources(sources);
  problem_builder->SetOutputs(outputs);
  problem_builder->SetSolverOptions(solver_options);

  problem_builder->FinalizeProblem();

  auto problem = problem_builder->ReturnProblem();

  hephaestus::InputParameters exec_params;
  exec_params.SetParam("Problem", static_cast<hephaestus::SteadyStateProblem *>(problem.get()));

  auto executioner = std::make_unique<hephaestus::SteadyExecutioner>(exec_params);

  executioner->Execute();
}
