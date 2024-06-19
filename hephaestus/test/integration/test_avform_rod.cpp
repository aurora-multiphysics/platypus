#include "hephaestus.hpp"
#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

class TestAVFormRod
{
protected:
  static double PotentialHigh(const mfem::Vector & x, double t)
  {
    double wj(2.0 * M_PI / 60.0);
    return 2 * cos(wj * t);
  }
  static double PotentialGround(const mfem::Vector & x, double t) { return 0.0; }

  static void AdotBc(const mfem::Vector & x, double t, mfem::Vector & dAdt) { dAdt = 0.0; }

  static void SourceCurrent(const mfem::Vector & x, double t, mfem::Vector & J) { J = 0.0; }

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

    coefficients._scalars.Register("magnetic_permeability",
                                   std::make_shared<mfem::ConstantCoefficient>(1.0));

    hephaestus::BCMap bc_map;
    auto adot_vec_coef = std::make_shared<mfem::VectorFunctionCoefficient>(3, AdotBc);
    coefficients._vectors.Register("surface_tangential_dAdt", adot_vec_coef);

    bc_map.Register("tangential_dAdt",
                    std::make_shared<hephaestus::VectorDirichletBC>(
                        std::string("dmagnetic_vector_potential_dt"),
                        mfem::Array<int>({1, 2, 3}),
                        adot_vec_coef.get()));

    coefficients._scalars.Register("high_potential_func",
                                   std::make_shared<mfem::FunctionCoefficient>(PotentialHigh));

    mfem::Array<int> high_terminal(1);
    high_terminal[0] = 1;
    bc_map.Register("high_potential",
                    std::make_shared<hephaestus::ScalarDirichletBC>(
                        std::string("electric_potential"),
                        high_terminal,
                        coefficients._scalars.Get("high_potential_func")));

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 2;

    coefficients._scalars.Register("ground_potential_func",
                                   std::make_shared<mfem::FunctionCoefficient>(PotentialGround));

    bc_map.Register("ground_potential",
                    std::make_shared<hephaestus::ScalarDirichletBC>(
                        std::string("electric_potential"),
                        ground_terminal,
                        coefficients._scalars.Get("ground_potential_func")));

    mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./cylinder-hex-q2.gen")).c_str(), 1, 1);

    hephaestus::Outputs outputs;
    outputs.Register("VisItDataCollection",
                     std::make_shared<mfem::VisItDataCollection>("AVFormVisIt"));
    outputs.Register("ParaViewDataCollection",
                     std::make_shared<mfem::ParaViewDataCollection>("AVFormParaView"));

    hephaestus::FESpaces fespaces;
    hephaestus::GridFunctions gridfunctions;
    hephaestus::AuxSolvers postprocessors;
    hephaestus::AuxSolvers preprocessors;
    hephaestus::Sources sources;

    hephaestus::InputParameters params;
    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("FESpaces", fespaces);
    params.SetParam("GridFunctions", gridfunctions);
    params.SetParam("PreProcessors", preprocessors);
    params.SetParam("PostProcessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);

    return params;
  }
};

TEST_CASE_METHOD(TestAVFormRod, "TestAVFormRod", "[CheckRun]")
{
  hephaestus::InputParameters params(TestParams());
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(params.GetParam<mfem::ParMesh>("Mesh"));

  auto problem_builder = std::make_unique<hephaestus::AVFormulation>("magnetic_reluctivity",
                                                                     "magnetic_permeability",
                                                                     "electrical_conductivity",
                                                                     "magnetic_vector_potential",
                                                                     "electric_potential");

  auto bc_map(params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
  auto coefficients(params.GetParam<hephaestus::Coefficients>("Coefficients"));
  auto preprocessors(params.GetParam<hephaestus::AuxSolvers>("PreProcessors"));
  auto postprocessors(params.GetParam<hephaestus::AuxSolvers>("PostProcessors"));
  auto sources(params.GetParam<hephaestus::Sources>("Sources"));
  auto outputs(params.GetParam<hephaestus::Outputs>("Outputs"));
  auto solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
      "SolverOptions", hephaestus::InputParameters()));

  problem_builder->SetMesh(pmesh);
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
