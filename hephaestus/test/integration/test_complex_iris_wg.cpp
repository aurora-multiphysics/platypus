#include "hephaestus.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

extern const char * DATA_DIR;

class TestComplexIrisWaveguide
{
protected:
  static void EBCR(const mfem::Vector & x, mfem::Vector & E)
  {
    E.SetSize(3);
    E = 0.0;
  }

  static void EBCI(const mfem::Vector & x, mfem::Vector & E)
  {
    E.SetSize(3);
    E = 0.0;
  }

  inline static const double epsilon0_ = 8.8541878176e-12; // F/m;
  inline static const double mu0_ = 4.0e-7 * M_PI;         // H/m;
  inline static const double freq_ = 9.3e9;                // 10/2pi
  double _port_length_vector[3] = {0.0, 22.86e-3, 0.0};
  double _port_width_vector[3] = {0.0, 0.0, 10.16e-3};

  hephaestus::InputParameters TestParams()
  {
    hephaestus::Subdomain air("air", 1);

    air._scalar_coefficients.Register("real_electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(0.0));
    air._scalar_coefficients.Register("imag_electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(0.0));
    air._scalar_coefficients.Register("real_rel_permittivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1.0));
    air._scalar_coefficients.Register("imag_rel_permittivity",
                                      std::make_shared<mfem::ConstantCoefficient>(0.0));
    air._scalar_coefficients.Register("real_rel_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(1.0));
    air._scalar_coefficients.Register("imag_rel_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(0.0));

    hephaestus::Coefficients coefficients(std::vector<hephaestus::Subdomain>({air}));

    coefficients._scalars.Register("frequency", std::make_shared<mfem::ConstantCoefficient>(freq_));
    coefficients._scalars.Register("magnetic_permeability",
                                   std::make_shared<mfem::ConstantCoefficient>(mu0_));
    coefficients._scalars.Register("dielectric_permittivity",
                                   std::make_shared<mfem::ConstantCoefficient>(epsilon0_));
    coefficients._scalars.Register("electrical_conductivity",
                                   std::make_shared<mfem::ConstantCoefficient>(0.0));

    coefficients._vectors.Register("EBCR",
                                   std::make_shared<mfem::VectorFunctionCoefficient>(3, EBCR));
    coefficients._vectors.Register("EBCI",
                                   std::make_shared<mfem::VectorFunctionCoefficient>(3, EBCI));

    hephaestus::BCMap bc_map;
    mfem::Array<int> dirichlet_attr(1);
    dirichlet_attr[0] = 1;
    bc_map.Register(
        "tangential_E",
        std::make_shared<hephaestus::VectorDirichletBC>(std::string("electric_field"),
                                                        dirichlet_attr,
                                                        coefficients._vectors.Get("EBCR"),
                                                        coefficients._vectors.Get("EBCI")));

    mfem::Array<int> wgi_in_attr(1);
    wgi_in_attr[0] = 2;
    bc_map.Register("WaveguidePortIn",
                    std::make_shared<hephaestus::RWTE10PortRBC>(std::string("electric_field"),
                                                                wgi_in_attr,
                                                                freq_,
                                                                _port_length_vector,
                                                                _port_width_vector,
                                                                true));

    mfem::Array<int> wgi_out_attr(1);
    wgi_out_attr[0] = 3;
    bc_map.Register("WaveguidePortOut",
                    std::make_shared<hephaestus::RWTE10PortRBC>(std::string("electric_field"),
                                                                wgi_out_attr,
                                                                freq_,
                                                                _port_length_vector,
                                                                _port_width_vector,
                                                                false));

    mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./irises.g")).c_str(), 1, 1);

    hephaestus::Outputs outputs;
    outputs.Register("VisItDataCollection",
                     std::make_shared<mfem::VisItDataCollection>("Hertz-AMR-Parallel-VisIt"));

    hephaestus::FESpaces fespaces;
    hephaestus::GridFunctions gridfunctions;
    hephaestus::AuxSolvers postprocessors;
    hephaestus::AuxSolvers preprocessors;
    hephaestus::Sources sources;

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-16));
    solver_options.SetParam("MaxIter", (unsigned int)1000);

    hephaestus::InputParameters params;
    params.SetParam("UseGLVis", true);

    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("FESpaces", fespaces);
    params.SetParam("GridFunctions", gridfunctions);
    params.SetParam("PreProcessors", preprocessors);
    params.SetParam("PostProcessors", postprocessors);
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);
    params.SetParam("SolverOptions", solver_options);

    return params;
  }
};

TEST_CASE_METHOD(TestComplexIrisWaveguide, "TestComplexIrisWaveguide", "[CheckRun]")
{
  hephaestus::InputParameters params(TestParams());
  std::shared_ptr<mfem::ParMesh> pmesh =
      std::make_shared<mfem::ParMesh>(params.GetParam<mfem::ParMesh>("Mesh"));

  auto problem_builder =
      std::make_unique<hephaestus::ComplexEFormulation>("magnetic_reluctivity",
                                                        "electrical_conductivity",
                                                        "dielectric_permittivity",
                                                        "frequency",
                                                        "electric_field",
                                                        "electric_field_real",
                                                        "electric_field_imag");

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

  problem_builder->AddFESpace("HDiv", "RT_3D_P0");
  problem_builder->AddFESpace("Scalar_L2", "L2_3D_P0");

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

  problem_builder->SetCoefficients(coefficients);
  problem_builder->SetPostprocessors(postprocessors);
  problem_builder->SetSources(sources);
  problem_builder->SetOutputs(outputs);
  problem_builder->SetSolverOptions(solver_options);

  problem_builder->FinalizeProblem();

  auto problem = problem_builder->ReturnProblem();

  hephaestus::InputParameters exec_params;
  exec_params.SetParam("Problem", static_cast<hephaestus::SteadyStateProblem *>(problem.get()));

  auto executioner = std::make_unique<hephaestus::SteadyExecutioner>(exec_params);

  executioner->Execute();

  mfem::Vector zero_vec(3);
  zero_vec = 0.0;
  mfem::VectorConstantCoefficient zero_coef(zero_vec);

  double norm_r = problem->_gridfunctions.Get("electric_field_real")->ComputeMaxError(zero_coef);
  double norm_i = problem->_gridfunctions.Get("electric_field_imag")->ComputeMaxError(zero_coef);
  REQUIRE_THAT(norm_r, Catch::Matchers::WithinAbs(4896.771, 0.001));
  REQUIRE_THAT(norm_i, Catch::Matchers::WithinAbs(5357.650, 0.001));
}
