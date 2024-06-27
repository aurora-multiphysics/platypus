// Based on an H form MMS test provided by Joseph Dean
#include "hephaestus.hpp"
#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

class TestHFormSource
{
protected:
  const int _var_order{2};
  static double EstimateConvergenceRate(
      HYPRE_BigInt n_i, HYPRE_BigInt n_imo, double error_i, double error_imo, int dim)
  {
    return std::log(error_i / error_imo) /
           std::log(std::pow(n_imo / static_cast<double>(n_i), 1.0 / dim));
  }

  static double PotentialGround(const mfem::Vector & x, double t) { return 0.0; }

  static void HdotBc(const mfem::Vector & x, double t, mfem::Vector & H)
  {
    H(0) = sin(x(1) * M_PI) * sin(x(2) * M_PI);
    H(1) = 0;
    H(2) = 0;
  }

  static void HExactExpr(const mfem::Vector & x, double t, mfem::Vector & H_exact)
  {
    H_exact(0) = sin(x(1) * M_PI) * sin(x(2) * M_PI) * t;
    H_exact(1) = 0;
    H_exact(2) = 0;
  }
  static double SigmaExpr(const mfem::Vector & x)
  {
    double variation_scale = 0.5;
    double sigma = 1.0 / (1.0 + variation_scale * cos(M_PI * x(0)) * cos(M_PI * x(1)));
    return sigma;
  }
  // Source field
  static void SourceField(const mfem::Vector & x, double t, mfem::Vector & f)
  {
    double variation_scale = 0.5;
    f(0) = t * M_PI * M_PI * sin(M_PI * x(1)) * sin(M_PI * x(2)) *
               (3 * variation_scale * cos(M_PI * x(0)) * cos(M_PI * x(1)) + 2) +
           sin(M_PI * x(1)) * sin(M_PI * x(2));
    f(1) = -variation_scale * M_PI * M_PI * t * sin(M_PI * x(0)) * cos(M_PI * x(1)) *
           cos(M_PI * x(1)) * sin(M_PI * x(2));
    f(2) = -0.5 * variation_scale * M_PI * M_PI * t * sin(M_PI * x(0)) * sin(2 * M_PI * x(1)) *
           cos(M_PI * x(2));
  }

  hephaestus::InputParameters TestParams()
  {
    hephaestus::Subdomain wire("wire", 1);
    wire._scalar_coefficients.Register("electrical_conductivity",
                                       std::make_shared<mfem::ConstantCoefficient>(1.0));
    hephaestus::Subdomain air("air", 2);
    air._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1.0));

    hephaestus::Coefficients coefficients(std::vector<hephaestus::Subdomain>({wire, air}));

    coefficients._scalars.Register("electrical_conductivity",
                                   std::make_shared<mfem::FunctionCoefficient>(SigmaExpr));

    hephaestus::BCMap bc_map;

    auto hdot_vec_coef = std::make_shared<mfem::VectorFunctionCoefficient>(3, HdotBc);
    coefficients._vectors.Register("surface_tangential_dHdt", hdot_vec_coef);

    bc_map.Register(
        "tangential_dHdt",
        std::make_shared<hephaestus::VectorDirichletBC>(
            std::string("dmagnetic_field_dt"), mfem::Array<int>({1, 2, 3}), hdot_vec_coef.get()));
    coefficients._scalars.Register("magnetic_permeability",
                                   std::make_shared<mfem::ConstantCoefficient>(1.0));

    auto d_bdt_src_coef = std::make_shared<mfem::VectorFunctionCoefficient>(3, SourceField);
    coefficients._vectors.Register("source", d_bdt_src_coef);

    auto h_exact = std::make_shared<mfem::VectorFunctionCoefficient>(3, HExactExpr);
    coefficients._vectors.Register("h_exact_coeff", h_exact);

    mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./beam-tet.mesh")).c_str(), 1, 1);

    hephaestus::Outputs outputs;
    outputs.Register("VisItDataCollection",
                     std::make_shared<mfem::VisItDataCollection>("HFormVisIt"));
    outputs.Register("ParaViewDataCollection",
                     std::make_shared<mfem::ParaViewDataCollection>("HFormParaView"));

    hephaestus::InputParameters l2errpostprocparams;
    l2errpostprocparams.SetParam("VariableName", std::string("magnetic_field"));
    l2errpostprocparams.SetParam("VectorCoefficientName", std::string("h_exact_coeff"));
    hephaestus::AuxSolvers postprocessors;
    postprocessors.Register(
        "L2ErrorPostprocessor",
        std::make_shared<hephaestus::L2ErrorVectorPostprocessor>(l2errpostprocparams));

    hephaestus::AuxSolvers preprocessors;
    preprocessors.Register("VectorCoefficientAux",
                           std::make_shared<hephaestus::VectorCoefficientAux>(
                               "analytic_magnetic_field", "h_exact_coeff"));

    hephaestus::Sources sources;
    hephaestus::InputParameters current_solver_options;
    current_solver_options.SetParam("Tolerance", float(1.0e-12));
    current_solver_options.SetParam("MaxIter", (unsigned int)200);

    sources.Register("source",
                     std::make_shared<hephaestus::DivFreeSource>("source",
                                                                 "source",
                                                                 "_HCurlFESpace",
                                                                 "H1",
                                                                 "magnetic_potential",
                                                                 current_solver_options,
                                                                 false));

    hephaestus::InputParameters solver_options;
    solver_options.SetParam("Tolerance", float(1.0e-16));
    solver_options.SetParam("MaxIter", (unsigned int)1000);

    hephaestus::InputParameters params;
    params.SetParam("TimeStep", float(0.05));
    params.SetParam("StartTime", float(0.00));
    params.SetParam("EndTime", float(0.05));
    params.SetParam("VisualisationSteps", int(1));
    params.SetParam("UseGLVis", true);

    params.SetParam("Mesh", mfem::ParMesh(MPI_COMM_WORLD, mesh));
    params.SetParam("BoundaryConditions", bc_map);
    params.SetParam("Coefficients", coefficients);
    params.SetParam("PreProcessors", preprocessors);
    params.SetParam("PostProcessors", postprocessors);
    params.SetParam("Sources", sources);
    params.SetParam("Outputs", outputs);
    params.SetParam("SolverOptions", solver_options);

    return params;
  }
};

TEST_CASE_METHOD(TestHFormSource, "TestHFormSource", "[CheckRun]")
{
  hephaestus::InputParameters params(TestParams());
  auto unrefined_pmesh(params.GetParam<mfem::ParMesh>("Mesh"));

  int num_conv_refinements = 3;
  for (int par_ref_levels = 0; par_ref_levels < num_conv_refinements; ++par_ref_levels)
  {

    std::shared_ptr<mfem::ParMesh> pmesh = std::make_shared<mfem::ParMesh>(unrefined_pmesh);

    for (int l = 0; l < par_ref_levels; l++)
    {
      pmesh->UniformRefinement();
    }

    auto problem_builder = std::make_unique<hephaestus::HFormulation>("electrical_resistivity",
                                                                      "electrical_conductivity",
                                                                      "magnetic_permeability",
                                                                      "magnetic_field");

    auto bc_map(params.GetParam<hephaestus::BCMap>("BoundaryConditions"));
    auto coefficients(params.GetParam<hephaestus::Coefficients>("Coefficients"));
    //   hephaestus::FESpaces fespaces(
    //       params.GetParam<hephaestus::FESpaces>("FESpaces"));
    //   hephaestus::GridFunctions gridfunctions(
    //       params.GetParam<hephaestus::GridFunctions>("GridFunctions"));
    auto preprocessors(params.GetParam<hephaestus::AuxSolvers>("PreProcessors"));
    auto postprocessors(params.GetParam<hephaestus::AuxSolvers>("PostProcessors"));
    auto sources(params.GetParam<hephaestus::Sources>("Sources"));
    auto outputs(params.GetParam<hephaestus::Outputs>("Outputs"));
    auto solver_options(params.GetOptionalParam<hephaestus::InputParameters>(
        "SolverOptions", hephaestus::InputParameters()));

    problem_builder->SetMesh(pmesh);
    problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P2"));
    problem_builder->AddFESpace(std::string("H1"), std::string("H1_3D_P2"));
    problem_builder->AddGridFunction(std::string("analytic_magnetic_field"), std::string("HCurl"));
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
    exec_params.SetParam("TimeStep", float(0.05));
    exec_params.SetParam("StartTime", float(0.00));
    exec_params.SetParam("EndTime", float(0.05));
    exec_params.SetParam("VisualisationSteps", int(1));
    exec_params.SetParam("Problem", static_cast<hephaestus::TimeDomainProblem *>(problem.get()));

    auto executioner = std::make_unique<hephaestus::TransientExecutioner>(exec_params);

    executioner->Execute();
  }

  auto l2errpostprocessor =
      params.GetParam<hephaestus::AuxSolvers>("PostProcessors")
          .Get<hephaestus::L2ErrorVectorPostprocessor>("L2ErrorPostprocessor");

  double r;
  for (std::size_t i = 1; i < l2errpostprocessor->_ndofs.Size(); ++i)
  {
    r = EstimateConvergenceRate(l2errpostprocessor->_ndofs[i],
                                l2errpostprocessor->_ndofs[i - 1],
                                l2errpostprocessor->_l2_errs[i],
                                l2errpostprocessor->_l2_errs[i - 1],
                                3);
    hephaestus::logger.info("{}", r);
    REQUIRE(r > _var_order - 0.15);
    REQUIRE(r < _var_order + 1.0);
  }
}
