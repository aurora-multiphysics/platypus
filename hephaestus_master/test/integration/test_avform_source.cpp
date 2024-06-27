// Based on an H form MMS test provided by Joseph Dean
#include "hephaestus.hpp"
#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

class TestAVFormSource
{
protected:
  const int _var_order{2};
  static double EstimateConvergenceRate(
      HYPRE_BigInt n_i, HYPRE_BigInt n_imo, double error_i, double error_imo, int dim)
  {
    return std::log(error_i / error_imo) /
           std::log(std::pow(n_imo / static_cast<double>(n_i), 1.0 / dim));
  }

  static double PotentialGround(const mfem::Vector & x, double t) { return -x(0); }

  static void AdotBC(const mfem::Vector & x, double t, mfem::Vector & H)
  {
    H(0) = 1 + sin(x(1) * M_PI) * sin(x(2) * M_PI);
    H(1) = 0;
    H(2) = 0;
  }

  static void AExactExpr(const mfem::Vector & x, double t, mfem::Vector & A_exact)
  {
    A_exact(0) = (1 + sin(x(1) * M_PI) * sin(x(2) * M_PI)) * t;
    A_exact(1) = 0;
    A_exact(2) = 0;
  }
  static double MuExpr(const mfem::Vector & x)
  {
    double variation_scale = 0.0;
    double mu = 1.0 / (1.0 + variation_scale * cos(M_PI * x(0)) * cos(M_PI * x(1)));
    return mu;
  }

  static void SourceField(const mfem::Vector & x, double t, mfem::Vector & f)
  {
    f(0) = sin(M_PI * x(1)) * sin(M_PI * x(2)) * (t * 2 * M_PI * M_PI + 1);
    f(1) = 0.0;
    f(2) = 0.0;
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

    coefficients._scalars.Register("magnetic_permeability",
                                   std::make_shared<mfem::FunctionCoefficient>(MuExpr));

    hephaestus::BCMap bc_map;

    auto adot_vec_coef = std::make_shared<mfem::VectorFunctionCoefficient>(3, AdotBC);
    coefficients._vectors.Register("surface_tangential_dAdt", adot_vec_coef);

    bc_map.Register("tangential_dAdt",
                    std::make_shared<hephaestus::VectorDirichletBC>(
                        std::string("dmagnetic_vector_potential_dt"),
                        mfem::Array<int>({1, 2, 3}),
                        adot_vec_coef.get()));

    coefficients._scalars.Register("electrical_conductivity",
                                   std::make_shared<mfem::ConstantCoefficient>(1.0));

    mfem::Array<int> ground_terminal(1);
    ground_terminal[0] = 1;

    auto ground_coeff = std::make_shared<mfem::FunctionCoefficient>(PotentialGround);
    coefficients._scalars.Register("ground_potential", ground_coeff);

    bc_map.Register(
        "ground_potential",
        std::make_shared<hephaestus::ScalarDirichletBC>(
            std::string("electric_potential"), mfem::Array<int>({1, 2, 3}), ground_coeff.get()));

    auto a_exact = std::make_shared<mfem::VectorFunctionCoefficient>(3, AExactExpr);
    coefficients._vectors.Register("a_exact_coeff", a_exact);

    mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./beam-tet.mesh")).c_str(), 1, 1);

    hephaestus::Outputs outputs;
    outputs.Register("VisItDataCollection",
                     std::make_shared<mfem::VisItDataCollection>("AVFormVisIt"));

    hephaestus::InputParameters l2errpostprocparams;
    l2errpostprocparams.SetParam("VariableName", std::string("magnetic_vector_potential"));
    l2errpostprocparams.SetParam("VectorCoefficientName", std::string("a_exact_coeff"));
    hephaestus::AuxSolvers postprocessors;
    postprocessors.Register(
        "L2ErrorPostprocessor",
        std::make_shared<hephaestus::L2ErrorVectorPostprocessor>(l2errpostprocparams));

    auto vec_coef_aux = std::make_shared<hephaestus::VectorCoefficientAux>(
        "analytic_vector_potential", "a_exact_coeff");
    vec_coef_aux->SetPriority(-1);
    postprocessors.Register("VectorCoefficientAux", vec_coef_aux);

    hephaestus::AuxSolvers preprocessors;

    hephaestus::Sources sources;
    auto j_src_coef = std::make_shared<mfem::VectorFunctionCoefficient>(3, SourceField);
    coefficients._vectors.Register("source", j_src_coef);

    sources.Register("source",
                     std::make_shared<hephaestus::DivFreeSource>("source",
                                                                 "source",
                                                                 "_HCurlFESpace",
                                                                 "_H1FESpace",
                                                                 "_source_potential",
                                                                 hephaestus::InputParameters(),
                                                                 false));

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
    params.SetParam("Outputs", outputs);
    params.SetParam("Sources", sources);

    return params;
  }
};

TEST_CASE_METHOD(TestAVFormSource, "TestAVFormSource", "[CheckRun]")
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

    auto problem_builder = std::make_unique<hephaestus::AVFormulation>("magnetic_reluctivity",
                                                                       "magnetic_permeability",
                                                                       "electrical_conductivity",
                                                                       "magnetic_vector_potential",
                                                                       "electric_potential");

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
    problem_builder->AddGridFunction(std::string("analytic_vector_potential"),
                                     std::string("HCurl"));
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
