#include "problem_builder.hpp"
#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

TEST_CASE("VariablesTest", "[CheckSetup]")
{
  mfem::Mesh mesh((std::string(DATA_DIR) + std::string("./beam-tet.mesh")).c_str(), 1, 1);

  auto problem_builder = std::make_unique<hephaestus::TimeDomainProblemBuilder>();

  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);

  problem_builder->SetMesh(pmesh);
  problem_builder->AddFESpace(std::string("HCurl"), std::string("ND_3D_P2"));
  problem_builder->AddGridFunction(std::string("vector_potential"), std::string("HCurl"));

  auto problem = problem_builder->ReturnProblem();

  auto stored_gf = problem->_gridfunctions.Get("vector_potential");
  auto stored_fespace = problem->_fespaces.Get("HCurl");

  REQUIRE(stored_fespace->GetVSize() == stored_gf->ParFESpace()->GetVSize());
  REQUIRE(stored_fespace->GetVSize() > 0);
}
