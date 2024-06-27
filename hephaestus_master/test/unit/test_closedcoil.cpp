#include "closed_coil.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

extern const char * DATA_DIR;

TEST_CASE("ClosedCoilTest", "[CheckData]")
{

  // Error tolerance
  const double eps{1e-2};

  int order = 1;

  mfem::Mesh mesh((std::string(DATA_DIR) + "team7.g").c_str(), 1, 1);
  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);

  mfem::ND_FECollection h_curl_collection(order, pmesh.get()->Dimension());

  auto h_curl_fe_space =
      std::make_shared<mfem::ParFiniteElementSpace>(pmesh.get(), &h_curl_collection);
  auto grad_phi = std::make_shared<mfem::ParGridFunction>(h_curl_fe_space.get());

  const double ival = 10.0;
  const double cond_val = 1e6;

  auto itot = std::make_shared<mfem::ConstantCoefficient>(ival);
  auto conductivity = std::make_shared<mfem::ConstantCoefficient>(cond_val);

  hephaestus::BCMap bc_maps;

  hephaestus::Coefficients coefficients;
  coefficients._scalars.Register(std::string("Itotal"), itot);
  coefficients._scalars.Register(std::string("Conductivity"), conductivity);

  hephaestus::FESpaces fespaces;
  fespaces.Register(std::string("HCurl"), h_curl_fe_space);

  hephaestus::GridFunctions gridfunctions;
  gridfunctions.Register(std::string("GradPhi"), grad_phi);

  int elec_attr = 7;
  int test_attr = 8;
  mfem::Array<int> submesh_domains({3, 4, 5, 6});

  hephaestus::ClosedCoilSolver closedcoil(
      "GradPhi", "HCurl", "H1", "Itotal", "Conductivity", submesh_domains, elec_attr, true);
  closedcoil.Init(gridfunctions, fespaces, bc_maps, coefficients);
  mfem::ParLinearForm dummy(h_curl_fe_space.get());
  closedcoil.Apply(&dummy);

  double flux = hephaestus::calcFlux(grad_phi.get(), test_attr, *conductivity);

  REQUIRE_THAT(-flux, Catch::Matchers::WithinAbs(ival, eps));
}
