#include "open_coil.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

extern const char * DATA_DIR;

TEST_CASE("OpenCoilTest", "[CheckData]")
{

  // Floating point error tolerance
  const double eps{1e-10};

  int par_ref_lvl = -1;
  int order = 1;

  mfem::Mesh mesh((std::string(DATA_DIR) + "coil.gen").c_str(), 1, 1);
  auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, mesh);

  for (int l = 0; l < par_ref_lvl; ++l)
    pmesh->UniformRefinement();

  mfem::ND_FECollection h_curl_collection(order, pmesh.get()->Dimension());
  auto h_curl_fe_space =
      std::make_shared<mfem::ParFiniteElementSpace>(pmesh.get(), &h_curl_collection);
  auto e = std::make_shared<mfem::ParGridFunction>(h_curl_fe_space.get());

  mfem::H1_FECollection h1_collection(order, pmesh.get()->Dimension());
  auto h1_fe_space = std::make_shared<mfem::ParFiniteElementSpace>(pmesh.get(), &h1_collection);
  auto v = std::make_shared<mfem::ParGridFunction>(h1_fe_space.get());

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
  fespaces.Register(std::string("H1"), h1_fe_space);

  hephaestus::GridFunctions gridfunctions;
  gridfunctions.Register(std::string("E"), e);
  gridfunctions.Register(std::string("V"), v);

  std::pair<int, int> elec_attrs{1, 2};
  mfem::Array<int> submesh_domains;
  submesh_domains.Append(1);

  hephaestus::OpenCoilSolver opencoil(
      "E", "V", "Itotal", "Conductivity", submesh_domains, elec_attrs);
  opencoil.Init(gridfunctions, fespaces, bc_maps, coefficients);
  mfem::ParLinearForm dummy(h_curl_fe_space.get());
  opencoil.Apply(&dummy);

  //- sign comes from the direction of the outward facing normal relative to elec_attrs order
  double flux = -hephaestus::calcFlux(e.get(), elec_attrs.first, *conductivity);

  REQUIRE_THAT(flux, Catch::Matchers::WithinAbs(ival, eps));
}
