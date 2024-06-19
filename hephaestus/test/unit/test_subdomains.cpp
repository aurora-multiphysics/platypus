#include "coefficients.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

extern const char * DATA_DIR;

// Floating point error tolerance
const double eps{1e-10};

double
scalar_f(const mfem::Vector & x)
{
  return x.Sum();
}

void
vector_f(const mfem::Vector & x, mfem::Vector & f)
{
  for (int i = 0; i < x.Size(); ++i)
  {
    f[i] = sin(x[i]);
  }
}

TEST_CASE("PWCoefficientTest", "[CheckSetup]")
{
  mfem::Vector one_vec(3);
  one_vec = 1.0;

  hephaestus::Subdomain air("air", 1);
  air._scalar_coefficients.Register("electrical_conductivity",
                                    std::make_shared<mfem::ConstantCoefficient>(1e-4));
  air._scalar_coefficients.Register("magnetic_permeability",
                                    std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));
  air._vector_coefficients.Register("vector_function",
                                    std::make_shared<mfem::VectorConstantCoefficient>(one_vec));

  hephaestus::Subdomain plate("plate", 2);
  plate._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::FunctionCoefficient>(scalar_f));
  plate._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));
  plate._vector_coefficients.Register(
      "vector_function", std::make_shared<mfem::VectorFunctionCoefficient>(3, vector_f));

  hephaestus::Subdomain coil1("coil1", 3);
  coil1._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1e4));
  coil1._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));
  coil1._vector_coefficients.Register(
      "vector_function", std::make_shared<mfem::VectorFunctionCoefficient>(3, vector_f));

  hephaestus::Subdomain coil2("coil2", 4);
  coil2._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1e4));
  coil2._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));
  coil2._vector_coefficients.Register(
      "vector_function", std::make_shared<mfem::VectorFunctionCoefficient>(3, vector_f));

  hephaestus::Subdomain coil3("coil3", 5);
  coil3._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1e4));
  coil3._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));
  coil3._vector_coefficients.Register(
      "vector_function", std::make_shared<mfem::VectorFunctionCoefficient>(3, vector_f));

  hephaestus::Subdomain coil4("coil4", 6);
  coil4._scalar_coefficients.Register("electrical_conductivity",
                                      std::make_shared<mfem::ConstantCoefficient>(1e4));
  coil4._scalar_coefficients.Register("magnetic_permeability",
                                      std::make_shared<mfem::ConstantCoefficient>(M_PI * 4.0e-7));
  coil4._vector_coefficients.Register(
      "vector_function", std::make_shared<mfem::VectorFunctionCoefficient>(3, vector_f));

  hephaestus::Coefficients coefficients(
      std::vector<hephaestus::Subdomain>({air, plate, coil1, coil2, coil3, coil4}));

  coefficients._scalars.Register("I", std::make_shared<mfem::ConstantCoefficient>(1.0));
  coefficients._vectors.Register("const_vector",
                                 std::make_shared<mfem::VectorConstantCoefficient>(one_vec));

  REQUIRE(coefficients._subdomains.size() == 6);

  REQUIRE(coefficients._scalars.Has("electrical_conductivity"));
  REQUIRE(coefficients._scalars.Has("magnetic_permeability"));
  REQUIRE(coefficients._scalars.Has("I"));

  REQUIRE(coefficients._vectors.Has("vector_function"));
  REQUIRE(coefficients._vectors.Has("const_vector"));
}

TEST_CASE("CoefficientsTest", "[CheckData]")
{

  hephaestus::Subdomain wire("wire", 1);
  wire._scalar_coefficients.Register("property_one",
                                     std::make_shared<mfem::ConstantCoefficient>(1.0));
  wire._scalar_coefficients.Register("property_two",
                                     std::make_shared<mfem::ConstantCoefficient>(150.0));

  hephaestus::Subdomain air("air", 2);
  air._scalar_coefficients.Register("property_one",
                                    std::make_shared<mfem::ConstantCoefficient>(26.0));
  air._scalar_coefficients.Register("property_two",
                                    std::make_shared<mfem::ConstantCoefficient>(152.0));

  hephaestus::Coefficients coefficients(std::vector<hephaestus::Subdomain>({wire, air}));

  // Verify predefined values
  mfem::IsoparametricTransformation t;
  mfem::IntegrationPoint ip;

  mfem::Coefficient * pw = coefficients._scalars.Get("property_one");
  t.Attribute = 1;
  REQUIRE_THAT(pw->Eval(t, ip), Catch::Matchers::WithinAbs(1.0, eps));
  t.Attribute = 2;
  REQUIRE_THAT(pw->Eval(t, ip), Catch::Matchers::WithinAbs(26.0, eps));

  pw = coefficients._scalars.Get("property_two");
  t.Attribute = 1;
  REQUIRE_THAT(pw->Eval(t, ip), Catch::Matchers::WithinAbs(150.0, eps));
  t.Attribute = 2;
  REQUIRE_THAT(pw->Eval(t, ip), Catch::Matchers::WithinAbs(152.0, eps));
}