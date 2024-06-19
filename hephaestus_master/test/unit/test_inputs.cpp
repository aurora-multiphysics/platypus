#include "inputs.hpp"
#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

TEST_CASE("InputParametersTest", "[CheckData]")
{
  hephaestus::InputParameters params;
  int example_int = 5;
  params.SetParam("IntegerParam", example_int);

  std::string example_string("ExampleString");
  params.SetParam("StringParam", example_string);

  mfem::Array<int> example_array({1, 2, 3});
  params.SetParam("ArrayParam", example_array);

  REQUIRE(params.GetParam<int>("IntegerParam") == example_int);

  REQUIRE(params.GetParam<std::string>("StringParam") == example_string);

  auto stored_array = params.GetParam<mfem::Array<int>>("ArrayParam");

  for (int i = 0; i < example_array.Size(); ++i)
    REQUIRE(example_array[i] == stored_array[i]);
}
