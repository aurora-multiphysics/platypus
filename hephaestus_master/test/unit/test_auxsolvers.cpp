#include "auxsolvers.hpp"
#include "gridfunctions.hpp"
#include <catch2/catch_test_macros.hpp>

extern const char * DATA_DIR;

// DummyAuxSolver that appends to its string every time it is
class DummyAuxSolver : public hephaestus::AuxSolver
{
public:
  DummyAuxSolver(std::string & _string_to_append, std::string & _string_modified)
    : _string_to_append(_string_to_append), _string_modified(_string_modified){};
  void Init(const hephaestus::GridFunctions & gridfunctions,
            hephaestus::Coefficients & coefficients) override{};

  void Solve(double t = 0.0) override { _string_modified = _string_modified + _string_to_append; };
  std::string & _string_to_append;
  std::string & _string_modified;
};

TEST_CASE("AuxSolverTest", "[CheckQueue]")
{
  std::string string_base("");
  std::string string_arg1("The order");
  std::string string_arg2(" should");
  std::string string_arg3(" be correct");
  std::string string_arg4("!");
  std::string string_arg5("!");

  hephaestus::AuxSolvers auxsolvers;

  auto auxsolver1 = std::make_shared<DummyAuxSolver>(string_arg1, string_base);
  auxsolver1->SetPriority(-1);
  auxsolvers.Register("OneAuxSolver", auxsolver1);

  auto auxsolver2 = std::make_shared<DummyAuxSolver>(string_arg2, string_base);
  auxsolvers.Register("AnotherAuxSolver", auxsolver2);

  auto auxsolver3 = std::make_shared<DummyAuxSolver>(string_arg3, string_base);
  auxsolver3->SetPriority(3);
  auxsolvers.Register("AThirdAuxSolver", auxsolver3);

  auto auxsolver4 = std::make_shared<DummyAuxSolver>(string_arg4, string_base);
  auxsolver4->SetPriority(4);
  auxsolvers.Register("YetAnotherAuxSolver", auxsolver4);

  auto auxsolver5 = std::make_shared<DummyAuxSolver>(string_arg5, string_base);
  auxsolver5->SetPriority(4);
  auxsolvers.Register("AFinalAuxSolver", auxsolver5);

  hephaestus::GridFunctions gridfunctions;
  hephaestus::Coefficients coefficients;
  auxsolvers.Init(gridfunctions, coefficients);
  auxsolvers.Solve();

  REQUIRE(string_base == "The order should be correct!!");
  auxsolvers.Solve();
  REQUIRE(string_base == "The order should be correct!!The order should be correct!!");
}
