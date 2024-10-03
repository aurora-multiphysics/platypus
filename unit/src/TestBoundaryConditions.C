#include "MFEMObjectUnitTest.h"
#include "MFEMVectorBoundaryIntegratedBC.h"
#include "boundary_conditions.h"
#include "mfem.hpp"
#include "gtest/gtest.h"
#include <memory>

TEST(CheckData, boundaryConditions)
{
  platypus::BCMap bc_map;
  mfem::Array<int> bdr_attrs({1, 2, 3});

  bc_map.Register(
      "tangential_dEdt",
      std::make_shared<platypus::BoundaryCondition>(std::string("boundary_1"), bdr_attrs));

  mfem::Array<int> ess_bdr = bc_map.Get("tangential_dEdt")->_bdr_attributes;

  for (int i = 0; i < bdr_attrs.Size(); ++i)
  {
    EXPECT_EQ(bdr_attrs[i], ess_bdr[i]);
  }
}

class MFEMBCTest : public MFEMObjectUnitTest
{
public:
  MFEMBCTest() : MFEMObjectUnitTest("PlatypusApp") {}
};

TEST_F(MFEMBCTest, MFEMVectorBoundaryIntegratedBC)
{
  // Build required kernel inputs
  InputParameters vec_coef_params = _factory.getValidParams("MFEMVectorConstantCoefficient");
  vec_coef_params.set<double>("value_x") = 1.0;
  vec_coef_params.set<double>("value_y") = 2.0;
  vec_coef_params.set<double>("value_z") = 3.0;
  _mfem_problem->addVectorCoefficient(
      "MFEMVectorConstantCoefficient", "vec_coef1", vec_coef_params);

  // Construct kernel
  InputParameters bc_params = _factory.getValidParams("MFEMVectorBoundaryIntegratedBC");
  bc_params.set<std::string>("variable") = "test_variable_name";
  bc_params.set<std::vector<BoundaryName>>("boundary") = {"1"};
  bc_params.set<std::string>("vector_coefficient") = "vec_coef1";
  auto & bc =
      addObject<MFEMVectorBoundaryIntegratedBC>("MFEMVectorBoundaryIntegratedBC", "bc1", bc_params);

  // Test MFEMVectorBoundaryIntegratedBC returns bc of the expected type
  std::shared_ptr<platypus::BoundaryCondition> generic_bc = bc.getBC();
  auto integrated_bc = std::dynamic_pointer_cast<platypus::IntegratedBC>(generic_bc);
  ASSERT_NE(integrated_bc, nullptr);

  auto integrator = dynamic_cast<mfem::VectorBoundaryLFIntegrator &>(*integrated_bc->_lfi_re);
}
