#include <algorithm>

#include "gtest/gtest.h"
#include "coefficient_map.h"

#include "mfem.hpp"

class CheckCoefficientMap : public testing::Test
{
protected:
  mfem::IsoparametricTransformation fe_transform;
  mfem::IntegrationPoint point1, point2;
  CheckCoefficientMap()
  {
    point1.Init(2);
    point1.Set2(0., 0.);
    point2.Init(2);
    point2.Set2(0.5, 0.1);
    fe_transform.SetIdentityTransformation(mfem::Geometry::SQUARE);
  }
};

TEST_F(CheckCoefficientMap, HasCoefficient)
{
  platypus::ScalarMap coeff_map;
  coeff_map.addCoefficient("test", coeff_map.make<mfem::ConstantCoefficient>(1.));
  EXPECT_TRUE(coeff_map.hasCoefficient("test"));
  EXPECT_FALSE(coeff_map.hasCoefficient("missing"));
}

TEST_F(CheckCoefficientMap, GetCoefficient)
{
  platypus::ScalarMap coeff_map;
  coeff_map.addCoefficient("resistivity", coeff_map.make<mfem::ConstantCoefficient>(2.));
  mfem::ConstantCoefficient * c =
      dynamic_cast<mfem::ConstantCoefficient *>(&coeff_map.getCoefficient("resistivity"));
  EXPECT_NE(c, nullptr);
  EXPECT_EQ(c->Eval(fe_transform, point1), 2.0);
  EXPECT_EQ(c->Eval(fe_transform, point2), 2.0);
}

TEST_F(CheckCoefficientMap, GetPWCoefficient)
{
  platypus::ScalarMap coeff_map;
  coeff_map.addPiecewiseBlocks("test", coeff_map.make<mfem::ConstantCoefficient>(2.), {"1", "2"});
  coeff_map.addPiecewiseBlocks("test", coeff_map.make<mfem::ConstantCoefficient>(1.), {"3"});
  mfem::PWCoefficient * c = dynamic_cast<mfem::PWCoefficient *>(&coeff_map.getCoefficient("test"));
  EXPECT_NE(c, nullptr);
  fe_transform.Attribute = 1;
  EXPECT_EQ(c->Eval(fe_transform, point1), 2.0);
  EXPECT_EQ(c->Eval(fe_transform, point2), 2.0);
  fe_transform.Attribute = 2;
  EXPECT_EQ(c->Eval(fe_transform, point1), 2.0);
  fe_transform.Attribute = 3;
  EXPECT_EQ(c->Eval(fe_transform, point1), 1.0);
  fe_transform.Attribute = 10;
  EXPECT_EQ(c->Eval(fe_transform, point1), 0.0);
}

TEST_F(CheckCoefficientMap, HasVecCoefficient)
{
  platypus::VectorMap coeff_map;
  coeff_map.addCoefficient(
      "test", coeff_map.make<mfem::VectorConstantCoefficient>(mfem::Vector({0., 1., 0.})));
  EXPECT_TRUE(coeff_map.hasCoefficient("test"));
  EXPECT_FALSE(coeff_map.hasCoefficient("missing"));
}

TEST_F(CheckCoefficientMap, GetVecCoefficient)
{
  platypus::VectorMap coeff_map;
  coeff_map.addCoefficient(
      "resistivity", coeff_map.make<mfem::VectorConstantCoefficient>(mfem::Vector({0., 1., 0.})));
  mfem::VectorConstantCoefficient * c =
      dynamic_cast<mfem::VectorConstantCoefficient *>(&coeff_map.getCoefficient("resistivity"));
  EXPECT_NE(c, nullptr);
  mfem::Vector vec;
  c->Eval(vec, fe_transform, point1);
  EXPECT_EQ(vec[0], 0.);
  EXPECT_EQ(vec[1], 1.);
  EXPECT_EQ(vec[2], 0.);
}

TEST_F(CheckCoefficientMap, GetPWVecCoefficient)
{
  platypus::VectorMap coeff_map;
  coeff_map.addPiecewiseBlocks(
      "test",
      coeff_map.make<mfem::VectorConstantCoefficient>(mfem::Vector({0., 1., 0.})),
      {"1", "2"});
  coeff_map.addPiecewiseBlocks(
      "test", coeff_map.make<mfem::VectorConstantCoefficient>(mfem::Vector({1., 0., 0.})), {"3"});
  mfem::PWVectorCoefficient * c =
      dynamic_cast<mfem::PWVectorCoefficient *>(&coeff_map.getCoefficient("test"));
  EXPECT_NE(c, nullptr);
  mfem::Vector vec;
  fe_transform.Attribute = 1;
  c->Eval(vec, fe_transform, point1);
  EXPECT_EQ(vec[0], 0.);
  EXPECT_EQ(vec[1], 1.);
  EXPECT_EQ(vec[2], 0.);
  fe_transform.Attribute = 2;
  c->Eval(vec, fe_transform, point1);
  EXPECT_EQ(vec[0], 0.);
  EXPECT_EQ(vec[1], 1.);
  EXPECT_EQ(vec[2], 0.);
  fe_transform.Attribute = 3;
  c->Eval(vec, fe_transform, point1);
  EXPECT_EQ(vec[0], 1.);
  EXPECT_EQ(vec[1], 0.);
  EXPECT_EQ(vec[2], 0.);
  fe_transform.Attribute = 10;
  c->Eval(vec, fe_transform, point1);
  EXPECT_EQ(vec[0], 0.);
  EXPECT_EQ(vec[1], 0.);
  EXPECT_EQ(vec[2], 0.);
}

TEST_F(CheckCoefficientMap, HasMatCoefficient)
{
  platypus::MatrixMap coeff_map;
  coeff_map.addCoefficient(
      "test",
      coeff_map.make<mfem::MatrixConstantCoefficient>(mfem::DenseMatrix({{0., 1.}, {1., 0.}})));
  EXPECT_TRUE(coeff_map.hasCoefficient("test"));
  EXPECT_FALSE(coeff_map.hasCoefficient("missing"));
}

TEST_F(CheckCoefficientMap, GetMatCoefficient)
{
  platypus::MatrixMap coeff_map;
  coeff_map.addCoefficient(
      "resistivity",
      coeff_map.make<mfem::MatrixConstantCoefficient>(mfem::DenseMatrix({{0., 1.}, {1., 0.}})));
  mfem::MatrixConstantCoefficient * c =
      dynamic_cast<mfem::MatrixConstantCoefficient *>(&coeff_map.getCoefficient("resistivity"));
  EXPECT_NE(c, nullptr);
  mfem::DenseMatrix mat;
  c->Eval(mat, fe_transform, point1);
  EXPECT_EQ(mat(0, 0), 0.);
  EXPECT_EQ(mat(0, 1), 1.);
  EXPECT_EQ(mat(1, 0), 1.);
  EXPECT_EQ(mat(1, 1), 0.);
}

TEST_F(CheckCoefficientMap, GetPWMatCoefficient)
{
  platypus::MatrixMap coeff_map;
  coeff_map.addPiecewiseBlocks(
      "test",
      coeff_map.make<mfem::MatrixConstantCoefficient>(mfem::DenseMatrix({{0., 1.}, {1., 0.}})),
      {"1", "2"});
  coeff_map.addPiecewiseBlocks(
      "test",
      coeff_map.make<mfem::MatrixConstantCoefficient>(mfem::DenseMatrix({{11., 0.}, {0., -1.}})),
      {"3"});
  mfem::PWMatrixCoefficient * c =
      dynamic_cast<mfem::PWMatrixCoefficient *>(&coeff_map.getCoefficient("test"));
  EXPECT_NE(c, nullptr);
  mfem::DenseMatrix mat;
  fe_transform.Attribute = 1;
  c->Eval(mat, fe_transform, point1);
  EXPECT_EQ(mat(0, 0), 0.);
  EXPECT_EQ(mat(0, 1), 1.);
  EXPECT_EQ(mat(1, 0), 1.);
  EXPECT_EQ(mat(1, 1), 0.);
  fe_transform.Attribute = 2;
  c->Eval(mat, fe_transform, point1);
  EXPECT_EQ(mat(0, 0), 0.);
  EXPECT_EQ(mat(0, 1), 1.);
  EXPECT_EQ(mat(1, 0), 1.);
  EXPECT_EQ(mat(1, 1), 0.);
  fe_transform.Attribute = 3;
  c->Eval(mat, fe_transform, point1);
  EXPECT_EQ(mat(0, 0), 11.);
  EXPECT_EQ(mat(0, 1), 0.);
  EXPECT_EQ(mat(1, 0), 0.);
  EXPECT_EQ(mat(1, 1), -1.);
  fe_transform.Attribute = 10;
  c->Eval(mat, fe_transform, point1);
  EXPECT_EQ(mat(0, 0), 0.);
  EXPECT_EQ(mat(0, 1), 0.);
  EXPECT_EQ(mat(1, 0), 0.);
  EXPECT_EQ(mat(1, 1), 0.);
}

TEST_F(CheckCoefficientMap, CoefficientDefinedOnBlock)
{
  platypus::ScalarMap coeff_map;
  coeff_map.addPiecewiseBlocks("a", coeff_map.make<mfem::ConstantCoefficient>(2.), {"1", "2"});
  coeff_map.addCoefficient("b", coeff_map.make<mfem::ConstantCoefficient>(5.));
  EXPECT_TRUE(coeff_map.propertyDefinedOnBlock("a", "1"));
  EXPECT_TRUE(coeff_map.propertyDefinedOnBlock("a", "2"));
  EXPECT_FALSE(coeff_map.propertyDefinedOnBlock("a", "3"));
  EXPECT_FALSE(coeff_map.propertyDefinedOnBlock("a", "0"));
  EXPECT_TRUE(coeff_map.propertyDefinedOnBlock("b", "0"));
  EXPECT_TRUE(coeff_map.propertyDefinedOnBlock("b", "1"));
  EXPECT_TRUE(coeff_map.propertyDefinedOnBlock("b", "2"));
  EXPECT_FALSE(coeff_map.propertyDefinedOnBlock("c", "1"));
  EXPECT_FALSE(coeff_map.propertyDefinedOnBlock("A", "1"));
  EXPECT_FALSE(coeff_map.propertyDefinedOnBlock("B", "1"));
}

TEST_F(CheckCoefficientMap, OverwriteCoefficient)
{
  platypus::ScalarMap coeff_map;
  coeff_map.addCoefficient("resistivity", coeff_map.make<mfem::ConstantCoefficient>(2.));
  EXPECT_THROW(
      coeff_map.addCoefficient("resistivity", coeff_map.make<mfem::ConstantCoefficient>(4.)),
      MooseException);
}

TEST_F(CheckCoefficientMap, CoefficientNotDeclared)
{
  platypus::ScalarMap coeff_map;
  EXPECT_THROW(coeff_map.getCoefficient("NotDeclared"), MooseException);
}

TEST_F(CheckCoefficientMap, AddBlocksForGlobalCoefficient)
{
  platypus::ScalarMap coeff_map;
  coeff_map.addCoefficient("resistivity", coeff_map.make<mfem::ConstantCoefficient>(2.));
  EXPECT_THROW(coeff_map.addPiecewiseBlocks(
                   "resistivity", coeff_map.make<mfem::ConstantCoefficient>(1.), {"1", "2"}),
               MooseException);
}

TEST_F(CheckCoefficientMap, OverwriteBlocks)
{
  platypus::ScalarMap coeff_map;
  coeff_map.addPiecewiseBlocks(
      "resistivity", coeff_map.make<mfem::ConstantCoefficient>(2.), {"2", "3"});
  EXPECT_THROW(coeff_map.addPiecewiseBlocks(
                   "resistivity", coeff_map.make<mfem::ConstantCoefficient>(1.), {"1", "2"}),
               MooseException);
}

TEST_F(CheckCoefficientMap, DifferentVecSize)
{
  platypus::VectorMap coeff_map;
  coeff_map.addPiecewiseBlocks(
      "test",
      coeff_map.make<mfem::VectorConstantCoefficient>(mfem::Vector({0., 1., 0.})),
      {"1", "2"});
  EXPECT_THROW(
      coeff_map.addPiecewiseBlocks(
          "test", coeff_map.make<mfem::VectorConstantCoefficient>(mfem::Vector({1., 0.})), {"3"}),
      MooseException);
}

TEST_F(CheckCoefficientMap, DifferentMatSize)
{
  platypus::MatrixMap coeff_map;
  coeff_map.addPiecewiseBlocks(
      "test",
      coeff_map.make<mfem::MatrixConstantCoefficient>(mfem::DenseMatrix({{0., 1.}, {1., 0.}})),
      {"1", "2"});
  EXPECT_THROW(
      coeff_map.addPiecewiseBlocks("test",
                                   coeff_map.make<mfem::MatrixConstantCoefficient>(
                                       mfem::DenseMatrix({{11., 0., 42.}, {0., -1., 10.}})),
                                   {"3"}),
      MooseException);
}
