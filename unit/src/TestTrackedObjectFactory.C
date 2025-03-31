#include <algorithm>
#include <memory>

#include "gtest/gtest.h"
#include "TrackedObjectFactory.h"

#include "mfem.hpp"

class CheckTrackedObjectFactory : public testing::Test
{
protected:
  platypus::TrackedObjectFactory<mfem::Coefficient> factory;
  std::shared_ptr<mfem::ConstantCoefficient> c1, c2, c3;

  CheckTrackedObjectFactory()
  {
    c1 = factory.make<mfem::ConstantCoefficient>(1.);
    c2 = factory.make<mfem::ConstantCoefficient>(2.);
    c3 = factory.make<mfem::ConstantCoefficient>(3.);
  }
};

TEST_F(CheckTrackedObjectFactory, Iter)
{
  int i = -1;
  for (auto coef : factory)
  {
    auto coef_const = std::dynamic_pointer_cast<mfem::ConstantCoefficient>(coef);
    ASSERT_NE(coef_const.get(), nullptr);
    coef_const->constant = i--;
  }
  EXPECT_EQ(c1->constant, -1);
  EXPECT_EQ(c2->constant, -2);
  EXPECT_EQ(c3->constant, -3);
}

TEST_F(CheckTrackedObjectFactory, IterMutable)
{
  int i = -1;
  for (auto & coef : factory)
  {
    EXPECT_EQ(coef.use_count(), 2);
    coef.reset();
  }
  EXPECT_EQ(c1.use_count(), 1);
  EXPECT_EQ(c2.use_count(), 1);
  EXPECT_EQ(c3.use_count(), 1);
}

TEST_F(CheckTrackedObjectFactory, ConstIter)
{
  int i = -1;
  for (auto it = factory.cbegin(); it != factory.cend(); ++it)
  {
    auto coef_const = std::dynamic_pointer_cast<mfem::ConstantCoefficient>(*it);
    ASSERT_NE(coef_const.get(), nullptr);
    coef_const->constant = i--;
  }
  EXPECT_EQ(c1->constant, -1);
  EXPECT_EQ(c2->constant, -2);
  EXPECT_EQ(c3->constant, -3);
}

TEST_F(CheckTrackedObjectFactory, ReverseIter)
{
  int i = -1;
  for (auto it = factory.rbegin(); it != factory.rend(); ++it)
  {
    auto coef_const = std::dynamic_pointer_cast<mfem::ConstantCoefficient>(*it);
    ASSERT_NE(coef_const.get(), nullptr);
    coef_const->constant = i--;
  }
  EXPECT_EQ(c1->constant, -3);
  EXPECT_EQ(c2->constant, -2);
  EXPECT_EQ(c3->constant, -1);
}

TEST_F(CheckTrackedObjectFactory, ReverseIterMutable)
{
  int i = -1;
  for (auto it = factory.rbegin(); it != factory.rend(); ++it)
  {
    EXPECT_EQ(it->use_count(), 2);
    it->reset();
  }
  EXPECT_EQ(c1.use_count(), 1);
  EXPECT_EQ(c2.use_count(), 1);
  EXPECT_EQ(c3.use_count(), 1);
}

TEST_F(CheckTrackedObjectFactory, ConstReverseIter)
{
  int i = -1;
  for (auto it = factory.crbegin(); it != factory.crend(); ++it)
  {
    auto coef_const = std::dynamic_pointer_cast<mfem::ConstantCoefficient>(*it);
    ASSERT_NE(coef_const.get(), nullptr);
    coef_const->constant = i--;
  }
  EXPECT_EQ(c1->constant, -3);
  EXPECT_EQ(c2->constant, -2);
  EXPECT_EQ(c3->constant, -1);
}
