/*! \file HyperplaneTest.cpp
 *  \brief Unit test for the Hyperplane class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <vector>

#include "gtest/gtest.h"
#include "../Hyperplane.h"
#include "../Point.h"


using namespace pareto_approximator;


namespace {


// The fixture for testing class Hyperplane.
class HyperplaneTest : public ::testing::Test {
  protected:
    HyperplaneTest() { }

    ~HyperplaneTest() { }
};


// Test that Hyperplane's constructors and accessors work as expected.
TEST_F(HyperplaneTest, HyperplaneConstructorsAndAccessorsWork) {
  Hyperplane h1(5, 10, 15);
  EXPECT_EQ(h1[0], 5.0);
  EXPECT_EQ(h1[1], 10.0);
  EXPECT_EQ(h1.b(), 15.0);
  EXPECT_EQ(h1.space_dimension(), 2);

  Hyperplane h2(0.0, 2.0, 5.0, -4.5);
  EXPECT_EQ(h2[0], 0.0);
  EXPECT_EQ(h2[1], 2.0);
  EXPECT_EQ(h2[2], 5.0);
  EXPECT_EQ(h2.b(), -4.5);
  EXPECT_EQ(h2.space_dimension(), 3);

  std::vector<int> coefficientsA(4);
  coefficientsA[0] = -2;
  coefficientsA[1] = 1;
  coefficientsA[2] = 0;
  coefficientsA[3] = 0;
  Hyperplane h3(coefficientsA.begin(), coefficientsA.end(), 12);
  EXPECT_EQ(h3[0], -2.0);
  EXPECT_EQ(h3[1], 1.0);
  EXPECT_EQ(h3[2], 0.0);
  EXPECT_EQ(h3[3], 0.0);
  EXPECT_EQ(h3.b(), 12.0);
  EXPECT_EQ(h3.space_dimension(), 4);

  std::vector<double> coefficientsB(4);
  coefficientsB[0] = -2.0;
  coefficientsB[1] = 1.0;
  coefficientsB[2] = 0.0;
  coefficientsB[3] = 0.0;
  Hyperplane h4(coefficientsB.begin(), coefficientsB.end(), 12.0);
  EXPECT_EQ(h4[0], -2.0);
  EXPECT_EQ(h4[1], 1.0);
  EXPECT_EQ(h4[2], 0.0);
  EXPECT_EQ(h4[3], 0.0);
  EXPECT_EQ(h4.b(), 12.0);
  EXPECT_EQ(h4.space_dimension(), 4);

  int coefficientsC[6] = {-1, 0, 1, 2, 3, 4};
  Hyperplane h5(coefficientsC, coefficientsC + 6, 5);
  EXPECT_EQ(h5[0], -1.0);
  EXPECT_EQ(h5[1], 0.0);
  EXPECT_EQ(h5[2], 1.0);
  EXPECT_EQ(h5[3], 2.0);
  EXPECT_EQ(h5[4], 3.0);
  EXPECT_EQ(h5[5], 4.0);
  EXPECT_EQ(h5.b(), 5.0);
  EXPECT_EQ(h5.space_dimension(), 6);

  double coefficientsD[6] = {-1.0, 0.0, 1.0, 2.0, 3.0, 4.0};
  Hyperplane h6(coefficientsD, coefficientsD + 6, 5.0);
  EXPECT_EQ(h6[0], -1.0);
  EXPECT_EQ(h6[1], 0.0);
  EXPECT_EQ(h6[2], 1.0);
  EXPECT_EQ(h6[3], 2.0);
  EXPECT_EQ(h6[4], 3.0);
  EXPECT_EQ(h6[5], 4.0);
  EXPECT_EQ(h6.b(), 5.0);
  EXPECT_EQ(h6.space_dimension(), 6);
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

