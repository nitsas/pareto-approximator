/*! \file HyperplaneTest.cpp
 *  \brief Unit test for the Hyperplane class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <vector>
#include <algorithm>

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
  // test that Hyperplane::getCoefficients() works
  EXPECT_EQ(h4.space_dimension(), 4);
  std::vector<double> coeffic = h4.getCoefficients();
  EXPECT_EQ(coefficientsB.size(), coeffic.size());
  EXPECT_EQ(std::equal(coefficientsB.begin(), 
                       coefficientsB.end(), coeffic.begin()), true);
  // test that it returns a copy of the actual vector of coefficients
  coeffic[0] = 8.0;
  coeffic[1] = 8.0;
  coeffic[2] = 8.0;
  coeffic[3] = 8.0;
  EXPECT_EQ(h4[0], -2.0);
  EXPECT_EQ(h4[1], 1.0);
  EXPECT_EQ(h4[2], 0.0);
  EXPECT_EQ(h4[3], 0.0);

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


TEST_F(HyperplaneTest, HyperplaneOperatorsAndStrWork) {
  Hyperplane h1(4.1, -2.2, 0.15, -2.1);
  EXPECT_EQ(h1.str(), "( 4.1 * x1 - 2.2 * x2 + 0.15 * x3 = -2.1 )");
  Hyperplane h2(-1.0, 0.0, 0.0);
  EXPECT_EQ(h2.str(), "( -1 * x1 + 0 * x2 = 0 )");

  EXPECT_EQ(h1 == h2, false);
  EXPECT_EQ(h1 != h2, true);

  Hyperplane h3(2 * 4.1, 2 * (-2.2), 2 * 0.15, 2 * (-2.1));
  EXPECT_EQ(h1 == h3, true);
  EXPECT_EQ(h1 != h3, false);
}


TEST_F(HyperplaneTest, HyperplaneParallelThroughAndIsParallelWork) {
  Hyperplane h1(4.0, -2.0, 0.0, -2.0);
  Hyperplane h2 = h1.parallelThrough(Point(1, 1, 1));
  EXPECT_EQ(h2[0], 4.0);
  EXPECT_EQ(h2[1], -2.0);
  EXPECT_EQ(h2[2], 0.0);
  EXPECT_EQ(h2.b(), 4.0 * 1 - 2.0 * 1 + 0.0 * 1);
  EXPECT_EQ(h1.isParallel(h2), true);

  Hyperplane h3 = h1.parallelThrough(Point(-1, -1, 1));
  EXPECT_EQ(h3, h1);
  EXPECT_EQ(h3.isParallel(h1), true);

  Hyperplane h4(1.0, 1.0, 1.0);
  EXPECT_EQ(h4.isParallel(h1), false);

  Hyperplane h5(2.0, 2.0, 5.0);
  EXPECT_EQ(h5.isParallel(h4), true);
}


TEST_F(HyperplaneTest, HyperplaneIntersectionWorksFor2Hyperplanes) {
  Hyperplane h1(1.0, -1.0, 0.0);
  Hyperplane h2(5.0, 2.0, 0.0);
  EXPECT_EQ(h1.intersection(h2), Point(0.0, 0.0));

  Hyperplane h3(-2.0, 1.0, -1.0);
  EXPECT_EQ(h3.intersection(h1), Point(1.0, 1.0));

  Hyperplane h4(0.0, 1.0, 3.3);
  EXPECT_EQ(h4.intersection(h1), Point(3.3, 3.3));
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

