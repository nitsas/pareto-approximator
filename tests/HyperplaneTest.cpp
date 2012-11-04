/*! \file HyperplaneTest.cpp
 *  \brief Unit test for the Hyperplane class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <set>
#include <vector>

#include "gtest/gtest.h"
#include "../Hyperplane.h"
#include "../Point.h"


using namespace pareto_approximator;


namespace {


// The fixture for testing class Hyperplane.
class HyperplaneTest : public ::testing::Test 
{
  protected:
    HyperplaneTest() { }

    ~HyperplaneTest() { }
};


// Test that Hyperplane's constructors and accessors work as expected.
TEST_F(HyperplaneTest, HyperplaneConstructorsAndAccessorsWork) 
{
  // test Hyperplane(double, double, double)
  Hyperplane h1(5, 10, 15);
  EXPECT_EQ(h1.a(0), 5.0);
  EXPECT_EQ(h1.a(1), 10.0);
  EXPECT_EQ(h1.b(), 15.0);
  EXPECT_EQ(h1.spaceDimension(), 2);

  // test Hyperplane(double, double, double, double)
  Hyperplane h2(0.0, 2.0, 5.0, -4.5);
  EXPECT_EQ(h2.a(0), 0.0);
  EXPECT_EQ(h2.a(1), 2.0);
  EXPECT_EQ(h2.a(2), 5.0);
  EXPECT_EQ(h2.b(), -4.5);
  EXPECT_EQ(h2.spaceDimension(), 3);

  // test Hyperplane(std::vector<int>::const_iterator, 
  //                 std::vector<int>::const_iterator)
  std::vector<int> coefficientsA(4);
  coefficientsA[0] = -2;
  coefficientsA[1] = 1;
  coefficientsA[2] = 0;
  coefficientsA[3] = 0;
  Hyperplane h3(coefficientsA.begin(), coefficientsA.end(), 12);
  EXPECT_EQ(h3.a(0), -2.0);
  EXPECT_EQ(h3.a(1), 1.0);
  EXPECT_EQ(h3.a(2), 0.0);
  EXPECT_EQ(h3.a(3), 0.0);
  EXPECT_EQ(h3.b(), 12.0);
  EXPECT_EQ(h3.spaceDimension(), 4);

  // test Hyperplane(std::vector<double>::const_iterator, 
  //                 std::vector<double>::const_iterator)
  std::vector<double> coefficientsB(4);
  coefficientsB[0] = -2.0;
  coefficientsB[1] = 1.0;
  coefficientsB[2] = 0.0;
  coefficientsB[3] = 0.0;
  Hyperplane h4(coefficientsB.begin(), coefficientsB.end(), 12.0);
  EXPECT_EQ(h4.a(0), -2.0);
  EXPECT_EQ(h4.a(1), 1.0);
  EXPECT_EQ(h4.a(2), 0.0);
  EXPECT_EQ(h4.a(3), 0.0);
  EXPECT_EQ(h4.b(), 12.0);
  EXPECT_EQ(h4.spaceDimension(), 4);

  // test Hyperplane(const int *, const int *, int)
  int coefficientsC[6] = {-1, 0, 1, 2, 3, 4};
  Hyperplane h5(coefficientsC, coefficientsC + 6, 5);
  EXPECT_EQ(h5.a(0), -1.0);
  EXPECT_EQ(h5.a(1), 0.0);
  EXPECT_EQ(h5.a(2), 1.0);
  EXPECT_EQ(h5.a(3), 2.0);
  EXPECT_EQ(h5.a(4), 3.0);
  EXPECT_EQ(h5.a(5), 4.0);
  EXPECT_EQ(h5.b(), 5.0);
  EXPECT_EQ(h5.spaceDimension(), 6);

  // test Hyperplane(const double *, const double *, double)
  double coefficientsD[6] = {-1.0, 0.0, 1.0, 2.0, 3.0, 4.0};
  Hyperplane h6(coefficientsD, coefficientsD + 6, 5.0);
  EXPECT_EQ(h6.a(0), -1.0);
  EXPECT_EQ(h6.a(1), 0.0);
  EXPECT_EQ(h6.a(2), 1.0);
  EXPECT_EQ(h6.a(3), 2.0);
  EXPECT_EQ(h6.a(4), 3.0);
  EXPECT_EQ(h6.a(5), 4.0);
  EXPECT_EQ(h6.b(), 5.0);
  EXPECT_EQ(h6.spaceDimension(), 6);

  // test Hyperplane(const Point &, const Point &)
  Hyperplane h7(Point(0, 1), Point(1, 2));
  EXPECT_EQ(h7.spaceDimension(), 2);
  EXPECT_EQ(h7.a(0), -h7.a(1));
  EXPECT_EQ(h7.a(1), h7.b());

  // test Hyperplane(std::vector<Point>::const_iterator, 
  //                 std::vector<Point>::const_iterator)
  std::vector<Point> pointsA;
  pointsA.push_back(Point(1, 0, 0));
  pointsA.push_back(Point(0, 1, 0));
  pointsA.push_back(Point(0, 0, 1));
  Hyperplane h8(pointsA.begin(), pointsA.end());
  EXPECT_EQ(h8.spaceDimension(), 3);
  EXPECT_EQ(h8.a(0), h8.a(1));
  EXPECT_EQ(h8.a(1), h8.a(2));
  EXPECT_EQ(h8.a(2), h8.b());

  // test Hyperplane(const Point *, const Point *)
  Point pointsB[4];
  pointsB[0] = Point(1.0, 0.0, 0.0, 0.0);
  pointsB[1] = Point(0.0, 1.0, 0.0, 0.0);
  pointsB[2] = Point(0.0, 0.0, 1.0, 0.0);
  pointsB[3] = Point(0.0, 0.0, 0.0, 1.0);
  Hyperplane h9(pointsB, pointsB + 4);
  EXPECT_EQ(h9.spaceDimension(), 4);
  EXPECT_EQ(h9.a(0), h9.a(1));
  EXPECT_EQ(h9.a(1), h9.a(2));
  EXPECT_EQ(h9.a(2), h9.a(3));
  EXPECT_EQ(h9.a(3), h9.b());

  // test Hyperplane constructors (special cases)
  // 1) Hyperplane parallel to the y-z (axes) plane:
  //    a_{1} * x_{1} = b
  //    We expect a_{1} to be equal to b and non-zero. 
  //    a_{2} and a_{3} should be zero.
  std::vector<Point> pointsC;
  pointsC.push_back(Point(1.0, 0.0, 0.0));
  pointsC.push_back(Point(1.0, 1.0, 0.0));
  pointsC.push_back(Point(1.0, 1.0, 1.0));
  Hyperplane h10(pointsC.begin(), pointsC.end());
  EXPECT_EQ(h10.spaceDimension(), 3);
  EXPECT_NE(h10.a(0), 0.0);
  EXPECT_EQ(h10.a(0), h10.b());
  EXPECT_EQ(h10.a(1), 0.0);
  EXPECT_EQ(h10.a(2), 0.0);

  // 2) 3 points on the same 3D-line. Line:
  //    x = 2t
  //    y = 3t
  //    z = 4t
  //    where t is a real number.
  //    An infinite number of different hyperplanes pass through points on 
  //    that line. We expect the all-zero 3D-hyperplane (0 = 0) to be returned.
  std::vector<Point> pointsD;
  pointsD.push_back(Point(0.0, 0.0, 0.0));
  pointsD.push_back(Point(2.0, 3.0, 4.0));
  pointsD.push_back(Point(4.0, 6.0, 8.0));
  Hyperplane h11(pointsD.begin(), pointsD.end());
  EXPECT_EQ(h11.spaceDimension(), 3);
  EXPECT_EQ(h11.a(0), 0.0);
  EXPECT_EQ(h11.a(1), 0.0);
  EXPECT_EQ(h11.a(2), 0.0);
  EXPECT_EQ(h11.b(), 0.0);

  // 3) Positive a_{i} and b coefficients even though the determinants used 
  //    in Hyperplane::init() would give negative ones.
  //    We want the hyperplane (line) that passes through (0, 1) and (2, 0).
  //    Both ( x_{1} + 2 * x_{2} = 2 ) and ( -x_{1} - 2 * x_{2} = - 2 ) are 
  //    valid answers (and represent the same hyperplane) but we expect the 
  //    first one.
  Hyperplane h12(Point(0.0, 1.0), Point(2.0, 0.0));
  EXPECT_EQ(h12.spaceDimension(), 2);
  EXPECT_GT(h12.a(0), 0.0);
  EXPECT_EQ(h12.a(1), 2 * h12.a(0));
  EXPECT_EQ(h12.b(), h12.a(1));
}


// Test that Hyperplane's operators and Hyperplane::str() work as expected.
TEST_F(HyperplaneTest, HyperplaneOperatorsAndStrWork) 
{
  Hyperplane h1(4.1, -2.2, 0.15, -2.1);
  EXPECT_EQ(h1.str(), "( 4.1 * x1 - 2.2 * x2 + 0.15 * x3 = -2.1 )");
  Hyperplane h2(-1.0, 0.0, 0.0);
  EXPECT_EQ(h2.str(), "( -1 * x1 + 0 * x2 = 0 )");

  EXPECT_FALSE(h1 == h2);
  EXPECT_TRUE(h1 != h2);

  Hyperplane h3(2 * 4.1, 2 * (-2.2), 2 * 0.15, 2 * (-2.1));
  EXPECT_TRUE(h1 == h3);
  EXPECT_FALSE(h1 != h3);
}


// Test that Hyperplane's iterators, Hyperplane::begin() and 
// Hyperplane::end() work as expected.
TEST_F(HyperplaneTest, HyperplaneIteratorsWork) 
{
  double coefficients[6] = {-1.0, 0.0, 1.0, 2.0, 3.0, 4.0};
  Hyperplane h1(coefficients, coefficients + 6, 5.0);
  EXPECT_TRUE(std::equal(h1.begin(), h1.end(), coefficients));

  Hyperplane::iterator h1_iter = h1.begin();
  *h1_iter = 10.0;
  *(h1_iter + 1) = 20.0;
  EXPECT_EQ(h1.a(0), 10.0);
  EXPECT_EQ(h1.a(1), 20.0);
}


// Test that Hyperplane::parallelThrough() and Hyperplane::isParallel() 
// work as expected.
TEST_F(HyperplaneTest, HyperplaneParallelThroughAndIsParallelWork) 
{
  Hyperplane h1(4.0, -2.0, 0.0, -2.0);
  Hyperplane h2 = h1.parallelThrough(Point(1, 1, 1));
  EXPECT_EQ(h2.a(0), 4.0);
  EXPECT_EQ(h2.a(1), -2.0);
  EXPECT_EQ(h2.a(2), 0.0);
  EXPECT_EQ(h2.b(), 4.0 * 1 - 2.0 * 1 + 0.0 * 1);
  EXPECT_TRUE(h1.isParallel(h2));

  Hyperplane h3 = h1.parallelThrough(Point(-1, -1, 1));
  EXPECT_EQ(h3, h1);
  EXPECT_TRUE(h3.isParallel(h1));

  Hyperplane h4(1.0, 1.0, 1.0);
  EXPECT_FALSE(h4.isParallel(h1));

  Hyperplane h5(2.0, 2.0, 5.0);
  EXPECT_TRUE(h5.isParallel(h4));
}


// Test that Hyperplane::hasAllAiCoefficientsNonPositive() works as expected.
TEST_F(HyperplaneTest, HyperplaneHasAllAiCoefficientsNonPositiveWorks)
{
  // Careful when testing Hyperplane::hasAllAiCoefficientsNonPositive(), 
  // because some hyperplane constructors call Hyperplane::init() and it will 
  // return hyperplanes with all positive a_{i} coefficients if possible. 
  // We can safely use constructors that directly set the coefficients.

  Hyperplane h1(1.0, -1.0, 1.0);
  EXPECT_FALSE(h1.hasAllAiCoefficientsNonPositive());

  Hyperplane h2(-1.0, -1.0, -1.0, 3.0);
  EXPECT_TRUE(h2.hasAllAiCoefficientsNonPositive());

  Hyperplane h3(-1.0, 0.0, 0.0, 2.0);
  EXPECT_TRUE(h3.hasAllAiCoefficientsNonPositive());
}


// Test that Hyperplane::hasAllAiCoefficientsNonNegative() works as expected.
TEST_F(HyperplaneTest, HyperplaneHasAllAiCoefficientsNonNegativeWorks)
{
  // We can safely use constructors that directly set the coefficients.

  Hyperplane h1(1.0, -1.0, 1.0);
  EXPECT_FALSE(h1.hasAllAiCoefficientsNonNegative());

  Hyperplane h2(1.0, 1.0, 1.0, -3.0);
  EXPECT_TRUE(h2.hasAllAiCoefficientsNonNegative());

  Hyperplane h3(1.0, 0.0, 0.0, -2.0);
  EXPECT_TRUE(h3.hasAllAiCoefficientsNonNegative());
}


// Test that Hyperplane::reverseCoefficientSigns() works as expected.
TEST_F(HyperplaneTest, HyperplaneReverseCoefficientSignsWorks) 
{
  Hyperplane h1(1.0, -1.0, 1.0);
  h1.reverseCoefficientSigns();
  EXPECT_EQ(h1.a(0), -1.0);
  EXPECT_EQ(h1.a(1), 1.0);
  EXPECT_EQ(h1.b(), -1.0);

  Hyperplane h2(1.0, 0.0, -4.0, -4.0);
  h2.reverseCoefficientSigns();
  EXPECT_EQ(h2.a(0), -1.0);
  EXPECT_EQ(h2.a(1), 0.0);
  EXPECT_EQ(h2.a(2), 4.0);
  EXPECT_EQ(h2.b(), 4.0);
}


// Test that Hyperplane::normalizeAiCoefficients() works as expected.
TEST_F(HyperplaneTest, HyperplaneNormalizeAiCoefficients) 
{
  std::vector<double> a(4, 4.0);
  double b = 32.0;
  Hyperplane h1(a.begin(), a.end(), b);

  h1.normalizeAiCoefficients();
  // All a_{i}'s should now be equal to 1/2 and b should be 4.
  ASSERT_EQ(h1.spaceDimension(), 4);
  EXPECT_EQ(h1.a(0), 0.5);
  EXPECT_EQ(h1.a(1), 0.5);
  EXPECT_EQ(h1.a(2), 0.5);
  EXPECT_EQ(h1.a(3), 0.5);
  EXPECT_EQ(h1.b(), 4.0);
}


// Test that Hyperplane::toVec() and Hyperplane::toRowVec() work.
TEST_F(HyperplaneTest, HyperplaneToArmadilloVectorMethodsWork) 
{
  Hyperplane h1(1.0, 2.0, 4.0);
  arma::vec h1v = h1.toVec();
  EXPECT_EQ(h1v.size(), 2);
  EXPECT_EQ(h1v(0), 1.0);
  EXPECT_EQ(h1v(1), 2.0);
  arma::rowvec h1rv = h1.toRowVec();
  EXPECT_EQ(h1rv.size(), 2);
  EXPECT_EQ(h1rv(0), 1.0);
  EXPECT_EQ(h1rv(1), 2.0);

  Hyperplane h2(-1.0, 0.0, 1.0, 4.0);
  arma::vec h2v = h2.toVec();
  EXPECT_EQ(h2v.size(), 3);
  EXPECT_EQ(h2v(0), -1.0);
  EXPECT_EQ(h2v(1), 0.0);
  EXPECT_EQ(h2v(2), 1.0);
  arma::rowvec h2rv = h2.toRowVec();
  EXPECT_EQ(h2rv.size(), 3);
  EXPECT_EQ(h2v(0), -1.0);
  EXPECT_EQ(h2v(1), 0.0);
  EXPECT_EQ(h2v(2), 1.0);
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

