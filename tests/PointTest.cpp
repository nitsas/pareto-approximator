/*! \file PointTest.cpp
 *  \brief Unit test for the Point class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <vector>
#include <armadillo>

#include "gtest/gtest.h"
#include "../Point.h"


using pareto_approximator::Point;
using pareto_approximator::exception_classes::NullObjectException;
using pareto_approximator::exception_classes::DifferentDimensionsException;
using pareto_approximator::exception_classes::NonExistentCoordinateException;
using pareto_approximator::exception_classes::NegativeApproximationRatioException;
using pareto_approximator::exception_classes::NotPositivePointException;
using pareto_approximator::exception_classes::NotStrictlyPositivePointException;


namespace {


// The fixture for testing class Point.
class PointTest : public ::testing::Test 
{
  protected:
    PointTest() 
    {
      nullPoint = Point();
    } 
    
    ~PointTest() { }

    Point nullPoint;
};


// Test that Point's empty constructor (Point::Point()) works. 
TEST_F(PointTest, PointsNullConstructorWorks)
{
  // test Point()
  // - Should create a null Point, i.e. a Point instance with 
  //   Point::isNull() (and Point::isNull_) true.
  // - We have called Point() to already (in PointTest's constructor) 
  //   to make nullPoint.
  // - All the methods of a null Point except Point::isNull() should throw 
  //   a NullObjectException when called.
  EXPECT_TRUE(nullPoint.isNull());
  EXPECT_THROW(nullPoint[0], NullObjectException);
  EXPECT_THROW(nullPoint.dimension(), NullObjectException);
}


// Test that Point's constructors work for ints.
TEST_F(PointTest, PointConstructorsWorkForInts) 
{
  // test Point(int)
  Point p1(5);
  EXPECT_FALSE(p1.isNull());
  EXPECT_EQ(p1[0], 5);
  EXPECT_EQ(p1.dimension(), 1);

  // test Point(int, int)
  Point p2(4, -1);
  EXPECT_FALSE(p2.isNull());
  EXPECT_EQ(p2[0], 4);
  EXPECT_EQ(p2[1], -1);
  EXPECT_EQ(p2.dimension(), 2);

  // test Point(int, int, int)
  Point p3(-10, 3, 7);
  EXPECT_FALSE(p3.isNull());
  EXPECT_EQ(p3[0], -10);
  EXPECT_EQ(p3[1], 3);
  EXPECT_EQ(p3[2], 7);
  EXPECT_EQ(p3.dimension(), 3);

  // test Point(int, int, int, int)
  Point p4(-1, 0, 1, 2);
  EXPECT_FALSE(p4.isNull());
  EXPECT_EQ(p4[0], -1);
  EXPECT_EQ(p4[1], 0);
  EXPECT_EQ(p4[2], 1);
  EXPECT_EQ(p4[3], 2);
  EXPECT_EQ(p4.dimension(), 4);

  // test Point(std::vector<int>::const_iterator, 
  //            std::vector<int>::const_iterator)
  std::vector<int> coordinatesA(5);
  coordinatesA[0] = 3;
  coordinatesA[1] = -2;
  coordinatesA[2] = 7;
  coordinatesA[3] = 0;
  coordinatesA[4] = -8;
  Point p5(coordinatesA.begin(), coordinatesA.end());
  EXPECT_FALSE(p5.isNull());
  EXPECT_EQ(p5[0], 3);
  EXPECT_EQ(p5[1], -2);
  EXPECT_EQ(p5[2], 7);
  EXPECT_EQ(p5[3], 0);
  EXPECT_EQ(p5[4], -8);
  EXPECT_EQ(p5.dimension(), 5);

  // test Point(const int *, const int *)
  int coordinatesB[6] = {-2, -1, 0, 1, 2, 3};
  Point p6(coordinatesB, coordinatesB + 6);
  EXPECT_FALSE(p6.isNull());
  EXPECT_EQ(p6[0], -2);
  EXPECT_EQ(p6[1], -1);
  EXPECT_EQ(p6[2], 0);
  EXPECT_EQ(p6[3], 1);
  EXPECT_EQ(p6[4], 2);
  EXPECT_EQ(p6[5], 3);
}


// Test that Point's constructors work for doubles.
TEST_F(PointTest, PointConstructorsWorkForDoubles) 
{
  // test Point(double)
  Point p1(7.4);
  EXPECT_FALSE(p1.isNull());
  EXPECT_EQ(p1[0], 7.4);
  EXPECT_EQ(p1.dimension(), 1);

  // test Point(double, double)
  Point p2(-4.73, 12.497);
  EXPECT_FALSE(p2.isNull());
  EXPECT_EQ(p2[0], -4.73);
  EXPECT_EQ(p2[1], 12.497);
  EXPECT_EQ(p2.dimension(), 2);

  // test Point(double, double, double)
  Point p3(8.888802, -0.0001, 29.3);
  EXPECT_FALSE(p3.isNull());
  EXPECT_EQ(p3[0], 8.888802);
  EXPECT_EQ(p3[1], -0.0001);
  EXPECT_EQ(p3[2], 29.3);
  EXPECT_EQ(p3.dimension(), 3);

  // test Point(double, double, double, double)
  Point p4(-1.1, 0.0, 1.1, 2.2);
  EXPECT_FALSE(p4.isNull());
  EXPECT_EQ(p4[0], -1.1);
  EXPECT_EQ(p4[1], 0.0);
  EXPECT_EQ(p4[2], 1.1);
  EXPECT_EQ(p4[3], 2.2);
  EXPECT_EQ(p4.dimension(), 4);

  // test Point(std::vector<double>::const_iterator, 
  //            std::vector<double>::const_iterator)
  std::vector<double> coordinatesA(5);
  coordinatesA[0] = -2.1;
  coordinatesA[1] = 5.49;
  coordinatesA[2] = 0.0;
  coordinatesA[3] = -3.0;
  coordinatesA[4] = 8.5;
  Point p5(coordinatesA.begin(), coordinatesA.end());
  EXPECT_FALSE(p5.isNull());
  EXPECT_EQ(p5[0], -2.1);
  EXPECT_EQ(p5[1], 5.49);
  EXPECT_EQ(p5[2], 0.0);
  EXPECT_EQ(p5[3], -3.0);
  EXPECT_EQ(p5[4], 8.5);
  EXPECT_EQ(p5.dimension(), 5);

  // test Point(const double *, const double *)
  double coordinatesB[6] = {-2.5, -1.5, -0.5, 0.5, 1.5, 2.5};
  Point p6(coordinatesB, coordinatesB + 6);
  EXPECT_FALSE(p6.isNull());
  EXPECT_EQ(p6[0], -2.5);
  EXPECT_EQ(p6[1], -1.5);
  EXPECT_EQ(p6[2], -0.5);
  EXPECT_EQ(p6[3], 0.5);
  EXPECT_EQ(p6[4], 1.5);
  EXPECT_EQ(p6[5], 2.5);
  EXPECT_EQ(p6.dimension(), 6);

  // test Point(arma::vec::const_iterator, arma::vec::const_iterator)
  arma::vec coordinatesC;
  coordinatesC << -1.0 << arma::endr << 0.0 << arma::endr << 1.0 << arma::endr;
  Point p7(coordinatesC.begin(), coordinatesC.end());
  EXPECT_FALSE(p7.isNull());
  EXPECT_EQ(p7[0], -1.0);
  EXPECT_EQ(p7[1], 0.0);
  EXPECT_EQ(p7[2], 1.0);
  EXPECT_EQ(p7.dimension(), 3);
}


// Test that Point::operator=() works as expected.
TEST_F(PointTest, PointAssignmentOperatorWorks) 
{
  Point p1(4.0, 3.5, -2.7);
  Point p2 = p1;

  EXPECT_EQ(p1.isNull(), p2.isNull());
  EXPECT_EQ(p1[0], p2[0]);
  EXPECT_EQ(p1[1], p2[1]);
  EXPECT_EQ(p1[2], p2[2]);
  EXPECT_EQ(p1.dimension(), p2.dimension());

  double coordinates[5] = {10.0, 9.0, 8.0, 7.0, 6.0};
  Point p3(coordinates, coordinates + 5);
  Point p4 = p3;
  EXPECT_EQ(p3.isNull(), p4.isNull());
  EXPECT_EQ(p3[0], p4[0]);
  EXPECT_EQ(p3[1], p4[1]);
  EXPECT_EQ(p3[2], p4[2]);
  EXPECT_EQ(p3[3], p4[3]);
  EXPECT_EQ(p3[4], p4[4]);
  EXPECT_EQ(p3.dimension(), p4.dimension());

  Point p5 = nullPoint;
  EXPECT_EQ(p5.isNull(), nullPoint.isNull());
}

// Test that Point::operator==() works as expected.
TEST_F(PointTest, PointEqualityOperatorWorks) 
{
  Point p1(4.0, 3.5, -2.7);
  Point p2(1.8, 2.1,  8.2);
  Point p3(4.0, 3.5, -2.7);

  EXPECT_TRUE(p1 == p3);
  EXPECT_FALSE(p1 == p2);

  double coordinatesA[5] = {10.0, 9.0, 8.0, 7.0, 6.0};
  double coordinatesB[5] = {-1.0, 0.0, 1.0, 2.0, 3.0};
  Point p4(coordinatesA, coordinatesA + 5);
  Point p5(coordinatesA, coordinatesA + 5);
  Point p6(coordinatesB, coordinatesB + 5);

  EXPECT_TRUE(p4 == p5);
  EXPECT_FALSE(p4 == p6);

  EXPECT_FALSE(p1 == nullPoint);
  EXPECT_FALSE(p2 == nullPoint);
  EXPECT_FALSE(p3 == nullPoint);
  EXPECT_FALSE(p4 == nullPoint);
  EXPECT_FALSE(p5 == nullPoint);
  EXPECT_FALSE(p6 == nullPoint);
  EXPECT_TRUE(nullPoint == Point());
}

// Test that Point::operator!=() works as expected.
TEST_F(PointTest, PointInequalityOperatorWorks) 
{
  Point p1(4.0, 3.5, -2.7);
  Point p2(1.8, 2.1,  8.2);
  Point p3(4.0, 3.5, -2.7);

  EXPECT_EQ(p1 != p3, false);
  EXPECT_EQ(p1 != p2, true);

  double coordinates[5] = {10.0, 9.0, 8.0, 7.0, 6.0};
  Point p4(coordinates, coordinates + 5);
  Point p5(coordinates, coordinates + 5);

  EXPECT_EQ(p4 != p5, false);
  EXPECT_EQ(p4 != p3, true);

  EXPECT_TRUE(p1 != nullPoint);
  EXPECT_TRUE(p2 != nullPoint);
  EXPECT_TRUE(p3 != nullPoint);
  EXPECT_TRUE(p4 != nullPoint);
  EXPECT_TRUE(p5 != nullPoint);
  EXPECT_FALSE(nullPoint != Point());
}

// Test that Point::operator<() works as expected.
TEST_F(PointTest, PointLessThanOperatorWorks) 
{
  Point p1(4.0, 3.5, -2.7);
  Point p2(1.8, 2.1,  8.2);
  Point p3(4.0, 3.5, -2.8);
  Point p4(1.8, 2.0, 17.5);
  Point p5(17.1, 15.4);
  Point p6(17.1, 13.1);
  
  EXPECT_EQ(p2 < p1, true);
  EXPECT_EQ(p1 < p2, false);
  EXPECT_EQ(p1 < p1, false);
  EXPECT_EQ(p3 < p1, true);
  EXPECT_EQ(p4 < p2, true);
  EXPECT_EQ(p6 < p5, true);
  EXPECT_EQ(p5 < p6, false);

  double coordinatesA[5] = {10.0, 9.0, 8.0, 7.0, 6.0};
  double coordinatesB[5] = {-1.0, 0.0, 1.0, 2.0, 3.0};
  Point p7(coordinatesA, coordinatesA + 5);
  Point p8(coordinatesB, coordinatesB + 5);
  
  EXPECT_EQ(p7 < p8, false);
  EXPECT_EQ(p8 < p7, true);

  EXPECT_THROW(p1 < p5, DifferentDimensionsException);

  EXPECT_THROW(p1 < nullPoint, NullObjectException);
  EXPECT_THROW(nullPoint < p1, NullObjectException);
}

// Test that Point::operator[]() works as expected.
TEST_F(PointTest, PointAccessOperatorWorks) 
{
  Point p1(4.0, 3.5, -2.7);

  EXPECT_EQ(p1[0], 4.0);
  EXPECT_EQ(p1[1], 3.5);
  EXPECT_EQ(p1[2], -2.7);

  double coordinates[5] = {10.0, 9.0, 8.0, 7.0, 6.0};
  Point p2(coordinates, coordinates + 5);
  EXPECT_EQ(p2[0], 10.0);
  EXPECT_EQ(p2[1], 9.0);
  EXPECT_EQ(p2[2], 8.0);
  EXPECT_EQ(p2[3], 7.0);
  EXPECT_EQ(p2[4], 6.0);

  EXPECT_THROW(p1[3], NonExistentCoordinateException);

  EXPECT_THROW(nullPoint[0], NullObjectException);
}

// Test that Point::operator+() works as expected.
TEST_F(PointTest, PointMinusOperatorWorks) 
{
  Point p1(4.0, 3.5, -4.0);
  Point p2(2.0, 2.0, 2.0);
  Point p3 = p1 + p2;
  EXPECT_EQ(p3[0], 6.0);
  EXPECT_EQ(p3[1], 5.5);
  EXPECT_EQ(p3[2], -2.0);

  Point p4(1.0, 1.0);
  EXPECT_THROW(p1 + p4, DifferentDimensionsException);

  EXPECT_THROW(p1 + nullPoint, NullObjectException);
  EXPECT_THROW(nullPoint + p1, NullObjectException);
}

// Test that Point::operator-() works as expected.
TEST_F(PointTest, PointPlusOperatorWorks) 
{
  Point p1(4.0, 3.5, -4.0);
  Point p2(2.0, 2.0, 2.0);
  Point p3 = p1 - p2;
  EXPECT_EQ(p3[0], 2.0);
  EXPECT_EQ(p3[1], 1.5);
  EXPECT_EQ(p3[2], -6.0);

  Point p4(1.0, 1.0);
  EXPECT_THROW(p1 - p4, DifferentDimensionsException);

  EXPECT_THROW(p1 - nullPoint, NullObjectException);
  EXPECT_THROW(nullPoint - p1, NullObjectException);
}

// Test that Point::toVec() and Point::toRowVec() work.
TEST_F(PointTest, PointToArmadilloVectorMethodsWork) 
{
  Point p1(1.4);
  arma::vec p1v = p1.toVec();
  EXPECT_EQ(p1v.size(), 1);
  EXPECT_EQ(p1v(0), 1.4);
  arma::rowvec p1rv = p1.toRowVec();
  EXPECT_EQ(p1rv.size(), 1);
  EXPECT_EQ(p1rv(0), 1.4);

  Point p2(-1.0, 0.0, 1.0, 2.0);
  arma::vec p2v = p2.toVec();
  EXPECT_EQ(p2v.size(), 4);
  EXPECT_EQ(p2v(0), -1.0);
  EXPECT_EQ(p2v(1), 0.0);
  EXPECT_EQ(p2v(2), 1.0);
  EXPECT_EQ(p2v(3), 2.0);
  arma::rowvec p2rv = p2.toRowVec();
  EXPECT_EQ(p2rv.size(), 4);
  EXPECT_EQ(p2v(0), -1.0);
  EXPECT_EQ(p2v(1), 0.0);
  EXPECT_EQ(p2v(2), 1.0);
  EXPECT_EQ(p2v(3), 2.0);

  EXPECT_THROW(nullPoint.toVec(), NullObjectException);
  EXPECT_THROW(nullPoint.toRowVec(), NullObjectException);
}

// Test that Point::isZero() works.
TEST_F(PointTest, PointIsZeroWorks) 
{
  Point p1(0.0);
  EXPECT_TRUE(p1.isZero());

  Point p2(0.0, 0.0);
  EXPECT_TRUE(p2.isZero());
  
  Point p3(0.0, 0.0, 0.0);
  EXPECT_TRUE(p3.isZero());

  Point p4(1.0, 0.0, 2.4);
  EXPECT_FALSE(p4.isZero());

  double coordinates[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  Point p5(coordinates, coordinates + 5);
  EXPECT_TRUE(p5.isZero());

  EXPECT_THROW(nullPoint.isZero(), NullObjectException);
}
  
// Test that Point::dimension() works as expected.
TEST_F(PointTest, PointGetDimensionWorks) 
{
  Point p1(3.9);
  Point p2(0.3, 1.0);
  double coordinates[4] = {10.0, 9.0, 8.0, 7.0};
  Point p3(coordinates, coordinates + 4);
  EXPECT_EQ(p1.dimension(), 1);
  EXPECT_EQ(p2.dimension(), 2);
  EXPECT_EQ(p3.dimension(), 4);

  EXPECT_THROW(nullPoint.dimension(), NullObjectException);
}


// Test that Point::euclideanDistance() works as expected.
TEST_F(PointTest, PointEuclideanDistanceWorks)
{
  Point p1(0, 0), p2(3, 4);
  EXPECT_EQ(p1.euclideanDistance(p2), 5.0);
  EXPECT_EQ(p2.euclideanDistance(p1), p1.euclideanDistance(p2));

  Point p3(-3, -4);
  EXPECT_EQ(p1.euclideanDistance(p3), 5.0);
  EXPECT_EQ(p3.euclideanDistance(p1), 5.0);

  Point p4(0, 0, 0);
  EXPECT_THROW(p4.euclideanDistance(p1), DifferentDimensionsException);
  EXPECT_THROW(p1.euclideanDistance(p4), DifferentDimensionsException);

  Point nullPoint;
  EXPECT_THROW(nullPoint.euclideanDistance(p1), NullObjectException);
  EXPECT_THROW(p1.euclideanDistance(nullPoint), NullObjectException);
}


// Test that Point::ratioDistance() works as expected.
TEST_F(PointTest, PointRatioDistanceWorks) 
{
  Point p1(2, 100);
  Point p2(4, 900);
  EXPECT_EQ(p1.ratioDistance(p2), 8);
  Point p3(4, 110);
  EXPECT_EQ(p1.ratioDistance(p3), 1);
  Point p4(1, 100);
  EXPECT_EQ(p1.ratioDistance(p4), 0);
  EXPECT_EQ(p1.ratioDistance(p1), 0);
  Point p5(1, 10, 100);
  Point p6(2, 30, 400);
  EXPECT_EQ(p5.ratioDistance(p6), 3);

  double coordinatesA[5] = {1.0, 10.0, 100.0, 1000.0, 10000.0};
  double coordinatesB[5] = {1.0, 20.0, 300.0, 4000.0, 50000.0};
  Point p7(coordinatesA, coordinatesA + 5);
  Point p8(coordinatesB, coordinatesB + 5);
  EXPECT_EQ(p7.ratioDistance(p8), 4);

  EXPECT_THROW(p1.ratioDistance(p5), DifferentDimensionsException);
  EXPECT_THROW(p1.ratioDistance(nullPoint), NullObjectException);
  EXPECT_THROW(nullPoint.ratioDistance(p1), NullObjectException);
}


// Test that Point::additiveDistance() works as expected.
TEST_F(PointTest, PointAdditiveDistanceWorks)
{
  Point p1(1, 1, 1), p2(2, 3, 4);
  EXPECT_EQ(p1.additiveDistance(p2), 3.0);
  EXPECT_EQ(p2.additiveDistance(p1), 0.0);

  Point p3(1, 5, 3);
  EXPECT_EQ(p3.additiveDistance(p2), 1.0);
  EXPECT_EQ(p2.additiveDistance(p3), 2.0);

  Point nullPoint = Point();
  EXPECT_THROW(nullPoint.additiveDistance(p1), NullObjectException);
  EXPECT_THROW(p1.additiveDistance(nullPoint), NullObjectException);

  Point pd(1, 1);
  EXPECT_THROW(pd.additiveDistance(p1), DifferentDimensionsException);
  EXPECT_THROW(p1.additiveDistance(pd), DifferentDimensionsException);
}


// Test that Point::dominatesAdditive() works as expected.
TEST_F(PointTest, PointDominatesAdditiveWorks)
{
  Point p1(1.0, 5.0, 2.0), p2(1.5, 7.0, 3.0);
  EXPECT_TRUE(p1.dominatesAdditive(p2));
  EXPECT_FALSE(p2.dominatesAdditive(p1));
  EXPECT_TRUE(p2.dominatesAdditive(p1, 2.0));
  EXPECT_TRUE(p1.dominatesAdditive(p1));

  Point nullPoint = Point();
  EXPECT_THROW(p1.dominatesAdditive(nullPoint), NullObjectException);
  EXPECT_THROW(nullPoint.dominatesAdditive(p1), NullObjectException);

  Point p3(1, 1);
  EXPECT_THROW(p1.dominatesAdditive(p3), DifferentDimensionsException);

  EXPECT_THROW(p1.dominatesAdditive(p2, -2.0), 
               NegativeApproximationRatioException);

  Point p4(-1, 2, 2);
  EXPECT_THROW(p4.dominatesAdditive(p1), NotPositivePointException);
  Point p5(0, 0, 0);
  EXPECT_NO_THROW(p5.dominatesAdditive(p1));
}


// Test that Point::dominatesMultiplicative() works as expected.
TEST_F(PointTest, PointDominatesMultiplicativeWorks) 
{
  Point p1(1.0, 5.0);
  Point p2(1.5, 7.0);
  EXPECT_TRUE(p1.dominatesMultiplicative(p2));
  EXPECT_FALSE(p2.dominatesMultiplicative(p1));
  EXPECT_TRUE(p2.dominatesMultiplicative(p1, 0.5));
  Point p3(1.6, 6.0);
  EXPECT_FALSE(p3.dominatesMultiplicative(p1, 0.5));
  Point p4(1.0, 1.0, 1.0);
  Point p5(2.0, 2.0, 2.0);
  EXPECT_TRUE(p4.dominatesMultiplicative(p5));
  EXPECT_FALSE(p5.dominatesMultiplicative(p4));

  double coordinatesA[5] = {1.0, 10.0, 100.0, 1000.0, 10000.0};
  double coordinatesB[5] = {1.0, 20.0, 300.0, 4000.0, 50000.0};
  Point p7(coordinatesA, coordinatesA + 5);
  Point p8(coordinatesB, coordinatesB + 5);
  EXPECT_TRUE(p7.dominatesMultiplicative(p8));
  EXPECT_FALSE(p8.dominatesMultiplicative(p7));
  EXPECT_TRUE(p8.dominatesMultiplicative(p7, 4));

  EXPECT_THROW(p1.dominatesMultiplicative(p2, -0.5), 
               NegativeApproximationRatioException);

  Point p9(-1.3, 8.7);
  EXPECT_THROW(p1.dominatesMultiplicative(p9), 
               NotStrictlyPositivePointException);
  Point p10(0.0, 0.0);
  EXPECT_THROW(p10.dominatesMultiplicative(p1), 
               NotStrictlyPositivePointException);
  EXPECT_THROW(p1.dominatesMultiplicative(p10), 
               NotStrictlyPositivePointException);

  Point p11(2.4, 8.97, 1.42);
  EXPECT_THROW(p1.dominatesMultiplicative(p11), DifferentDimensionsException);

  EXPECT_THROW(p1.dominatesMultiplicative(nullPoint), NullObjectException);
  EXPECT_THROW(nullPoint.dominatesMultiplicative(p1), NullObjectException);
}


// Test that Point::str() works as expected.
TEST_F(PointTest, PointStrWorks) 
{
  Point p1(1, 1000);
  EXPECT_EQ(p1.str(), "(1, 1000)");
  EXPECT_EQ(p1.str(true), "1.0000000000000e+00 1.0000000000000e+03");
  Point p2(49.75, 5000000.2);
  EXPECT_EQ(p2.str(), "(49.75, 5e+06)");
  EXPECT_EQ(p2.str(true), "4.9750000000000e+01 5.0000002000000e+06");
  Point p3(-4.9, 0.0);
  EXPECT_EQ(p3.str(), "(-4.9, 0)");
  EXPECT_EQ(p3.str(true), "-4.9000000000000e+00 0.0000000000000e+00");
  double coordinates[4] = {2.2, 4.2, 8.2, 16.2};
  Point p4(coordinates, coordinates + 4);
  EXPECT_EQ(p4.str(), "(2.2, 4.2, 8.2, 16.2)");
  EXPECT_EQ(p4.str(true), "2.2000000000000e+00 4.2000000000000e+00 8.2000000000000e+00 1.6200000000000e+01");

  EXPECT_EQ(nullPoint.str(), "()");
  EXPECT_EQ(nullPoint.str(true), "");
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

