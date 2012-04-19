/*! \file PointTest.cpp
 *  \brief Unit test for the Point class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include "gtest/gtest.h"
#include "../Point.h"


using namespace pareto_approximator;


namespace {


// The fixture for testing class Point.
class PointTest : public ::testing::Test {
  protected:
    PointTest() { }
    
    ~PointTest() { }

};


// Test that Point's constructors work for ints.
TEST_F(PointTest, PointConstructorsWorkForInts) {
  Point p1(5);
  EXPECT_EQ(p1[0], 5);
  EXPECT_EQ(p1.dimension(), 1);

  Point p2(4, -1);
  EXPECT_EQ(p2[0], 4);
  EXPECT_EQ(p2[1], -1);
  EXPECT_EQ(p2.dimension(), 2);

  Point p3(-10, 3, 7);
  EXPECT_EQ(p3[0], -10);
  EXPECT_EQ(p3[1], 3);
  EXPECT_EQ(p3[2], 7);
  EXPECT_EQ(p3.dimension(), 3);
}


// Test that Point's constructors work for doubles.
TEST_F(PointTest, PointConstructorsWorkForDoubles) {
  Point p1(7.4);
  EXPECT_EQ(p1[0], 7.4);
  EXPECT_EQ(p1.dimension(), 1);

  Point p2(-4.73, 12.497);
  EXPECT_EQ(p2[0], -4.73);
  EXPECT_EQ(p2[1], 12.497);
  EXPECT_EQ(p2.dimension(), 2);

  Point p3(8.888802, -0.0001, 29.3);
  EXPECT_EQ(p3[0], 8.888802);
  EXPECT_EQ(p3[1], -0.0001);
  EXPECT_EQ(p3[2], 29.3);
  EXPECT_EQ(p3.dimension(), 3);
}


// Test that Point::operator=() works as expected.
TEST_F(PointTest, PointAssignmentOperatorWorks) {
  Point p1(4.0, 3.5, -2.7);
  Point p2(1.8, 2.1,  8.2);
  Point p3 = p1;

  EXPECT_EQ(p1[0], p3[0]);
  EXPECT_EQ(p1[1], p3[1]);
  EXPECT_EQ(p1[2], p3[2]);
  EXPECT_EQ(p1.dimension(), p3.dimension());
}

// Test that Point::operator==() works as expected.
TEST_F(PointTest, PointEqualityOperatorWorks) {
  Point p1(4.0, 3.5, -2.7);
  Point p2(1.8, 2.1,  8.2);
  Point p3(4.0, 3.5, -2.7);

  EXPECT_EQ(p1 == p3, true);
  EXPECT_EQ(p1 == p2, false);
}

// Test that Point::operator!=() works as expected.
TEST_F(PointTest, PointInequalityOperatorWorks) {
  Point p1(4.0, 3.5, -2.7);
  Point p2(1.8, 2.1,  8.2);
  Point p3(4.0, 3.5, -2.7);

  EXPECT_EQ(p1 != p3, false);
  EXPECT_EQ(p1 != p2, true);
}

// Test that Point::operator<() works as expected.
TEST_F(PointTest, PointLessThanOperatorWorks) {
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

  EXPECT_THROW(p1 < p5, DifferentDimensionsException);
}

// Test that Point::operator[]() works as expected.
TEST_F(PointTest, PointAccessOperatorWorks) {
  Point p1(4.0, 3.5, -2.7);

  EXPECT_EQ(p1[0], 4.0);
  EXPECT_EQ(p1[1], 3.5);
  EXPECT_EQ(p1[2], -2.7);

  EXPECT_THROW(p1[3], NonExistentCoordinateException);
}

// Test that Point::dimension() and Point::dimension(int) work as expected.
TEST_F(PointTest, PointGetDimensionAndSetDimensionWork) {
  Point p1(3.9);
  Point p2(0.3, 1.0);
  Point p3(4.0, 3.5, -2.7);
  EXPECT_EQ(p1.dimension(), 1);
  EXPECT_EQ(p2.dimension(), 2);
  EXPECT_EQ(p3.dimension(), 3);

  p3.dimension(2);
  EXPECT_EQ(p3.dimension(), 2);
  EXPECT_EQ(p3[0], 4.0);
  EXPECT_EQ(p3[1], 3.5);

  p3.dimension(1);
  EXPECT_EQ(p3.dimension(), 1);
  EXPECT_EQ(p3[0], 4.0);
}


// Test that Point::ratioDistance() works as expected.
TEST_F(PointTest, PointRatioDistanceWorks) {
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

  EXPECT_THROW(p1.ratioDistance(p5), DifferentDimensionsException);
}


// Test that Point::dominates() works as expected.
TEST_F(PointTest, PointDominatesWorks) {
  Point p1(1.0, 5.0);
  Point p2(1.5, 7.0);
  EXPECT_EQ(p1.dominates(p2), true);
  EXPECT_EQ(p2.dominates(p1), false);
  EXPECT_EQ(p2.dominates(p1, 0.5), true);
  Point p3(1.6, 6.0);
  EXPECT_EQ(p3.dominates(p1, 0.5), false);
  Point p4(1.0, 1.0, 1.0);
  Point p5(2.0, 2.0, 2.0);
  EXPECT_EQ(p4.dominates(p5), true);
  EXPECT_EQ(p5.dominates(p4), false);

  EXPECT_THROW(p1.dominates(p2, -0.5), NegativeApproximationRatioException);
  Point p6(-1.3, 8.7);
  EXPECT_THROW(p1.dominates(p6), NotPositivePointException);
  Point p7(2.4, 8.97, 1.42);
  EXPECT_THROW(p1.dominates(p7), DifferentDimensionsException);
}


// Test that Point::str() works as expected.
TEST_F(PointTest, PointStrWorks) {
  Point p1(1, 1000);
  EXPECT_EQ(p1.str(), "(1, 1000)");
  Point p2(49.75, 5000000.2);
  EXPECT_EQ(p2.str(), "(49.75, 5e+06)");
  Point p3(-4.9, 0.0);
  EXPECT_EQ(p3.str(), "(-4.9, 0)");
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

