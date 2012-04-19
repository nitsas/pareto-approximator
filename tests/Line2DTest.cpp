/*! \file Line2DTest.cpp
 *  \brief Unit test for the Line2D class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <string>

#include "gtest/gtest.h"
#include "../Line2D.h"
#include "../Point.h"


using std::string;

using namespace pareto_approximator;


namespace {


// The fixture for testing class Line2D.
class Line2DTest : public ::testing::Test {
  protected:
    Line2DTest() { }
    
    ~Line2DTest() { }

};


// Test that Line2D's constructors (and attribute getters) work as expected.
TEST_F(Line2DTest, Line2DConstructorsAndAccessorsWork) {
  // Create line   y = 3.7 x - 1.15
  Line2D l1(3.7, -1.15);
  EXPECT_EQ(3.7, l1.m());
  EXPECT_EQ(-1.15, l1.b());

  // Create line   x = 5
  Line2D l2(5);
  // l2.m() should throw a VerticalLineException
  EXPECT_EQ(-5, l2.b());
  EXPECT_EQ(true, l2.isVertical());

  // Create the line that passes through points (1, 1) and (4, 4).
  Point p1(1, 1);
  Point p2(4, 4);
  Line2D l3(p1, p2);
  EXPECT_EQ(1, l3.m());
  EXPECT_EQ(0, l3.b());

  // Check possible exceptions.
  Point p3(2, 2, 2);
  EXPECT_THROW(Line2D l4(p1, p3), Not2DPointsException);
  Point p4(1, 1);
  EXPECT_THROW(Line2D l5(p1, p1), SamePointsException);
  EXPECT_THROW(l2.m(), VerticalLineException);
}


// Test that Line2D::isVertical() works correctly (and that we correctly 
// identify vertical lines).
TEST_F(Line2DTest, Line2DIsVerticalWorks) {
  Line2D l1(5);
  EXPECT_EQ(true, l1.isVertical());
  string willChange = "this will change";
  try {
    // this should throw a VerticalLineException
    l1.m();
  }
  catch (VerticalLineException) {
    // this should be executed
    willChange = "VerticalLineException thrown and caught!";
  }
  EXPECT_EQ("VerticalLineException thrown and caught!", willChange);

  Line2D l2(4.7);
  EXPECT_EQ(true, l2.isVertical());
  willChange = "this will change";
  try {
    // this should throw a VerticalLineException
    l2.m();
  }
  catch (VerticalLineException) {
    // this should be executed
    willChange = "VerticalLineException thrown and caught!";
  }

  Point p1(4.0, 3.0);
  Point p2(4.0, 2.0);
  Line2D l3(p1, p2);
  EXPECT_EQ(true, l3.isVertical());
  willChange = "this will change";
  try {
    // this should throw a VerticalLineException
    l3.m();
  }
  catch (VerticalLineException) {
    // this should be executed
    willChange = "VerticalLineException thrown and caught!";
  }

  // now a non-vertical line
  Line2D l4(2.0, 4.1);
  EXPECT_EQ(false, l4.isVertical());
  string willNotChange;
  willNotChange = "no VerticalLineException, as expected!";
  try {
    // this shouldn't throw a VerticalLineException
    l4.m();
  }
  catch (VerticalLineException) {
    // this shouldn't be executed
    willNotChange = "VerticalLineException was thrown when it shouldn't have!";
  }
  EXPECT_EQ("no VerticalLineException, as expected!", willNotChange);

  // and another non-vertical one
  Point p3(5.0, 5.0);
  Line2D l5(p1, p3);
  EXPECT_EQ(false, l5.isVertical());
  willNotChange = "no VerticalLineException, as expected!";
  try {
    // this shouldn't throw a VerticalLineException
    l5.m();
  }
  catch (VerticalLineException) {
    // this shouldn't be executed
    willNotChange = "VerticalLineException was thrown when it shouldn't have!";
  }
  EXPECT_EQ("no VerticalLineException, as expected!", willNotChange);
}


// Test that Line2D::intersection() works as expected.
TEST_F(Line2DTest, Line2DIntersectionWorks) {
  Point p1(4.0, 3.0);
  Point p2(2.0, 8.0);
  Point p3(3.0, 4.0);
  Line2D l1(p1, p2);
  Line2D l2(p2, p3);
  EXPECT_EQ(p2, l1.intersection(l2));

  Point p4(4.0, 2.0);
  Point p5(4.0, 0.0);
  Line2D l3(p1, p4);        // vertical line   x = 4.0
  EXPECT_EQ(p5, l3.intersection(l2));

  Line2D l4(-10);           // vertical line   x = -10
  EXPECT_THROW(l4.intersection(l4), ParallelLinesException);
  Line2D l5(1.0, 1.0);
  Line2D l6(1.0, 2.0);
  EXPECT_THROW(l5.intersection(l6), ParallelLinesException);
}


// Test that Line2D::ratioDistance() works as expected.
TEST_F(Line2DTest, Line2DRatioDistanceWorks) {
  Line2D l1(-1, 3);
  Point p1(1, 1);
  EXPECT_EQ(0.5, l1.ratioDistance(p1));
  Point p2(1.5, 1.5);
  EXPECT_EQ(0.0, l1.ratioDistance(p2));
  Point p3(2, 0);
  EXPECT_EQ(0.5, l1.ratioDistance(p3));
  Point p4(10, 10);
  EXPECT_EQ(0.0, l1.ratioDistance(p4));

  // Test against a vertical line   x = 5
  Line2D l2(5);
  EXPECT_EQ(1.5, l2.ratioDistance(p3));
  EXPECT_EQ(0.0, l2.ratioDistance(p4));
  Point p5(5, 7);
  EXPECT_EQ(0.0, l2.ratioDistance(p5));
}


// Test that Line2D::parallelThrough() works as expected.
TEST_F(Line2DTest, Line2DParallelThroughWorks) {
  Line2D l1(-1, 3);
  Point p1(0, 4);
  Point p2(0, 3);
  Line2D l2(-1, 4);
  EXPECT_EQ(l2, l1.parallelThrough(p1));
  EXPECT_EQ(l1, l1.parallelThrough(p2));

  // now on a vertical line
  Line2D l3(4);
  Line2D l4(0);
  EXPECT_EQ(l4, l3.parallelThrough(p1));
  Point p3(4, 5);
  EXPECT_EQ(l3, l3.parallelThrough(p3));
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

