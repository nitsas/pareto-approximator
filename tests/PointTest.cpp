/* PointTest.cpp */


#include "gtest/gtest.h"
#include "../Point.h"


using namespace pareto_approximator;


namespace {


// The fixture for testing class Point.
class PointTest : public ::testing::Test {
  protected:
    PointTest() {
      
    }
    
    ~PointTest() {}

};


// Test that Point's constructors work for ints.
TEST_F(PointTest, PointConstructorsWorkForInts) {
  Point p1(5);
  EXPECT_EQ(5, p1.x);
  EXPECT_EQ(0, p1.y);
  EXPECT_EQ(0, p1.z);
  EXPECT_EQ(1, p1.dimension());

  Point p2(4, -1);
  EXPECT_EQ(4, p2.x);
  EXPECT_EQ(-1, p2.y);
  EXPECT_EQ(0, p2.z);
  EXPECT_EQ(2, p2.dimension());

  Point p3(-10, 3, 7);
  EXPECT_EQ(-10, p3.x);
  EXPECT_EQ(3, p3.y);
  EXPECT_EQ(7, p3.z);
  EXPECT_EQ(3, p3.dimension());
}


// Test that Point's constructors work for doubles.
TEST_F(PointTest, PointConstructorsWorkForDoubles) {
  Point p1(7.4);
  EXPECT_EQ(7.4, p1.x);
  EXPECT_EQ(0, p1.y);
  EXPECT_EQ(0, p1.z);
  EXPECT_EQ(1, p1.dimension());

  Point p2(-4.73, 12.497);
  EXPECT_EQ(-4.73, p2.x);
  EXPECT_EQ(12.497, p2.y);
  EXPECT_EQ(0, p2.z);
  EXPECT_EQ(2, p2.dimension());

  Point p3(8.888802, -0.0001, 29.3);
  EXPECT_EQ(8.888802, p3.x);
  EXPECT_EQ(-0.0001, p3.y);
  EXPECT_EQ(29.3, p3.z);
  EXPECT_EQ(3, p3.dimension());
}


// Test that Point's operators work as expected.
TEST_F(PointTest, PointOperatorsWorkCorrectly) {
  Point p1( 4.0, 3.5, -2.7);
  Point p2(-1.8, 2.1,  8.2);
  
  // test operator=
  Point p3 = p1;
  EXPECT_EQ(p1.x, p3.x);
  EXPECT_EQ(p1.y, p3.y);
  EXPECT_EQ(p1.z, p3.z);
  EXPECT_EQ(p1.dimension(), p3.dimension());

  // test operator==
  Point p4(4.0, 3.5, -2.7);
  EXPECT_EQ(true,  p1 == p4);
  EXPECT_EQ(false, p1 == p2);

  // test operator!=
  EXPECT_EQ(false, p1 != p4);
  EXPECT_EQ(true,  p1 != p2);
}


//Test that Point::dimension() and Point::dimension(int) work as expected.
TEST_F(PointTest, PointGetDimensionAndSetDimensionWork) {
  Point p1(3.9);
  Point p2(0.3, 1.0);
  Point p3(4.0, 3.5, -2.7);
  EXPECT_EQ(1, p1.dimension());
  EXPECT_EQ(2, p2.dimension());
  EXPECT_EQ(3, p3.dimension());

  p3.dimension(2);
  EXPECT_EQ(2, p3.dimension());
  EXPECT_EQ(0, p3.z);

  p3.dimension(1);
  EXPECT_EQ(1, p3.dimension());
  EXPECT_EQ(0, p3.y);
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

