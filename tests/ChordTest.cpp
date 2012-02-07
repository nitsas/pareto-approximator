/* ChordTest.cpp */


#include <string>
#include <list>
#include <tr1/functional>

#include "gtest/gtest.h"
#include "../Point.h"
#include "../Line2D.h"
#include "../PointAndSolution.h"
#include "../chord.h"


using std::string;

using namespace pareto_approximator;


namespace {


// A very simple comb routine we'll pass to chord as a callback.
// xWeight is objective x's weight and yWeight objective y's weight
PointAndSolution<string> simpleComb(double xWeight, double yWeight);


// The fixture for testing our chord implementation.
class ChordTest : public ::testing::Test {
  protected:
    ChordTest() 
    {
      smallEps = 0.1;
      bigEps = 0.5;
    }
    
    ~ChordTest() { }

    double smallEps;
    double bigEps;
};


// A very simple comb routine we'll pass to chord as a callback.
// xWeight is objective x's weight and yWeight objective y's weight
PointAndSolution<string> 
verySimpleComb(double xWeight, double yWeight)
{
  if (xWeight > 2 * yWeight) {
    return PointAndSolution<string>(Point(1.0, 4.0), "west");
  }
  else if (yWeight > 2 * xWeight) {
    return PointAndSolution<string>(Point(4.0, 1.0), "south");
  }
  else {
    return PointAndSolution<string>(Point(2.0, 2.0), "southwest");
  }
}


// Test that doChord does what expected if we pass it verySimpleComb() 
// and the west and south points.
TEST_F(ChordTest, DoChordTestWithVerySimpleComb) {
  list< PointAndSolution<string> > paretoSet;
  PointAndSolution<string> west(Point(1.0, 4.0), "west");
  PointAndSolution<string> south(Point(4.0, 1.0), "south");
  Point ideal(1.0, 1.0);
  paretoSet = doChord<string>(verySimpleComb, west, south, ideal, smallEps);
  ASSERT_EQ(3, paretoSet.size());
  paretoSet.sort();
  list< PointAndSolution<string> >::iterator it = paretoSet.begin();
  ++it;
  EXPECT_EQ(Point(2.0, 2.0), it->point);
  EXPECT_EQ("southwest", it->solution);
}


// Test that chordAlgorithm does what expected if we pass it verySimpleComb().
TEST_F(ChordTest, ChordTestWithVerySimpleComb) {
  list< PointAndSolution<string> > paretoSet;
  paretoSet = chordAlgorithm<string>(verySimpleComb, smallEps);
  ASSERT_EQ(3, paretoSet.size());
  paretoSet.sort();
  list< PointAndSolution<string> >::iterator it = paretoSet.begin();
  EXPECT_EQ(Point(1.0, 4.0), it->point);
  EXPECT_EQ("west", it->solution);
  ++it;
  EXPECT_EQ(Point(2.0, 2.0), it->point);
  EXPECT_EQ("southwest", it->solution);
  ++it;
  EXPECT_EQ(Point(4.0, 1.0), it->point);
  EXPECT_EQ("south", it->solution);
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

