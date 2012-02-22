/* BaseProblemTest.cpp */


#include <string>
#include <list>

#include "gtest/gtest.h"
#include "../Point.h"
#include "../PointAndSolution.h"
#include "VerySmallProblem.h"
#include "SmallGraphProblem.h"


using std::string;

using pareto_approximator::Point;
using pareto_approximator::PointAndSolution;
using pareto_approximator::BaseProblem;


namespace {


// The fixture for testing the BaseProblem wrapper class template.
// (and our implementation of the chord algorithm)
class BaseProblemTest : public ::testing::Test {
  protected:
    BaseProblemTest() { }
    
    ~BaseProblemTest() { }

    static const double bigEps       = 0.5;
    static const double smallEps     = 0.1;
    static const double verySmallEps = 0.001;
};


/*
// Test that doChord does what expected if we give it verySimpleComb() 
// and the west and south points.
TEST_F(BaseProblemTest, DoChordWithVerySimpleComb) 
{
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
*/


// Test that the approximateParetoSet() method does what expected for 
// the VerySmallProblem class (child of BaseProblem).
TEST_F(BaseProblemTest, ChordForVerySmallProblem) 
{
  VerySmallProblem vsp;
  list< PointAndSolution<string> > paretoSet;
  paretoSet = vsp.approximateParetoSet(smallEps);
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


TEST_F(BaseProblemTest, ChordForSmallGraphProblem)
{
  SmallGraphProblem sgp;
  list< PointAndSolution<PredecessorMap> > paretoSet;
  paretoSet = sgp.approximateParetoSet(verySmallEps);
  EXPECT_EQ(4, paretoSet.size());
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

