/*! \file BaseProblemTest.cpp
 *  \brief Unit test for the BaseProblem class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <string>
#include <list>

#include "gtest/gtest.h"
#include "../Point.h"
#include "../PointAndSolution.h"
#include "NonOptimalStartingPointsProblem.h"
#include "SmallGraphProblem.h"


using std::string;

using pareto_approximator::Point;
using pareto_approximator::PointAndSolution;


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


// Test that the computeConvexParetoSet() method does what's expected 
// for the SmallGraphProblem class (child of BaseProblem). 
// SmallGraphProblem is a multiobjective shortest path problem on a 
// small boost graph.
TEST_F(BaseProblemTest, ChordForSmallGraphProblem)
{
  SmallGraphProblem sgp;
  list< PointAndSolution<PredecessorMap> > paretoSet;
  paretoSet = sgp.computeConvexParetoSet(verySmallEps);
  EXPECT_EQ(4, paretoSet.size());
}


// Test computeConvexParetoSet().
// Test that non-optimal starting points are correctly deleted when 
// points that dominate them are found.
TEST_F(BaseProblemTest, ChordForNonOptimalStartingPointsProblem) 
{
  NonOptimalStartingPointsProblem nospp;
  list< PointAndSolution<string> > paretoSet;
  paretoSet = nospp.computeConvexParetoSet(verySmallEps);
  EXPECT_EQ(3, paretoSet.size());
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

