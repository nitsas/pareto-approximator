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
#include "SmallBiobjectiveSPProblem.h"
#include "SmallTripleobjectiveSPProblem.h"


using std::string;

using pareto_approximator::Point;
using pareto_approximator::PointAndSolution;


namespace {


// The fixture for testing the BaseProblem wrapper class template.
// (and our implementation of the chord algorithm)
class BaseProblemTest : public ::testing::Test 
{
  protected:
    BaseProblemTest() { }
    
    ~BaseProblemTest() { }

    static const double bigEps       = 0.5;
    static const double smallEps     = 0.1;
    static const double verySmallEps = 0.001;
    static const double zeroEps      = 0.0;
};


// Test computeConvexParetoSet().
// Test that non-optimal starting points are correctly deleted when 
// points that dominate them are found.
TEST_F(BaseProblemTest, ChordForNonOptimalStartingPointsProblem) 
{
  using non_optimal_starting_points_problem::NonOptimalStartingPointsProblem;

  NonOptimalStartingPointsProblem nospp;
  std::list< PointAndSolution<string> > paretoSet;
  unsigned int numObjectives = 3;
  paretoSet = nospp.computeConvexParetoSet(numObjectives, zeroEps);

  ASSERT_EQ(6, paretoSet.size());
  std::list< PointAndSolution<string> >::iterator psi = paretoSet.begin();
  EXPECT_EQ(Point(1.0, 4.0, 4.0), psi->point);
  EXPECT_EQ("best-on-x", psi->solution);
  ++psi;
  EXPECT_EQ(Point(2.0, 2.0, 3.0), psi->point);
  EXPECT_EQ("other", psi->solution);
  ++psi;
  EXPECT_EQ(Point(2.0, 3.0, 2.0), psi->point);
  EXPECT_EQ("other", psi->solution);
  ++psi;
  /*
  EXPECT_EQ(Point(2.5, 2.5, 2.5), psi->point);
  EXPECT_EQ("other", psi->solution);
  ++psi;
  */
  EXPECT_EQ(Point(3.0, 2.0, 2.0), psi->point);
  EXPECT_EQ("other", psi->solution);
  ++psi;
  EXPECT_EQ(Point(5.0, 1.0, 6.0), psi->point);
  EXPECT_EQ("best-on-y", psi->solution);
  ++psi;
  EXPECT_EQ(Point(5.0, 3.0, 1.0), psi->point);
  EXPECT_EQ("best-on-z", psi->solution);
  ++psi;
  EXPECT_TRUE(psi == paretoSet.end());
}


// Test that the computeConvexParetoSet() method does what's expected for 
// the SmallBiobjectiveSPProblem class (child of BaseProblem). 
// SmallBiobjectiveSPProblem is a biobjective shortest path problem on a 
// small boost graph.
TEST_F(BaseProblemTest, ChordForSmallBiobjectiveSPProblem)
{
  using small_biobjective_sp_problem::SmallBiobjectiveSPProblem;
  using small_biobjective_sp_problem::PredecessorMap;

  SmallBiobjectiveSPProblem sbspp;
  std::list< PointAndSolution<PredecessorMap> > paretoSet;
  unsigned int numObjectives = 2;
  paretoSet = sbspp.computeConvexParetoSet(numObjectives, 0.0);

  EXPECT_EQ(4, paretoSet.size());
  std::list< PointAndSolution<PredecessorMap> >::iterator psi = paretoSet.begin();
  EXPECT_EQ(Point(2, 16), psi->point);
  EXPECT_EQ(Point(3, 12), (++psi)->point);
  EXPECT_EQ(Point(4, 10), (++psi)->point);
  EXPECT_EQ(Point(14, 2), (++psi)->point);
}


// Test that the computeConvexParetoSet() method does what's expected for 
// the SmallTripleobjectiveSPProblem class (child of BaseProblem). 
// SmallTripleobjectiveSPProblem is a biobjective shortest path problem 
// on a small boost graph.
TEST_F(BaseProblemTest, ChordForSmallTripleobjectiveSPProblem)
{
  using small_tripleobjective_sp_problem::SmallTripleobjectiveSPProblem;
  using small_tripleobjective_sp_problem::PredecessorMap;

  SmallTripleobjectiveSPProblem stspp;
  std::list< PointAndSolution<PredecessorMap> > paretoSet;
  unsigned int numObjectives = 3;
  paretoSet = stspp.computeConvexParetoSet(numObjectives, 0.0);

  EXPECT_EQ(4, paretoSet.size());
  std::list< PointAndSolution<PredecessorMap> >::iterator psi = paretoSet.begin();
  EXPECT_EQ(Point(1, 9, 9), psi->point);
  EXPECT_EQ(Point(5, 5, 5), (++psi)->point);
  EXPECT_EQ(Point(9, 1, 9), (++psi)->point);
  EXPECT_EQ(Point(9, 9, 1), (++psi)->point);
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

