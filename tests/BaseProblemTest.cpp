/*! \file BaseProblemTest.cpp
 *  \brief Unit test for the BaseProblem class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <string>
#include <list>
#include <algorithm>

#include "gtest/gtest.h"
#include "../Point.h"
#include "../PointAndSolution.h"
#include "NonOptimalStartingPointsProblem.h"
#include "SmallBiobjectiveSPProblem.h"
#include "SmallTripleobjectiveSPProblem.h"
#include "TripleobjectiveWithNegativeWeightsProblem.h"


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

    static const double smallEpsilon = 0.001;
    static const double verySmallEpsilon = 1e-12;
};


// Test that the computeConvexParetoSet() method finds the correct 
// (approximate) convex Pareto set for the SmallBiobjectiveSPProblem problem 
// class (child of BaseProblem). 
// SmallBiobjectiveSPProblem is a biobjective shortest path problem on a 
// small boost graph.
TEST_F(BaseProblemTest, SmallBiobjectiveSPProblem)
{
  using small_biobjective_sp_problem::SmallBiobjectiveSPProblem;
  using small_biobjective_sp_problem::PredecessorMap;

  SmallBiobjectiveSPProblem sbspp;
  std::vector< PointAndSolution<PredecessorMap> > paretoSet;
  unsigned int numObjectives = 2;
  paretoSet = sbspp.computeConvexParetoSet(numObjectives, verySmallEpsilon);
  std::sort(paretoSet.begin(), paretoSet.end());

  EXPECT_EQ(4, paretoSet.size());
  std::vector< PointAndSolution<PredecessorMap> >::iterator psi = paretoSet.begin();
  EXPECT_EQ(Point(2, 16), psi->point);
  EXPECT_EQ(Point(3, 12), (++psi)->point);
  EXPECT_EQ(Point(4, 10), (++psi)->point);
  EXPECT_EQ(Point(14, 2), (++psi)->point);
}


// Test computeConvexParetoSet().
// Test that non-optimal starting points are correctly deleted when 
// points that dominate them are found.
TEST_F(BaseProblemTest, BiobjectiveNonOptimalStartingPointsProblem) 
{
  using non_optimal_starting_points_problem::NonOptimalStartingPointsProblem;

  unsigned int numObjectives = 2;
  NonOptimalStartingPointsProblem nospp(numObjectives);
  std::vector< PointAndSolution<string> > paretoSet;
  paretoSet = nospp.computeConvexParetoSet(numObjectives, verySmallEpsilon);
  std::sort(paretoSet.begin(), paretoSet.end());

  ASSERT_EQ(4, paretoSet.size());
  std::vector< PointAndSolution<string> >::iterator psi = paretoSet.begin();
  EXPECT_EQ(Point(1.0, 5.0), psi->point);
  EXPECT_EQ("best-on-x", psi->solution);
  ++psi;
  EXPECT_EQ(Point(2.0, 3.0), psi->point);
  EXPECT_EQ("other", psi->solution);
  ++psi;
  EXPECT_EQ(Point(3.0, 2.0), psi->point);
  EXPECT_EQ("other", psi->solution);
  ++psi;
  EXPECT_EQ(Point(5.0, 1.0), psi->point);
  EXPECT_EQ("best-on-y", psi->solution);
  ++psi;
  EXPECT_TRUE(psi == paretoSet.end());
}


// Test computeConvexParetoSet().
// Test that non-optimal starting points are correctly deleted when 
// points that dominate them are found.
TEST_F(BaseProblemTest, TripleObjectiveNonOptimalStartingPointsProblem) 
{
  using non_optimal_starting_points_problem::NonOptimalStartingPointsProblem;

  unsigned int numObjectives = 3;
  NonOptimalStartingPointsProblem nospp(numObjectives);
  std::vector< PointAndSolution<string> > paretoSet;
  paretoSet = nospp.computeConvexParetoSet(numObjectives, verySmallEpsilon);
  std::sort(paretoSet.begin(), paretoSet.end());

  ASSERT_EQ(6, paretoSet.size());
  std::vector< PointAndSolution<string> >::iterator psi = paretoSet.begin();
  EXPECT_EQ(Point(1.0, 4.0, 4.0), psi->point);
  EXPECT_EQ("best-on-x", psi->solution);
  ++psi;
  EXPECT_EQ(Point(2.0, 2.0, 3.0), psi->point);
  EXPECT_EQ("other", psi->solution);
  ++psi;
  EXPECT_EQ(Point(2.0, 3.0, 2.0), psi->point);
  EXPECT_EQ("other", psi->solution);
  ++psi;
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


// Test that the computeConvexParetoSet() method finds the correct 
// (approximate) convex Pareto set for the SmallTripleobjectiveSPProblem 
// problem class (child of BaseProblem). 
// SmallTripleobjectiveSPProblem is a triple-objective shortest path problem 
// on a small boost graph.
TEST_F(BaseProblemTest, SmallTripleobjectiveSPProblem)
{
  using small_tripleobjective_sp_problem::SmallTripleobjectiveSPProblem;
  using small_tripleobjective_sp_problem::PredecessorMap;

  SmallTripleobjectiveSPProblem stspp;
  std::vector< PointAndSolution<PredecessorMap> > paretoSet;
  unsigned int numObjectives = 3;
  paretoSet = stspp.computeConvexParetoSet(numObjectives, verySmallEpsilon);
  std::sort(paretoSet.begin(), paretoSet.end());

  EXPECT_EQ(4, paretoSet.size());
  std::vector< PointAndSolution<PredecessorMap> >::iterator psi = paretoSet.begin();
  EXPECT_EQ(Point(1, 9, 9), psi->point);
  EXPECT_EQ(Point(5, 5, 5), (++psi)->point);
  EXPECT_EQ(Point(9, 1, 9), (++psi)->point);
  EXPECT_EQ(Point(9, 9, 1), (++psi)->point);
}


// Test that computeConvexParetoSet() can detect and avoid using negative 
// weights.
// We do that by testing if the computeConvexParetoSet() doesn't return 
// the specially placed dominated point in the 
// TripleobjectiveWithNegativeWeightsProblem problem (child class of 
// BaseProblem). 
// TripleobjectiveWithNegativeWeightsProblem is a very simple triple-objective 
// problem (child of BaseProblem) that leads to negative weights (and they 
// to a specially placed dominated point) if we use the simple version of the 
// Chord algorithm. 
TEST_F(BaseProblemTest, TripleobjectiveWithNegativeWeightsProblem)
{
  using tripleobjective_with_negative_weights_problem::TripleobjectiveWithNegativeWeightsProblem;

  TripleobjectiveWithNegativeWeightsProblem twnwp;
  std::vector< PointAndSolution<string> > paretoSet;
  unsigned int numObjectives = 3;

  ASSERT_NO_THROW(paretoSet = twnwp.computeConvexParetoSet(numObjectives, verySmallEpsilon));
  std::sort(paretoSet.begin(), paretoSet.end());

  ASSERT_EQ(4, paretoSet.size());
  std::vector< PointAndSolution<string> >::iterator psi = paretoSet.begin();
  EXPECT_EQ(Point(36.0, 349.0, 280.0), psi->point);
  EXPECT_EQ("best-on-x", psi->solution);
  ++psi;
  EXPECT_EQ(Point(92.0, 142.0, 42.0), psi->point);
  EXPECT_EQ("best-other", psi->solution);
  ++psi;
  EXPECT_EQ(Point(221.0, 54.0, 261.0), psi->point);
  EXPECT_EQ("best-on-y", psi->solution);
  ++psi;
  EXPECT_EQ(Point(342.0, 351.0, 37.0), psi->point);
  EXPECT_EQ("best-on-z", psi->solution);
  ++psi;
  EXPECT_TRUE(psi == paretoSet.end());
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

