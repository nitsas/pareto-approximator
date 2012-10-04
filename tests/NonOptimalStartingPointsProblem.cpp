/*! \file NonOptimalStartingPointsProblem.cpp
 *  \brief Implementation of the NonOptimalStartingPointsProblem class, 
 *         a simple problem class used in BaseProblemTest.cpp.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <assert.h>
#include <iterator>
#include <algorithm>
#include <iostream>

#include "../Point.h"
#include "NonOptimalStartingPointsProblem.h"


using pareto_approximator::Point;


namespace non_optimal_starting_points_problem {


NonOptimalStartingPointsProblem::NonOptimalStartingPointsProblem() 
{
  // Make some weakly Pareto optimal points (PointAndSolution instances).
  PointAndSolution<string> pas1(Point(1.0, 8.0, 9.0), "weakly-best-on-x");
  PointAndSolution<string> pas2(Point(7.0, 1.0, 8.0), "weakly-best-on-y");
  PointAndSolution<string> pas3(Point(9.0, 7.0, 1.0), "weakly-best-on-z");
  weaklyOptimalPoints_.push_back(pas1);
  weaklyOptimalPoints_.push_back(pas2);
  weaklyOptimalPoints_.push_back(pas3);

  // Make some Pareto optimal points (PointAndSolution instances). 
  PointAndSolution<string> pas4(Point(1.0, 4.0, 4.0), "best-on-x");
  PointAndSolution<string> pas5(Point(5.0, 1.0, 6.0), "best-on-y");
  PointAndSolution<string> pas6(Point(5.0, 3.0, 1.0), "best-on-z");
  
  PointAndSolution<string> pas7(Point(2.0, 2.0, 3.0), "other");
  PointAndSolution<string> pas8(Point(2.0, 3.0, 2.0), "other");
  PointAndSolution<string> pas9(Point(3.0, 2.0, 2.0), "other");
  
  // This one is Pareto optimal but inside the convex hull of the Pareto set.
  // Algorithms that find Pareto optimal points by minimizing weighted linear 
  // combinations of the objectives should not be able to find it.
  PointAndSolution<string> pas10(Point(2.5, 2.5, 2.5), "other");

  optimalPoints_.push_back(pas4);
  optimalPoints_.push_back(pas5);
  optimalPoints_.push_back(pas6);
  
  optimalPoints_.push_back(pas7);
  optimalPoints_.push_back(pas8);
  optimalPoints_.push_back(pas9);

  optimalPoints_.push_back(pas10);
}


NonOptimalStartingPointsProblem::~NonOptimalStartingPointsProblem() { }


PointAndSolution<string> 
NonOptimalStartingPointsProblem::comb(std::vector<double>::const_iterator first,
                                      std::vector<double>::const_iterator last)
{
  assert(std::distance(first, last) == 3);
  double xWeight = *first;
  double yWeight = *(first + 1);
  double zWeight = *(first + 2);

  if (yWeight == 0.0 && xWeight == 0.0 && zWeight == 0.0)
    return PointAndSolution<string>(Point(0.0, 0.0, 0.0), "error: all weights were 0");
  
  PointAndSolution<string> result;
  if (xWeight == 0.0 && yWeight == 0.0)
    // return a weakly optimal (starting) point with respect to z
    result = weaklyOptimalPoints_[2];
  else if (xWeight == 0.0 && zWeight == 0.0)
    // return a weakly optimal (starting) point with respect to y
    result = weaklyOptimalPoints_[1];
  else if (yWeight == 0.0 && zWeight == 0.0)
    // return a weakly optimal (starting) point with respect to x
    result = weaklyOptimalPoints_[0];
  else {
    // return the best point from the optimalPoints_ vector
    // use a CompareUsingWeights object to compare PointAndSolution instances
    CompareUsingWeightsFunctor comp(xWeight, yWeight, zWeight);
    std::vector< PointAndSolution<string> >::iterator min;
    min = std::min_element(optimalPoints_.begin(), optimalPoints_.end(), comp);
    result = *min;
  }

  return result;
}


CompareUsingWeightsFunctor::CompareUsingWeightsFunctor(double xWeight, double yWeight, double zWeight) : xWeight_(xWeight), yWeight_(yWeight), zWeight_(zWeight) { }


bool 
CompareUsingWeightsFunctor::operator() (PointAndSolution<string> a, PointAndSolution<string> b) 
{
  double valueA, valueB;
  valueA = xWeight_ * a.point[0] + yWeight_ * a.point[1] + zWeight_ * a.point[2];
  valueB = xWeight_ * b.point[0] + yWeight_ * b.point[1] + zWeight_ * b.point[2];

  return valueA < valueB;
}


}  // namespace non_optimal_starting_points_problem
