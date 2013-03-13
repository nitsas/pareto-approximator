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


NonOptimalStartingPointsProblem::NonOptimalStartingPointsProblem(unsigned int dimension) :
                                                      dimension_(dimension)
{
  // We will fill in weaklyOptimalPoints_ and optimalPoints_.
  // Different contents for the biobjective and tripleobjective problem.
  if (dimension_ == 2) {
    // Make some weakly Pareto optimal points (PointAndSolution instances).
    PointAndSolution<string> pas1(Point(1.0, 6.0), "weakly-best-on-x");
    PointAndSolution<string> pas2(Point(6.0, 1.0), "weakly-best-on-y");

    weaklyOptimalPoints_.push_back(pas1);
    weaklyOptimalPoints_.push_back(pas2);

    // Make some Pareto optimal points (PointAndSolution instances).
    PointAndSolution<string> pas3(Point(1.0, 5.0), "best-on-x");
    PointAndSolution<string> pas4(Point(2.0, 3.0), "other");
    PointAndSolution<string> pas5(Point(3.0, 2.0), "other");
    PointAndSolution<string> pas6(Point(5.0, 1.0), "best-on-y");
      // The next one is Pareto optimal but inside the convex hull of the 
      // Pareto set. Algorithms that find Pareto optimal points by 
      // minimizing weighted linear combinations of the objectives should 
      // not be able to find it.
    PointAndSolution<string> pas7(Point(2.8, 2.8), "should-not-be-found");

    optimalPoints_.push_back(pas3);
    optimalPoints_.push_back(pas4);
    optimalPoints_.push_back(pas5);
    optimalPoints_.push_back(pas6);
    optimalPoints_.push_back(pas7);
  }
  else {
    assert(dimension_ == 3);

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
    
      // The next one is Pareto optimal but inside the convex hull of the 
      // Pareto set. Algorithms that find Pareto optimal points by 
      // minimizing weighted linear combinations of the objectives should 
      // not be able to find it.
    PointAndSolution<string> pas10(Point(2.5, 2.5, 2.5), "other");

    optimalPoints_.push_back(pas4);
    optimalPoints_.push_back(pas5);
    optimalPoints_.push_back(pas6);
    
    optimalPoints_.push_back(pas7);
    optimalPoints_.push_back(pas8);
    optimalPoints_.push_back(pas9);

    optimalPoints_.push_back(pas10);
  }
}


NonOptimalStartingPointsProblem::~NonOptimalStartingPointsProblem() { }


PointAndSolution<string> 
NonOptimalStartingPointsProblem::comb(
                          std::vector<double>::const_iterator first,
                          std::vector<double>::const_iterator last) 
{
  PointAndSolution<string> result;

  if (dimension_ == 2) {
    assert(std::distance(first, last) == 2);
    double xWeight = *first;
    double yWeight = *(first + 1);

    if (xWeight == 0.0 && yWeight == 0.0)
      return PointAndSolution<string>(Point(0.0, 0.0), "error: all weights were 0");

    if (xWeight == 0.0)
      return weaklyOptimalPoints_[1];
    else if (yWeight == 0.0)
      return weaklyOptimalPoints_[0];
    else {
      // Return the best point from the optimalPoints_ vector.
      // - Use a CompareUsingWeightsFunctor object to compare 
      //   PointAndSolution instances.
      CompareUsingWeightsFunctor comp(xWeight, yWeight);
      std::vector< PointAndSolution<string> >::const_iterator min;
      min = std::min_element(optimalPoints_.begin(), 
                             optimalPoints_.end(), comp);
      result = *min;
    }
  }
  else {
    assert(dimension_ == 3);

    assert(std::distance(first, last) == 3);
    double xWeight = *first;
    double yWeight = *(first + 1);
    double zWeight = *(first + 2);

    if (xWeight == 0.0 && yWeight == 0.0 && zWeight == 0.0)
      return PointAndSolution<string>(Point(0.0, 0.0, 0.0), "error: all weights were 0");
    
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
      // Return the best point from the optimalPoints_ vector.
      // - Use a CompareUsingWeightsFunctor object to compare 
      //   PointAndSolution instances.
      CompareUsingWeightsFunctor comp(xWeight, yWeight, zWeight);
      std::vector< PointAndSolution<string> >::const_iterator min;
      min = std::min_element(optimalPoints_.begin(), 
                             optimalPoints_.end(), comp);
      result = *min;
    }
  }

  return result;
}


CompareUsingWeightsFunctor::CompareUsingWeightsFunctor(double xWeight, 
                                                       double yWeight) : 
                    xWeight_(xWeight), yWeight_(yWeight), dimension_(2) { }


CompareUsingWeightsFunctor::CompareUsingWeightsFunctor(double xWeight, 
                                                       double yWeight, 
                                                       double zWeight) : 
                    xWeight_(xWeight), yWeight_(yWeight), 
                    zWeight_(zWeight), dimension_(3) { }


bool 
CompareUsingWeightsFunctor::operator() (const PointAndSolution<string> a, 
                                        PointAndSolution<string> b) 
{
  assert(a.point.dimension() == dimension_);
  assert(b.point.dimension() == dimension_);

  double valueA, valueB;

  valueA = xWeight_ * a.point[0] + yWeight_ * a.point[1];
  valueB = xWeight_ * b.point[0] + yWeight_ * b.point[1];

  if (dimension_ == 3) {
    valueA += zWeight_ * a.point[2];
    valueB += zWeight_ * b.point[2];
  }

  return valueA < valueB;
}


}  // namespace non_optimal_starting_points_problem
