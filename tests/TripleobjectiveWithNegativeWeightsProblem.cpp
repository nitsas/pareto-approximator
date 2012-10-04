/*! \file TripleobjectiveWithNegativeWeightsProblem.cpp
 *  \brief Implementation of the TripleobjectiveWithNegativeWeightsProblem 
 *         class, a simple problem class used in BaseProblemTest.cpp.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <assert.h>
#include <iterator>
#include <algorithm>

#include "../Point.h"
#include "TripleobjectiveWithNegativeWeightsProblem.h"


using pareto_approximator::Point;


namespace tripleobjective_with_negative_weights_problem {


TripleobjectiveWithNegativeWeightsProblem::TripleobjectiveWithNegativeWeightsProblem() 
{
  // Make some non-dominated points and a dominated one.
  PointAndSolution<string> pas1(Point(36.0, 349.0, 280.0), "best-on-x");
  PointAndSolution<string> pas2(Point(221.0, 54.0, 261.0), "best-on-y");
  PointAndSolution<string> pas3(Point(342.0, 351.0, 37.0), "best-on-z");

  PointAndSolution<string> pas4(Point(92, 142, 42), "best-other");
  PointAndSolution<string> pas5(Point(100, 143, 43), "dominated");

  points_.push_back(pas1);
  points_.push_back(pas2);
  points_.push_back(pas3);
  points_.push_back(pas4);
  points_.push_back(pas5);
}


TripleobjectiveWithNegativeWeightsProblem::~TripleobjectiveWithNegativeWeightsProblem() { }


PointAndSolution<string> 
TripleobjectiveWithNegativeWeightsProblem::comb(std::vector<double>::const_iterator first,
                                      std::vector<double>::const_iterator last)
{
  assert(std::distance(first, last) == 3);
  double xWeight = *first;
  double yWeight = *(first + 1);
  double zWeight = *(first + 2);

  if (xWeight < 0.0 || yWeight < 0.0 || zWeight < 0.0)
    throw NegativeWeightsException();

  if (yWeight == 0.0 && xWeight == 0.0 && zWeight == 0.0)
    return PointAndSolution<string>(Point(0.0, 0.0, 0.0), "error: all weights were 0");
  
  // return the best point from the points_ vector
  // use a CompareUsingWeights object to compare PointAndSolution instances
  CompareUsingWeightsFunctor comp(xWeight, yWeight, zWeight);
  std::vector< PointAndSolution<string> >::iterator min;
  min = std::min_element(points_.begin(), points_.end(), comp);
  return *min;
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


}  // namespace tripleobjective_with_negative_weights_problem
