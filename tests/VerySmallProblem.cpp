/* VerySmallProblem.cpp */


#include "../Point.h"
#include "VerySmallProblem.h"


using pareto_approximator::Point;


VerySmallProblem::VerySmallProblem() { }


VerySmallProblem::~VerySmallProblem() { }


PointAndSolution<string> 
VerySmallProblem::comb(double xWeight, double yWeight)
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


