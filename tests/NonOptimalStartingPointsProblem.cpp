/* NonOptimalStartingPointsProblem.cpp */


#include "../Point.h"
#include "NonOptimalStartingPointsProblem.h"


using pareto_approximator::Point;


NonOptimalStartingPointsProblem::NonOptimalStartingPointsProblem() { }


NonOptimalStartingPointsProblem::~NonOptimalStartingPointsProblem() { }


PointAndSolution<string> 
NonOptimalStartingPointsProblem::comb(double xWeight, double yWeight)
{
  if (yWeight == 0.0 && xWeight == 0.0)
    return PointAndSolution<string>(Point(0.0, 0.0), "error: both weights were 0");
  
  if (yWeight == 0.0) 
    return PointAndSolution<string>(Point(1.0, 8.0), "west-dominated");
  
  if (xWeight == 0.0)
    return PointAndSolution<string>(Point(8.0, 1.0), "south-dominated");
  
  if (yWeight < xWeight / 2.0) 
    return PointAndSolution<string>(Point(1.0, 4.0), "west");
  else if (yWeight > 2 * xWeight) 
    return PointAndSolution<string>(Point(4.0, 1.0), "south");
  else 
    return PointAndSolution<string>(Point(2.0, 2.0), "southwest");
}


