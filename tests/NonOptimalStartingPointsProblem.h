/*! \file NonOptimalStartingPointsProblem.h
 *  \brief Declaration of the NonOptimalStartingPointsProblem class, 
 *         a simple problem class used in BaseProblemTest.cpp.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef EXAMPLE_CLASS_NON_OPTIMAL_STARTING_POINTS_PROBLEM_H
#define EXAMPLE_CLASS_NON_OPTIMAL_STARTING_POINTS_PROBLEM_H


#include <string>
#include <vector>

#include "../PointAndSolution.h"
#include "../BaseProblem.h"


using std::string;

using pareto_approximator::PointAndSolution;
using pareto_approximator::BaseProblem;


class NonOptimalStartingPointsProblem : public BaseProblem<string>
{
  public:
    NonOptimalStartingPointsProblem();
    ~NonOptimalStartingPointsProblem();

    PointAndSolution<string> comb(std::vector<double>::const_iterator first, 
                                  std::vector<double>::const_iterator last);
};


#endif  // EXAMPLE_CLASS_NON_OPTIMAL_STARTING_POINTS_PROBLEM_H
