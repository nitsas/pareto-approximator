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


namespace non_optimal_starting_points_problem {


class NonOptimalStartingPointsProblem : public BaseProblem<string>
{
  public:
    NonOptimalStartingPointsProblem(unsigned int dimension=2);
    ~NonOptimalStartingPointsProblem();

    PointAndSolution<string> comb(
                        std::vector<double>::const_iterator first, 
                        std::vector<double>::const_iterator last) const;

  private:
    std::vector< PointAndSolution<string> > optimalPoints_;
    std::vector< PointAndSolution<string> > weaklyOptimalPoints_;
    // problem dimension:
    unsigned int dimension_;
};


// A simple class we'll use as a comparison function (functor).
class CompareUsingWeightsFunctor {
  public:
    CompareUsingWeightsFunctor(double xWeight, double yWeight);
    CompareUsingWeightsFunctor(double xWeight, double yWeight, double zWeight);
    
    bool operator() (PointAndSolution<string> a, PointAndSolution<string> b);
  
  private:
    double xWeight_, yWeight_, zWeight_;
    // problem dimension:
    unsigned int dimension_;
};


}  // namespace non_optimal_starting_points_problem


#endif  // EXAMPLE_CLASS_NON_OPTIMAL_STARTING_POINTS_PROBLEM_H
