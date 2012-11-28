/*! \file TripleobjectiveWithNegativeWeightsProblem.h
 *  \brief Declaration of the TripleobjectiveWithNegativeWeightsProblem class, 
 *         a simple problem class used in BaseProblemTest.cpp.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef EXAMPLE_CLASS_TRIPLEOBJECTIVE_WITH_NEGATIVE_WEIGHTS_PROBLEM_H
#define EXAMPLE_CLASS_TRIPLEOBJECTIVE_WITH_NEGATIVE_WEIGHTS_PROBLEM_H


#include <string>
#include <vector>

#include "../PointAndSolution.h"
#include "../BaseProblem.h"


using std::string;

using pareto_approximator::PointAndSolution;
using pareto_approximator::BaseProblem;


namespace tripleobjective_with_negative_weights_problem {


class NegativeWeightsException : public std::exception
{
  public:
    // Return a simple char* message.
    const char * what() const throw()
    {
      return "Comb was called with some negative weights.";
    }
};


class TripleobjectiveWithNegativeWeightsProblem : public BaseProblem<string>
{
  public:
    TripleobjectiveWithNegativeWeightsProblem();
    ~TripleobjectiveWithNegativeWeightsProblem();

    PointAndSolution<string> comb(
                        std::vector<double>::const_iterator first, 
                        std::vector<double>::const_iterator last) const;

  private:
    std::vector< PointAndSolution<string> > points_;
};


// A simple class we'll use as a comparison function (functor).
class CompareUsingWeightsFunctor {
  public:
    CompareUsingWeightsFunctor(double xWeight, double yWeight, double zWeight);
    
    bool operator() (PointAndSolution<string> a, PointAndSolution<string> b);
  
  private:
    double xWeight_, yWeight_, zWeight_;
};


}  // namespace tripleobjective_with_negative_weights_problem


#endif  // EXAMPLE_CLASS_TRIPLEOBJECTIVE_WITH_NEGATIVE_WEIGHTS_PROBLEM_H
