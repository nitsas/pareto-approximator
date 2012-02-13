#ifndef PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H
#define PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H


#include "PointAndSolution.h"


using pareto_approximator::PointAndSolution;


namespace pareto_approximator {


// A problem wrapper. Users should inherit from this and implement comb.
template <class S>
class SimpleProblem 
{
  public:
    SimpleProblem() { }
    virtual ~SimpleProblem() { }

    virtual PointAndSolution<S> 
    comb(double xWeight, double yWeight)
    {
      return PointAndSolution<S>();
    }

    PointAndSolution<S> 
    operator() (double xWeight, double yWeight)
    { 
      return comb(xWeight, yWeight); 
    }
};


}  // namespace pareto_approximator


#endif  // PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H
