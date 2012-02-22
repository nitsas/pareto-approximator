/* VerySmallProblem.h */


#ifndef EXAMPLE_CLASS_VERY_SMALL_PROBLEM_H
#define EXAMPLE_CLASS_VERY_SMALL_PROBLEM_H


#include <string>

#include "../PointAndSolution.h"
#include "../BaseProblem.h"


using std::string;

using pareto_approximator::PointAndSolution;
using pareto_approximator::BaseProblem;


class VerySmallProblem : public BaseProblem<string>
{
  public:
    VerySmallProblem();
    ~VerySmallProblem();

    PointAndSolution<string> comb(double xWeight, double yWeight);
};


#endif  // EXAMPLE_CLASS_VERY_SMALL_PROBLEM_H
