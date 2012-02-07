/* PointAndSolution.h */


#ifndef POINT_AND_SOLUTION_H
#define POINT_AND_SOLUTION_H

#include "Point.h"


namespace pareto_approximator {


template <class S> 
class PointAndSolution
{
  public:
    PointAndSolution() {}
    PointAndSolution(const Point& p, const S& s) : point(p), solution(s) {}
    ~PointAndSolution() {}

    bool operator< (const PointAndSolution<S>& pas) const 
            throw(DifferentDimensionsException) { return point < pas.point; }

    Point point;
    S solution;
};


}  // namespace pareto_approximator


#endif  // POINT_AND_SOLUTION_H
