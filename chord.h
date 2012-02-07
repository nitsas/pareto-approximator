/* chord.h */


#ifndef CHORD_ALGORITHM_H
#define CHORD_ALGORITHM_H

#include <list>
#include <tr1/functional>

#include "Point.h"
#include "Line2D.h"
#include "PointAndSolution.h"


using std::list;
using std::tr1::function;


namespace pareto_approximator {


template <class S>
list< PointAndSolution<S> > 
chordAlgorithm(function<PointAndSolution<S> (double, double)> comb, double eps);


template <class S> 
list< PointAndSolution<S> > 
doChord(function<PointAndSolution<S> (double, double)> comb, 
        const PointAndSolution<S>& west, const PointAndSolution<S>& south, 
        const Point& tip, double eps);


}  // namespace pareto_approximator


// We've got to include the implementation here because we are describing 
// function templates, not simple functions.
#include "chord.cpp"


#endif  // CHORD_ALGORITHM_H
