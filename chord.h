/*! \file chord.h
 *  \brief A file containing the declaration of the chord algorithm 
 *         function templates.
 *  
 *  chord.h contains function template declarations but the problem with 
 *  templates is that declaration and definition should be in the same 
 *  file. To overcome this without making the file harder to read we made 
 *  a separate file, chord.cpp, which is #included after all the function 
 *  declarations.
 */


#ifndef CHORD_ALGORITHM_H
#define CHORD_ALGORITHM_H

#include <list>
#include <tr1/functional>

#include "Point.h"
#include "Line2D.h"
#include "PointAndSolution.h"


using std::list;
using std::tr1::function;


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! The main chord algorithm function.
/*! 
 *  \param comb A std::tr1::function object with arguments (double, double) 
 *              that returns a PointAndSolution<S> instance (where S is a 
 *              type of an instance representing a problem solution).  
 *  \param eps The degree of approximation. chordAlgorithm() will find an 
 *             (1+eps)-convex Pareto set of the problem.
 *  \return An (1+eps)-convex Pareto set of the problem whose linear 
 *          combinations of objectives comb optimizes.
 *  
 *  How to use:
 *  Users should create a class (let's call it Problem), deriving from 
 *  BaseProblem<S> with all the data needed for the user's problem and 
 *  they should implement its comb() method. All the comb() method needs to 
 *  do is optimize linear combinations of the problem's objectives and return 
 *  the resulting problem solution and the corresponding point in objective 
 *  space (see BaseProblem for more). The Problem class's instances will 
 *  be functors since we have defined BaseProblem<S>'s operator()() (it 
 *  just calls comb() which is declared virtual in BaseProblem<S>). 
 *  After all the above the user can call chordAlgorithm() with a Problem 
 *  instance and the eps they want.
 *
 *  \sa BaseProblem, PointAndSolution and Point
 */
template <class S>
list< PointAndSolution<S> > 
chordAlgorithm(function<PointAndSolution<S> (double, double)> comb, double eps);


//! A recursive function called by chordAlgorithm() to do the bulk of the work.
/*! 
 *  \param comb A std::tr1::function object with arguments (double, double) 
 *              that returns a PointAndSolution<S> instance (where S is a 
 *              type of an instance representing a problem solution).  
 *  \param west A PointAndSolution<S> instance (where S is the type of the 
 *              problem solutions). 
 *  \param south A PointAndSolution<S> instance (where S is the type of the 
 *               problem solutions).
 *  \param tip A Point instance which, together with the points in "west" 
 *             and "south" forms the triangle inside which doChord() will 
 *             search for points of the (1+eps)-convex Pareto set.
 *  \param eps The degree of approximation. doChord() will find a subset 
 *             of an (1+eps)-convex Pareto set of the problem.
 *  \return The part of the problem's (1+eps)-convex Pareto set between the 
 *          point "tip", the point in "west" and the one in "south".
 *  
 *  Users don't need to use doChord(). It's just a recursive routine the 
 *  chordAlgorithm() function uses to do the bulk of the work.
 *  
 *  Each time it's called doChord() finds at most one new (1+eps)-convex 
 *  Pareto set point, splits the problem into two subproblems and calls 
 *  itself recursivelly on the subproblems until the requested degree of 
 *  approximation is met.
 *  
 *  Please read "How good is the Chord lgorithm?" by Constantinos Daskalakis, 
 *  Ilias Diakonikolas and Mihalis Yannakakis for in-depth info on how the 
 *  chord algorithm works.
 *  
 *  \sa chordAlgorithm(), BaseProblem, PointAndSolution and Point
 */
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
