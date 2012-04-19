/*! \file BaseProblem.cpp
 *  \brief The definition of the BaseProblem<S> class template.
 *  \author Christos Nitsas
 *  \date 2012
 *  
 *  Won't `#include` "BaseProblem.h". In fact "BaseProblem.h" will 
 *  `#include` "BaseProblem.cpp" because it describes a class template 
 *  (which doesn't allow us to split declaration from definition).
 */


using std::min;


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! BaseProblem's default constructor. (empty)
template <class S> 
BaseProblem<S>::BaseProblem() { }


//! BaseProblem's default destructor. (virtual and empty)
template <class S>
BaseProblem<S>::~BaseProblem() { }


//! Compute an (1+eps)-convex Pareto set of the problem.
/*! 
 *  \param eps The degree of approximation. computeConvexParetoSet() will 
 *             find an (1+eps)-convex Pareto set of the problem.
 *  \return An (1+eps)-convex Pareto set of the problem whose linear 
 *          combinations of objectives comb optimizes.
 *  
 *  How to use:
 *  Users should create a class (let's call it Problem), deriving from 
 *  BaseProblem<S> with all the data needed for the user's problem and 
 *  they should implement its comb() method. All the comb() method needs 
 *  to do is optimize linear combinations of the problem's objectives and 
 *  return the resulting problem solution and the corresponding point in 
 *  objective space. After all the above users can make a Problem 
 *  instance and call its computeConvexParetoSet() method with the eps 
 *  they want.
 *
 *  computeConvexParetoSet() will use the comb() method the user 
 *  implemented. That is why comb() is declared virtual.
 *
 *  \sa BaseProblem, PointAndSolution and Point
 */
template <class S> 
list< PointAndSolution<S> > 
BaseProblem<S>::computeConvexParetoSet(double eps)
{
  // reminder: comb's arguments are x objective's weight and y's weight

  // find the westmost (best on x objective) and the southmost (best on y
  // objective) solutions (and the corresponding points in objective space)
  PointAndSolution<S> west = comb(1.0, 0.0);
  PointAndSolution<S> south = comb(0.0, 1.0);

  // find that point in objective space which corresponds to the best possible 
  // solution we could expect (we'll pass it to doChord as a sort of lower 
  // limit on the search space)
  Point tip(min(west.point[0], south.point[0]), min(west.point[1], south.point[1]));

  // let doChord do all the work (it's recursive)
  return doChord(west, south, tip, eps);
}


/*! \brief A recursive function called by computeConvexParetoSet() to do 
 *         the bulk of the work.
 * 
 *  \param west A PointAndSolution<S> instance (where S is the type of 
 *              the problem solutions). 
 *  \param south A PointAndSolution<S> instance (where S is the type of 
 *               the problem solutions).
 *  \param tip A Point instance which, together with the points in "west" 
 *             and "south" forms the triangle inside which doChord() will 
 *             search for points of the (1+eps)-convex Pareto set.
 *  \param eps The degree of approximation. doChord() will find a subset 
 *             of an (1+eps)-convex Pareto set of the problem.
 *  \return The part of the problem's (1+eps)-convex Pareto set between 
 *          the point "tip", the point in "west" and the one in "south".
 *  
 *  Users don't need to use doChord() - that is why it's declared private. 
 *  It's just a recursive routine the computeConvexParetoSet() method uses 
 *  to do the bulk of the work.
 *  
 *  Each time it's called doChord() finds at most one new (1+eps)-convex 
 *  Pareto set point, splits the problem into two subproblems and calls 
 *  itself recursivelly on the subproblems until the requested degree of 
 *  approximation is met.
 *  
 *  Please read "How good is the Chord lgorithm?" by Constantinos 
 *  Daskalakis, Ilias Diakonikolas and Mihalis Yannakakis for in-depth 
 *  info on how the chord algorithm works.
 *  
 *  \sa computeConvexParetoSet(), BaseProblem, PointAndSolution and Point
 */
template <class S> 
list< PointAndSolution<S> > 
BaseProblem<S>::doChord(const PointAndSolution<S>& west, 
                        const PointAndSolution<S>& south, 
                        const Point& tip, double eps)
{
  // reminder: comb's arguments are x objective's weight and y's weight

  // check if the westmost or the southmost point given is the same as
  // the tip (ideal point) and if either is stop searching and return 
  // the west and south objects
  if (west.point.dominates(south.point)) {
    list< PointAndSolution<S> > resultList;
    // if west.point is better than (dominates) south.point return 
    // only the west object
    resultList.push_back(west);
    return resultList;
  }
  else if (south.point.dominates(west.point)) {
    list< PointAndSolution<S> > resultList;
    // if south.point is better than (dominates) west.point return 
    // only the south object
    resultList.push_back(south);
    return resultList;
  }

  // check if the best possible point is approximately dominated by the 
  // points we have so far
  Line2D ws(west.point, south.point);
  if (ws.ratioDistance(tip) <= eps) {
    list< PointAndSolution<S> > resultList;
    resultList.push_back(west);
    resultList.push_back(south);
    return resultList;
  }
  // else

  // check if ws is a vertical line (slope is infinite) and call comb 
  // using its slope
  PointAndSolution<S> southwest;
  if (!ws.isVertical())
    southwest = comb(-ws.m(), 1.0);
  else
    southwest = comb(1.0, 0.0);

  // check if the point we just found is approximately dominated by the 
  // points we have so far
  if (ws.ratioDistance(southwest.point) <= eps) {
    list< PointAndSolution<S> > resultList;
    resultList.push_back(west);
    resultList.push_back(south);
    return resultList;
  }
  // else

  // split the problem into two subproblems, the west one and the south one
  Line2D parallel = ws.parallelThrough(southwest.point);
  Line2D wt(west.point, tip);
  Line2D ts(tip, south.point);
  Point westTip  = parallel.intersection(wt);
  Point southTip = parallel.intersection(ts);
  list< PointAndSolution<S> > westList, southList;
  // call comb on the two subproblems
  westList  = doChord(west,  southwest, westTip,  eps);
  southList = doChord(southwest, south, southTip, eps);
  // remove southList's first element (it's the same as westList's last one)
  // and merge the lists
  southList.pop_front();
  westList.merge(southList);

  return westList;
}


}  // namespace pareto_approximator


