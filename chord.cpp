/*! \file chord.cpp
 *  \brief A file containing the definition of the chord algorithm 
 *         function templates.
 *
 *  Won't #include "chord.h". In fact "chord.h" will #include "chord.cpp" 
 *  because it describes function templates (which don't allow us to split 
 *  declaration from definition).
 */


using std::min;


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
chordAlgorithm(function<PointAndSolution<S> (double, double)> comb, double eps)
{
  // reminder: comb's arguments are x objective's weight and y's weight

  // find the westmost (best on x objective) and the southmost (best on y
  // objective) solutions (and the corresponding points in objective space)
  PointAndSolution<S> west = comb(1.0, 0.0);
  PointAndSolution<S> south = comb(0.0, 1.0);

  // find that point in objective space which corresponds to the best possible 
  // solution we could expect (we'll pass it to doChord as a sort of lower 
  // limit on the search space)
  Point tip(min(west.point.x, south.point.x), min(west.point.y, south.point.y));

  // let doChord do all the work (it's recursive)
  return doChord<S>(comb, west, south, tip, eps);
}


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
        const Point& tip, double eps)
{
  // reminder: comb's arguments are x objective's weight and y's weight

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
  westList  = doChord<S>(comb, west,  southwest, westTip,  eps);
  southList = doChord<S>(comb, southwest, south, southTip, eps);
  // remove southList's first element (it's the same as westList's last one)
  // and merge the lists
  southList.pop_front();
  westList.merge(southList);

  return westList;
}


}  // namespace pareto_approximator


