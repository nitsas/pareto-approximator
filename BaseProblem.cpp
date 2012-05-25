/*! \file BaseProblem.cpp
 *  \brief The definition of the BaseProblem<S> class template.
 *  \author Christos Nitsas
 *  \date 2012
 *  
 *  Won't `#include` "BaseProblem.h". In fact "BaseProblem.h" will 
 *  `#include` "BaseProblem.cpp" because it describes a class template 
 *  (which doesn't allow us to split declaration from definition).
 */


#include <assert.h>
#include <algorithm>


using std::min;


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


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
 *  \param numObjectives The number of objectives to minimize. Note: The 
 *                       user's comb() routine should be able to handle a 
 *                       std::vector<double> of \#numObjectives weights.
 *  \param eps The degree of approximation. computeConvexParetoSet() will 
 *             find an (1+eps)-convex Pareto set of the problem.
 *  \return A NonDominatedSet< PointAndSolution<S> > instance containing 
 *          a (1+eps)-convex Pareto set of the problem whose linear 
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
 *  Possible exceptions:
 *  - May throw a NotEnoughBasePointsException exception if at some step
 *    the number of points for the new base are less than \#numObjectives.
 *  
 *  \sa BaseProblem, PointAndSolution and Point
 */
template <class S> 
std::list< PointAndSolution<S> > 
BaseProblem<S>::computeConvexParetoSet(unsigned int numObjectives, double eps)
{
  // reminder: comb's arguments are a set of iterators over a 
  // std::vector<double> of weights (one for each objective)

  assert(numObjectives <= 3);

  // Find a best solution for each objective. We'll end up with up to  
  // \#numObjectives (possibly less) different solutions.
  std::vector<double> weights(numObjectives, 0.0);
  std::vector< PointAndSolution<S> > base;
  for (unsigned int i = 0; i != numObjectives; ++i) {
    weights[i] = 1.0;
    base.push_back(comb(weights.begin(), weights.end()));
    weights[i] = 0.0;
  }

  NonDominatedSet< PointAndSolution<S> > nds(base.begin(), base.end());
  if (nds.size() < numObjectives)
    throw NotEnoughBasePointsException();

  // make a point with the worst (biggest) of the base points' coordinates
  // (will need it as a direction to avoid when minimizing the weighted sum) 
  std::vector<double> maxCoords(numObjectives, 0.0);
  typename std::vector< PointAndSolution<S> >::iterator bi;
  for (bi = base.begin(); bi != base.end(); ++bi)
    for (unsigned int i = 0; i != numObjectives; ++i)
      if (bi->point[i] > maxCoords[i])
        maxCoords[i] = bi->point[i];
  Point pointToMoveAwayFrom(maxCoords.begin(), maxCoords.end());
  
  // let doChord do all the work (it's recursive)
  NonDominatedSet< PointAndSolution<S> > resultSet(base.begin(), base.end());
  std::list< PointAndSolution<S> > resultsFromDoChord;
  resultsFromDoChord = doChord(numObjectives, base, pointToMoveAwayFrom, eps);
  resultSet.insert(resultsFromDoChord.begin(), resultsFromDoChord.end());
  std::list< PointAndSolution<S> > resultList(resultSet.begin(), resultSet.end());

  return resultList;
}


/*! \brief A recursive function called by computeConvexParetoSet() to do 
 *         the bulk of the work.
 * 
 *  \param numObjectives The number of objectives to minimize. Note: The 
 *                       user's comb() routine should be able to handle a 
 *                       std::vector<double> of \#numObjectives weights.
 *  \param base A std::vector of PointAndSolution<S> instances (where S is 
 *              the type of the problem solutions).
 *  \param pointToMoveAwayFrom A Point instance. Helps us determine the 
 *                             right direction to minimize towards.
 *  \param eps The degree of approximation. doChord() will find a subset 
 *             of an (1+eps)-convex Pareto set of the problem.
 *  \return The part of the problem's (1+eps)-convex Pareto set between 
 *          the points in base and 0.
 *  
 *  Users don't need to use doChord() - that is why it's declared private. 
 *  It's just a recursive routine the computeConvexParetoSet() method uses 
 *  to do the bulk of the work.
 *  
 *  Each time it's called doChord() finds at most one new (1+eps)-convex 
 *  Pareto set point (inside the convex polytope defined by the given 
 *  points: tip and the points in base), splits the problem into 
 *  subproblems and calls itself recursivelly on the subproblems until 
 *  the requested degree of approximation is met.
 *  
 *  Please read "How good is the Chord Algorithm?" by Constantinos 
 *  Daskalakis, Ilias Diakonikolas and Mihalis Yannakakis for in-depth 
 *  info on how the chord algorithm works.
 *
 *  Possible exceptions:
 *  - May throw a NotEnoughBasePointsException exception if at some step
 *    the number of points for the new base are less than \#numObjectives.
 *  
 *  \sa computeConvexParetoSet(), BaseProblem, PointAndSolution and Point
 */
template <class S> 
std::list< PointAndSolution<S> > 
BaseProblem<S>::doChord(unsigned int numObjectives, 
                        std::vector< PointAndSolution<S> > base, 
                        const Point & pointToMoveAwayFrom, double eps)
{
  // reminder: comb accepts a set of iterators to the objectives' weights

  assert(base.size() <= 3);
  assert(numObjectives <= 3);
  assert(base.size() == numObjectives);

  // make a \#numObjectives-dimensional hyperplane that passes 
  // through all the points in base
  std::vector<Point> basePoints;
  typename std::vector< PointAndSolution<S> >::iterator bi;
  for (bi = base.begin(); bi != base.end(); ++bi)
    basePoints.push_back(bi->point);
  Hyperplane baseHyperplane(basePoints.begin(), basePoints.end());

  // Make sure the hyperplane faces in the right direction (away from 
  // pointToMoveAwayFrom). That is, when we'll call comb to minimize a 
  // weighted sum of the objectives (with the hyperplane's coefficients 
  // as the weights) minimizing will move us away from pointToMoveAwayFrom.
  baseHyperplane.faceAwayFrom(pointToMoveAwayFrom);

  // call comb using baseHyperplane's coefficients (essentially 
  // baseHyperplane's slope) as weights 
  PointAndSolution<S> opt = comb(baseHyperplane.begin(), baseHyperplane.end());

  // check if the point we just found is approximately dominated by the 
  // points we have so far
  if (baseHyperplane.ratioDistance(opt.point) <= eps)
    return std::list< PointAndSolution<S> >();
  // else

  // keep (for the new base) only those base points that opt doesn't dominate
  std::vector< PointAndSolution<S> > newBase;
  Point newPointToMoveAwayFrom;
  for (bi = base.begin(); bi != base.end(); ++bi)
    if (!opt.dominates(*bi))
      newBase.push_back(*bi);
    else
      // won't use it for the new base since opt dominates it but we can 
      // use it as a point to move away from (when "combing")
      newPointToMoveAwayFrom = bi->point;

  // split the problem into subproblems and call doChord() on each of them
  std::list< PointAndSolution<S> > resultList;
  resultList.push_back(opt);
  if (newBase.size() < numObjectives - 1)
    // not enough points to form a new #numObjectives-hyperplane
    throw NotEnoughBasePointsException();
  else if (newBase.size() == numObjectives - 1) {
    // enough points (with opt) for exactly one subproblem
    newBase.push_back(opt);
    std::list< PointAndSolution<S> > subproblemResults;
    subproblemResults = doChord(numObjectives, newBase, newPointToMoveAwayFrom, eps);
    resultList.splice(resultList.end(), subproblemResults);
  }
  else {
    // enough points (with opt) for #numObjectives subproblems
    assert(newBase.size() == numObjectives);
    // opt is not yet included in newBase - newBase is the same as base
    for (unsigned int i = 0; i != numObjectives; ++i) {
      // use the point that will not be in the new base as a 
      // point to move away from
      newPointToMoveAwayFrom = newBase[i].point;
      // temporarily replace one of the old base points with opt 
      newBase[i] = opt;
      // call doChord on the subproblem 
      std::list< PointAndSolution<S> > subproblemResults;
      subproblemResults = doChord(numObjectives, newBase, newPointToMoveAwayFrom, eps);
      // append the subproblem results to resultList
      resultList.splice(resultList.end(), subproblemResults);
      // restore newBase (replace opt with the old base point)
      newBase[i] = base[i];
    }
  }

  return resultList;
}


}  // namespace pareto_approximator


/* @} */
