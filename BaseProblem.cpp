/*! \file BaseProblem.cpp
 *  \brief The definition of the BaseProblem<S> class template.
 *  \author Christos Nitsas
 *  \date 2012
 *  
 *  Won't `include` BaseProblem.h. In fact BaseProblem.h will 
 *  `include` BaseProblem.cpp because it describes a class template 
 *  (which doesn't allow us to split declaration from definition).
 */


#include <assert.h>
#include <algorithm>
#include <stack>

#include "Point.h"
#include "NonDominatedSet.h"
#include "utility.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
namespace pareto_approximator {


//! BaseProblem's default constructor. (empty)
template <class S> 
BaseProblem<S>::BaseProblem() 
{
  // not really needed (it is done automatically)
  // - added this line just to make the code more explicit
  usedWeightVectors_ = std::list< std::vector<double> >();
}


//! BaseProblem's default destructor. (virtual and empty)
template <class S>
BaseProblem<S>::~BaseProblem() { }


//! Compute an (1+eps)-approximate convex Pareto set of the problem.
/*! 
 *  \param numObjectives The number of objectives to minimize. Note: The 
 *                       user's comb() routine should be able to handle a 
 *                       std::vector<double> of \#numObjectives weights.
 *  \param eps The degree of approximation. computeConvexParetoSet() will 
 *             find an (1+eps)-approximate convex Pareto set of the problem.
 *  \return An (1+eps)-approximate convex Pareto set of the problem whose 
 *          linear combinations of objectives comb optimizes.
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
 *  computeConvexParetoSet() will use the comb() method that the user 
 *  implemented. That is why comb() is declared virtual.
 *
 *  computeConvexParetoSet() initializes the usedWeightVectors_ 
 *  attribute to an empty list every time it is called (before it 
 *  calls any other method).
 *
 *  \sa BaseProblem, PointAndSolution and Point
 */
template <class S> 
std::vector< PointAndSolution<S> > 
BaseProblem<S>::computeConvexParetoSet(unsigned int numObjectives, 
                                       double eps) 
{
  // reminder: comb's arguments are a set of iterators over a 
  // std::vector<double> of weights (one for each objective)

  assert(eps >= 0.0);
  assert(numObjectives >= 2);
  assert(numObjectives <= 3);

  // Clear the list of used weight vectors.
  // - In case computeConvexParetoSet() was called earlier.
  usedWeightVectors_.clear();

  // Find a best solution for each objective. 
  // - We'll end up with up to \#numObjectives (possibly less) different 
  //   solutions. 
  // - Let's call the corresponding points (in objective space) anchor 
  //   points and the PointAndSolution instances that contain them anchors.
  // - Using 0 as weights in the linear combination of objective functions 
  //   may get us a weakly Pareto optimal point (i.e. a point that can be 
  //   dominated but not strongly dominated). In case we do not want this, 
  //   we can use a very small positive number (e.g. eps/2) where we would 
  //   use 0.
  // CHANGE temporary
  std::vector<double> weights(numObjectives, eps/2);
//  std::vector<double> weights(numObjectives, 0.0);
  std::vector< PointAndSolution<S> > anchors;
  for (unsigned int i = 0; i != numObjectives; ++i) {
    // make only the i'th element of the weight vector non-zero
    weights[i] = 1.0;
    // generate an anchor
    PointAndSolution<S> anchor = generateNewParetoPoint(weights);
    assert(not anchor.isNull());
    anchors.push_back(anchor);
    // restore the weight vector's i'th element (all zero again)
    // CHANGE temporary
    weights[i] = eps/2;
//    weights[i] = 0.0;
  }

  // Filter the anchor points (some might be weakly-dominated by others).
  // - We might even have 1 anchor point that dominates all the others. In 
  //   that case just return the single anchor point as the result.
  // - We use a NonDominatedSet for the filtering.
  std::vector< PointAndSolution<S> > results;
  NonDominatedSet< PointAndSolution<S> > nds(anchors.begin(), anchors.end());

  assert( (nds.size() > 0) && (nds.size() <= numObjectives) );
  if (nds.size() == 1) 
    // We are very lucky, we got a single solution that is optimum in 
    // every objective!!!
    results.assign(nds.begin(), nds.end());
  else if ( (nds.size() > 1) && (nds.size() < numObjectives) ) 
    // Not enough anchor points to continue.
    // - Return the anchor points we have so far.
    results.assign(nds.begin(), nds.end());
  else {
    // Exactly \#numObjectives anchor points - enough to continue.
    // (no anchor point was dominated by any other)

    // make the convex hull of the anchor points (it is just a single facet)
    Facet<S> anchorFacet(anchors.begin(), anchors.end());

    // Call doChord() for biobjective problems or doPgen() for 
    // more than two objectives.
    std::vector< PointAndSolution<S> > unfilteredResults;
    if (numObjectives == 2) {
      // anchorFacet is just a single line segment for numObjectives == 2

      // Let doChord do all the work.
      unfilteredResults = doChord(anchorFacet, eps);
    }
    else {
      assert(numObjectives > 2);

      // Let doPgen do all the work.
      unfilteredResults = doPgen(numObjectives, anchorFacet, eps);
    }

    // Filter the results.
    // - Some of the anchor points might be weakly Pareto optimal, so 
    //   some of the points computed by doChord might dominate them. 
    results = pareto_approximator::utility::
                      filterDominatedPoints<S>(unfilteredResults.begin(), 
                                               unfilteredResults.end());
  }

  return results;
}


/*! \brief A function called by computeConvexParetoSet() to do most of 
 *         the work. (for biobjective optimization problems)
 * 
 *  \param anchorFacet The Facet defined by the anchor points.
 *  \param eps The degree of approximation. 
 *  \return A vector of Pareto optimal points (PointAndSolution instances).
 *          It might contain weakly-dominated points (some of the anchor 
 *          points). 
 *  
 *  Note: doChord() is only called for problems with exactly 2 criteria.
 *  
 *  Users don't need to use doChord() - that is why it's declared private. 
 *  It's just a routine that BaseProblem::computeConvexParetoSet() uses 
 *  to do most of the work.
 *  
 *  doChord() has a big while loop that processes facets (from a stack). 
 *  On each iteration doChord() finds at most one new Pareto optimal point, 
 *  makes new facets using that point and pushes the new facets on the 
 *  stack (iff they improve the approximation; if the requested degree of 
 *  approximation has been met, no new facets are pushed onto the stack).
 *  
 *  Please read "How good is the Chord Algorithm?" by Constantinos 
 *  Daskalakis, Ilias Diakonikolas and Mihalis Yannakakis for in-depth 
 *  info on how the chord algorithm works.
 *  
 *  \sa computeConvexParetoSet(), BaseProblem, PointAndSolution and Point
 */
template <class S> 
std::vector< PointAndSolution<S> > 
BaseProblem<S>::doChord(Facet<S> anchorFacet, double eps) 
{
  // reminder: comb accepts a set of iterators to the objectives' weights

  assert(anchorFacet.spaceDimension() == 2);

  // a vector that will hold all the approximation points:
  std::vector< PointAndSolution<S> > results;
  results.assign(anchorFacet.beginVertex(), anchorFacet.endVertex());

  // a stack of Facets to try (for generating new Pareto optimal points):
  std::stack< Facet<S> > facetsToTry;
  facetsToTry.push(anchorFacet);

  while (not facetsToTry.empty()) {
    // Get a facet from the stack and try to generate a new Pareto 
    // optimal point using that facet.
    Facet<S> generatingFacet = facetsToTry.top();
    facetsToTry.pop();

    // if the facet's local approximation error upper bound is less 
    // than the tolerance move on to the next facet
    if (generatingFacet.getLocalApproximationErrorUpperBound() <= eps)
      continue;

    PointAndSolution<S> opt = 
                        generateNewParetoPointUsingFacet(generatingFacet);

    // We will never encounter the same weight vector twice in 
    // biobjective problems. 
    //assert(not opt.isNull());
    if (opt.isNull())
      continue;

    // Note: the generatingFacet will always have an all-positive normal 
    //       vector in biobjective problems

    // Check if the point we just found is approximately dominated by the 
    // facet that made it (i.e. dominated by some convex combination of 
    // the facet's two vertices). 
    // - If it is dominated ignore it. (try the next facet)
    // - It is dominated if it is one of the facet's two vertices.
    if (generatingFacet.dominates(opt.point, eps))
      continue;
    // else

    // Add opt to the list of approximation points.
    results.push_back(opt);

    // Keep (for the new facet/facets) only those vertices of 
    // generatingFacet that opt doesn't dominate. 
    std::vector< PointAndSolution<S> > newFacetVertices;
    typename Facet<S>::ConstVertexIterator fvi;
    for (fvi = generatingFacet.beginVertex(); 
         fvi != generatingFacet.endVertex(); ++fvi) 
      if (!opt.dominates(*fvi))
        newFacetVertices.push_back(*fvi);

    // Make the new facets (using generatingFacet and opt) and push them 
    // onto the stack.
    if (newFacetVertices.size() == 0)
      // opt dominated both of generatingFacet's vertices.
      // - This can only happen if generatingFacet's vertices were both 
      //   anchor points. 
      // - Don't make any new facets, opt is the utopia point - we will 
      //   not find any more points. (don't need any)
      continue;
    else if (newFacetVertices.size() == 1) {
      // enough points (including opt) for exactly one new facet
      newFacetVertices.push_back(opt);
      Facet<S> newFacet(newFacetVertices.begin(), newFacetVertices.end());
      facetsToTry.push(newFacet);
    }
    else {
      assert(newFacetVertices.size() == 2);
      // Enough points (including opt) for two new facets.
      // opt is not yet included in newFacetVertices - newFacetVertices 
      // currently contains the vertices of generatingFacet
      PointAndSolution<S> tempVertex;
      for (unsigned int i = 0; i != 2; ++i) {
        // Temporarily replace one of the old facet vertices with opt.
        tempVertex = newFacetVertices[i];
        newFacetVertices[i] = opt;
        // Push the new facet on top of the stack.
        Facet<S> newFacet(newFacetVertices.begin(), newFacetVertices.end());
        facetsToTry.push(newFacet);
        // Restore newFacetVertices (replace opt with the old facet vertex).
        newFacetVertices[i] = tempVertex;
      }
    }
  }   // while (not facetsToTry.empty())

  return results;
}


/*! \brief A function that uses the PGEN algorithm (Craft et al.) to 
 *         approximate the Pareto set.
 *  
 *  \param numObjectives The number of objectives to minimize. 
 *                       Note: The user's comb() routine should be able to 
 *                       handle a std::vector<double> of \#numObjectives 
 *                       weights.
 *  \param anchors The Facet defined by the anchor points.
 *  \param eps The degree of approximation. 
 *  \return A vector of Pareto optimal points (PointAndSolution instances).
 *          BaseProblem::computeConvexParetoSet() will filter them to make 
 *          the (1+eps)-approximate convex Pareto set.
 *
 *  Please read "Approximating convex Pareto surfaces in multiobjective 
 *  radiotherapy planning" by David L. Craft et al. (2006) for more 
 *  info on the algorithm.
 *  
 *  \sa computeConvexParetoSet(), BaseProblem, PointAndSolution and Point
 */
template <class S> 
std::vector< PointAndSolution<S> > 
BaseProblem<S>::doPgen(unsigned int numObjectives, Facet<S> anchorFacet, 
                        double eps) 
{
  // reminder: comb accepts a set of iterators to the objectives' weights

  assert(numObjectives <= 3);
  assert(anchorFacet.spaceDimension() <= 3);
  assert(anchorFacet.spaceDimension() == numObjectives);

  unsigned int spaceDimension = numObjectives;

  std::vector< PointAndSolution<S> > 
                            approximationPoints(anchorFacet.beginVertex(), 
                                                anchorFacet.endVertex());

  // We need to add one more point (interior point) before we can call qconvex.

  // Make a Pareto point using anchorFacet as a generating facet.
  PointAndSolution<S> interiorPoint = 
                      generateNewParetoPointUsingFacet(anchorFacet);

  // Is interiorPoint either an existing point or coplanar with the facet?
  if ( anchorFacet.isCoplanarWith(interiorPoint.point) || 
       (std::find(approximationPoints.begin(), 
                 approximationPoints.end(), 
                 interiorPoint) != approximationPoints.end()) ) {
    // InteriorPoint is either an existing point (no new Pareto points 
    // found) or is coplanar with the anchor facet (no new Pareto points 
    // found beneath the anchor facet).
    // - No facets to make except for anchorFacet, which we already tried.
    //   We cannot generate any new points.
    // - Return approximation points found so far.
    // - If anchorFacet did not have an all-positive normal vector there 
    //   might be other Pareto optimal points we couldn't find (on the convex 
    //   hull of the Pareto set of course; we can't find Pareto points 
    //   inside the convex hull either way). The facet's normal vector not 
    //   being all-positive  would make us use the mean of its vertices' 
    //   weightsUsed attributes as weights (which is kind of an arbitrary 
    //   choice) and apparently they did not produce a new Pareto point. Is 
    //   is not completely unlikely that some other weight vector might 
    //   produce one but we have no systematic way to try every one of the 
    //   infinite possible weight vectors.
    return approximationPoints;
  }
  // else 

  // Add the new point to the existing set of approximation points.
  approximationPoints.push_back(interiorPoint);

  // Compute the convex hull of the approximation points. 
  // - Each facet has a local approximation error upper bound built-in.
  // - That local approximation error upper bound is computed during the 
  //   construction of the Facet instance.
  std::list< Facet<S> > facets;
  facets = pareto_approximator::utility::
                       computeConvexHullFacets<S>(approximationPoints, 
                                                  spaceDimension);
  // Discard facets with all-negative normal vectors.
  pareto_approximator::utility::discardUselessFacets<S>(facets);

  while (not facets.empty()) {
    // Choose the facet with largest local approximation error upper bound.
    typename std::list< Facet<S> >::iterator generatingFacet;
    generatingFacet = pareto_approximator::utility::
                    chooseFacetWithLargestLocalApproximationErrorUpperBound<S>(
                                          facets.begin(), facets.end());

    // Were there any facets (except boundary facets)? 
    if (generatingFacet != facets.end()) {
      // generatingFacet is not a boundary facet

      // Have we reached the required approximation factor?
      // - the facet is surely not a boundary facet, it surely has a 
      //   local approximation error upper bound
      // - if we have reached the required approximation factor stop 
      //   the algorithm
      if (generatingFacet->getLocalApproximationErrorUpperBound() <= eps) 
        break;
      // else 
    }
    else {
      // Choose the boundary facet with the smallest angle.
      // The angle for a boundary facet is defined as the angle 
      // between the facet normal vector n and the average of the 
      // vertex weight vectors w.
      generatingFacet = pareto_approximator::utility::
                        chooseBoundaryFacetWithSmallestAngle<S>(
                                           facets.begin(), facets.end());

      if (generatingFacet == facets.end()) {
        // There are no more facets. Exit
        break;
      }

      assert(generatingFacet->isBoundaryFacet());
    }

    // Make a new Pareto point using generatingFacet as a generating facet.
    // Reminder: generatingFacet is actually an iterator pointing to 
    //           the actual facet - that is why we use the * operator
    PointAndSolution<S> opt = 
                        generateNewParetoPointUsingFacet(*generatingFacet);

    // Is opt a null instance or an existing point?
    if ( opt.isNull() or 
         std::find(approximationPoints.begin(), 
                   approximationPoints.end(), 
                   opt) != approximationPoints.end() ) {
      // Either we have already tried this set of weights
      // or opt has already been found (using a different weight vector).
      // - discard generatingFacet and go to the next iteration (i.e. 
      //   choose another facet)
      // Reminder: generatingFacet is actually an iterator pointing to one 
      //           of "facets"'s elements
      facets.erase(generatingFacet);
      continue;
    }
    // else 

    // opt is a new point - we will add it to the set of approximation 
    // points and calculate the new convex hull of the set
    approximationPoints.push_back(opt);
    facets = pareto_approximator::utility::
                              computeConvexHullFacets<S>(approximationPoints, 
                                                         spaceDimension);
    pareto_approximator::utility::discardUselessFacets<S>(facets);
  }

  return approximationPoints;
}


/*! 
 *  \brief Generate a new Pareto optimal point using the given Facet 
 *         instance as a generating facet.
 *
 *  \param facet A Facet instance. (Its vertices' weightsUsed 
 *               attributes will be needed if the facet's normal vector 
 *               is not all-positive.)
 *  \return A Pareto optimal point (inside a PointAndSolution<S>  
 *          object) generated using the given facet, i.e. the weights 
 *          generated from the facet, if the weights were not used 
 *          before; a null PointAndSolution<S> object otherwise.
 *          
 *  This method generates a weight vector using the given facet and 
 *  delegates the jobs of making a Pareto point and updating the 
 *  usedWeightVectors_ attribute to generateNewParetoPoint().
 *  
 *  This method will call:
 *  - pareto_approximator::generateNewWeightVector() (using the given 
 *    facet instance as a parameter) to get a weight vector 
 *  - BaseProblem::generateNewParetoPoint() (using the weights it got 
 *    in the previous step) which will in turn call comb() to make a 
 *    Pareto point (if the weights were not used before)
 *  
 *  \sa BaseProblem, comb(), generateNewParetoPoint(), 
 *      pareto_approximator::generateNewWeightVector(), 
 *      PointAndSolution and Point
 */
template <class S> 
PointAndSolution<S> 
BaseProblem<S>::generateNewParetoPointUsingFacet(const Facet<S> & facet) 
{
  // Get a weight vector (using the given facet as a generating facet).
  std::vector<double> weights = pareto_approximator::utility::
                                generateNewWeightVector<S>(facet);

  return generateNewParetoPoint(weights);
}


/*!
 *  \brief Generate a new Pareto optimal point using the given weights
 *         to call comb().
 *
 *  \param weights A vector of weights for comb().
 *  \return A Pareto optimal point (inside a PointAndSolution<S> object) 
 *          generated using the given weights if the weights were not 
 *          used before; a null PointAndSolution<S> object otherwise.
 *          
 *  This method will call the user-implemented comb() method (using the 
 *  given weight vector) to make a Pareto point.
 *  
 *  If the user returns a point that is not strictly positive (i.e. not 
 *  every coordinate is greater than zero) a 
 *  NotStrictlyPositivePointException exception will be thrown.
 *  
 *  Every time the method is called with a weight vector W it checks 
 *  if W has been used before (using the usedWeightVectors_ attribute):
 *  - If they have, it returns a null PointAndSolution instance without 
 *    calling comb().
 *  - If they have not, it calls comb() using the given weights (W) 
 *    and adds W to the usedWeightsVectors_ list. 
 *  
 *  Possible exceptions:
 *  - May throw a NotStrictlyPositivePointException exception if the 
 *    point returned by comb() is not strictly positive. (i.e. if one
 *    or more of its coordinates is not greater than zero)
 *  
 *  \sa BaseProblem, comb(), generateNewParetoPointUsingFacet(), 
 *      PointAndSolution and Point
 */
template <class S> 
PointAndSolution<S> 
BaseProblem<S>::generateNewParetoPoint(const std::vector<double> & weights)
{
  // Check if the given weights have been used before.
  bool haveUsedTheseWeightsBefore = false;
  std::list< std::vector<double> >::iterator it;
  for (it = usedWeightVectors_.begin(); it != usedWeightVectors_.end(); ++it)
    if ( std::equal(weights.begin(), weights.end(), it->begin()) ) {
      haveUsedTheseWeightsBefore = true;
      break;
    }

  // Have the given weights been used before?
  if (haveUsedTheseWeightsBefore)
    // Yes, return a null PointAndSolution<S> instance.
    return PointAndSolution<S>();
  // else

  // Call comb() with the given weights.
  PointAndSolution<S> newPoint = comb(weights.begin(), weights.end());

  // Make sure the user didn't return an invalid point:
  // - We are talking about the Point instance contained inside the 
  //   PointAndSolution<S> instance. The PointAndSolution<S> instance's 
  //   _isNull attribute does not have to have been initialized already
  //   (we'll do it below).
  assert(not newPoint.point.isNull());
  // Is the point returned strictly positive? (it should)
  if (not newPoint.point.isStrictlyPositive()) 
    throw exception_classes::NotStrictlyPositivePointException();
  // else

  // Initialize newPoint's weightsUsed and _isNull attributes.
  // - So that the user doesn't have to do it inside comb().
  newPoint.weightsUsed.assign(weights.begin(), weights.end());
  newPoint._isNull = false;
  // Add the newly used weight vector to the usedWeightVectors_ list.
  usedWeightVectors_.push_front(weights);

  return newPoint;
}


}  // namespace pareto_approximator


/* @} */
