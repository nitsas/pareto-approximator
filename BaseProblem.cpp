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
#include <stack>


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
 *  computeConvexParetoSet() will use the comb() method that the user 
 *  implemented. That is why comb() is declared virtual.
 *
 *  Possible exceptions:
 *  - May throw a NotEnoughAnchorPointsException exception if at some step
 *    the number of points for a newFacet are less than \#numObjectives.
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
  // Let's call the corresponding points (in objective space) anchor points
  // and the PointAndSolution instances that contain them anchors.
  std::vector<double> weights(numObjectives, 0.0);
  Facet anchors;
  for (unsigned int i = 0; i != numObjectives; ++i) {
    weights[i] = 1.0;
    PointAndSolution<S> pas = comb(weights.begin(), weights.end());
    pas.weightsUsed.assign(weights.begin(), weights.end());
    anchors.push_back(pas);
    weights[i] = 0.0;
  }

  NonDominatedSet< PointAndSolution<S> > nds(anchors.begin(), anchors.end());
  if (nds.size() < numObjectives)
    throw NotEnoughAnchorPointsException();

  // Let doChord do all the work.
  NonDominatedSet< PointAndSolution<S> > filter(anchors.begin(), anchors.end());
  std::list< PointAndSolution<S> > resultsFromDoChord;
  resultsFromDoChord = doChord(numObjectives, anchors, eps);
  filter.insert(resultsFromDoChord.begin(), resultsFromDoChord.end());
  std::list< PointAndSolution<S> > resultList(filter.begin(), filter.end());

  return resultList;
}


/*! \brief A function called by computeConvexParetoSet() to do most of 
 *         the work.
 * 
 *  \param numObjectives The number of objectives to minimize. 
 *                       Note: The user's comb() routine should be able to 
 *                       handle a std::vector<double> of \#numObjectives 
 *                       weights.
 *  \param anchors The Facet defined by the anchor points.
 *  \param eps The degree of approximation. 
 *  \return A list of Pareto optimal points (PointAndSolution instances).
 *          BaseProblem::computeConvexParetoSet() will filter them to make 
 *          the (1+eps)-approximate convex Pareto set.
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
 *  You can also read "Approximating convex Pareto surfaces in 
 *  multiobjective radiotherapy planning" by David L. Craft et al. for 
 *  info on how we handle facets whose normal vector has both positive 
 *  and negative components.
 *
 *  Possible exceptions:
 *  - May throw a NotEnoughAnchorPointsException exception if at some step
 *    the number of points for the new base are less than \#numObjectives.
 *  
 *  \sa computeConvexParetoSet(), BaseProblem, PointAndSolution and Point
 */
template <class S> 
std::list< PointAndSolution<S> > 
BaseProblem<S>::doChord(unsigned int numObjectives, Facet anchors, double eps)
{
  // reminder: comb accepts a set of iterators to the objectives' weights

  assert(anchors.size() <= 3);
  assert(numObjectives <= 3);
  assert(anchors.size() == numObjectives);

  // a stack of Facets to try (for generating new Pareto optimal points)
  std::stack< Facet > facetsToTry;
  facetsToTry.push(anchors);

  std::list< PointAndSolution<S> > resultList;
  while (not facetsToTry.empty()) {
    // Get a facet from the stack and try to generate a new Pareto 
    // optimal point using that facet.
    Facet generatingFacet = facetsToTry.top();
    facetsToTry.pop();

    // Make a \#numObjectives-dimensional hyperplane that passes 
    // through all the points in "anchors".
    std::vector<Point> facetPoints;
    typename Facet::iterator fi;
    for (fi = generatingFacet.begin(); fi != generatingFacet.end(); ++fi)
      facetPoints.push_back(fi->point);
    Hyperplane generatingHyperplane(facetPoints.begin(), facetPoints.end());

    PointAndSolution<S> opt;
    if (generatingHyperplane.hasAllAiCoefficientsNonNegative()) {
      // Call comb using generatingHyperplane's coefficients (that is, 
      // generatingHyperplane's slope) as weights.
      generatingHyperplane.normalizeAiCoefficients();
      opt = comb(generatingHyperplane.begin(), generatingHyperplane.end());
      opt.weightsUsed.assign(generatingHyperplane.begin(), 
                             generatingHyperplane.end());
    }
    else {
      // "meanWeights" is a std::vector<double> of weights mw_{i}, where:
      // /f$ mw_{i} = sum_{j=1}^{j=generatingFacet.size()} ( w_{ij} ) /f$,
      // where w_{ij} is the i'th of the weights used to obtain the j'th 
      // point in "generatingFacet".

      std::vector<double> meanWeights = computeMeanWeights(generatingFacet);
      Hyperplane meanHyperplane(meanWeights.begin(), meanWeights.end(), 0.0);
      meanHyperplane.normalizeAiCoefficients();

      // Now call comb() using meanHyperplane's coefficients.
      opt = comb(meanHyperplane.begin(), meanHyperplane.end());
      opt.weightsUsed.assign(meanHyperplane.begin(), meanHyperplane.end());
    }

    // Check if the point we just found is approximately dominated by the 
    // points we have so far. If it is ignore it. (try the next facet)
    if (generatingHyperplane.ratioDistance(opt.point) <= eps)
      continue;
    // else
    resultList.push_back(opt);

    // Keep (for the new facet) only those points in generatingFacet 
    // that opt doesn't dominate.
    Facet newFacet;
    for (fi = generatingFacet.begin(); fi != generatingFacet.end(); ++fi)
      if (!opt.dominates(*fi))
        newFacet.push_back(*fi);

    // Make the new facets (using generatingFacet and opt) and push them 
    // onto the stack.
    if (newFacet.size() < numObjectives - 1)
      // not enough points to form a new facet
      throw NotEnoughAnchorPointsException();
    else if (newFacet.size() == numObjectives - 1) {
      // enough points (including opt) for exactly one new facet
      newFacet.push_back(opt);
      facetsToTry.push(newFacet);
    }
    else {
      // enough points (including opt) for #numObjectives new facets
      assert(newFacet.size() == numObjectives);
      // opt is not yet included in newFacet - newFacet is currently 
      // the same as generatingFacet.
      for (unsigned int i = 0; i != numObjectives; ++i) {
        // Temporarily replace one of the old facet points with opt.
        newFacet[i] = opt;
        // Push the new facet on top of the stack.
        facetsToTry.push(newFacet);
        // Restore newFacet (replace opt with the old facet point).
        newFacet[i] = generatingFacet[i];
      }
    }
  }   // while (not facetsToTry.empty())

  return resultList;
}


//! Computes the mean of all the weight vectors in "base".
/*!
 *  \param base A std::vector of PointAndSolution<S> instances (where S is 
 *              the type of the problem solutions).
 *  \return A weight vector W of size base.size(). Each element W_{j} is
 *          the mean of all w_{ij}'s, where w_{i} is the weight vector 
 *          inside the i'th element of "base".
 *  
 *  \sa BaseProblem and PointAndSolution
 */
template <class S>
std::vector<double>
BaseProblem<S>::computeMeanWeights(Facet base)
{
  typename Facet::iterator bi;

  std::vector<double> meanWeights(base.size(), 0.0);
  for (unsigned int i = 0; i != base.size(); ++i) {
    for (bi = base.begin(); bi != base.end(); ++bi)
      meanWeights[i] += bi->weightsUsed[i];
    meanWeights[i] = meanWeights[i] / base.size();
  }
  
  return meanWeights;
}


}  // namespace pareto_approximator


/* @} */
