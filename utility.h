/*! \file utility.h
 *  \brief The declaration of some utility functions.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef PARETO_APPROXIMATOR_UTILITIES_H
#define PARETO_APPROXIMATOR_UTILITIES_H

#include <vector>
#include <list>

#include "Point.h"
#include "PointAndSolution.h"
#include "Facet.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
namespace pareto_approximator {


//! The namespace containing utility functions. 
/*! 
 *  (i.e. functions called by class methods e.t.c.)
 */
namespace utility {


//! Compute the facets of the convex hull of the given set of points.
/*!
 *  \param points A (const reference to a) std::vector of points. 
 *                (PointAndSolution<S> instances)
 *  \param spaceDimension The dimension of the space that the points live in.
 *  \return A list containing all the facets (Facet<S> instances) of the 
 *          convex hull.
 *  
 *  This function requires that the external program/tool qconvex, 
 *  distributed with qhull (see www.qhull.org) be installed on the system 
 *  (and be on the PATH).
 *  
 *  This function currently only works for Unix-like systems, it won't 
 *  work on Windows. It has only been tested on Mac OS X Mountain Lion 
 *  but should work on other Unix-like systems as well.
 *  
 *  \sa BaseProblem::doCraft()
 */
template <class S> 
typename std::list< Facet<S> > 
computeConvexHullFacets(const std::vector< PointAndSolution<S> > & points, 
                        unsigned int spaceDimension);


//! Compute the convex hull of the given set of points.
/*!
 *  \param points A (const reference to a) std::vector of points. 
 *                (PointAndSolution<S> instances)
 *  \param spaceDimension The dimension of the space that the points live in.
 *  \return A vector containing all the extreme points (PointAndSolution<S> 
 *          instances) of the convex hull.
 *  
 *  This function requires that the external program/tool qconvex, 
 *  distributed with qhull (see www.qhull.org) be installed on the system 
 *  (and be on the PATH).
 *  
 *  This function currently only works for Unix-like systems, it won't 
 *  work on Windows. It has only been tested on Mac OS X Mountain Lion 
 *  but should work on other Unix-like systems as well.
 *  
 *  \sa BaseProblem::doCraft()
 */
template <class S> 
typename std::list< PointAndSolution<S> > 
computeConvexHull(const std::vector< PointAndSolution<S> > & points, 
                  unsigned int spaceDimension);


/*!
 *  \brief Filter a sequence of PointAndSolution instances and return 
 *         only the non-dominated ones.
 *
 *  \param first An iterator to the first element in the sequence.
 *  \param last An iterator to the past-the-end element in the sequence.
 *  
 *  We will use a pareto_approximator::NonDominatedSet to discard dominated 
 *  points.
 *  
 *  \sa NonDominatedSet
 */
template <class S> 
std::vector< PointAndSolution<S> > 
filterDominatedPoints(
        typename std::vector< PointAndSolution<S> >::const_iterator first, 
        typename std::vector< PointAndSolution<S> >::const_iterator last);


//! Discard facets not useful for generating new Pareto points.
/*!
 *  \param facets A (reference to a) list of facets.
 *  
 *  Discard facets with all normal vector coefficients non-positive (<= 0).
 *  Facets with no positive normal vector coefficient are not useful for 
 *  generating new Pareto optimal points.
 *  
 *  Only facets with all-positive or mixed (i.e. containing at least some 
 *  positive coefficients) normal vectors can be used to generate new 
 *  Pareto optimal points.
 *  
 *  \sa Facet, BaseProblem::doChord() and BaseProblem::doCraft()
 */
template <class S> 
void 
discardUselessFacets(std::list< Facet<S> > & facets);


/*! \brief Choose the Facet instance with the largest local approximation 
 *         error upper bound from sequence of Facet instances.
 *  
 *  \param first A const_iterator to the first element in the sequence.
 *  \param last A const_iterator to the past-the-end element in the sequence.
 *  \return A const_iterator to the first element in the range that has the  
 *          largest local approximation error upper bound. If no element 
 *          is a non-boundary facet the function returns "last".
 *  
 *  The Facet::getLocalApproximationErrorUpperBound() method is used (of 
 *  course) for the facet's local approximation error upper bound.
 *  
 *  Boundary facets (i.e. those with isBoundaryFacet() == true) are ignored.
 *  
 *  If all the facets in the sequence are boundary facets the iterator 
 *  "last" is returned.
 *  
 *  \sa Facet
 */
template <class S> 
typename std::list< Facet<S> >::iterator 
chooseFacetWithLargestLocalApproximationErrorUpperBound(
                    typename std::list< Facet<S> >::iterator first, 
                    typename std::list< Facet<S> >::iterator last);


//! Generate a weight vector (for comb()) using the given facet.
/*!
 *  \param facet A Facet instance. (Its vertices' weightsUsed attributes 
 *               will be needed if the facet's normal vector is not 
 *               all-positive.)
 *  \return A std::vector of objective weights (for BaseProblem::comb()).
 *  
 *  The resulting weights will be:
 *  - Either the facet's normal vector. 
 *    (if it has no negative elements)
 *  - Or the mean of the weights used to obtain the facet's vertices.
 *  
 *  \sa BaseProblem, BaseProblem::comb() and 
 *      BaseProblem::generateNewParetoPoint()
 */
template <class S> 
std::vector<double> 
generateNewWeightVector(const Facet<S> & facet);


}  // namespace utility


}  // namespace pareto_approximator


/* @} */


// We've got to #include the implementation here because this file contains 
// function templates, not only simple functions.
#include "utility.cpp"


#endif  // PARETO_APPROXIMATOR_UTILITIES_H
