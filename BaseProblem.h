/*! \file BaseProblem.h
 *  \brief The declaration of the BaseProblem<S> class template.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H
#define PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H


#include <list>
#include <vector>

#include "Point.h"
#include "Hyperplane.h"
#include "PointAndSolution.h"
#include "NonDominatedSet.h"
#include "NotEnoughAnchorPointsException.h"


using pareto_approximator::PointAndSolution;


/*!
 *  \defgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! A template for a problem wrapper class.
/*!
 *  BaseProblem is a wrapper base class for user problems. Users should 
 *  create a derived class containing problem-specific methods and attributes 
 *  and implement its comb method. 
 *  
 *  Instances of the user-created derived class will represent bi-objective 
 *  problems. Examples of problems: 
 *  - A graph with two "costs" for each edge and the problem of finding a 
 *    minimum spanning tree.
 *  - A similar graph with two "capacities" for each edge and the problem 
 *    of finding a maximum flow.
 *  - A linear optimization problem with two objective functions.
 *  - etc
 *  
 *  BaseProblem is a class template. The template argument is the 
 *  representation of a problem solution (it depends on the problem). 
 *  e.g. a list of edges forming a minimum spanning tree of a graph.
 *
 *  Instances of derived classes will still have access to BaseProblem's 
 *  computeConvexParetoSet(). It will use the user-implemented comb() method 
 *  to find an approximation to the problem's Pareto set (or, depending on 
 *  the problem and the approximation parameter, the actual Pareto set).
 *
 *  \sa BaseProblem(), ~BaseProblem(), comb() and operator()()
 */

template <class S>
class BaseProblem 
{
  public:
    //! A set of points representing a facet of the current approximation.
    typedef typename std::vector< PointAndSolution<S> > Facet;

    //! BaseProblem's default constructor. (empty)
    BaseProblem();

    //! BaseProblem's default destructor. (virtual and empty)
    virtual ~BaseProblem();

    //! Optimize a linear combination of the objectives.
    /*! 
     *  \param first Iterator to the initial position in a std::vector<double> 
     *               containing the weights w_{i}.
     *  \param last Iterator to the final position in the std::vector<double> 
     *              containing the weights w_{i}. (the position right after 
     *              the last weight we want (w_{n}))
     *  \return A PointAndSolution<S> object containing: 
     *          - An optimal solution of the problem with respect to the 
     *            linear combination:
     *            \f$ w_{1} * f_{1} + w_{2} * f_{2} + ... + w_{n} * f_{n} \f$ 
     *            of the objectives (f_{i}) and the weights (w_{i}).
     *          - The corresponding point in objective space. Points returned 
     *            by comb() must have positive coordinates.
     *          - The weights w_{i} used in the linear combination.
     *
     *  computeConvexParetoSet() uses the instance's comb() to optimize linear 
     *  combinations of the objectives in order to come up with an 
     *  approximation of the Pareto curve.
     *  
     *  BaseProblem's comb() is declared pure virtual (doesn't do anything).
     *  The user is expected to derive a class from BaseProblem and 
     *  implement its comb() method.
     *  
     *  /sa BaseProblem(), ~BaseProblem() and operator()()
     */
    virtual PointAndSolution<S> 
    comb(std::vector<double>::const_iterator first, 
         std::vector<double>::const_iterator last) = 0;

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
     *  computeConvexParetoSet() will use the comb() method the user 
     *  implemented. That is why comb() is declared virtual.
     *
     *  Possible exceptions:
     *  - May throw a NotEnoughAnchorPointsException exception if at some step
     *    the number of points for the new base are less than \#numObjectives.
     *  
     *  \sa BaseProblem, PointAndSolution and Point
     */
    std::vector< PointAndSolution<S> > 
    computeConvexParetoSet(unsigned int numObjectives, double eps=1e-12);

  private:
    /*! \brief A function called by computeConvexParetoSet() to do most of 
     *         the work.
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
    std::vector< PointAndSolution<S> > 
    doChord(unsigned int numObjectives, Facet anchors, double eps);

    /*! \brief Generate a new Pareto optimal point using the given facet 
     *         as a generating facet.
     *
     *  \sa comb(), BaseProblem, PointAndSolution and Point
     */
    PointAndSolution<S> 
    generateNewParetoPoint(Facet facet);

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
    std::vector<double>
    computeMeanWeights(Facet base);
};


}  // namespace pareto_approximator


/* @} */


// We have got to #include the implementation here because we are 
// describing a class template, not a simple class.
#include "BaseProblem.cpp"


#endif  // PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H
