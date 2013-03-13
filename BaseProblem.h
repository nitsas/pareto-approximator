/*! \file BaseProblem.h
 *  \brief The declaration of the BaseProblem<S> class template.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H
#define PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H


#include <vector>
#include <list>

#include "Facet.h"
#include "PointAndSolution.h"


using pareto_approximator::Facet;
using pareto_approximator::PointAndSolution;


/*!
 *  \defgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
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
 *    minimum spanning tree. (biobjective problem)
 *  - A similar graph with three "capacities" for each edge and the problem 
 *    of finding a maximum flow. (triple objective problem)
 *  - A linear optimization problem with four objective functions.
 *    (four objectives problem)
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
    //! BaseProblem's default constructor. (empty)
    BaseProblem();

    //! BaseProblem's default destructor. (virtual and empty)
    virtual ~BaseProblem();

    //! Optimize a linear combination of the objectives.
    /*! 
     *  \param first Iterator to the first element in a std::vector<double> 
     *               containing the weights w_{i}.
     *  \param last Iterator to the past-the-end element in the 
     *              std::vector<double> containing the weights w_{i}. (the 
     *              position right after the last weight we want (w_{n}))
     *  \return A PointAndSolution<S> object containing: 
     *          - An optimal solution of the problem with respect to the 
     *            linear combination:
     *            \f$ w_{1} * f_{1} + w_{2} * f_{2} + ... + w_{n} * f_{n} \f$ 
     *            of the objectives (f_{i}) and the weights (w_{i}).
     *          - The corresponding point in objective space. Points returned 
     *            by comb() must be strictly positive (i.e. all their 
     *            coordinates must be strictly greater than zero).
     *          - The weights w_{i} used in the linear combination. The 
     *            comb() does not need to set the weights, it will be 
     *            done automatically after comb() returns.
     *          - (The PointAndSolution instance's _isNull attribute does 
     *            not need to be set inside comb(), it too will be set 
     *            automatically after comb() returns.)
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

    //! Compute an (1+eps)-approximate convex Pareto set of the problem.
    /*! 
     *  \param numObjectives The number of objectives to minimize. Note: The 
     *                       user's comb() routine should be able to handle a 
     *                       std::vector<double> of \#numObjectives weights.
     *  \param eps The degree of approximation. computeConvexParetoSet() will 
     *             find an (1+eps)-approximate convex Pareto set of the 
     *             problem.
     *  \return An (1+eps)-approximate convex Pareto set of the problem whose 
     *          linear combinations of objectives comb optimizes.
     *  
     *  How to use:
     *  Users should create a class (let's call it Problem), deriving from 
     *  BaseProblem<S> with all the data needed for the user's problem and 
     *  they should implement its comb() method. All the comb() method needs 
     *  to do is optimize linear combinations of the problem's objectives and 
     *  return the resulting problem solution and the corresponding point in 
     *  objective space. After all the above, users can make a Problem 
     *  instance and call its computeConvexParetoSet() method with the eps 
     *  they want.
     *
     *  computeConvexParetoSet() will use the comb() method the user 
     *  implemented. That is why comb() is declared virtual.
     *  
     *  computeConvexParetoSet() initializes the usedWeightVectors_ 
     *  attribute to an empty list every time it is called (before it 
     *  calls any other method).
     *
     *  \sa BaseProblem, PointAndSolution and Point
     */
    std::vector< PointAndSolution<S> > 
    computeConvexParetoSet(unsigned int numObjectives, double eps=1e-12);

  private:
    /*! \brief A function called by computeConvexParetoSet() to do most of 
     *         the work. (for biobjective optimization problems)
     * 
     *  \param anchors The Facet defined by the anchor points.
     *  \param eps The degree of approximation. 
     *  \return A vector of Pareto optimal points (PointAndSolution instances).
     *          BaseProblem::computeConvexParetoSet() will filter them to make 
     *          the (1+eps)-approximate convex Pareto set.
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
    std::vector< PointAndSolution<S> > 
    doChord(Facet<S> anchors, double eps);

    /*! \brief A function that uses Craft's (et al.) algorithm to 
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
    std::vector< PointAndSolution<S> > 
    doCraft(unsigned int numObjectives, Facet<S> anchors, double eps);

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
    PointAndSolution<S> 
    generateNewParetoPointUsingFacet(const Facet<S> & facet);

    /*!
     *  \brief Generate a new Pareto optimal point using the given weights
     *         to call comb().
     *
     *  \param weights A vector of weights for comb().
     *  \return A Pareto optimal point (inside a PointAndSolution<S>  
     *          object) generated using the given weights if the weights 
     *          were not used before; a null PointAndSolution<S> object 
     *          otherwise.
     *          
     *  This method will call the user-implemented comb() method (using 
     *  the given weight vector) to make a Pareto point.
     *  
     *  Every time the method is called with a weight vector W it 
     *  checks if W has been used before (using the usedWeightVectors_
     *  attribute):
     *  - If they have, it returns a null PointAndSolution instance 
     *    without calling comb().
     *  - If they have not, it calls comb() using the given weights (W) 
     *    and adds W to the usedWeightsVectors_ list. 
     *  
     *  \sa BaseProblem, comb(), generateNewParetoPointUsingFacet(), 
     *      PointAndSolution and Point
     */
    PointAndSolution<S> 
    generateNewParetoPoint(const std::vector<double> & weights);

    /*! 
     *  \brief A list of already used weight vectors so that we never
     *         call comb() with the same weights a second time.
     *  
     *  Every time generateNewParetoPoint() is called with a weight vector 
     *  W it checks usedWeightVectors_ for W.
     *  - If W has been used before, it will not call comb().
     *  - If it has not, it calls comb() using W as weights and adds W
     *    to this list. 
     *  
     *  Places where it is used (and how it is used):
     *  - Initialized to an empty list inside (the constructor and)
     *    BaseProblem::computeConvexParetoSet(). 
     *  - Maintained inside the BaseProblem::generateNewParetoPoint() 
     *    method.
     *  
     *  \sa BaseProblem, computeConvexParetoSet() and 
     *      generateNewParetoPoint()
     */
    std::list< std::vector<double> > usedWeightVectors_;
};


}  // namespace pareto_approximator


/* @} */


// We have got to #include the implementation here because we are 
// describing a class template, not a simple class.
#include "BaseProblem.cpp"


#endif  // PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H
