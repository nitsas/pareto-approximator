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


using std::list;

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
     *          - the corresponding point in objective space. Points returned 
     *            by comb() must have positive coordinates.
     *
     *  computeConvexParetoSet() uses the instance's comb() to optimize linear 
     *  combinations of the objectives in order to come up with an 
     *  approximation of the Pareto curve.
     *  
     *  BaseProblem's comb() is declared pure virtual (doesn't do anything).
     *  The user is expected to derive a class from the BaseProblem class 
     *  and implement its comb() method.
     *  
     *  /sa BaseProblem(), ~BaseProblem() and operator()()
     */
    virtual PointAndSolution<S> 
    comb(std::vector<double>::const_iterator first, 
         std::vector<double>::const_iterator last) = 0;

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
    list< PointAndSolution<S> > 
    computeConvexParetoSet(unsigned int numObjectives, double eps=0.0);

  private:
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
     *  Please read "How good is the Chord Algorithm?" by Constantinos 
     *  Daskalakis, Ilias Diakonikolas and Mihalis Yannakakis for in-depth 
     *  info on how the chord algorithm works.
     *  
     *  \sa computeConvexParetoSet(), BaseProblem, PointAndSolution and Point
     */
    list< PointAndSolution<S> > 
    doChord(const PointAndSolution<S>& west, 
            const PointAndSolution<S>& south, const Point& tip, double eps);
};


}  // namespace pareto_approximator


/* @} */


// We've got to #include the implementation here because we are describing 
// a class template, not a simple class.
#include "BaseProblem.cpp"


#endif  // PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H
