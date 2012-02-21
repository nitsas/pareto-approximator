/*! \file SimpleProblem.h
 *  \brief A file containing the declaration and definition of the 
 *         SimpleProblem<S> class template.
 */


#ifndef PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H
#define PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H


#include "PointAndSolution.h"


using pareto_approximator::PointAndSolution;


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! A template for a problem wrapper class.
/*!
 *  SimpleProblem is a wrapper parent class for user problems. Users should 
 *  create a child class containing problem-specific methods and attributes 
 *  and implement its comb method. 
 *  
 *  Instances of the user-created child class will represent bi-objective 
 *  problems. Examples of problems: 
 *  - A graph with two "costs" for each edge and the problem of finding a 
 *    minimum spanning tree.
 *  - A similar graph with two "capacities" for each edge and the problem 
 *    of finding a maximum flow.
 *  - A linear optimization problem with two objective functions.
 *  - etc
 *  
 *  SimpleProblem is a class template. The template argument is the 
 *  representation of a problem solution (it depends on the problem). 
 *  e.g. a list of edges forming a minimum spanning tree of a graph.
 *
 *  Users can then pass an instance of the child class to chordAlgorithm() 
 *  to find an approximation to the problem's Pareto set (or, depending on 
 *  the problem and the approximation parameter, the actual Pareto set).
 *
 *  \sa SimpleProblem(), ~SimpleProblem(), comb() and operator()()
 */

template <class S>
class SimpleProblem 
{
  public:
    //! SimpleProblem's default constructor. (empty)
    SimpleProblem() { }
    //! SimpleProblem's default destructor. (empty)
    virtual ~SimpleProblem() { }

    //! Optimize a linear combination of the objectives.
    /*! 
     *  \param xWeight The weight of the x objective.
     *  \param yWeight The weight of the y objective.
     *  \return A PointAndSolution<S> object containing: 
     *          - An optimum solution of the problem with respect to the 
     *            linear combination \f$ xWeight * x + yWeight * y \f$ of the 
     *            objectives and
     *          - the corresponding point in objective space.
     *
     *  chordAlgorithm() uses an instance's comb() to optimize linear 
     *  combinations of the objectives in order to come up with an 
     *  approximation of the Pareto curve.
     *  
     *  SimpleProblem's comb() is declared virtual (and doesn't do anything).
     *  The user is supposed to make a child of the SimpleProblem class and 
     *  implement its comb() method.
     *  
     *  /sa SimpleProblem(), ~SimpleProblem() and operator()()
     */
    virtual PointAndSolution<S> 
    comb(double xWeight, double yWeight)
    {
      return PointAndSolution<S>();
    }

    //! Call the object's comb() method with the given arguments.
    /*! 
     *  \param xWeight The weight of the x objective.
     *  \param yWeight The weight of the y objective.
     *  \return A PointAndSolution<S> object containing: 
     *          - An optimum solution of the problem with respect to the 
     *            linear combination \f$ xWeight * x + yWeight * y \f$ of the 
     *            objectives and
     *          - the corresponding point in objective space.
     *
     *  /sa SimpleProblem(), ~SimpleProblem, comb()
     */
    PointAndSolution<S> 
    operator() (double xWeight, double yWeight)
    { 
      return comb(xWeight, yWeight); 
    }
};


}  // namespace pareto_approximator


#endif  // PARETO_APPROXIMATOR_SIMPLE_PROBLEM_H
