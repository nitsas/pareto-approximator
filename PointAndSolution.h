/*! \file PointAndSolution.h
 *  \brief The declaration of the PointAndSolution<S> class template.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef POINT_AND_SOLUTION_H
#define POINT_AND_SOLUTION_H

#include "Point.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


/*! 
 *  \brief A template for a class containing a point in objective space 
 *         and a problem solution. 
 *  
 *  PointAndSolution instances are a simple way of wrapping together a 
 *  problem solution and the corresponding point in objective space.
 *  
 *  Examples of instances:
 *  - Say the problem is a linear programming problem with two objective 
 *    functions. A solution will be a set of values V for the problem's 
 *    variables (which optimizes some linear combination of the objective 
 *    functions). The corresponding PointAndSolution instance will contain 
 *    V and a 2-dimensional Point with coordinates the values of the two 
 *    objective functions for the values in V.
 *  - Say the problem is finding a minimum spanning tree of a graph with 
 *    two functions C1 and C2 mapping each edge to a cost. A solution will 
 *    be a set of edges K forming a minimum spanning tree of the same graph 
 *    for some linear combination of C1 and C2. The corresponding 
 *    PointAndSolution instance will contain K and a 2-dimensional Point 
 *    with coordinates: 
 *      - the sum of the costs of the edges in K with respect to C1 and 
 *      - the sum of the costs of the edges in K with respect to C2
 *  - etc
 *  
 *  PointAndSolution is a class template. The template argument is the 
 *  representation of a problem solution (it depends on the problem).
 *  e.g. a list of edges (for the second example above)
 *
 *  \sa PointAndSolution(), ~PointAndSolution() and operator<()
 */
template <class S> 
class PointAndSolution
{
  public:
    //! PointAndSolution's default constructor. (empty)
    PointAndSolution();
    //! A constructor initializing PointAndSolution's attributes.
    PointAndSolution(const Point& p, const S& s);
    //! PointAndSolution's default destructor. (empty)
    ~PointAndSolution();

    //! PointAndSolution equality operator.
    /*!
     *  \param pas A PointAndSolution instance to compare with the current 
     *             instance.
     *  \return true if the point in pas is the same as the one in the 
     *          current instance; false otherwise.
     *
     *  Compare the points in the two PointAndSolutionInstances.
     *  
     *  \sa PointAndSolution
     */
    inline bool operator==(const PointAndSolution& pas) const;

    //! The PointAndSolution less-than operator. 
    /*! 
     *  Let's call L the instance on the left of the operator and R the one 
     *  on the right. Returns true if L.point is less than R.point using 
     *  Point::operator<().
     *  
     *  Compare the two PointAndSolution instances according to the Point 
     *  instances they contain.
     *  
     *  Possible exceptions:
     *  - May throw a DifferentDimensionsException if the two Points are of 
     *    different dimensions (can't be compared).
     *  
     *  \sa PointAndSolution and Point::operator<()
     */
    inline bool operator< (const PointAndSolution<S>& pas) const;

    //! A point in objective space.
    Point point;
    //! A problem solution (its type is the template argument).
    S solution;
};


}  // namespace pareto_approximator


/* @} */


// We've got to #include the implementation here because we are describing 
// a class template, not a simple class.
#include "PointAndSolution.cpp"


#endif  // POINT_AND_SOLUTION_H
