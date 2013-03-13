/*! \file PointAndSolution.h
 *  \brief The declaration of the PointAndSolution<S> class template.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef POINT_AND_SOLUTION_H
#define POINT_AND_SOLUTION_H

#include <vector>

#include "Point.h"
#include "NullObjectException.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
namespace pareto_approximator {


/*! 
 *  \brief A template for a class containing a point in objective space, 
 *         a problem solution and the weights used to obtain that solution.
 *  
 *  PointAndSolution instances are a simple way of wrapping together a 
 *  problem solution, the corresponding point in objective space and the 
 *  weights used (inside comb()) to obtain that solution.
 *  
 *  We use one more bool attribute, namely _isNull, to distinguish a 
 *  special, null, instance. The PointAndSolution::isNull() method will 
 *  return true for null instances and false otherwise. isNull() does 
 *  not check the contained point (it could be a null Point instance 
 *  even the PointAndSolution instance is not null, i.e. 
 *  PointAndSolution::isNull() returns false). The user should never 
 *  assign null Point instances to a PointAndSolution instance though.
 *  
 *  Instances with different template parameters cannot interact. In some 
 *  cases this is obvious but some of its implications are:
 *  - PointAndSolution instances containing the same Point instance but 
 *    with different template parameters (e.g. a PointAndSolution<int> 
 *    a PointAndSolution<std::string> both containing the Point(1.0, 1.0))
 *    are not the same and are not even comparable.
 *  - Null instances with different template parameters (e.g. a null 
 *    PointAndSolution<int> and a null PointAndSolution<double>) are not 
 *    comparable nor interchangeable even though they are both null.
 *
 *  Examples of instances:
 *  - Say the problem is a linear programming problem with two objective 
 *    functions. A solution will be a set of values V for the problem's 
 *    variables (which optimizes some linear combination of the objective 
 *    functions L). The corresponding PointAndSolution instance will contain 
 *    V, a 2-dimensional Point with coordinates the values of the two 
 *    objective functions for the values in V and the weights used in the 
 *    linear combination L. (and the `_isNull = false` attribute)
 *  - Say the problem is finding a minimum spanning tree of a graph with 
 *    two functions C1 and C2 mapping each edge to a cost. A solution will 
 *    be a set of edges K forming a minimum spanning tree of the same graph 
 *    for some linear combination L of C1 and C2. The corresponding 
 *    PointAndSolution instance will contain K, the weights used in the 
 *    linear combination L and a 2-dimensional Point with coordinates:
 *      - the sum of the costs of the edges in K with respect to C1 and 
 *      - the sum of the costs of the edges in K with respect to C2
 *    (and the `_isNull = false` attribute)
 *  - A null PointAndSolution<std::string> instance. It will have 
 *    `_isNull = true` and we should not care about any of it's other 
 *    attributes (they can be anything and can be set to anything). 
 *    Its methods (except for PointAndSolution::isNull()) will throw 
 *    NullObjectException exceptions when called. 
 *    Note that a null instance with a different template parameter (e.g. 
 *    a null PointAndSolution<int> instance) will not be comparable to 
 *    nor interchangeable with the null PointAndSolution<std::string> 
 *    instance even though they are both null. Instances with different 
 *    template parameters cannot interact.
 *  - etc
 *  
 *  PointAndSolution is a class template. The template argument is the 
 *  representation of a problem solution (it depends on the problem).
 *  e.g. a list of edges (for the second example above)
 *
 *  Users do not have to set the weights themselves (inside comb()) - it 
 *  will be done automatically after comb() returns.
 *  
 *  \sa PointAndSolution(), ~PointAndSolution() and operator<()
 */
template <class S> 
class PointAndSolution
{
  public:
    //! The empty constructor. Creates a null PointAndSolution<S> instance.
    /*!
     *  Instances with different template parameters cannot interact so, 
     *  for example, a null PointAndSolution<std::string> and a null 
     *  PointAndSolution<int> will not be comparable nor interchangeable 
     *  even though they are both null.
     *  
     *  \sa PointAndSolution
     */
    PointAndSolution();

    //! A constructor initializing all attributes except weightsUsed.
    PointAndSolution(const Point & p, const S & s);

    //! A constructor initializing PointAndSolution's attributes.
    /*! 
     *  \param p A Point object.
     *  \param s An S object. A problem solution.
     *  \param first Iterator to the initial position in a std::vector<double> 
     *               containing the weights used to obtain p and s.
     *  \param last Iterator to the final (past-the-end) position in a 
     *              std::vector<double> containing the weights used to 
     *              obtain p and s.
     */
    PointAndSolution(const Point & p, const S & s, 
                     std::vector<double>::const_iterator first,
                     std::vector<double>::const_iterator last);

    //! \brief PointAndSolution's default destructor. (empty)
    ~PointAndSolution();

    //! \brief Sets the object's "point" attribute. (also makes it non-null)
    void setPoint(const Point & p);

    //! \brief Sets the object's attributes. (also makes it non-null)
    void setAttributes(const Point & p, const S & s);

    //! Is the instance null?
    /*!
     *  \return true if the instance is a null instance; false otherwise.
     *  
     *  \sa PointAndSolution and Point
     */
    bool isNull() const;

    /*! 
     *  \brief Check if the contained Point instance is strictly positive 
     *         (i.e. all coordinates strictly greater than zero).
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if either this or the 
     *    contained Point instance is null.
     *  
     *  \sa Point
     */
    bool isStrictlyPositive() const;

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
    inline bool operator== (const PointAndSolution & pas) const;

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
     *  - May throw a NullObjectException exception if either 
     *    PointAndSolution instance (or either of the contained Point 
     *    instances) is null.
     *  - May throw a DifferentDimensionsException if the two Points are of 
     *    different dimensions (can't be compared).
     *  
     *  \sa PointAndSolution and Point::operator<()
     */
    inline bool operator< (const PointAndSolution<S> & pas) const;

    //! Check if this instance's point eps-covers the given instance's point.
    /*!
     *  \param pas A PointAndSolution<S> instance whose point (q) has 
     *             \f$ q_{i} \ge 0 \f$ for all i.
     *  \param eps An approximation factor.
     *  \return true if this instance's point (p) eps-covers the given 
     *          instance's point (q); false otherwise.
     *  
     *  This method is just a pass-through to Point::dominates().
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if either 
     *    PointAndSolution instance (or either of the contained Point 
     *    instances) is null.
     *  - May throw a NotPositivePointException (or 
     *    NotStrictlyPositivePointException if Point::dominates() is using the 
     *    multiplicative error measure) exception if either p or q is not 
     *    positive (strictly positive, respectively), i.e. some coordinate is 
     *    less than 0.0 (less than or equal to 0.0, respectively).
     *  - May throw a NegativeApproximationRatioException if \f$ eps < 0 \f$.
     *  - May throw a DifferentDimensionsException if p and q are of 
     *    different dimensions.
     *  
     *  \sa PointAndSolution and Point::dominates()
     */
    inline bool dominates(const PointAndSolution<S> & pas, 
                          double eps=0.0) const;

    //! The dimension of the space that the contained point lives in.
    /*!
     *  Just a shortcut for point.dimension().
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the instance is null.
     *  
     *  \sa PointAndSolution
     */
    unsigned int dimension() const;

    //! A point in objective space.
    Point point;
    //! A problem solution (its type is the template argument).
    S solution;
    //! The weights used (inside comb()) to obtain the solution.
    std::vector<double> weightsUsed;
    //! Is the PointAndSolution instance null?
    bool _isNull;
};


}  // namespace pareto_approximator


/* @} */


// We've got to #include the implementation here because we are describing 
// a class template, not a simple class.
#include "PointAndSolution.cpp"


#endif  // POINT_AND_SOLUTION_H
