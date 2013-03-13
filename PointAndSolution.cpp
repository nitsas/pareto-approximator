/*! \file PointAndSolution.cpp
 *  \brief The definition of the PointAndSolution<S> class template.
 *  \author Christos Nitsas
 *  \date 2012
 *  
 *  Won't `include` PointAndSolution.h. In fact, PointAndSolution.h 
 *  will `include` PointAndSolution.cpp because it describes a class 
 *  template (which doesn't allow us to split declaration from definition).
 */


/*!
 *  \weakgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
namespace pareto_approximator {


template <class S> 
//! The empty constructor. Creates a null PointAndSolution<S> instance.
PointAndSolution<S>::PointAndSolution() : _isNull(true) { }


//! A constructor initializing all attributes except weightsUsed.
template <class S> 
PointAndSolution<S>::PointAndSolution(const Point & p, const S & s) : 
                              point(p), solution(s), _isNull(false) { }


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
template <class S> 
PointAndSolution<S>::PointAndSolution(const Point & p, const S & s, 
                             std::vector<double>::const_iterator first,
                             std::vector<double>::const_iterator last) :
                                point(p), solution(s), _isNull(false)
{
  weightsUsed.assign(first, last);
}


//! PointAndSolution's default destructor. (empty)
template <class S> 
PointAndSolution<S>::~PointAndSolution() { }


//! \brief Sets the object's "point" attribute. (also makes it non-null)
template <class S> 
void 
PointAndSolution<S>::setPoint(const Point & p)
{
  point = p;
  _isNull = false;
}


//! \brief Sets the object's attributes. (also makes it non-null)
template <class S> 
void 
PointAndSolution<S>::setAttributes(const Point & p, const S & s)
{
  point    = p;
  solution = s;
  _isNull  = false;
}


//! Is the instance null?
/*!
 *  \return true if the instance is a null instance; false otherwise.
 *  
 *  \sa PointAndSolution and Point
 */
template <class S> 
bool 
PointAndSolution<S>::isNull() const
{
  return _isNull;
}


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
template <class S> 
bool 
PointAndSolution<S>::isStrictlyPositive() const
{
  if (isNull())
    throw exception_classes::NullObjectException();

  return point.isStrictlyPositive();
}


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
template <class S> 
bool 
PointAndSolution<S>::operator== (const PointAndSolution & pas) const
{
  if (isNull() and pas.isNull())
    return true;
  else if (isNull() or pas.isNull())
    return false;
  // else

  return (this->point == pas.point);
}


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
 *  - May throw a DifferentDimensionsException exception if the two Points 
 *    are of different dimensions (can't be compared).
 *  
 *  \sa PointAndSolution and Point::operator<()
 */
template <class S> 
bool 
PointAndSolution<S>::operator< (const PointAndSolution<S> & pas) const 
{ 
  if (isNull() or pas.isNull()) 
    throw exception_classes::NullObjectException();
  // else 

  return (this->point < pas.point); 
}


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
template <class S> 
bool 
PointAndSolution<S>::dominates(const PointAndSolution<S> & pas, 
                               double eps) const
{
  if (isNull() or pas.isNull())
    throw exception_classes::NullObjectException();

  return this->point.dominates(pas.point, eps);
}


//! The dimension of the space that the contained point lives in.
/*!
 *  Just a shortcut for point.dimension().
 *  
 *  Possible exceptions:
 *  - May throw a NullObjectException exception if the instance is null.
 *  
 *  \sa PointAndSolution
 */
template <class S> 
unsigned int 
PointAndSolution<S>::dimension() const
{
  if (isNull())
    throw exception_classes::NullObjectException();

  return point.dimension();
}


}  // namespace pareto_approximator


/*! @} */
