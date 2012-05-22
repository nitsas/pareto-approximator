/*! \file PointAndSolution.cpp
 *  \brief The definition of the PointAndSolution<S> class template.
 *  \author Christos Nitsas
 *  \date 2012
 *  
 *  Won't `#include` "PointAndSolution.h". In fact, "PointAndSolution.h" 
 *  will `#include` "PointAndSolution.cpp" because it describes a class 
 *  template (which doesn't allow us to split declaration from definition).
 */


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


template<class S> 
//! PointAndSolution's default constructor. (empty)
PointAndSolution<S>::PointAndSolution() { }


//! A constructor initializing PointAndSolution's attributes.
template<class S> 
PointAndSolution<S>::PointAndSolution(const Point& p, const S& s) : point(p), solution(s) { }


//! PointAndSolution's default destructor. (empty)
template<class S> 
PointAndSolution<S>::~PointAndSolution() { }


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
template<class S> 
bool 
PointAndSolution<S>::operator== (const PointAndSolution& pas) const
{
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
 *  - May throw a DifferentDimensionsException if the two Points are of 
 *    different dimensions (can't be compared).
 *  
 *  \sa PointAndSolution and Point::operator<()
 */
template<class S> 
bool 
PointAndSolution<S>::operator< (const PointAndSolution<S>& pas) const 
{ 
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
 *  From here on down, we'll call this instance's point p and the given 
 *  instance's point q for convenience.
 *  
 *  Note that both p and q must be greater than zero (dominated by 0);
 *  that is both \f$ p_{i} \ge 0 \f$ and \f$ q_{i} \ge 0 \f$ must hold
 *  for all i.
 *  
 *  We say that p \f$ \epsilon \f$-covers q (\f$\epsilon \ge 0 \f$) if 
 *  \f$ p_{i} \le (1 + \epsilon) q_{i} \f$, for all i. Both p and 
 *  q must be of the same dimension.
 *  
 *  If eps=0.0 the method simply checks whether or not p dominates 
 *  q and that is how it got its name.
 *  
 *  Possible exceptions:
 *  - May throw a NotPositivePointException if either p or q is not 
 *    greater than 0 (dominated by 0).
 *  - May throw a NegativeApproximationRatioException if \f$ eps < 0 \f$.
 *  - May throw a DifferentDimensionsException if p and q are of 
 *    different dimensions.
 *  
 *  \sa PointAndSolution and Point::dominates()
 */
template <class S> 
bool 
PointAndSolution<S>::dominates(const PointAndSolution<S>& pas, 
                               double eps) const
{
  return this->point.dominates(pas.point, eps);
}


}  // namespace pareto_approximator


/*! @} */
