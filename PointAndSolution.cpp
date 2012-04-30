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
  return (point == pas.point);
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
  return point < pas.point; 
}


}  // namespace pareto_approximator


/*! @} */
