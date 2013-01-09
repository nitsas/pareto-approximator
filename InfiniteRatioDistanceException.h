/*! \file InfiniteRatioDistanceException.h
 *  \brief The declaration and definition of the InfiniteRatioDistanceException 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef INFINITE_RATIO_DISTANCE_EXCEPTION_H
#define INFINITE_RATIO_DISTANCE_EXCEPTION_H

#include <exception>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
namespace pareto_approximator {


//! The namespace containing all the exception classes.
namespace exception_classes {


/*! 
 *  \brief Exception thrown when we were asked to compute a ratio distance 
 *         and it turns out that it is infinite.
 *
 *  An exception thrown when:
 *  - We were trying to compute the ratio distance between a point p 
 *    and a hyperplane h but h's normal vector was perpendicular to p's 
 *    coordinate vector (and p was not on h), so multiplying p by a 
 *    constant moves p in a direction parallel to h. (multiplying p will 
 *    never land it on h)
 *  - We were trying to compute the ratio distance between a point p 
 *    and a facet f but f's normal vector was perpendicular to p's 
 *    coordinate vector (and p was not on f's supporting hyperplane), so 
 *    multiplying p by a constant moves p in a direction parallel to f. 
 *    (multiplying p will never land it on f's supporting hyperplane)
 */
class InfiniteRatioDistanceException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw()
    {
      return "The requested ratio distance is infinite.";
    }
};


}  // namespace exception_classes


}  // namespace pareto_approximator


/*! @} */


#endif  // INFINITE_RATIO_DISTANCE_EXCEPTION_H
