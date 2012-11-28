/*! \file NullObjectException.h
 *  \brief The declaration and definition of the NullObjectException 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef PARETO_APPROXIMATOR_NULL_OBJECT_EXCEPTION_H
#define PARETO_APPROXIMATOR_NULL_OBJECT_EXCEPTION_H

#include <exception>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


/*! 
 *  \brief Exception thrown by most of Point's methods and other 
 *         functions/methods that use Point instances and cannot handle
 *         null Points.
 *
 *  An exception thrown when:
 *  - A null Point's method (except for Point::isNull()) was called on a 
 *    null Point instance (i.e. instance created using the Point::Point()
 *    constructor).
 *  - A function/other-class's-method that operates on Point instances was 
 *    called on a null Point instance and can't handle it.
 */
class NullObjectException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw()
    {
      return "Cannot operate on a null object.";
    }
};


}  // namespace pareto_approximator


/*! @} */


#endif  // PARETO_APPROXIMATOR_NULL_OBJECT_EXCEPTION_H
