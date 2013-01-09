/*! \file NonExistentCoordinateException.h
 *  \brief The declaration and definition of the NonExistentCoordinateException 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef NON_EXISTENT_COORDINATE_EXCEPTION_H
#define NON_EXISTENT_COORDINATE_EXCEPTION_H

#include <exception>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
namespace pareto_approximator {


//! The namespace containing all the exception classes.
namespace exception_classes {


//! Exception thrown by Point::operator[]().
/*! 
 *  An exception thrown when the requested Point coordinate does not exist. 
 */
class NonExistentCoordinateException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw()
    {
      return "The requested coordinate does not exist. (out of bounds)";
    }
};


}  // namespace exception_classes


}  // namespace pareto_approximator


/* @} */


#endif  // NON_EXISTENT_COORDINATE_EXCEPTION_H
