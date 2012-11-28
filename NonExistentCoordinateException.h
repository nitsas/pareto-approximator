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
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


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


}  // namespace pareto_approximator


/* @} */


#endif  // NON_EXISTENT_COORDINATE_EXCEPTION_H
