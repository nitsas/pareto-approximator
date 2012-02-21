/*! \file DifferentDimensionsException.h
 *  \brief A file containing the declaration and definition of the 
 *         DifferentDimensionsException exception class.
 */


#ifndef DIFFERENT_DIMENSIONS_EXCEPTION_H
#define DIFFERENT_DIMENSIONS_EXCEPTION_H

#include <exception>


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Exception thrown by Point::operator<() and Point::ratioDistance().
/*! 
 *  An exception thrown when two Point instances were supposed to be 
 *  compared but they have different dimensions.
 */
class DifferentDimensionsException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const 
    {
      return "The points have different dimensions.";
    }
};


}  // namespace pareto_approximator


#endif  // DIFFERENT_DIMENSIONS_EXCEPTION_H
