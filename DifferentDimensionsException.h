/*! \file DifferentDimensionsException.h
 *  \brief The declaration and definition of the DifferentDimensionsException 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef DIFFERENT_DIMENSIONS_EXCEPTION_H
#define DIFFERENT_DIMENSIONS_EXCEPTION_H

#include <exception>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Exception thrown by Point::operator<() and Point::ratioDistance().
/*! 
 *  An exception thrown when two Point instances or two Hyperplane instances 
 *  were supposed to be compared but they have different dimensions.
 *  (they belong to spaces of different dimensions)
 */
class DifferentDimensionsException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw() 
    {
      return "The instances have different dimensions.";
    }
};


}  // namespace pareto_approximator


/*! @} */


#endif  // DIFFERENT_DIMENSIONS_EXCEPTION_H
