/*! \file NegativeApproximationRatioException.h
 *  \brief The declaration and definition of the 
 *         NegativeApproximationRatioException exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef NEGATIVE_APPROXIMATION_RATIO_EXCEPTION_H
#define NEGATIVE_APPROXIMATION_RATIO_EXCEPTION_H

#include <exception>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Exception thrown by Point::dominates().
/*! 
 *  An exception thrown when the given approximation ratio is negative.
 *  (approximation ratios must always be positive!)
 */
class NegativeApproximationRatioException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw()
    {
      return "The given approximation ratio is negative.";
    }
};


}  // namespace pareto_approximator


/* @} */


#endif  // NEGATIVE_APPROXIMATION_RATIO_EXCEPTION_H
