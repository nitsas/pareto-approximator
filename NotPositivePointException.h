/*! \file NotPositivePointException.h
 *  \brief A file containing the declaration and definition of the 
 *         NotPositivePointException exception class.
 */


#ifndef NOT_POSITIVE_POINT_EXCEPTION_H
#define NOT_POSITIVE_POINT_EXCEPTION_H

#include <exception>


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Exception thrown by Point::dominates().
/*! 
 *  An exception thrown when some given Point instance (p) is not positive, 
 *  i.e. \f$ p_{i} \ge 0 \f$ does not hold for all i.
 *
 *  When a point we are working with is not positive, concepts like ratio 
 *  distance and approximation break down.
 */
class NotPositivePointException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw() 
    {
      return "Some given Point instance is not positive.";
    }
};


}  // namespace pareto_approximator


#endif  // NOT_POSITIVE_POINT_EXCEPTION_H
