/*! \file NotPositivePointException.h
 *  \brief A file containing the declaration and definition of the 
 *         NotPositivePointException exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef NOT_POSITIVE_POINT_EXCEPTION_H
#define NOT_POSITIVE_POINT_EXCEPTION_H

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
 *  \brief Exception thrown when we expected a positive point (i.e. all 
 *         coordinates >= 0.0) but the given point was not positive.
 *
 *  An exception thrown when some given Point instance (p) is not positive, 
 *  i.e. \f$ p_{i} >= 0 \f$ does not hold for all i.
 *
 *  \sa Point and Facet
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


}  // namespace exception_classes


}  // namespace pareto_approximator


/* @} */


#endif  // NOT_POSITIVE_POINT_EXCEPTION_H
