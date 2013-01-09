/*! \file NotStrictlyPositivePointException.h
 *  \brief A file containing the declaration and definition of the 
 *         NotStrictlyPositivePointException exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef NOT_STRICTLY_POSITIVE_POINT_EXCEPTION_H
#define NOT_STRICTLY_POSITIVE_POINT_EXCEPTION_H

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
 *  \brief Exception thrown when we expected a strictly positive point 
 *         (i.e. all coordinates > 0.0) but the given point was not 
 *         strictly positive.
 *
 *  An exception thrown when some given Point instance (p) is not positive, 
 *  i.e. \f$ p_{i} > 0 \f$ does not hold for all i.
 *
 *  When a point we are working with is not strictly positive, concepts 
 *  like ratio distance and approximation (in its multiplicative sense) 
 *  break down.
 *  
 *  \sa Point and Facet
 */
class NotStrictlyPositivePointException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw()
    {
      return "Some given Point instance is not strictly positive.";
    }
};


}  // namespace exception_classes


}  // namespace pareto_approximator


/* @} */


#endif  // NOT_STRICTLY_POSITIVE_POINT_EXCEPTION_H
