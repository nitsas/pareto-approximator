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
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


/*! 
 *  \brief Exception thrown when we expected a strictly positive point 
 *         (i.e. all coordinates > 0.0) but the given point was not 
 *         strictly positive.by Point::dominates() and Point::ratioDistance()
 *
 *  An exception thrown when some given Point instance (p) is not positive, 
 *  i.e. \f$ p_{i} > 0 \f$ does not hold for all i.
 *
 *  When a point we are working with is not strictly positive, concepts 
 *  like ratio distance and approximation break down.
 *  
 *  Thrown by:
 *  - Point::dominates()
 *  - Point::ratioDistance()
 *  - PointAndSolution::dominates()
 *  - Facet::ratioDistance()
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


}  // namespace pareto_approximator


/* @} */


#endif  // NOT_STRICTLY_POSITIVE_POINT_EXCEPTION_H
