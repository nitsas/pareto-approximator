/*! \file PointIsOnTheHyperplaneException.h
 *  \brief The declaration and definition of the PointIsOnTheHyperplaneException 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef POINT_IS_ON_THE_HYPERPLANE_EXCEPTION_H
#define POINT_IS_ON_THE_HYPERPLANE_EXCEPTION_H

#include <exception>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


/*! 
 *  \brief Exception thrown by Hyperplane::faceAwayFrom().
 *
 *  Thrown when we were trying to make a hyperplane face away from a 
 *  point but that point was on the hyperplane.
 */
class PointIsOnTheHyperplaneException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw() 
    {
      return "The given point is on the hyperplane. (satisfies the hyperplane equation)";
    }
};


}  // namespace pareto_approximator


/*! @} */


#endif  // POINT_IS_ON_THE_HYPERPLANE_EXCEPTION_H
