/*! \file NotEnoughBasePointsException.h
 *  \brief The declaration and definition of the NotEnoughBasePointsException 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef NOT_ENOUGH_BASE_POINTS_EXCEPTION_H
#define NOT_ENOUGH_BASE_POINTS_EXCEPTION_H

#include <exception>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


/*! 
 *  \brief Exception thrown by BaseProblem<S>::doChord() and 
 *         BaseProblem<S>::computeConvexParetoSet().
 *
 *  Thrown when we don't have enough points to form a base and an 
 *  n-hyperplane in BaseProblem.
 */
class NotEnoughBasePointsException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw() 
    {
      return "Not enough base points to form a hyperplane.";
    }
};


}  // namespace pareto_approximator


/*! @} */


#endif  // NOT_ENOUGH_BASE_POINTS_EXCEPTION_H
