/*! \file SamePointsException.h
 *  \brief The declaration and definition of the SamePointsException 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef SAME_POINTS_EXCEPTION_H
#define SAME_POINTS_EXCEPTION_H

#include <exception>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Exception thrown by some Hyperplane constructors.
/*! 
 *  An exception thrown when some of the given Point instances represent 
 *  the same point. (unique Point instances are expected when creating a 
 *  Hyperplane instance)
 */
class SamePointsException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw() 
    {
      return "Encountered two Point instances that represent the same point. (expected unique points!)";
    }
};


}  // namespace pareto_approximator


/* @} */


#endif  // SAME_POINTS_EXCEPTION_H
