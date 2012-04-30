/*! \file Not2DHyperplanesException.h
 *  \brief The declaration and definition of the Not2DHyperplanesException 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef NOT_2D_HYPERPLANES_EXCEPTION_H
#define NOT_2D_HYPERPLANES_EXCEPTION_H

#include <exception>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Exception thrown by Hyperplane::intersection().
/*! 
 *  An exception thrown by Hyperplane::intersection() when the given 
 *  hyperplanes were not 2-hyperplanes (lines).
 */
class Not2DHyperplanesException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw() 
    {
      return "Expected 2-hyperplanes (lines).";
    }
};


}  // namespace pareto_approximator


/* @} */


#endif  // NOT_2D_HYPERPLANES_EXCEPTION_H
