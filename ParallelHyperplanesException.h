/*! \file ParallelHyperplanesException.h
 *  \brief The declaration and definition of the ParallelHyperplanesException 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef PARALLEL_HYPERPLANES_EXCEPTION_H
#define PARALLEL_HYPERPLANES_EXCEPTION_H

#include <exception>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Exception thrown by Hyperplane::intersection().
/*! 
 *  An exception thrown when the intersection of two Hyperplane instances was 
 *  requested but the two hyperplanes are parallel (they don't intersect) or 
 *  the same hyperplane (infinite intersection points).
 */
class ParallelHyperplanesException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw() 
    {
      return "The hyperplanes are parallel or the same hyperplane.";
    }
};


}  // namespace pareto_approximator


/* @} */


#endif  // PARALLEL_HYPERPLANES_EXCEPTION_H
