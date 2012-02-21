/*! \file ParallelLinesException.h
 *  \brief A file containing the declaration and definition of the 
 *         ParallelLinesException exception class.
 */


#ifndef PARALLEL_LINES_EXCEPTION_H
#define PARALLEL_LINES_EXCEPTION_H

#include <exception>


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Exception thrown by Line2D::intersection().
/*! 
 *  An exception thrown when the intersection of two Line2D instances was 
 *  requested but the two lines are parallel (they don't intersect) or the 
 *  same line (infinite intersection points).
 */
class ParallelLinesException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const 
    {
      return "The lines are parallel or the same line.";
    }
};


}  // namespace pareto_approximator


#endif  // PARALLEL_LINES_EXCEPTION_H
