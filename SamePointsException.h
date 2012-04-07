/*! \file SamePointsException.h
 *  \brief The declaration and definition of the SamePointsException 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef SAME_POINTS_EXCEPTION_H
#define SAME_POINTS_EXCEPTION_H

#include <exception>


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Exception thrown by the Line2D(const Point&, const Point&) constructor.
/*! 
 *  An exception thrown when two different Point instances were expected 
 *  but they turn out to be the same point.
 */
class SamePointsException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw() 
    {
      return "The points given are the same point.";
    }
};


}  // namespace pareto_approximator


#endif  // SAME_POINTS_EXCEPTION_H
