/*! \file VerticalLineException.h
 *  \brief The declaration and definition of the VerticalLineException 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef VERTICAL_LINE_EXCEPTION_H
#define VERTICAL_LINE_EXCEPTION_H

#include <exception>


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Exception thrown by Line2D::m().
/*! 
 *  An exception thrown when a Line2D instance's slope was requested but 
 *  the line is vertical (infinite slope).
 */
class VerticalLineException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw() 
    {
      return "The line is vertical (slope is infinite).";
    }
};


}  // namespace pareto_approximator


#endif  // VERTICAL_LINE_EXCEPTION_H
