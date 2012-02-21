/*! \file Not2DPointsException.h
 *  \brief A file containing the declaration and definition of the 
 *         Not2DPointsException exception class.
 */


#ifndef NOT_2D_POINTS_EXCEPTION_H
#define NOT_2D_POINTS_EXCEPTION_H

#include <exception>


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Exception thrown by the Line2D(const Point&, const Point&) constructor.
/*! 
 *  An exception thrown from a Line2D constructor when 2-dimensional Point 
 *  instances were expected (to determine the line passing through them) 
 *  but at least one of the given points was not 2-dimensional.
 */
class Not2DPointsException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw() 
    {
      return "Expected 2-dimensional points.";
    }
};


}  // namespace pareto_approximator


#endif  // NOT_2D_POINTS_EXCEPTION_H
