/* SamePointsException.h */


#ifndef SAME_POINTS_EXCEPTION_H
#define SAME_POINTS_EXCEPTION_H

#include <exception>


namespace pareto_approximator {


class SamePointsException : public std::exception
{
  public:
    const char* what() const throw()
    {
      return "The points given are the same point.";
    }
};


}  // namespace pareto_approximator


#endif  // SAME_POINTS_EXCEPTION_H
