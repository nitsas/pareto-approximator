/* VerticalLineException.h */


#ifndef VERTICAL_LINE_EXCEPTION_H
#define VERTICAL_LINE_EXCEPTION_H

#include <exception>


namespace pareto_approximator {


class VerticalLineException : public std::exception
{
  public:
    const char* what() const throw()
    {
      return "The line is vertical (slope is infinite).";
    }
};


}  // namespace pareto_approximator


#endif  // VERTICAL_LINE_EXCEPTION_H
