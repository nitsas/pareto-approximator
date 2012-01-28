/* ParallelLinesException.h */


#ifndef PARALLEL_LINES_EXCEPTION_H
#define PARALLEL_LINES_EXCEPTION_H

#include <exception>


namespace pareto_approximator {


class ParallelLinesException : public std::exception
{
  public:
    const char* what() const throw()
    {
      return "The lines are parallel or the same line.";
    }
};


}  // namespace pareto_approximator


#endif  // PARALLEL_LINES_EXCEPTION_H
