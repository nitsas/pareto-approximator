/* DifferentDimensionsException.h */


#ifndef DIFFERENT_DIMENSIONS_EXCEPTION_H
#define DIFFERENT_DIMENSIONS_EXCEPTION_H

#include <exception>


namespace pareto_approximator {


class DifferentDimensionsException : public std::exception
{
  public:
    DifferentDimensionsException();
    DifferentDimensionsException(int p1_dim, int p2_dim);

    const char* what() const throw();

  private:
    const char* message_;
};


}  // namespace pareto_approximator


#endif  // DIFFERENT_DIMENSIONS_EXCEPTION_H
