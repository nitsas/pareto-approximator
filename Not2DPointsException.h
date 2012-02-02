/* Not2DPointsException.h */


#ifndef NOT_2D_POINTS_EXCEPTION_H
#define NOT_2D_POINTS_EXCEPTION_H

#include <exception>


namespace pareto_approximator {


class Not2DPointsException : public std::exception
{
  public:
    Not2DPointsException();
    Not2DPointsException(int p1_dim, int p2_dim);

    const char* what() const throw();

  private:
    const char* message_;
};


}  // namespace pareto_approximator


#endif  // NOT_2D_POINTS_EXCEPTION_H
