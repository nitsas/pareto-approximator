/* DifferentDimensionsException.cpp */


#include <sstream>

#include "DifferentDimensionsException.h"


namespace pareto_approximator {


// --- Constructors ---
// Default constructor
DifferentDimensionsException::DifferentDimensionsException()
{
  message_ = "The points have different dimensions.";
}


// More detailed constructor.
DifferentDimensionsException::DifferentDimensionsException(int p1_dim,
                                                           int p2_dim)
{
  std::stringstream ss;
  ss << "The points have different dimensions. " << p1_dim << " vs " << p2_dim;
  message_ = ss.str().c_str();
}


// --- Show message ---
const char* 
DifferentDimensionsException::what() const throw()
{
  return message_;
}


}  // namespace pareto_approximator
