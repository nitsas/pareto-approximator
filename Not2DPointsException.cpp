/* Not2DPointsException.cpp */


#include <sstream>

#include "Not2DPointsException.h"


namespace pareto_approximator {


// --- Constructors ---
// Default constructor
Not2DPointsException::Not2DPointsException()
{
  message_ = "Expected two-dimensional points.";
}


// More detailed constructor.
Not2DPointsException::Not2DPointsException(int p1_dim, int p2_dim)
{
  std::stringstream ss;
  ss << "Expected two-dimensional points. Got: " << 
        p1_dim << "D and " << p2_dim << "D";
  message_ = ss.str().c_str();
}


// --- Show message ---
const char* 
Not2DPointsException::what() const throw()
{
  return message_;
}


}  // namespace pareto_approximator
