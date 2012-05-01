/*! \file Point.cpp
 *  \brief The declaration of Point. (a simple point class)
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <assert.h>
#include <string>
#include <sstream>
#include <vector>

#include "Point.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! The empty constructor. Creates a zero-dimensional Point.
Point::Point() { }


//! An 1-dimensional Point constructor.
Point::Point(int x)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
}


//! An 1-dimensional Point constructor.
Point::Point(double x)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
}


//! A 2-dimensional Point constructor.
Point::Point(int x, int y)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
}


//! A 2-dimensional Point constructor.
Point::Point(double x, double y)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
}


//! A 3-dimensional Point constructor.
Point::Point(int x, int y, int z)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
  coordinates_.push_back(z);
}


//! A 3-dimensional Point constructor.
Point::Point(double x, double y, double z)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
  coordinates_.push_back(z);
}


//! An n-dimensional Point constructor.
/*! 
 *  \param first Iterator to the initial position in a std::vector<double>.
 *  \param last Iterator to the final position in a std::vector<double>.
 *  
 *  The range used is [first, last), which includes all the elements between 
 *  first and last, including the element pointed by first but not the 
 *  element pointed by last.
 *  
 *  The resulting point's coordinates will be doubles, not ints. 
 */
Point::Point(std::vector<int>::iterator first, 
             std::vector<int>::iterator last)
{
  assert(coordinates_.size() == 0);
  coordinates_.assign(first, last);
}


//! An n-dimensional Point constructor.
/*! 
 *  \param first Iterator to the initial position in a std::vector<double>.
 *  \param last Iterator to the final position in a std::vector<double>.
 *  
 *  The range used is [first, last), which includes all the elements between 
 *  first and last, including the element pointed by first but not the 
 *  element pointed by last.
 */
Point::Point(std::vector<double>::iterator first, 
             std::vector<double>::iterator last)
{
  assert(coordinates_.size() == 0);
  coordinates_.assign(first, last);
}


//! An n-dimensional Point constructor.
/*! 
 *  \param first Iterator (pointer) to the initial position in an array of 
 *               double.
 *  \param last Iterator (pointer) to the final position in an array of 
 *              double. (the position just beyond the last element we want)
 *  
 *  The range used is [first, last), which includes all the elements between 
 *  first and last, including the element pointed by first but not the 
 *  element pointed by last.
 *  
 *  The resulting point's coordinates will be doubles, not ints. 
 */
Point::Point(int* first, int* last)
{
  assert(coordinates_.size() == 0);
  coordinates_.assign(first, last);
}


//! An n-dimensional Point constructor.
/*! 
 *  \param first Iterator (pointer) to the initial position in an array of 
 *               double.
 *  \param last Iterator (pointer) to the final position in an array of 
 *              double. (the position just beyond the last element we want)
 *  
 *  The range used is [first, last), which includes all the elements between 
 *  first and last, including the element pointed by first but not the 
 *  element pointed by last.
 */
Point::Point(double* first, double* last)
{
  assert(coordinates_.size() == 0);
  coordinates_.assign(first, last);
}


//! A simple (and empty) destructor.
Point::~Point() { }


//! Get a Point instance's dimension. (1D, 2D or 3D point)
/*! 
 *  \return 1, 2 or 3 for 1-dimensional, 2-dimensional or 3-dimensional 
 *          Point instances respectively.
 *  
 *  \sa Point and bool Point::dimension(int)
 */
unsigned int 
Point::dimension() const
{
  return coordinates_.size();
}


//! Set the point's dimension. (1D, 2D, 3D, etc point)
/*! 
 *  \param dimension The dimension we want to change the Point instance to.
 *  
 *  If "dimension" is smaller than the current Point dimension only the 
 *  first "dimension" coordinates will be kept, the rest being dropped.
 *  
 *  \sa Point and int Point::dimension() const
 */
void 
Point::dimension(unsigned int dimension)
{
  // Resize the coordinates_ vector to "dimension" elements. 
  // Initialize any newly inserted elements to 0.0.
  coordinates_.resize(dimension, 0.0);
}


//! The Point access coordinate operator.
/*! 
 *  \param pos The position (coordinate) to access. 
 *             (0 <= pos < dimension())
 *  \return The Point's "pos" coordinate.
 *  
 *  Possible exceptions:
 *  - May throw a NonExistentCoordinateException exception if the 
 *    requested coordinate does not exist. ("pos" is greater than or 
 *    equal to dimension())
 *  
 *  \sa Point, dimension(), operator==(), operator!=(), operator<(), 
 *      std::ostream& operator<<(std::ostream&, Point&) and 
 *      std::istream& operator>>(std::istream&, Point&)
 */
double 
Point::operator[] (unsigned int pos) const
{
  if (pos >= dimension())
    throw NonExistentCoordinateException();
  else
    return coordinates_[pos];
}


//! The Point equality operator.
/*! 
 *  \param p A Point instance we want to compare with the current instance.
 *  \return true if the Points are equal, false otherwise.
 *  
 *  Checks if both Point instances are of the same dimension and have 
 *  equal coordinates. Returns true if all the above hold, false 
 *  otherwise. 
 *  
 *  \sa Point, operator!=(), operator<(), operator[](), 
 *      std::ostream& operator<<(std::ostream&, Point&) and 
 *      std::istream& operator>>(std::istream&, Point&)
 */
bool 
Point::operator== (const Point& p) const
{
  if (dimension() != p.dimension())
    return false;

  for (unsigned int i=0; i<dimension(); i++) 
    if (coordinates_[i] != p[i]) 
      return false;

  return true;
}


//! The Point inequality operator.
/*! 
 *  \param p A Point instance we want to compare with the current instance.
 *  \return true if the Points are not equal, false otherwise.
 *  
 *  Checks if the two Point instances are of different dimensions or differ
 *  in at least one of their coordinates. Returns true if at least one of 
 *  the above holds, false otherwise. 
 *  
 *  \sa Point, operator==(), operator<(), operator[](), 
 *      std::ostream& operator<<(std::ostream&, Point&) and 
 *      std::istream& operator>>(std::istream&, Point&)
 */
bool 
Point::operator!= (const Point& p) const
{
  if (dimension() != p.dimension())
    return true;

  for (unsigned int i=0; i<dimension(); i++) 
    if (coordinates_[i] != p[i]) 
      return true;

  return false;
}


//! The Point less-than operator.
/*! 
 *  \param p A Point instance we want to compare with the current instance.
 *  \return true if the Point instance on the left of the operator is 
 *          lexicographically smaller than p, false otherwise.
 *  
 *  Compare the current Point instance with the given instance (p) 
 *  lexicographically. A point p1 is lexicographically smaller than a 
 *  point p2 if   \f$ p1_{q} < p2_{q} \f$   where 
 *  \f$ q = \min{k : p1_{k} \ne p2_k} \f$.
 *  
 *  Possible exceptions:
 *  - May throw a DifferentDimensionsException exception if the two Point 
 *    instances are of different dimensions (can't be compared).
 *  
 *  \sa Point, operator==(), operator!=(), operator[](), 
 *      std::ostream& operator<<(std::ostream&, Point&) and 
 *      std::istream& operator>>(std::istream&, Point&)
 */
bool 
Point::operator< (const Point& p) const 
{
  if (dimension() != p.dimension()) 
    throw DifferentDimensionsException();
  // else
  for (unsigned int i=0; i<dimension(); i++) {
    if (coordinates_[i] < p[i])
      return true;
    else if (coordinates_[i] > p[i])
      return false;
    else 
      continue;
  }

  return false;
}


//! The Point output stream operator. A friend of the Point class.
/*! 
 *  The Point instance's dimensions will be output inside parentheses. 
 *  Examples:
 *  - (1.0, 4.27, 0.883)
 *  - (3.0)
 *  - (5, 1.99204e+09)
 *  - ()      <-- zero-dimensional point
 *  - etc
 *
 *  operator<<() uses the Point instance's Point::str() method to create 
 *  the output.
 *
 *  Usage:
 *  - std::cout << Point(4, 3, -8);
 *  - std::cout << "some text " << Point(2.7, -2.7) << std::endl;
 *
 *  \sa Point, Point::str(), Point::operator==(), Point::operator!=(), 
 *      Point::operator<(), Point::operator[]() and 
 *      std::istream& operator>>(std::istream&, Point&)
 */
std::ostream& 
operator<< (std::ostream& out, const Point& p)
{
  return out << p.str();
}


//! The Point input stream operator. A friend of the Point class.
/*! 
 *  The accepted format is similar to the one operator<<() uses for output.
 *  Zero-dimensional points of the form "()" are not accepted.
 *  
 *  \sa Point, Point::str(), Point::operator==(), Point::operator!=(), 
 *      Point::operator<(), Point::operator[]() and 
 *      std::ostream& operator<<(std::ostream&, Point&)
 */
std::istream& 
operator>> (std::istream& istr, Point& p)
{
  char c;
  double d;
  p.coordinates_.clear();
  istr >> c;          // skip '('
  while (c != ')') {
    istr >> d;
    p.coordinates_.push_back(d);
    istr >> c;        // get ',' or ')'
  }

  return istr;
}


//! Return the Point instance as a string.
/*! 
 *  Makes a string with the Point instance's dimensions inside parentheses. 
 *  Examples:
 *  - (1.0, 4.27, 0.883)
 *  - (3.0)
 *  - (5, 1.99204e+09)
 *  - ()      <-- zero-dimensional point
 *  - etc
 *
 *  \sa Point, Point::operator==(), Point::operator!=(), 
 *      Point::operator<(), Point::operator[]() and 
 *      std::ostream& operator<<(std::ostream&, const Point&)
 */
std::string 
Point::str() const
{
  if (dimension() == 0)
    return "()";
  else {
    std::stringstream ss;

    ss << "(" << coordinates_[0];
    for (unsigned int i=1; i<dimension(); i++) 
      ss << ", " << coordinates_[i];
    ss << ")";

    return ss.str();
  }
}


//! Return the ratio distance from the current to the given Point instance.
/*! 
 *  We define the ratio distance from point p to q as:
 *  \f$ RD(p, q) = \max{ x(q)/x(p) - 1, y(q)/y(p) - 1, 0 } \f$.
 *  Intuitively, it is the minimum value of \f$ \epsilon \ge 0 \f$ such 
 *  that  \f$ q \epsilon -covers p \f$.
 *  
 *  Possible exceptions:
 *  - May throw a DifferentDimensionsException exception if the two 
 *    Point instances are of different dimensions (can't be compared).
 *
 *  \sa Point, Point::dominates() and Line2D::ratioDistance()
 */
double 
Point::ratioDistance(const Point& q) const 
{
  if (dimension() != q.dimension())
    throw DifferentDimensionsException();
  // else
  double max = 0.0;
  for (unsigned int i=0; i<dimension(); i++) {
    double r = (q[i] - coordinates_[i]) / coordinates_[i];
    if (r > max) 
      max = r;
  }

  return max;
}


//! Check if the current point (p) eps-covers the given point (q).
/*! 
 *  \param q A Point instance with \f$ q_{i} \ge 0 \f$ for all i.
 *  \param eps An approximation factor. (default 0.0)
 *  \return true if p eps-covers q, false otherwise.
 *  
 *  Note that both p and q must be greater than zero (dominated by 0),
 *  that is both \f$ p_{i} \ge 0 \f$ and \f$ q_{i} \ge 0 \f$ must hold
 *  for all i.
 *  
 *  We say that p \f$ \epsilon \f$-covers q (\f$\epsilon \ge 0 \f$) if 
 *  \f$ p_{i} \le (1 + \epsilon) q_{i} \f$, for all i. Both p and 
 *  q must be of the same dimension.
 *  
 *  If eps=0.0 the method simply checks whether or not p dominates 
 *  q and that is how it got its name.
 *  
 *  Possible exceptions:
 *  - May throw a NotPositivePointException if either p or q is not 
 *    greater than 0 (dominated by 0).
 *  - May throw a NegativeApproximationRatioException if \f$ eps < 0 \f$.
 *  - May throw a DifferentDimensionsException if p and q are of different 
 *    dimensions.
 *  
 *  \sa Point and Point::ratioDistance()
 */
bool 
Point::dominates(const Point& q, double eps) const
{
  if (dimension() != q.dimension())
    throw DifferentDimensionsException();
  if (eps < 0.0)
    throw NegativeApproximationRatioException();
  for (unsigned int i=0; i<dimension(); i++) 
    if (coordinates_[i] < 0.0 || q[i] < 0.0) 
      throw NotPositivePointException();

  double r = 1+eps;
  for (unsigned int i=0; i<dimension(); i++) 
    if (coordinates_[i] > r * q[i])
      return false;

  return true;
}


}  // namespace pareto_approximator


/* @} */
