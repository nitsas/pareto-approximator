/*! \file Point.cpp
 *  \brief A file containing the definition of the Point class.
 */


#include <sstream>

#include "Point.h"


using std::max;


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! The empty constructor. Creates an all-zero 3-dimensional Point.
Point::Point()
{
  dimension_ = 3; 
  x = y = z = 0;
}


//! An 1-dimensional Point constructor. 
/*! The resulting Point's dimensions will be doubles, not ints. */
Point::Point(int xx)
{
  dimension_ = 1;
  x = xx;
  y = z = 0;
}


//! An 1-dimensional Point constructor.
Point::Point(double xx)
{
  dimension_ = 1;
  x = xx;
  y = z = 0;
}


//! A 2-dimensional Point constructor.
/*! The resulting Point's dimensions will be doubles, not ints. */
Point::Point(int xx, int yy)
{
  dimension_ = 2;
  x = xx;
  y = yy;
  z = 0;
}


//! A 2-dimensional Point constructor.
Point::Point(double xx, double yy)
{
  dimension_ = 2;
  x = xx;
  y = yy;
  z = 0;
}


//! A 3-dimensional Point constructor.
/*! The resulting Point's dimensions will be doubles, not ints. */
Point::Point(int xx, int yy, int zz)
{
  dimension_ = 3;
  x = xx;
  y = yy;
  z = zz;
}


//! A 3-dimensional Point constructor.
Point::Point(double xx, double yy, double zz)
{
  dimension_ = 3;
  x = xx;
  y = yy;
  z = zz;
}


//! A simple (and empty) Destructor.
Point::~Point() {}


//! Get a Point instance's dimension. (1D, 2D or 3D point)
/*! 
 *  \return 1, 2 or 3 for 1-dimensional, 2-dimensional or 3-dimensional 
 *          Point instances respectively.
 *  
 *  \sa Point and bool Point::dimension(int)
 */
int 
Point::dimension() const
{
  return dimension_;
}


//! Set the point's dimension. (1D, 2D or 3D point)
/*! 
 *  \param dim The dimension we want to change the Point instance to.
 *  \return true if everything went ok, false otherwise.
 *          (we get false only if dim was not 1, 2 or 3)
 *  
 *  Set higher dimension values to 0. (e.g. z = 0 if we want to make the 
 *  Point instance 2-dimensional or 1-dimensional)
 *  
 *  \sa Point and int Point::dimension() const
 */
bool 
Point::dimension(int dimension)
{
  switch (dimension) {
    case 1 :
      y = 0;
      // Fallthrough
      
    case 2 :
      z = 0;
      // Fallthrough
      
    case 3 :
      dimension_ = dimension;
      return true;
      
    default :
      return false;
  }
}


//! The Point equality operator.
/*! 
 *  \param p A Point instance we want to compare with the current instance.
 *  \return true if the Points are equal, false otherwise.
 *  
 *  Checks if both Point instances are of the same dimension and have 
 *  equal x, y and z dimensions. Returns true if all the above hold, false 
 *  otherwise. Undefined dimensions will not be checked (e.g. z for 
 *  2-dimensional or 1-dimensional Points).
 *  
 *  \sa Point, operator!=(), operator<(), 
 *      std::ostream& operator<<(std::ostream&, Point&) and 
 *      std::istream& operator>>(std::istream&, Point&)
 */
bool 
Point::operator== (const Point& p) const
{
  if (dimension_ != p.dimension())
    return false;

  switch (dimension_) {
    case 1 :
      return (x == p.x);

    case 2 :
      return (x == p.x && y == p.y);

    case 3 :
      // Fallthrough

    default :
      return (x == p.x && y == p.y && z == p.z);
  }
}


//! The Point inequality operator.
/*! 
 *  \param p A Point instance we want to compare with the current instance.
 *  \return true if the Points are not equal, false otherwise.
 *  
 *  Checks if the two Point instances are of different dimensions or differ
 *  in at least one of the x, y and z dimensions. Returns true if at least 
 *  one of the above holds, false otherwise. Undefined dimensions will not 
 *  be checked (e.g. z for 2-dimensional or 1-dimensional Points).
 *  
 *  \sa Point, operator==(), operator<(), 
 *      std::ostream& operator<<(std::ostream&, Point&) and 
 *      std::istream& operator>>(std::istream&, Point&)
 */
bool 
Point::operator!= (const Point& p) const
{
  if (dimension_ != p.dimension())
    return true;

  switch (dimension_) {
    case 1 :
      return (x != p.x);

    case 2 :
      return (x != p.x || y != p.y);

    case 3 :
      // Fallthrough

    default :
      return (x != p.x || y != p.y || z != p.z);
  }
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
 *  - May throw a DifferentDimensionException exception if the two Point 
 *    instances are of different dimensions (can't be compared).
 *  
 *  \sa Point, operator==(), operator!=(), 
 *      std::ostream& operator<<(std::ostream&, Point&) and 
 *      std::istream& operator>>(std::istream&, Point&)
 */
bool 
Point::operator< (const Point& p) const 
{
  if (dimension_ != p.dimension())
    throw DifferentDimensionsException();
  // else
  switch (dimension_) {
    case 1 :
      return x < p.x;

    case 2 :
      return (x < p.x) || (x == p.x && y < p.y);

    case 3 :
      return (x < p.x) || (x == p.x && ( y < p.y || (y == p.y && z < p.z) ) );

    default :
      return (x < p.x) || (x == p.x && ( y < p.y || (y == p.y && z < p.z) ) );
  }
}


//! The Point output stream operator. A friend of the Point class.
/*! 
 *  The Point instance's dimensions will be output inside parentheses. 
 *  Examples:
 *  - (1.0, 4.27, 0.883)
 *  - (3.0)
 *  - (5, 1.99204e+09)
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
 *      Point::operator<() and 
 *      std::istream& operator>>(std::istream&, Point&)
 */
std::ostream& 
operator<< (std::ostream& ostr, const Point& p)
{
  return ostr << p.str();
}


//! The Point input stream operator. A friend of the Point class.
/*! 
 *  The accepted format is similar to the one operator<<() uses for output.
 *  
 *  \sa Point, Point::operator==(), Point::operator!=(), 
 *      Point::operator<() and 
 *      std::ostream& operator<<(std::ostream&, Point&)
 */
std::istream& 
operator>> (std::istream& istr, Point& p)
{
  char c;
  istr >> c;          // skip '('
  istr >> p.x;
  istr >> c;
  if (c == ')') {
    p.dimension(1);   // 1-dimensional point
    return istr;
  }
  // else             // skip ','
  istr >> p.y;
  istr >> c;
  if (c == ')') {
    p.dimension(2);   // 2-dimensional point
    return istr;
  }
  // else             // skip ','
  istr >> p.z;
  p.dimension(3);     // 3-dimensional point
  istr >> c;
  return istr;
}


//! Return the Point instance as a string.
/*! 
 *  Makes a string with the Point instance's dimensions inside parentheses. 
 *  Examples:
 *  - (1.0, 4.27, 0.883)
 *  - (3.0)
 *  - (5, 1.99204e+09)
 *  - etc
 *
 *  /sa Point and std::ostream& operator<<(std::ostream&, const Point&)
 */
std::string 
Point::str() const
{
  std::stringstream ss;

  switch (dimension_) {
    case 1 :
      ss << "(" << x << ")";
      break;

    case 2 :
      ss << "(" << x << ", " << y << ")";
      break;
      
    case 3 :
      // Fallthrough

    default :
      ss << "(" << x << ", " << y << ", " << z << ")";
      break;
  }

  return ss.str();
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
 *  \sa Point and Line2D::ratioDistance()
 */
double 
Point::ratioDistance(const Point& q) const 
{
  if (q.dimension() != dimension_)
    throw DifferentDimensionsException();
  // else
  switch (dimension_) {
    case 1 :
      return max( (q.x - x)/x, 0.0 );
      break;

    case 2 :
      return max( (q.x - x)/x, max( (q.y - y)/y, 0.0 ) );
      break;
      
    case 3 :
      return max( (q.x - x)/x, max( (q.y - y)/y, max( (q.z - z)/z, 0.0) ) );
      break;

    default :
      return max( (q.x - x)/x, max( (q.y - y)/y, max( (q.z - z)/z, 0.0) ) );
  }
}


}  // namespace pareto_approximator
