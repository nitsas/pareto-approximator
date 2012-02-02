/* Point.cpp */


#include "Point.h"


using std::max;


namespace pareto_approximator {


// --- Constructors ---
Point::Point()
{
  dimension_ = 3; 
  x = y = z = 0;
}


// 1-dimensional point
Point::Point(int xx)
{
  dimension_ = 1;
  x = xx;
  y = z = 0;
}


Point::Point(double xx)
{
  dimension_ = 1;
  x = xx;
  y = z = 0;
}


// 2-dimensional point
Point::Point(int xx, int yy)
{
  dimension_ = 2;
  x = xx;
  y = yy;
  z = 0;
}


Point::Point(double xx, double yy)
{
  dimension_ = 2;
  x = xx;
  y = yy;
  z = 0;
}


// 3-dimensional point
Point::Point(int xx, int yy, int zz)
{
  dimension_ = 3;
  x = xx;
  y = yy;
  z = zz;
}


Point::Point(double xx, double yy, double zz)
{
  dimension_ = 3;
  x = xx;
  y = yy;
  z = zz;
}


// --- Destructor ---
Point::~Point() {}


// --- dimension_ setter and getter ---
int 
Point::dimension() const
{
  return dimension_;
}


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


// --- Operators ---
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


// Convert Point to a string with format: "(%f)", "(%f, %f)" or "(%f, %f, %f)".
std::ostream& 
operator<< (std::ostream& ostr, const Point& p)
{
  switch (p.dimension()) {
    case 1 :
      ostr << "(" << p.x << ")";
      break;

    case 2 :
      ostr << "(" << p.x << ", " << p.y << ")";
      break;
      
    case 3 :
      // Fallthrough

    default :
      ostr << "(" << p.x << ", " << p.y << ", " << p.z << ")";
      break;
  }
  return ostr;
}


// Read input Point with format: "(%f)", "(%f, %f)" or "(%f, %f, %f)".
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


// --- Methods ---
double 
Point::ratioDistance(const Point& q) const throw(DifferentDimensionsException)
{
  if (q.dimension() != dimension_)
    throw DifferentDimensionsException(dimension_, q.dimension());
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
