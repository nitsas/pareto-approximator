/* Line2D.cpp */


#include "Line2D.h"


using std::max;


namespace pareto_approximator {


// --- Constructors ---
// Construct default line   y = x
Line2D::Line2D() : m_(1), b_(0), isVertical_(false) { }


Line2D::Line2D(int m, int b) : m_(m), b_(b), isVertical_(false) { }


Line2D::Line2D(double m, double b) : m_(m), b_(b), isVertical_(false) { }


Line2D::Line2D(const Point& p1, const Point& p2) throw(Not2DPointsException, 
                                                       SamePointsException)
{
  if (p1 == p2)
    throw SamePointsException();
  if (p1.dimension() != 2 || p2.dimension() != 2)
    throw Not2DPointsException(p1.dimension(), p2.dimension());
  // else
  if (p1.x != p2.x) {
    m_ = (p2.y - p1.y) / (p2.x - p1.x);
    b_ = p1.y - m_ * p1.x;
    isVertical_ = false;
  }
  else {
    // the line passing through them is vertical
    m_ = 1;               // will not be used
    b_ = -p1.x;           // so that   0 = x + b   holds
    isVertical_ = true;
  }
}


// Construct a vertical line   x = c. (b_ = -c so that   0 = x + b   holds)
Line2D::Line2D(int c) : m_(1), b_(-c), isVertical_(true) { }


// Construct a vertical line   x = c. (b_ = -c so that   0 = x + b   holds)
Line2D::Line2D(double c) : m_(1), b_(-c), isVertical_(true) { }


// --- Attribute getters ---
double 
Line2D::m() const throw(VerticalLineException)
{
  if (isVertical_)
    throw VerticalLineException();
  // else
  return m_;
}


double 
Line2D::b() const
{
  return b_;
}


bool 
Line2D::isVertical() const
{
  return isVertical_;
}


// --- Operators ---
bool 
Line2D::operator== (const Line2D& line) const
{
  if (isVertical_) {
    if ( line.isVertical() && b_ == line.b() )
      return true;
    else
      return false;
  }
  else {
    if ( !line.isVertical() && m_ == line.m() && b_ == line.b() )
      return true;
    else
      return false;
  }
}


bool 
Line2D::operator!= (const Line2D& line) const
{
  if (isVertical_) {
    if ( !line.isVertical() || b_ != line.b() )
      return true;
    else
      return false;
  }
  else {
    if ( line.isVertical() || m_ != line.m() || b_ != line.b() )
      return true;
    else
      return false;
  }
}


// Output stream operator
std::ostream& 
operator<< (std::ostream& ostr, const Line2D& line)
{
  if (line.isVertical()) 
    ostr << "( x = " << -line.b() << " )";
  else {
    ostr << "( y = " << line.m() << " x ";
    if (line.b() >= 0.0)
      ostr << "+ " << line.b() << " )";
    else 
      ostr << "- " << -line.b() << " )";
  }
  return ostr;
}


// --- Methods ---
Point 
Line2D::intersection(const Line2D& line) const throw(ParallelLinesException)
{
  // If the lines are parallel:
  if ( (isVertical_ && line.isVertical()) || m_ == line.m() )
    throw ParallelLinesException();
  // else
  if (isVertical_) 
    // "this" is vertical
    return Point(-b_, line.m() * (-b_) + line.b());
  else if (line.isVertical()) 
    // "line" is vertical
    return Point(-line.b(), m_ * (-line.b()) + b_);
  else 
    // neither line is vertical
    return Point( (b_ - line.b()) / (line.m() - m_), 
                  (b_ - line.b()) / (line.m() - m_) * m_ + b_ );
}


double 
Line2D::ratioDistance(const Point& p) const
{
  if (isVertical_) 
    return max( (-b_ - p.x) / p.x, 0.0 );
  else
    return max( (p.y - m_ * p.x - b_) / (m_ * p.x - p.y), 0.0 );
}


Line2D 
Line2D::parallelThrough(const Point& p) const
{
  if (isVertical_)
    return Line2D(p.x);
  else
    return Line2D(m_, p.y - m_ * p.x);
}


}  // namespace pareto_approximator
