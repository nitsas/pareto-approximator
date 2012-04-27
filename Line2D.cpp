/*! \file Line2D.cpp
 *  \brief The definition of the Line2D class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <sstream>

#include "Line2D.h"


using std::max;


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! The empty constructor. Creates line \f$ y = x \f$.
Line2D::Line2D() : m_(1), b_(0), isVertical_(false) { }


//! A non-vertical line constructor. Creates line \f$ y = m x + b \f$.
Line2D::Line2D(double m, double b) : m_(m), b_(b), isVertical_(false) { }


//! Creates the line that passes through points p1 and p2 (might be 
//! vertical).
/*! 
 *  \param p1 A 2-dimensional Point instance.
 *  \param p2 A 2-dimensional Point instance.
 *  
 *  Possible exceptions:
 *  - May throw a SamePointsException exception if both the given Point 
 *    instances represent the same point.
 *  - May throw a Not2DPointsException exception if any of the given 
 *    Point instances doesn't represent a 2-dimensional point.
 */
Line2D::Line2D(const Point& p1, const Point& p2) 
{
  if (p1 == p2)
    throw SamePointsException();
  if (p1.dimension() != 2 || p2.dimension() != 2)
    throw Not2DPointsException();
  // else
  if (p1[0] != p2[0]) {
    m_ = (p2[1] - p1[1]) / (p2[0] - p1[0]);
    b_ = p1[1] - m_ * p1[0];
    isVertical_ = false;
  }
  else {
    // the line passing through them is vertical
    m_ = 1;                // will not be used
    b_ = -p1[0];           // so that   0 = x + b   holds
    isVertical_ = true;
  }
}


//! A vertical line constructor. Creates line \f$ x = c \f$.
/*! 
 *  \param c The resulting line's equation should be \f$ x = c \f$.
 *  
 *  Will set b_ to -c so that \f$ 0 = x + b \f$ holds.
 */
Line2D::Line2D(double c) : m_(1), b_(-c), isVertical_(true) { }


//! A simple (and empty) destructor.
Line2D::~Line2D() { }


//! Return the Line2D instance's slope (m_).
/*! 
 *  Possible exceptions:
 *  - May throw a VerticalLineException exception if the line is vertical.
 */
double 
Line2D::m() const 
{
  if (isVertical_)
    throw VerticalLineException();
  // else
  return m_;
}


//! Return the Line2D instance's y-intercept (b_).
double 
Line2D::b() const
{
  return b_;
}


//! Return true if the Line2D instance represents a vertical line, 
//! false otherwise.
bool 
Line2D::isVertical() const
{
  return isVertical_;
}


//! The Line2D equality operator.
/*! 
 *  \param line A Line2D instance we want to compare with the current
 *              instance.
 *  \return true if the two instances represent the same line, false 
 *          otherwise.
 *  
 *  \sa Line2D, operator!=() and operator<<()
 */
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


//! The Line2D inequality operator.
/*! 
 *  \param line A Line2D instance we want to compare with the current
 *              instance.
 *  \return true if the two instances represent different lines, false 
 *          otherwise.
 *  
 *  \sa Line2D, operator==() and operator<<()
 */
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


//! The Line2D output stream operator. A friend of the Line2D class.
/*! 
 *  The Line2D instance's equation will be output inside parentheses.
 *  Examples:
 *  - ( y = 5.3 x + 2.7 )
 *  - ( y = 0.0 x - 15.0 )
 *  - ( x = 9 )
 *  - etc
 *  
 *  operator<<() uses the Line2D instance's Line2D::str() method to 
 *  create the output.
 *
 *  Usage:
 *  - std::cout << Line2D(5.3, 2.7);
 *  - std::cout << Line2D(-90.5);
 *  - std::cout << "some text " << Line2D(-2.8, 1.5) << endl;
 *  
 *  \sa Line2D, str(), operator==(), operator!=(), m(), b() and
 *      isVertical()
 */
std::ostream& 
operator<< (std::ostream& ostr, const Line2D& line)
{
  return ostr << line.str();
}


//! Return the Line2D instance as a string. (its equation in parentheses)
/*! 
 *  Makes a string with the Line2D instance's equation in parentheses.
 *  Examples:
 *  - ( y = 5.3 x + 2.7 )
 *  - ( y = 0.0 x - 15.0 )
 *  - ( x = 9 )
 *  - etc
 *  
 *  \sa Line2D and std::ostream& operator<<(std::ostream&, const Line2D&)
 */
std::string 
Line2D::str() const
{
  std::stringstream ss;

  if (isVertical_) 
    ss << "( x = " << -b_ << " )";
  else {
    ss << "( y = " << m_ << " x ";
    if (b_ >= 0.0)
      ss << "+ " << b_ << " )";
    else 
      ss << "- " << -b_ << " )";
  }
  return ss.str();
}


//! Find the Point where two Line2D instances intersect.
/*! 
 *  \param line A Line2D instance.
 *  \return A Point instance representing the point where the given
 *          line intersects with the current one.
 *
 *  Possible exceptions:
 *  - May throw a ParallelLinesException exception if the two Line2D 
 *    instances represent parallel lines or the same line.
 *  
 *  \sa Line2D, m(), b(), isVertical(), ratioDistance(), 
 *      parallelThrough() and Point.
 */
Point 
Line2D::intersection(const Line2D& line) const
{
  if (isVertical_ && line.isVertical())
    // both lines are vertical (parallel lines)
    throw ParallelLinesException();
  // else
  if (isVertical_) 
    // "this" is vertical
    return Point(-b_, line.m() * (-b_) + line.b());
  else if (line.isVertical()) 
    // "line" is vertical
    return Point(-line.b(), m_ * (-line.b()) + b_);
  else { 
    // neither line is vertical
    if (m_ != line.m())
      // the lines intersect
      return Point( (b_ - line.b()) / (line.m() - m_), 
                    (b_ - line.b()) / (line.m() - m_) * m_ + b_ );
    else
      // the lines are parallel
      throw ParallelLinesException();
  }
}


//! Compute the ratio distance from the given Point instance to the line.
/*! 
 *  \param p A Point instance.
 *  \return The ratio distance from the given Point instance to the line.
 *
 *  The ratio distance from a point p to a line L is defined as:
 *  \f$ RD(p, L) = \min_{q \in L} RD(p, q) \f$, where q is a point on L.
 *  The ratio distance from a point p to a point q is defined as:
 *  \f$ RD(p, q) = \max\{ \max_{i}\{(q_{i} - p_{i}) / p_{i}\}, 0.0 \} \f$.
 *  
 *  Intuitively it is the minimum value of \f$ \epsilon \ge 0 \f$ such 
 *  that some point on L \f$ \epsilon -covers p \f$.
 *  
 *  \sa Line2D, m(), b(), isVertical(), intersection(), parallelThrough(), 
 *      Point and Point::ratioDistance()
 */
double 
Line2D::ratioDistance(const Point& p) const
{
  if (isVertical_) 
    return max( (-b_ - p[0]) / p[0], 0.0 );
  else
    return max( (p[1] - m_ * p[0] - b_) / (m_ * p[0] - p[1]), 0.0 );
}


//! Create a new Line2D instance parallel to the current instance, and 
//! passing through the given point.
/*! 
 *  \param p A Point instance.
 *  \return A Line2D instance parallel to the current instance and 
 *          passing through the given point.
 *  
 *  \sa Line2D, m(), b(), isVertical(), intersection(), ratioDistance() 
 *      and Point
 */
Line2D 
Line2D::parallelThrough(const Point& p) const
{
  if (isVertical_)
    return Line2D(p[0]);
  else
    return Line2D(m_, p[1] - m_ * p[0]);
}


}  // namespace pareto_approximator
