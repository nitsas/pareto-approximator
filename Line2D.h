/*! \file Line2D.h
 *  \brief The declaration of the Line2D class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef LINE_2D_H
#define LINE_2D_H

#include <iostream>
#include <string>

#include "Point.h"
#include "Not2DPointsException.h"
#include "SamePointsException.h"
#include "VerticalLineException.h"
#include "ParallelLinesException.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! A simple 2-dimensional line.
/*! 
 *  A Line2D instance represents a line using an equation of the form:
 *  \f$ y = m x + b \f$. A Line2D instance stores m and b (as private 
 *  attributes m_ and b_ respectively) and keeps track of vertical lines 
 *  using an extra private attribute called isVertical_.
 *
 *  A vertical line is stored using an equation of the form:
 *  \f$ x = -b \f$, so for the line \f$ x = 5 \f$ b (and b_) will be -5.
 *  
 *  As mentioned above, all attributes are private and can be set only by 
 *  the constructor (it wouldn't make sense to change a line after it was 
 *  created). m_, b_ and isVertical_ can be accessed through methods m(), 
 *  b() and isVertical() respectively.
 *
 *  Operators defined:
 *  - operator=(), (defined automatically by the compiler)
 *  - operator==(), 
 *  - operator!=() and 
 *  - operator<<() (as a friend of the class)
 *  
 *  \sa Line2D(), ~Line2D(), operator==(), operator!=(), operator<<(), 
 *      m(), b(), isVertical(), str(), intersection(), ratioDistance() 
 *      and parallelThrough()
 */
class Line2D
{
  public:
    //! The empty constructor. Creates line \f$ y = x \f$.
    Line2D();

    //! A non-vertical line constructor. Creates line \f$ y = m x + b \f$.
    Line2D(double m, double b);

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
    Line2D(const Point& p1, const Point& p2);

    //! A vertical line constructor. Creates line \f$ x = c \f$.
    /*! 
     *  \param c The resulting line's equation should be \f$ x = c \f$.
     *  
     *  Will set b_ to -c so that \f$ 0 = x + b \f$ holds.
     */
    Line2D(double c);

    //! A simple (and empty) destructor.
    ~Line2D();

    //! Return the Line2D instance's slope (m_).
    /*! 
     *  Possible exceptions:
     *  - May throw a VerticalLineException exception if the line is vertical.
     */
    double m() const;

    //! Return the Line2D instance's y-intercept (b_).
    double b() const;

    //! Return true if the Line2D instance represents a vertical line, 
    //! false otherwise.
    bool   isVertical() const;

    //! The Line2D equality operator.
    /*! 
     *  \param line A Line2D instance we want to compare with the current
     *              instance.
     *  \return true if the two instances represent the same line, false 
     *          otherwise.
     *  
     *  \sa Line2D, operator!=() and operator<<()
     */
    bool operator== (const Line2D& line) const;

    //! The Line2D inequality operator.
    /*! 
     *  \param line A Line2D instance we want to compare with the current
     *              instance.
     *  \return true if the two instances represent different lines, false 
     *          otherwise.
     *  
     *  \sa Line2D, operator==() and operator<<()
     */
    bool operator!= (const Line2D& line) const;

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
    friend std::ostream& operator<< (std::ostream& ostr, const Line2D& line);

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
    std::string str() const;

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
    Point  intersection(const Line2D& line) const;

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
    double ratioDistance(const Point& p) const;

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
    Line2D parallelThrough(const Point& p) const;

  private:
    //! The line's slope.
    double m_;
    //! The line's y-intercept.
    double b_;
    //! Is the line vertical (infinite slope)?
    bool isVertical_;
};


}  // namespace pareto_approximator


/* @} */


#endif  // LINE_2D_H
