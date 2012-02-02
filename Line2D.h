/* Line2D.h */


#ifndef LINE_2D_H
#define LINE_2D_H

#include <iostream>

#include "Point.h"
#include "Not2DPointsException.h"
#include "SamePointsException.h"
#include "VerticalLineException.h"
#include "ParallelLinesException.h"


namespace pareto_approximator {


// A simple 2-dimensional line.
// We use line equations of the form:   y = m*x + b
// and isVertical_ to keep track of vertical lines.
class Line2D
{
  public:
    // Constructors
    Line2D();
    Line2D(int m, int b);
    Line2D(double m, double b);
    Line2D(const Point& p1, const Point& p2) throw(Not2DPointsException, 
                                                   SamePointsException);
    // Vertical line constructors (line equation: x = c)
    Line2D(int c);
    Line2D(double c);

    // Attribute getters
    double m() const throw(VerticalLineException);
    double b() const;
    bool   isVertical() const;

    // Operators
    // operator= : use the default (copies all members)
    bool operator== (const Line2D& line) const;
    bool operator!= (const Line2D& line) const;
    // Output stream operator
    friend std::ostream& operator<< (std::ostream& ostr, const Line2D& line);

    // Methods
    Point  intersection(const Line2D& line) const throw(ParallelLinesException);
    double ratioDistance(const Point& p) const;
    Line2D parallelThrough(const Point& p) const;

  private:
    // slope
    double m_;
    // y-intercept
    double b_;
    // is the line vertical (slope = infinity) ?
    bool isVertical_;
};


}  // namespace pareto_approximator


#endif  // LINE_2D_H
