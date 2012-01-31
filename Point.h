/* Point.h */


#ifndef SIMPLE_POINT_H
#define SIMPLE_POINT_H

#include <iostream>

#include "DifferentDimensionsException.h"


namespace pareto_approximator {


// A simple 3-dimensional point class.
// It can represent 2-dimensional and 1-dimensional points too.
class Point {
  public:
    // attributes
    double x, y, z;

    // Constructors
    Point();
    // 1-dimensional point
    Point(int xx);
    Point(double xx);
    // 2-dimensional point
    Point(int xx, int yy);
    Point(double xx, double yy);
    // 3-dimensional point
    Point(int xx, int yy, int zz);
    Point(double xx, double yy, double zz);

    // Destructor
    ~Point();

    // Operators
    // operator= : use the default (copies all members)
    bool operator== (const Point& p);
    bool operator!= (const Point& p);

    // Input/Output stream operators
    friend std::ostream& operator<< (std::ostream& ostr, Point& p);
    friend std::istream& operator>> (std::istream& istr, Point& p);

    // Get the point's dimension. (1D, 2D or 3D point)
    int dimension() const;
    // Set the point's dimension. (1D, 2D or 3D point)
    bool dimension(int dim);

  private:
    // The point's dimension. (1D, 2D or 3D point)
    int dimension_;
};




}  // namespace pareto_approximator


#endif  // SIMPLE_POINT_H
