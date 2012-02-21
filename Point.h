/*! \file Point.h
 *  \brief A file containing the declaration of the Point class.
 */


#ifndef SIMPLE_POINT_H
#define SIMPLE_POINT_H

#include <iostream>
#include <string>

#include "DifferentDimensionsException.h"


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! A simple 3-dimensional point class. Can represent 2D and 1D points too.
/*! 
 *  A point with 3 dimensions called x, y and z. In the case of 2-dimensional 
 *  or 1-dimensional points higher dimensions contain undefined values. Users 
 *  can differentiate between points like (3.0, 10.0, 0.0) and (3.0, 10.0) by 
 *  checking each point's Point::dimension() method.
 *  
 *  Operators defined:
 *  - operator=(), (defined automatically by the compiler)
 *  - operator==(), 
 *  - operator!=(), 
 *  - operator<(), 
 *  - operator<<() (as a friend of the class) and
 *  - operator>>() (as a friend of the class).
 *  
 *  /sa Point(), ~Point(), operator==(), operator!=(), operator<(), 
 *      operator<<(), operator>>(), int dimension() const, bool dimension(int),
 *      str() and ratioDistance()
 */

class Point {
  public:
    //! A Point's 1st dimension.
    double x;
    //! A Point's 2nd dimension. (value is undefined for 1-dimensional points)
    double y;
    //! A Point's 3rd dimension. (value is undefined for 1-dimensional and 
    //! 2-dimensional points)
    double z;

    //! The empty constructor. Creates an all-zero 3-dimensional Point.
    Point();
    
    //! An 1-dimensional Point constructor. 
    /*! The resulting Point's dimensions will be doubles, not ints. */
    Point(int xx);
    
    //! An 1-dimensional Point constructor.
    Point(double xx);
    
    //! A 2-dimensional Point constructor.
    /*! The resulting Point's dimensions will be doubles, not ints. */
    Point(int xx, int yy);
    
    //! A 2-dimensional Point constructor.
    Point(double xx, double yy);
    
    //! A 3-dimensional Point constructor.
    /*! The resulting Point's dimensions will be doubles, not ints. */
    Point(int xx, int yy, int zz);
    
    //! A 3-dimensional Point constructor.
    Point(double xx, double yy, double zz);

    //! A simple (and empty) Destructor.
    ~Point();

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
    bool operator== (const Point& p) const;
    
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
    bool operator!= (const Point& p) const;

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
    bool operator<  (const Point& p) const;

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
    friend std::ostream& operator<< (std::ostream& ostr, const Point& p);

    //! The Point input stream operator. A friend of the Point class.
    /*! 
     *  The accepted format is similar to the one operator<<() uses for output.
     *  
     *  \sa Point, Point::operator==(), Point::operator!=(), 
     *      Point::operator<() and 
     *      std::ostream& operator<<(std::ostream&, Point&)
     */
    friend std::istream& operator>> (std::istream& istr, Point& p);

    //! Get a Point instance's dimension. (1D, 2D or 3D point)
    /*! 
     *  \return 1, 2 or 3 for 1-dimensional, 2-dimensional or 3-dimensional 
     *          Point instances respectively.
     *  
     *  \sa Point and bool Point::dimension(int)
     */
    int dimension() const;

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
    bool dimension(int dim);

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
    std::string str() const;

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
    double ratioDistance(const Point& q) const;

  private:
    //! The Point instance's dimension. (1D, 2D or 3D point)
    int dimension_;
};


}  // namespace pareto_approximator


#endif  // SIMPLE_POINT_H
