/*! \file Point.h
 *  \brief The definition of Point. (a simple point class)
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef SIMPLE_POINT_CLASS_H
#define SIMPLE_POINT_CLASS_H

#include <string>
#include <vector>
#include <armadillo>

#include "DifferentDimensionsException.h"
#include "NonExistentCoordinateException.h"
#include "NegativeApproximationRatioException.h"
#include "NotPositivePointException.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! A simple 3-dimensional point class. Can represent 2D and 1D points too.
/*! 
 *  Users can check (and set) the point's dimension using the dimension() 
 *  (and dimension(int)) method.
 *  
 *  Operators defined:
 *  - operator=(), (defined automatically by the compiler)
 *  - operator==(), 
 *  - operator!=(), 
 *  - operator[](), 
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
    //! The empty constructor. Creates a zero-dimensional Point.
    Point();
    
    //! An 1-dimensional Point constructor.
    explicit Point(int x);
    
    //! An 1-dimensional Point constructor.
    explicit Point(double x);
    
    //! A 2-dimensional Point constructor.
    Point(int x, int y);
    
    //! A 2-dimensional Point constructor.
    Point(double x, double y);
    
    //! A 3-dimensional Point constructor.
    Point(int x, int y, int z);

    //! A 3-dimensional Point constructor.
    Point(double x, double y, double z);

    //! A 4-dimensional Point constructor.
    Point(int x, int y, int z, int w);

    //! A 4-dimensional Point constructor.
    Point(double x, double y, double z, double w);

    //! An n-dimensional Point constructor.
    /*! 
     *  \param first Iterator to the initial position in a std::vector<int>.
     *  \param last Iterator to the final position in a std::vector<int>.
     *  
     *  The range used is [first, last), which includes all the elements between 
     *  first and last, including the element pointed by first but not the 
     *  element pointed by last.
     *  
     *  The resulting point's coordinates will be doubles, not ints. 
     *  
     *  \sa Point
     */
    Point(std::vector<int>::iterator first, std::vector<int>::iterator last);

    //! An n-dimensional Point constructor.
    /*! 
     *  \param first Iterator to the initial position in a std::vector<double>.
     *  \param last Iterator to the final position in a std::vector<double>.
     *  
     *  The range used is [first, last), which includes all the elements between 
     *  first and last, including the element pointed by first but not the 
     *  element pointed by last.
     *  
     *  \sa Point
     */
    Point(std::vector<double>::iterator first, 
          std::vector<double>::iterator last);

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
     *  
     *  \sa Point
     */
    Point(int* first, int* last);

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
     *  \sa Point
     */
    Point(double* first, double* last);

    //! A simple (and empty) destructor.
    ~Point();

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
    double operator[] (unsigned int pos) const;

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
    bool operator== (const Point & p) const;
    
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
    bool operator!= (const Point & p) const;

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
    bool operator< (const Point & p) const;

    //! The Point plus operator.
    /*!
     *  \param p A Point instance.
     *  \return A new Point instance, having the same dimensions as the 
     *          current instance and p.
     *
     *  The new point's i'th coordinate (for all \f$ 1 \le i \le dimension()\f$)
     *  will be the sum of the i'th coordinate of p and the i'th coordinate 
     *  of the current instance.
     *  
     *  \sa Point
     */
    Point operator+ (const Point & p) const;

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
    friend std::ostream & operator<< (std::ostream & out, const Point & p);

    //! The Point input stream operator. A friend of the Point class.
    /*! 
     *  The accepted format is similar to the one operator<<() uses for output.
     *  Zero-dimensional points of the form "()" are not accepted.
     *  
     *  \sa Point, Point::str(), Point::operator==(), Point::operator!=(), 
     *      Point::operator<(), Point::operator[]() and 
     *      std::ostream& operator<<(std::ostream&, Point&)
     */
    friend std::istream & operator>> (std::istream & istr, Point & p);

    //! Get a Point instance's dimension. (1D, 2D or 3D point)
    /*! 
     *  \return 1, 2 or 3 for 1-dimensional, 2-dimensional or 3-dimensional 
     *          Point instances respectively.
     *  
     *  \sa Point and bool Point::dimension(int)
     */
    unsigned int dimension() const;

    //! Set the point's dimension. (1D, 2D, 3D, etc. point)
    /*! 
     *  \param dimension The dimension we want to change the Point instance to.
     *  
     *  If "dimension" is smaller than the current Point dimension only the 
     *  first "dimension" coordinates will be kept, the rest being dropped.
     *  
     *  \sa Point and int Point::dimension() const
     */
    void dimension(unsigned int dimension);

    //! Check if the Point's coordinates are all zero.
    /*!
     *  \return true if all the Point's coordinates are zero; false otherwise.
     *          (false if the Point is dimensionless)
     */
    bool isZero() const;

    //! Return the Point instance as a string.
    /*! 
     *  Makes a string with the Point instance's dimensions inside parentheses. 
     *  Examples:
     *  - (1.0, 4.27, 0.883)
     *  - (3.0)
     *  - (5, 1.99204e+09)
     *  - ()      <-- dimensionless point
     *  - etc
     *
     *  /sa Point and std::ostream& operator<<(std::ostream&, const Point&)
     */
    std::string str() const;

    /*! 
     *  \brief Return the point's coordinates as an arma::vec (armadillo vector).
     */
    arma::vec toVec() const;

    /*! 
     *  \brief Return the point's coordinates as an arma::rowvec (armadillo row vector).
     */
    arma::rowvec toRowVec() const;

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
     *  \sa Point and Point::dominates()
     */
    double ratioDistance(const Point & q) const;

    //! Check if the current point (p) eps-covers the given point (q).
    /*! 
     *  \param q A Point instance with \f$ q_{i} \ge 0 \f$ for all i.
     *  \param eps An approximation factor.
     *  \return true if p eps-covers q; false otherwise.
     *  
     *  Note that both p and q must be greater than zero (dominated by 0);
     *  that is both \f$ p_{i} \ge 0 \f$ and \f$ q_{i} \ge 0 \f$ must hold
     *  for all i.
     *  
     *  We say that p \f$ \epsilon \f$-covers q (\f$\epsilon \ge 0 \f$) iff 
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
     *  - May throw a DifferentDimensionsException if p and q are of 
     *    different dimensions.
     *  
     *  \sa Point and Point::ratioDistance()
     */
    bool dominates(const Point & q, double eps=0.0) const;

  private:
    //! The Point instance's coordinates. The instance's dimension will 
    //! always be equal to the vector's size.
    std::vector<double> coordinates_;
};


}  // namespace pareto_approximator


/* @} */


#endif  // SIMPLE_POINT_CLASS_H
