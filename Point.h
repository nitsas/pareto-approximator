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

#include "NullObjectException.h"
#include "DifferentDimensionsException.h"
#include "NonExistentCoordinateException.h"
#include "NegativeApproximationRatioException.h"
#include "NotPositivePointException.h"
#include "NotStrictlyPositivePointException.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
namespace pareto_approximator {


//! A simple point class. 
/*! 
 *  Users can check the point's dimension using the dimension() method.
 *  
 *  The empty constructor will make a special kind of point, a null point.
 *  The method Point::isNull() will return true for null points and false 
 *  for any other point.
 *  
 *  Most of the rest of a null point's methods will throw a 
 *  NullObjectException since no operations on a null point would be valid.
 */
class Point {
  public:
    //! The empty constructor. Creates a null Point instance.
    Point();
    
    //! An 1-dimensional Point constructor.
    explicit Point(int x);
    
    //! An 1-dimensional Point constructor.
    explicit Point(double x);
    
    //! A 2-dimensional Point constructor.
    Point(int x, int y);

    //! A 2-dimensional Point constructor.
    Point(unsigned int x, unsigned int y);
    
    //! A 2-dimensional Point constructor.
    Point(double x, double y);
    
    //! A 3-dimensional Point constructor.
    Point(int x, int y, int z);

    //! A 3-dimensional Point constructor.
    Point(unsigned int x, unsigned int y, unsigned int z);
    
    //! A 3-dimensional Point constructor.
    Point(double x, double y, double z);

    //! A 4-dimensional Point constructor.
    Point(int x, int y, int z, int w);

    //! A 4-dimensional Point constructor.
    Point(unsigned int x, unsigned int y, unsigned int z, unsigned int w);
    
    //! A 4-dimensional Point constructor.
    Point(double x, double y, double z, double w);

    //! An n-dimensional Point constructor.
    /*! 
     *  \param first Iterator to the initial position in a std::vector<int>.
     *  \param last Iterator to the final position in a std::vector<int>.
     *  
     *  The range used is [first, last), which includes all the elements 
     *  between first and last, including the element pointed by first but 
     *  not the element pointed by last.
     *  
     *  The resulting point's coordinates will be doubles, not ints. 
     *  
     *  \sa Point
     */
    Point(std::vector<int>::const_iterator first, 
          std::vector<int>::const_iterator last);

    //! An n-dimensional Point constructor.
    /*! 
     *  \param first Iterator to the initial position in a 
     *               std::vector<unsigned int>.
     *  \param last Iterator to the final position in a 
     *              std::vector<unsigned int>.
     *  
     *  The range used is [first, last), which includes all the elements 
     *  between first and last, including the element pointed by first but 
     *  not the element pointed by last.
     *  
     *  The resulting point's coordinates will be doubles, not ints. 
     *  
     *  \sa Point
     */
    Point(std::vector<unsigned int>::const_iterator first, 
          std::vector<unsigned int>::const_iterator last);

    //! An n-dimensional Point constructor.
    /*! 
     *  \param first Iterator to the initial position in a std::vector<double>.
     *  \param last Iterator to the final position in a std::vector<double>.
     *  
     *  The range used is [first, last), which includes all the elements 
     *  between first and last, including the element pointed by first but 
     *  not the element pointed by last.
     *  
     *  \sa Point
     */
    Point(std::vector<double>::const_iterator first, 
          std::vector<double>::const_iterator last);

    //! An n-dimensional Point constructor.
    /*! 
     *  \param first Iterator (pointer) to the initial position in an array of 
     *               double.
     *  \param last Iterator (pointer) to the final position in an array of 
     *              double. (the position just beyond the last element we want)
     *  
     *  The range used is [first, last), which includes all the elements 
     *  between first and last, including the element pointed by first but 
     *  not the element pointed by last.
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
     *  The range used is [first, last), which includes all the elements 
     *  between first and last, including the element pointed by first but 
     *  not the element pointed by last.
     *  
     *  \sa Point
     */
    Point(double* first, double* last);

    //! An n-dimensional Point constructor.
    /*!
     *  \param first Iterator to the initial position in an armadillo column 
     *               vector of double (arma::vec).
     *  \param last Iterator to the final position in an armadillo column 
     *              vector of double (arma::vec). (the position just beyond 
     *              the last element we want)
     *  
     *  The range used is [first, last), which includes all the elements 
     *  between first and last, including the element pointed by first but 
     *  not the element pointed by last.
     *  
     *  \sa Point
     */
    Point(arma::vec::const_iterator first, arma::vec::const_iterator last);

    //! A simple (and empty) destructor.
    ~Point();

    //! Check if the Point instance is null. (i.e. isNull_ == true)
    /*! 
     *  \return true if the Point instance is null; false otherwise
     *
     *  \sa Point
     */
    bool isNull() const;

    //! Get a Point instance's dimension. (1D, 2D or 3D point)
    /*! 
     *  \return 1, 2 or 3 for 1-dimensional, 2-dimensional or 3-dimensional 
     *          Point instances respectively.
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  
     *  \sa Point
     */
    unsigned int dimension() const;

    //! Check if the Point's coordinates are all zero.
    /*!
     *  \return true if all the Point's coordinates are zero; false otherwise.
     *          (false if the Point is dimensionless)
     *
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  
     *  \sa Point
     */
    bool isZero() const;

    /*!
     *  \brief Check if the Point is positive (i.e. all coordinates greater 
     *         than zero).
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  
     *  \sa Point
     */
    bool isPositive() const;

    /*! 
     *  \brief Check if the Point is strictly positive (i.e. all coordinates 
     *         strictly greater than zero).
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  
     *  \sa Point
     */
    bool isStrictlyPositive() const;

    //! The Point access coordinate operator. 
    /*! 
     *  \param pos The position (coordinate) to access. 
     *             (0 <= pos < dimension())
     *  \return The Point's "pos" coordinate.
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  - May throw a NonExistentCoordinateException exception if the 
     *    requested coordinate does not exist. ("pos" is greater than or 
     *    equal to dimension())
     *  
     *  \sa Point, dimension(), operator==(), operator!=(), operator<(), 
     *      std::ostream & operator<<(std::ostream &, Point &) and 
     *      std::istream & operator>>(std::istream &, Point &)
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
     *      std::ostream & operator<<(std::ostream &, Point &) and 
     *      std::istream & operator>>(std::istream &, Point &)
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
     *      std::ostream & operator<<(std::ostream &, Point &) and 
     *      std::istream & operator>>(std::istream &, Point &)
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
     *  - May throw a NullObjectException exception if the point is null.
     *  - May throw a DifferentDimensionsException exception if the two Point 
     *    instances are of different dimensions (can't be compared).
     *  
     *  \sa Point, operator==(), operator!=(), operator[](), 
     *      std::ostream & operator<<(std::ostream &, Point &) and 
     *      std::istream & operator>>(std::istream &, Point &)
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
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  - May throw a DifferentDimensionsException exception if the two Point 
     *    instances are of different dimensions (can't be added).
     *  
     *  \sa Point
     */
    Point operator+ (const Point & p) const;

    //! The Point minus operator.
    /*!
     *  \param p A Point instance.
     *  \return A new Point instance, having the same dimensions as the 
     *          current instance and p.
     *
     *  The new point's i'th coordinate (for all \f$ 1 \le i \le dimension()\f$)
     *  will be the this instance's i'th coordinate minus p's i'th coordinate.
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  - May throw a DifferentDimensionsException exception if the two Point 
     *    instances are of different dimensions (can't be added).
     *  
     *  \sa Point
     */
    Point operator- (const Point & p) const;

    //! Return the Point instance as a string.
    /*! 
     *  \param rawCoordinates If true, output the coordinates seperated 
     *                        by spaces; otherwise output the coordinates 
     *                        separated by ", ", inside parentheses.
     *  
     *  Makes a string with the Point instance's coordinates. 
     *  - If `rawCoordinates==true` the result will be the point's 
     *    coordinates separated by spaces.
     *  - If `rawCoordinates==false` the result will be the point's 
     *    coordinates inside parentheses, separated by ", ". 
     *  
     *  Examples for rawCoordinates == true:
     *  - 1.0 4.27 0.883
     *  - 3.0
     *  - 5 1.99204e+09
     *  -       <-- dimensionless point (no spaces)
     *  -       <-- null point (no spaces)
     *  - etc
     *  
     *  Examples for rawCoordinates == false:
     *  - (1.0, 4.27, 0.883)
     *  - (3.0)
     *  - (5, 1.99204e+09)
     *  - ()      <-- dimensionless point
     *  - ()      <-- null point
     *  - etc
     *
     *  /sa Point and std::ostream & operator<<(std::ostream &, const Point &)
     */
    std::string str(bool rawCoordinates=false) const;

    /*! 
     *  \brief Return the point's coordinates as an arma::vec (armadillo 
     *         vector).
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  
     *  \sa Point
     */
    arma::vec toVec() const;

    /*! 
     *  \brief Return the point's coordinates as an arma::rowvec (armadillo 
     *         row vector).
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  
     *  \sa Point
     */
    arma::rowvec toRowVec() const;

    //! Return the distance from the current to the given Point instance.
    /*!
     *  \param q A Point instance.
     *  \return The distance from the current Point instance (*this) to 
     *          the given Point instance.
     *  
     *  There are different possible distance metrics we could use (e.g. 
     *  ratio distance, Euclidean distance, additive distance etc.). 
     *  
     *  Currently using the additive distance metric.
     *
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  - May throw a DifferentDimensionsException exception if the two 
     *    Point instances are of different dimensions (can't be compared).
     *  - May throw a NotStrictlyPositivePointException exception if either the 
     *    current instance or the given point is not strictly positive.
     *  
     *  \sa Point, Point::additiveDistance(), Point::euclideanDistance(), 
     *      Point::ratioDistance() and Point::dominates()
     */
    double distance(const Point & q) const;

    /*! 
     *  \brief Return the Euclidean distance from the current to the given 
     *         Point instance.
     *
     *  The formula for the Euclidean distance between two d-dimensional 
     *  points p and q is:
     *  \f$ ED(p, q) = \sqrt{(p_1 - q_1)^2 + (p_2 - q_2)^2 + ... + 
     *                       (p_d - q_d)^2} \f$
     *  
     *  In contrast to ratioDistance(), euclideanDistance() permits points 
     *  with negative coordinates.
     *  
     *  \sa Point, Point::distance() and Point::dominates()
     */
    double euclideanDistance(const Point & q) const;

    //! Return the ratio distance from the current to the given Point instance.
    /*! 
     *  \param q A Point instance.
     *  \return The ratio distance from the current Point instance (*this) 
     *          to the given Point instance.
     *  
     *  We define the ratio distance from point p to q as:
     *  \f$ RD(p, q) = \max\{ \max_{i}\{ q_{i}/p{i} - 1 \}, 0 \} \f$.
     *  
     *  Intuitively, it is the minimum value of \f$ \epsilon \ge 0 \f$ such 
     *  that q \f$\epsilon\f$ -dominates (\f$\epsilon\f$ -covers) p in the 
     *  multiplicative sense.
     *  
     *  The concept of ratio distance breaks down if one (or both) of the 
     *  two points is not strictly positive. A point p is strictly positive 
     *  iff \f$ p_{i} > 0.0 \f$ holds for all i, i.e. every coordinate is 
     *  strictly greater than zero.
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  - May throw a DifferentDimensionsException exception if the two 
     *    Point instances are of different dimensions (can't be compared).
     *  - May throw a NotStrictlyPositivePointException exception if either 
     *    the current instance or the given point is not strictly positive.
     *
     *  \sa Point, Point::distance() and Point::dominates()
     */
    double ratioDistance(const Point & q) const;

    /*!
     *  \brief Computes the minimum value of \f$\epsilon\f$ such that q 
     *         \f$\epsilon\f$ -dominates the current instance in the additive 
     *         sense.
     *  
     *  We say that a point p is \f$\epsilon\f$ -dominated (\f$\epsilon\f$ 
     *  -covered) by a point q in the additive sense if: 
     *  \f$ q_i \le p_i + \epsilon \f$ for all i.
     *  
     *  In other words, the additive distance from a point p to a point q 
     *  is defined as:
     *  \f$ AD(p, q) = \max\{ \max_{i}\{(q_{i} - p_{i})\}, 0.0 \} \f$.
     *  
     *  Intuitively it is the minimum value of \f$ \epsilon \ge 0 \f$ such 
     *  that q \f$\epsilon\f$ -dominates (\f$\epsilon\f$ -covers) p in the 
     *  additive sense.
     *  
     *  Note that this method does not require the points to be positive.
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  - May throw a DifferentDimensionsException exception if the two 
     *    Point instances are of different dimensions (can't be compared).
     *  
     *  \sa Point, Point::distance() and Point::dominatesAdditive()
     */
    double additiveDistance(const Point & q) const;

    //! Check if the current point (p) eps-dominates the given point (q).
    /*!
     *  \param q A Point instance.
     *  \param eps An approximation parameter.
     *  \return true if p eps-covers q; false otherwise.
     *  
     *  There are two different definitions of eps-dominance we could use:
     *  - Additive eps-dominance. Where a point q is \f$\epsilon\f$ 
     *    -dominated by a point p if: 
     *    \f$ p_{i} \le q_{i} + \epsilon \f$ for all i.
     *  - Multiplicative eps-dominance. Where a point q is \f$\epsilon\f$
     *    -dominated by a point p if: 
     *    \f$ p_{i} \le (1 + \epsilon) q_{i} \f$ for all i.
     *
     *  Note that for the additive variant both points must be positive (i.e. 
     *  both \f$ p_{i} \ge 0 \f$ and \f$ q_{i} \ge 0 \f$ must hold for all 
     *  i), whereas for the multiplicative variant both points must be 
     *  strictly positive (i.e. both \f$ p_{i} > 0 \f$ and \f$ q_{i} > 0 \f$ 
     *  must hold for all i).
     *  
     *  This method checks for errors and then calls either the 
     *  dominatesAdditive() or the dominatesMultiplicative() method. 
     *  (which one will be called is currently hardcoded in this method's 
     *  implementation) 
     *  
     *  Currently using the additive error measure.
     *  
     *  If eps=0.0 the method simply checks whether or not p dominates 
     *  q and that is how it got its name.
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  - May throw a NegativeApproximationRatioException if \f$ eps < 0 \f$.
     *  - May throw a DifferentDimensionsException if p and q are of different 
     *    dimensions.
     *  
     *  \sa Point, Point::dominatesAdditive() and 
     *      Point::dominatesMultiplicative()
     */
    bool dominates(const Point & q, double eps=0.0) const;

    /*!
     *  \brief Check if the current point (p) eps-dominates (in the additive 
     *         sense) the given point (q).
     *  
     *  We say that a point q is \f$\epsilon\f$ -dominated (\f$\epsilon\f$ 
     *  -covered) by a point p in the additive sense if: 
     *  \f$ p_i \le q_i + \epsilon \f$ for all i.
     *
     *  Note that both p and q must be greater than zero for the additive 
     *  variant, that is both \f$ p_{i} \ge 0 \f$ and \f$ q_{i} \ge 0 \f$ 
     *  must hold for all i.
     *  
     *  If eps=0.0 the method simply checks whether or not p dominates 
     *  q and that is how it got its name.
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  - May throw a NegativeApproximationRatioException if \f$ eps < 0 \f$.
     *  - May throw a DifferentDimensionsException if p and q are of different 
     *    dimensions.
     *  - May throw a NotPositivePointException exception if either p or q 
     *    is not positive (i.e. some coordinate is less than zero).
     *  
     *  \sa Point, Point::dominates() and Point::distance()
     */
    bool dominatesAdditive(const Point & q, double eps=0.0) const;

    /*! 
     *  \brief Check if the current point (p) eps-dominates (in the 
     *         multiplicative sense) the given point (q).
     *  
     *  \param q A Point instance with \f$ q_{i} \ge 0 \f$ for all i.
     *  \param eps An approximation factor.
     *  \return true if p eps-covers q; false otherwise.
     *  
     *  We say that a point q is \f$\epsilon\f$ -dominated (\f$\epsilon\f$ 
     *  -covered) by a point p in the multiplicative sense if: 
     *  \f$ p_i \le (1 + \epsilon) q_i \f$ for all i.
     *  
     *  Note that both p and q must be strictly greater than zero for the 
     *  multiplicative variant, that is both \f$ p_{i} > 0 \f$ and 
     *  \f$ q_{i} > 0 \f$ must hold for all i.
     *  
     *  If eps=0.0 the method simply checks whether or not p dominates 
     *  q and that is how it got its name.
     *  
     *  Possible exceptions:
     *  - May throw a NullObjectException exception if the point is null.
     *  - May throw a NegativeApproximationRatioException if \f$ eps < 0 \f$.
     *  - May throw a DifferentDimensionsException if p and q are of different 
     *    dimensions.
     *  - May throw a NotStrictlyPositivePointException if either p or q 
     *    is not strictly positive (i.e. some coordinate is not strictly 
     *    greater than zero).
     *  
     *  \sa Point and Point::ratioDistance()
     */
    bool dominatesMultiplicative(const Point & q, double eps=0.0) const;

    //! The Point output stream operator. A friend of the Point class.
    /*! 
     *  The Point instance's coordinates will be output separated by spaces. 
     *  Examples:
     *  - 1.0 4.27 0.883
     *  - 3.0
     *  - 5 1.99204e+09
     *  -       <-- dimensionless point (no spaces)
     *  -       <-- null point (no spaces)
     *  - etc
     *
     *  operator<<() uses the Point instance's Point::str(true) method to 
     *  create the output.
     *
     *  Usage:
     *  - std::cout << Point(4, 3, -8);
     *  - std::cout << "some text " << Point(2.7, -2.7) << " more text";
     *
     *  \sa Point, Point::str(), Point::operator==(), Point::operator!=(), 
     *      Point::operator<(), Point::operator[]() and 
     *      std::istream & operator>>(std::istream &, Point &)
     */
    friend std::ostream & operator<< (std::ostream & out, const Point & p);

    //! The Point input stream operator. A friend of the Point class.
    /*! 
     *  The accepted format is similar to the one operator<<() uses for output:
     *  Space-separated coordinates. One point per line.
     *
     *  Null points of the form "" are accepted - we encourage you to avoid 
     *  them though (a null point might also be returned in case of an error).
     *  
     *  \sa Point, Point::str(), Point::operator==(), Point::operator!=(), 
     *      Point::operator<(), Point::operator[]() and 
     *      std::ostream & operator<<(std::ostream &, Point &)
     */
    friend std::istream & operator>> (std::istream & istr, Point & p);

  private:
    /*! \brief The Point instance's coordinates. The instance's dimension 
     *         will always be equal to the vector's size (when the instance 
     *         is not null, i.e. isNull_ == false).
     */
    std::vector<double> coordinates_;

    /*! \brief Is the instance null? (i.e. Was it created using the empty 
     *         constructor?)
     */
    bool isNull_;
};


}  // namespace pareto_approximator


/* @} */


// We will #include the implementation here because we want to make a 
// header-only code base.
#include "Point.cpp"


#endif  // SIMPLE_POINT_CLASS_H
