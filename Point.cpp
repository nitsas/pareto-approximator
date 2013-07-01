/*! \file Point.cpp
 *  \brief The implementation of Point. (a simple point class)
 *  \author Christos Nitsas
 *  \date 2012
 *  
 *  Won't `include` Point.h. In fact Point.h will `include` 
 *  Point.cpp because we want a header-only code base.
 */


#include <assert.h>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
namespace pareto_approximator {


//! The empty constructor. Creates a null Point instance.
Point::Point() : isNull_(true) { }


//! An 1-dimensional Point constructor.
Point::Point(int x) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
}


//! An 1-dimensional Point constructor.
Point::Point(double x) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
}


//! A 2-dimensional Point constructor.
Point::Point(int x, int y) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
}


//! A 2-dimensional Point constructor.
Point::Point(unsigned int x, unsigned int y) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
}


//! A 2-dimensional Point constructor.
Point::Point(double x, double y) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
}


//! A 3-dimensional Point constructor.
Point::Point(int x, int y, int z) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
  coordinates_.push_back(z);
}


//! A 3-dimensional Point constructor.
Point::Point(unsigned int x, unsigned int y, unsigned int z) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
  coordinates_.push_back(z);
}


//! A 3-dimensional Point constructor.
Point::Point(double x, double y, double z) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
  coordinates_.push_back(z);
}


//! A 4-dimensional Point constructor.
Point::Point(int x, int y, int z, int w) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
  coordinates_.push_back(z);
  coordinates_.push_back(w);
}


//! A 4-dimensional Point constructor.
Point::Point(unsigned int x, unsigned int y, 
             unsigned int z, unsigned int w) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
  coordinates_.push_back(z);
  coordinates_.push_back(w);
}


//! A 4-dimensional Point constructor.
Point::Point(double x, double y, double z, double w) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.push_back(x);
  coordinates_.push_back(y);
  coordinates_.push_back(z);
  coordinates_.push_back(w);
}


//! An n-dimensional Point constructor.
/*! 
 *  \param first Iterator to the initial position in a std::vector<double>.
 *  \param last Iterator to the final position in a std::vector<double>.
 *  
 *  The range used is [first, last), which includes all the elements between 
 *  first and last, including the element pointed by first but not the 
 *  element pointed by last.
 *  
 *  The resulting point's coordinates will be doubles, not ints. 
 */
Point::Point(std::vector<int>::const_iterator first, 
             std::vector<int>::const_iterator last) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.assign(first, last);
}


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
Point::Point(std::vector<unsigned int>::const_iterator first, 
             std::vector<unsigned int>::const_iterator last)
{
  assert(coordinates_.size() == 0);
  coordinates_.assign(first, last);
}


//! An n-dimensional Point constructor.
/*! 
 *  \param first Iterator to the initial position in a std::vector<double>.
 *  \param last Iterator to the final position in a std::vector<double>.
 *  
 *  The range used is [first, last), which includes all the elements between 
 *  first and last, including the element pointed by first but not the 
 *  element pointed by last.
 */
Point::Point(std::vector<double>::const_iterator first, 
             std::vector<double>::const_iterator last) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.assign(first, last);
}


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
 */
Point::Point(int* first, int* last) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.assign(first, last);
}


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
 */
Point::Point(double* first, double* last) : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.assign(first, last);
}


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
Point::Point(arma::vec::const_iterator first, arma::vec::const_iterator last) 
                                          : isNull_(false)
{
  assert(coordinates_.size() == 0);
  coordinates_.assign(first, last);
}


//! A simple (and empty) destructor.
Point::~Point() { }


//! Check if the Point instance is null. (i.e. isNull_ == true)
/*! 
 *  \return true if the Point instance is null; false otherwise
 *  
 *  \sa Point
 */
bool 
Point::isNull() const
{
  return isNull_;
}


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
unsigned int 
Point::dimension() const
{
  if (isNull())
    throw exception_classes::NullObjectException();
  // else

  return coordinates_.size();
}


//! Check if the Point's coordinates are all zero.
/*!
 *  \return true if all the Point's coordinates are zero; false otherwise.
 *
 *  Possible exceptions:
 *  - May throw a NullObjectException exception if the point is null.
 *  
 *  \sa Point
 */
bool 
Point::isZero() const
{
  if (isNull())
    throw exception_classes::NullObjectException();
  // else

  assert(dimension() > 0);

  std::vector<double>::const_iterator it;
  for (it = coordinates_.begin(); it != coordinates_.end(); ++it)
    if (*it != 0.0)
      return false;

  return true;
}


/*!
 *  \brief Check if the Point is positive (i.e. all coordinates greater 
 *         than zero).
 *  
 *  Possible exceptions:
 *  - May throw a NullObjectException exception if the point is null.
 *  
 *  \sa Point
 */
bool 
Point::isPositive() const
{
  if (isNull())
    throw exception_classes::NullObjectException();
  // else

  assert(dimension() > 0);

  std::vector<double>::const_iterator it;
  for (it = coordinates_.begin(); it != coordinates_.end(); ++it)
    if (*it < 0.0)
      return false;

  return true;
}


/*! 
 *  \brief Check if the Point is strictly positive (i.e. all coordinates 
 *         strictly greater than zero).
 *  
 *  Possible exceptions:
 *  - May throw a NullObjectException exception if the point is null.
 *  
 *  \sa Point
 */
bool 
Point::isStrictlyPositive() const
{
  if (isNull())
    throw exception_classes::NullObjectException();
  // else

  assert(dimension() > 0);

  std::vector<double>::const_iterator it;
  for (it = coordinates_.begin(); it != coordinates_.end(); ++it)
    if (*it <= 0.0)
      return false;

  return true;
}


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
double 
Point::operator[] (unsigned int pos) const
{
  if (isNull())
    throw exception_classes::NullObjectException();
  // else

  if (pos >= dimension())
    throw exception_classes::NonExistentCoordinateException();
  else
    return coordinates_[pos];
}


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
bool 
Point::operator== (const Point & p) const
{
  if (isNull() and p.isNull())
    return true;
  else if (isNull() or p.isNull())
    return false;
  // else 
  
  // both points should be non-null now
  assert(not (isNull() or p.isNull()));

  if (dimension() != p.dimension())
    return false;

  for (unsigned int i=0; i<dimension(); ++i) 
    if (coordinates_[i] != p[i]) 
      return false;

  return true;
}


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
bool 
Point::operator!= (const Point & p) const
{
  if (isNull() and p.isNull())
    return false;
  else if (isNull() or p.isNull())
    return true;
  // else 

  // both points should be non-null now
  assert(not (isNull() or p.isNull()));

  if (dimension() != p.dimension())
    return true;

  for (unsigned int i=0; i<dimension(); ++i) 
    if (coordinates_[i] != p[i]) 
      return true;

  return false;
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
 *  - May throw a NullObjectException exception if either instance is null.
 *  - May throw a DifferentDimensionsException exception if the two Point 
 *    instances are of different dimensions (can't be compared).
 *  
 *  \sa Point, operator==(), operator!=(), operator[](), 
 *      std::ostream & operator<<(std::ostream &, Point &) and 
 *      std::istream & operator>>(std::istream &, Point &)
 */
bool 
Point::operator< (const Point & p) const 
{
  if (isNull() or p.isNull()) 
    throw exception_classes::NullObjectException();
  // else

  if (dimension() != p.dimension()) 
    throw exception_classes::DifferentDimensionsException();
  // else

  for (unsigned int i=0; i<dimension(); ++i) {
    if (coordinates_[i] < p[i]) 
      return true;
    else if (coordinates_[i] > p[i]) 
      return false;
    else 
      continue;
  }

  return false;
}


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
Point 
Point::operator+ (const Point & p) const
{
  if (isNull() or p.isNull())
    throw exception_classes::NullObjectException();
  // else

  if (dimension() != p.dimension()) 
    throw exception_classes::DifferentDimensionsException();
  // else
  std::vector<double> newCoords(dimension());
  for (unsigned int i = 0; i < dimension(); ++i) 
    newCoords[i] = coordinates_[i] + p[i];
  return Point(newCoords.begin(), newCoords.end());
}


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
Point 
Point::operator- (const Point & p) const
{
  if (isNull() or p.isNull())
    throw exception_classes::NullObjectException();
  // else

  if (dimension() != p.dimension()) 
    throw exception_classes::DifferentDimensionsException();
  // else
  std::vector<double> newCoords(dimension());
  for (unsigned int i = 0; i < dimension(); ++i) 
    newCoords[i] = coordinates_[i] - p[i];
  return Point(newCoords.begin(), newCoords.end());
}


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
std::string 
Point::str(bool rawCoordinates) const
{
  std::string separator, beginning, end;

  if (not rawCoordinates) {
    // rawCoordinates == false
    separator = ", ";
    beginning = "(";
    end = ")";
  }
  else {
    // rawCoordinates == true
    separator = " ";
    beginning = "";
    end = "";
  }

  if ( isNull() or (dimension() == 0) ) {
    std::stringstream ss;
    ss << beginning << end;
    return ss.str();
  }
  else {
    std::stringstream ss;

    if (rawCoordinates) {
      // coordinates in scientific notation, with 13 digits after the dot
      ss.precision(20);
      ss << std::scientific;
    }

    ss << beginning << coordinates_[0];
    for (unsigned int i=1; i<dimension(); ++i) 
      ss << separator << coordinates_[i];
    ss << end;

    return ss.str();
  }
}


/*! 
 *  \brief Return the point's coordinates as an arma::vec (armadillo 
 *         vector).
 *
 *  Possible exceptions:
 *  - May throw a NullObjectException exception if the point is null.
 *  
 *  \sa Point
 */
arma::vec 
Point::toVec() const
{
  if (isNull())
    throw exception_classes::NullObjectException();
  // else

  return arma::vec(coordinates_);
}


/*! 
 *  \brief Return the point's coordinates as an arma::rowvec (armadillo 
 *         row vector).
 *  
 *  Possible exceptions:
 *  - May throw a NullObjectException exception if the point is null.
 *  
 *  \sa Point
 */
arma::rowvec 
Point::toRowVec() const
{
  if (isNull())
    throw exception_classes::NullObjectException();
  // else

  return arma::rowvec(coordinates_);
}


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
double 
Point::distance(const Point & q) const
{
  return additiveDistance(q);
}


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
double 
Point::euclideanDistance(const Point & q) const
{
  if (isNull() or q.isNull())
    throw exception_classes::NullObjectException();
  if (dimension() != q.dimension())
    throw exception_classes::DifferentDimensionsException();
  // else

  double sumOfSquaredDifferences = 0.0;
  for (unsigned int i=0; i!=dimension(); ++i)
    sumOfSquaredDifferences += std::pow((coordinates_[i] - q[i]), 2);
  
  return std::sqrt(sumOfSquaredDifferences);
}


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
 *  - May throw a NotStrictlyPositivePointException exception if either the 
 *    current instance or the given point is not strictly positive.
 *
 *  \sa Point, Point::distance() and Point::dominates() 
 */
double 
Point::ratioDistance(const Point & q) const 
{
  if (isNull() or q.isNull())
    throw exception_classes::NullObjectException();
  if (dimension() != q.dimension())
    throw exception_classes::DifferentDimensionsException();
  if (not isStrictlyPositive() or not q.isStrictlyPositive()) 
    throw exception_classes::NotStrictlyPositivePointException();
  // else

  double max = 0.0;
  for (unsigned int i=0; i<dimension(); ++i) {
    double r = (q[i] - coordinates_[i]) / coordinates_[i];
    if (r > max) 
      max = r;
  }

  return max;
}


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
double 
Point::additiveDistance(const Point & q) const
{
  if (isNull() or q.isNull())
    throw exception_classes::NullObjectException();
  if (dimension() != q.dimension())
    throw exception_classes::DifferentDimensionsException();
  // else 

  double minEpsilon = 0.0;
  for (unsigned int i=0; i<dimension(); ++i) {
    double d = q[i] - coordinates_[i];
    if (d > minEpsilon)
      // in here only when d > minEpsilon and d > 0.0
      minEpsilon = d;
  }

  return minEpsilon;
}


//! Check if the current point (p) eps-dominates the given point (q).
/*!
 *  \param q A Point instance.
 *  \param eps An approximation parameter.
 *  \return true if p eps-covers q; false otherwise.
 *  
 *  There are two different definitions of eps-dominance we could use:
 *  - Additive eps-dominance. Where a point q is \f$\epsilon\f$ -dominated 
 *    by a point p if: 
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
bool 
Point::dominates(const Point & q, double eps) const
{
  return dominatesAdditive(q, eps);
}


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
bool 
Point::dominatesAdditive(const Point & q, double eps) const
{
  if (isNull() or q.isNull()) 
    throw exception_classes::NullObjectException();
  if (dimension() != q.dimension()) 
    throw exception_classes::DifferentDimensionsException();
  if (eps < 0.0) 
    throw exception_classes::NegativeApproximationRatioException();
  if (not isPositive() or not q.isPositive())
    throw exception_classes::NotPositivePointException();
  // else

  for (unsigned int i=0; i<dimension(); ++i)
    if (coordinates_[i] > q[i] + eps)
      return false;

  return true;
}


/*! 
 *  \brief Check if the current point (p) eps-dominates (in the 
 *         multiplicative sense) the given point (q).
 *  
 *  \param q A Point instance with \f$ q_{i} \ge 0 \f$ for all i.
 *  \param eps An approximation factor. (default 0.0)
 *  \return true if p eps-covers q, false otherwise.
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
 *  \sa Point, Point::dominates() and Point::distance()
 */
bool 
Point::dominatesMultiplicative(const Point & q, double eps) const
{
  if (isNull() or q.isNull()) 
    throw exception_classes::NullObjectException();
  if (dimension() != q.dimension()) 
    throw exception_classes::DifferentDimensionsException();
  if (eps < 0.0) 
    throw exception_classes::NegativeApproximationRatioException();
  if (not isStrictlyPositive() or not q.isStrictlyPositive()) 
    throw exception_classes::NotStrictlyPositivePointException();
  // else

  double r = 1+eps;
  for (unsigned int i=0; i<dimension(); ++i) 
    if (coordinates_[i] > r * q[i])
      return false;

  return true;
}


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
std::ostream & 
operator<< (std::ostream & out, const Point & p)
{
  return out << p.str(true);
}


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
std::istream & 
operator>> (std::istream & istr, Point & p)
{
  double coordinate = 0.0;
  std::string line;

  // Clear the point's coordinates.
  p.coordinates_.clear();

  // Read a line into a std::stringstream.
  std::getline(istr, line);
  std::stringstream lineStream(line);

  // Peek at the first character. 
  // - Sets the eofbit if the first character is EOF.
  lineStream.peek();
  // Read the coordinates from the std::stringstream.
  while (lineStream.good()) {
    lineStream >> coordinate;
    p.coordinates_.push_back(coordinate);
  }

  // If anything was read set the point to not null.
  // - it might have been null before
  if (p.coordinates_.size() > 0)
    p.isNull_ = false;

  return istr;
}


}  // namespace pareto_approximator


/* @} */
