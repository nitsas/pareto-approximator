/*! \file Hyperplane.cpp
 *  \brief The declaration of the Hyperplane class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <assert.h>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <cmath>
#include <armadillo>

#include "Hyperplane.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Constructor for a hyperplane on a 2D space. (a simple line)
/*!
 *  \param a1 The first dimension's (x_{1}) coefficient.
 *  \param a2 The second dimension's (x_{2}) coefficient.
 *  \param b The right hand side of the hyperplane equation.
 *  
 *  Constructs the hyperplane described by equation:
 *  \f$ a_{1} x_{1} + a_{2} x_{2} = b \f$
 *  
 *  \sa Hyperplane
 */
Hyperplane::Hyperplane(double a1, double a2, double b) : b_(b)
{
  assert(coefficients_.size() == 0);
  coefficients_.push_back(a1);
  coefficients_.push_back(a2);
}


//! Constructor for a hyperplane on a 3D space. (a simple plane)
/*!
 *  \param a1 The first dimension's (x_{1}) coefficient.
 *  \param a2 The second dimension's (x_{2}) coefficient.
 *  \param a3 The third dimension's (x_{3}) coefficient.
 *  \param b The right hand side of the hyperplane equation.
 *  
 *  Constructs the hyperplane described by equation:
 *  \f$ a_{1} x_{1} + a_{2} x_{2} + a_{3} x_{3} = b \f$
 *  
 *  \sa Hyperplane
 */
Hyperplane::Hyperplane(double a1, double a2, double a3, double b) : b_(b)
{
  assert(coefficients_.size() == 0);
  coefficients_.push_back(a1);
  coefficients_.push_back(a2);
  coefficients_.push_back(a3);
}


//! Constructor for a hyperplane on an n-dimensional space.
/*!
 *  \param first Iterator to the initial position in a std::vector<int>.
 *  \param last Iterator to the final position in a std::vector<int>.
 *  \param b The right hand side of the hyperplane equation.
 *  
 *  Constructs the hyperplane described by equation:
 *  \f$ a_{1} x_{1} + a_{2} x_{2} + ... + a_{n} x_{n} = b \f$
 *  
 *  The underlying vector should contain the hyperplane's coefficients, 
 *  ordered starting from a_{1} to a_{n}. The range used will be 
 *  [first, last), which includes all the elements between first and 
 *  last, including the element pointed by first but not the one pointed 
 *  by last.
 *  
 *  The resulting hyperplane's coefficients will be doubles, not ints.
 *  
 *  \sa Hyperplane
 */
Hyperplane::Hyperplane(std::vector<int>::const_iterator first, 
                       std::vector<int>::const_iterator last, int b) : b_(b)
{
  assert(coefficients_.size() == 0);
  coefficients_.assign(first, last);
}


//! Constructor for a hyperplane on an n-dimensional space.
/*!
 *  \param first Iterator to the initial position in a std::vector<double>.
 *  \param last Iterator to the final position in a std::vector<double>.
 *  \param b The right hand side of the hyperplane equation.
 *  
 *  Constructs the hyperplane described by equation:
 *  \f$ a_{1} x_{1} + a_{2} x_{2} + ... + a_{n} x_{n} = b \f$
 *  
 *  The underlying vector should contain the hyperplane's coefficients, 
 *  ordered starting from a_{1} to a_{n}. The range used will be 
 *  [first, last), which includes all the elements between first and 
 *  last, including the element pointed by first but not the one pointed 
 *  by last.
 *  
 *  \sa Hyperplane
 */
Hyperplane::Hyperplane(std::vector<double>::const_iterator first, 
                       std::vector<double>::const_iterator last, 
                       double b) : b_(b)
{
  assert(coefficients_.size() == 0);
  coefficients_.assign(first, last);
}


//! Constructor for a hyperplane on an n-dimensional space.
/*!
 *  \param first Iterator (pointer) to the initial position in an array 
 *               of int.
 *  \param last Iterator (pointer) to the final position in an array 
 *              of int. (the position just beyond the last element we want)
 *  \param b The right hand side of the hyperplane equation.
 *  
 *  Constructs the hyperplane described by equation:
 *  \f$ a_{1} x_{1} + a_{2} x_{2} + ... + a_{n} x_{n} = b \f$
 *  
 *  The underlying array should contain the hyperplane's coefficients, 
 *  ordered starting from a_{1} to a_{n}. The range used will be 
 *  [first, last), which includes all the elements between first and 
 *  last, including the element pointed by first but not the one pointed 
 *  by last.
 *  
 *  The resulting hyperplane's coefficients will be doubles, not ints.
 *  
 *  \sa Hyperplane
 */
Hyperplane::Hyperplane(const int * first, const int * last, int b) : b_(b)
{
  assert(coefficients_.size() == 0);
  coefficients_.assign(first, last);
}


//! Constructor for a hyperplane on an n-dimensional space.
/*!
 *  \param first Iterator (pointer) to the initial position in an array 
 *               of double.
 *  \param last Iterator (pointer) to the final position in an array of
 *              double. (the position just beyond the last element we want)
 *  \param b The right hand side of the hyperplane equation.
 *  
 *  Constructs the hyperplane described by equation:
 *  \f$ a_{1} x_{1} + a_{2} x_{2} + ... + a_{n} x_{n} = b \f$
 *  
 *  The underlying array should contain the hyperplane's coefficients, 
 *  ordered starting from a_{1} to a_{n}. The range used will be 
 *  [first, last), which includes all the elements between first and 
 *  last, including the element pointed by first but not the one pointed 
 *  by last.
 *  
 *  The resulting hyperplane's coefficients will be doubles, not ints.
 *  
 *  \sa Hyperplane
 */
Hyperplane::Hyperplane(const double * first, const double * last, 
                       double b) : b_(b)
{
  assert(coefficients_.size() == 0);
  coefficients_.assign(first, last);
}


//! Constructor for a hyperplane on a 2D space. (line)
/*!
 *  \param p1 A 2D Point instance.
 *  \param p2 A 2D Point instance.
 *  \return A 2-hyperplane (line) passing through both p1 and p2.
 *  
 *  Possible exceptions:
 *  - May throw a SamePointsException exception if p1 and p2 represent 
 *    the same point.
 *  - May throw a DifferentDimensionsException exception if a given 
 *    point's dimension is not 2.
 *  
 *  \sa Hyperplane, init() and Point
 */
Hyperplane::Hyperplane(const Point & p1, const Point & p2)
{
  if (p1 == p2)
    throw SamePointsException();

  std::set<Point> points;
  points.insert(p1);
  points.insert(p2);
  init(points);
}


//! Constructor for a hyperplane on a 3D space. (line)
/*!
 *  \param p1 A 3D Point instance.
 *  \param p2 A 3D Point instance.
 *  \param p3 A 3D Point instance.
 *  
 *  Constructs a 3-hyperplane (line) that passes through points p1, p2 
 *  and p3.
 *  
 *  Possible exceptions:
 *  - May throw a SamePointsException exception if two of the given Point 
 *    instances represent the same point.
 *  - May throw a DifferentDimensionsException exception if a given 
 *    point's dimension is not 3.
 *  
 *  \sa Hyperplane, init() and Point
 */
Hyperplane::Hyperplane(const Point & p1, const Point & p2, const Point & p3)
{
  std::set<Point> points;
  points.insert(p1);
  points.insert(p2);
  points.insert(p3);
  if (points.size() < 3)
    throw SamePointsException();

  init(points);
}

//! Constructor for a hyperplane on a 4D space. (line)
/*!
 *  \param p1 A 4D Point instance.
 *  \param p2 A 4D Point instance.
 *  \param p3 A 4D Point instance.
 *  \param p4 A 4D Point instance.
 *  
 *  Constructs a 4-hyperplane (line) that passes through points p1, p2, 
 *  p3 and p4.
 *  
 *  Possible exceptions:
 *  - May throw a SamePointsException exception if two of the given Point 
 *    instances represent the same point.
 *  - May throw a DifferentDimensionsException exception if a given 
 *    point's dimension is not 4.
 *  
 *  \sa Hyperplane, init() and Point
 */
Hyperplane::Hyperplane(const Point & p1, const Point & p2, 
                       const Point & p3, const Point & p4)
{
  std::set<Point> points;
  points.insert(p1);
  points.insert(p2);
  points.insert(p3);
  points.insert(p4);
  if (points.size() < 4)
    throw SamePointsException();

  init(points);
}


//! Constructor for a hyperplane on an n-dimensional space.
/*!
 *  \param first Iterator to the first element in a std::vector<Point>.
 *  \param last Iterator to the past-the-end element in a std::vector<Point>.
 *  
 *  Let n be the number of elements in the std::vector<Point> that first
 *  and last refer to.
 *  
 *  Constructs an n-hyperplane that passes through all the points in the 
 *  std::vector<Point> that first and last refer to.
 *  
 *  Possible exceptions:
 *  - May throw a SamePointsException exception if two of the given Point 
 *    instances represent the same point.
 *  - May throw a DifferentDimensionsException exception if a given 
 *    point's dimension is not n.
 *  
 *  \sa Hyperplane, init() and Point
 */
Hyperplane::Hyperplane(std::vector<Point>::const_iterator first, 
                       std::vector<Point>::const_iterator last)
{
  std::set<Point> points(first, last);
  if ( ((int) points.size()) < std::distance(first, last))
    throw SamePointsException();
  init(points);
}


//! Constructor for a hyperplane on an n-dimensional space.
/*!
 *  \param first Iterator (pointer) to the first element in an array 
 *               of Point instances.
 *  \param last Iterator (pointer) to the last element in an array of 
 *              Point instances.
 *  
 *  Let n be the number of elements in the array of Point instances 
 *  that first and last refer to.
 *  
 *  Constructs an n-hyperplane that passes through all the points in 
 *  array of Point instances that first and last refer to.
 *  
 *  Possible exceptions:
 *  - May throw a SamePointsException exception if two of the given Point 
 *    instances represent the same point.
 *  - May throw a DifferentDimensionsException exception if a given 
 *    point's dimension is not n.
 *  
 *  \sa Hyperplane, init() and Point
 */
Hyperplane::Hyperplane(const Point * first, const Point * last)
{
  std::set<Point> points(first, last);
  if ( ((int) points.size()) < std::distance(first, last))
    throw SamePointsException();
  init(points);
}


//! Initializer for a hyperplane on an n-dimensional space.
/*!
 *  \param points An std::set<Point> of Point instances.
 *  
 *  Let n be the number of elements in the std::set<Point> that "points"
 *  refers to.
 *  
 *  Initializes the current instance to an n-hyperplane that passes 
 *  through all the points in the std::set<Point> that "points" refers to.
 *  
 *  Possible exceptions:
 *  - May throw a DifferentDimensionsException exception if a given 
 *    point's dimension is not n.
 *  
 *  \sa Hyperplane and Point
 */
void 
Hyperplane::init(const std::set<Point> points)
{
  assert(coefficients_.size() == 0);

  unsigned int n = points.size();

  // check if all points have dimension n
  std::set<Point>::const_iterator pit;
  for (pit = points.begin(); pit != points.end(); ++pit)
    if (pit->dimension() != n)
      throw DifferentDimensionsException();
  
  // fill a matrix will each point's coordinates
  arma::mat M;
  for (pit = points.begin(); pit != points.end(); ++pit)
    M.insert_rows(M.n_rows, pit->toRowVec());
  // add a column of ones at the end (will make the following easier)
  M.insert_cols(M.n_cols, arma::ones<arma::vec>(n));
  
  // find the hyperplane's a_{i} coefficients
  for (unsigned int i = 0; i != n; ++i) {
    M.swap_cols(i, M.n_cols - 1);
    coefficients_.push_back(arma::det(M.cols(0, M.n_cols - 2)));
    M.swap_cols(i, M.n_cols - 1);
  }

  // find the hyperplane's b coefficient
  b_ = arma::det(M.cols(0, M.n_cols - 2));
}


//! A simple (and empty) destructor.
Hyperplane::~Hyperplane() { }


//! Access operator for the hyperplane's coefficients (except b).
/*!
 *  \param pos The position (in the coefficients vector) to access.
 *  \return The hyperplane's a_{pos+1} coefficient. 
 *  
 *  Remember that in the hyperplane equation we labeled the coefficients 
 *  a_{1}, a_{2}, ..., a_{n}. When accessing them though we'll be refering 
 *  to them starting from 0 to n-1 to comply with C++'s array notation. 
 *  So, for example to access the a_{1} coefficient of the myPlane Hyperplane 
 *  instance we'll have to use myPlane[0].
 *  
 *  - May throw a NonExistentCoefficientException if the requested 
 *    coefficient does not exist (pos is out of bounds).
 *  
 *  \sa Hyperplane
 */
double 
Hyperplane::operator[](unsigned int pos) const
{
  if (pos >= space_dimension())
    throw NonExistentCoefficientException();
  else
    return coefficients_[pos];
}


//! Access the instance's b coefficient. (equation's right hand side)
/*!
 *  Return the Hyperplane instance's b coefficient, the hyperplane 
 *  equation's right hand side.
 */
double 
Hyperplane::b() const
{
  return b_;
}


//! Return the dimension of the space the hyperplane lives in.
unsigned int 
Hyperplane::space_dimension() const
{
  return coefficients_.size();
}


//! Return iterator to beginning of the vector of a_{i} coefficients.
/*!
 *  \return An iterator referring to the first of the a_{i} coefficients.
 *  
 *  \sa Hyperplane
 */
Hyperplane::iterator 
Hyperplane::begin() 
{
  return coefficients_.begin();
}


//! Return const_iterator to beginning of the vector of a_{i} coefficients.
/*!
 *  \return A const_iterator referring to the first of the a_{i} 
 *          coefficients.
 *  
 *  \sa Hyperplane
 */
Hyperplane::const_iterator 
Hyperplane::begin() const 
{
  return coefficients_.begin();
}


//! Return iterator to end of the vector of a_{i} coefficients.
/*!
 *  \return An iterator pointing just after the last a_{i} coefficient.
 *  
 *  \sa Hyperplane
 */
Hyperplane::iterator 
Hyperplane::end() 
{
  return coefficients_.end();
}


//! Return const_iterator to end of the vector of a_{i} coefficients.
/*!
 *  \return A const_iterator pointing just after the last a_{i} 
 *          coefficient.
 *  
 *  \sa Hyperplane
 */
Hyperplane::const_iterator 
Hyperplane::end() const 
{
  return coefficients_.end();
}


//! Get the hyperplane's equation in a string.
/*!
 *  \return A std::string with the hyperplane's equation in parentheses.
 *  
 *  Examples:
 *  - ( 2.2 * x1 + 5 * x2 - 1.7 * x3 = 9.2 )
 *  - ( 1.3 * x1 - 6.7 * x2 + 0 * x3 + 0 * x4 - 1.0 * x5 = 10.1 )
 *  - ()         <--- null Hyperplane
 *  
 *  \sa Hyperplane
 */
std::string 
Hyperplane::str() const
{
  if (space_dimension() == 0)
    return "()";
  else {
    std::stringstream ss;

    ss << "( " << coefficients_[0] << " * x1";
    for (unsigned int i=1; i!=space_dimension(); ++i) {
      ss << ( coefficients_[i] >= 0 ? " + " : " - " );
      ss << fabs(coefficients_[i]) << " * x" << i+1;
    }
    ss << " = " << std::noshowpos << b_ << " )";

    return ss.str();
  }
}


/*!
 *  \brief Return the hyperplane's \f$ a_{i} \f$ coefficients as an 
 *         arma::vec (armadillo vector).
 */
arma::vec 
Hyperplane::toVec() const
{
  arma::vec coefficients(coefficients_.size());
  for (unsigned int i = 0; i != coefficients_.size(); ++i)
    coefficients(i) = coefficients_[i];

  return coefficients;
}


/*!
 *  \brief Return the hyperplane's \f$ a_{i} \f$ coefficients as an
 *         arma::rowvec (armadillo row vector).
 */
arma::rowvec 
Hyperplane::toRowVec() const
{
  arma::rowvec coefficients(coefficients_.size());
  for (unsigned int i = 0; i != coefficients_.size(); ++i)
    coefficients(i) = coefficients_[i];

  return coefficients;
}


//! The Hyperplane equality operator.
/*!
 *  \param hyperplane A Hyperplane instance we want to compare with 
 *                    the current instance.
 *  \return true if the two instances represent the same hyperplane, 
 *          false otherwise.
 *  
 *  We should check that the two hyperplanes' coefficients are the same 
 *  but the problem is that one hyperplane's coefficients might be scaled.
 *  We overcome this problem by multiplying each hyperplane's a_{i} 
 *  coefficients with the other hyperplane's b when comparing. If the 
 *  hyperplanes are different at least one pair of (scaled) coefficients 
 *  will not match.
 *  
 *  \sa Hyperplane and operator!=()
 */
bool 
Hyperplane::operator== (const Hyperplane & hyperplane) const
{
  if (space_dimension() != hyperplane.space_dimension())
    return false;
  // else
  for (unsigned int i=0; i!=space_dimension(); ++i)
    if (coefficients_[i] * hyperplane.b() != hyperplane[i] * b_)
      return false;

  return true;
}


//! The Hyperplane inequality operator.
/*!
 *  \param hyperplane A Hyperplane instance we want to compare with 
 *                    the current instance.
 *  \return true if the two instances represent different hyperplanes, 
 *          false otherwise.
 *  
 *  \sa Hyperplane and operator==()
 */
bool 
Hyperplane::operator!= (const Hyperplane & hyperplane) const
{
  if (space_dimension() != hyperplane.space_dimension())
    return true;
  // else
  for (unsigned int i=0; i!=space_dimension(); ++i)
    if (coefficients_[i] * hyperplane.b() != hyperplane[i] * b_)
      return true;

  return false;
}


//! The Hyperplane output stream operator.
/*!
 *  Output the hyperplane's equation in parentheses.
 *  
 *  Examples:
 *  - ( 2.2 x1 + 5.4 x2 - 1.7 x3 = 9.2 )
 *  - ( 1.3 x1 - 6.7 x2 + 0.0 x3 + 0.0 x4 - 1.0 x5 = 10.1 )
 *  
 *  \sa Hyperplane and str()
 */
std::ostream & 
operator<< (std::ostream & out, const Hyperplane & hyperplane)
{
  return out << hyperplane.str();
}


//! Compute the ratio distance from the given Point to the hyperplane.
/*!
 *  \param p A Point instance.
 *  \return The ratio distance from p to the hyperplane.
 *  
 *  The ratio distance from a point p to a hyperplane H is defined as:
 *  \f$ RD(p, H) = \min_{q \in H} RD(p, q) \f$, where q is a point on H.
 *  The ratio distance from a point p to a point q is defined as:
 *  \f$ RD(p, q) = \max\{ \max_{i}\{(q_{i} - p_{i}) / p_{i}\}, 0.0 \} \f$.
 *  
 *  Intuitively it is the minimum value of \f$ \epsilon \ge 0 \f$ such 
 *  that some point on H \f$ \epsilon -covers p \f$.
 *  
 *  Possible exceptions:
 *  - May throw a DifferentDimensionsException exception if the given point 
 *    and the hyperplane belong in spaces of different dimensions.
 *  
 *  \sa Hyperplane and Point
 */
double 
Hyperplane::ratioDistance(const Point & p) const
{
  if (space_dimension() != p.dimension())
    throw DifferentDimensionsException();
  // else
  double dotProduct = 0.0;
  for (unsigned int i=0; i!=space_dimension(); ++i) 
    dotProduct += coefficients_[i] * p[i];

  return std::max( (b_ - dotProduct) / dotProduct, 0.0 );
}


//! Create a new Hyperplane parallel to the current one (through a point).
/*!
 *  \param p A Point instance through which the new Hyperplane instance 
 *           must pass.
 *  \return A new Hyperplane instance, parallel to the current one and 
 *          passing through p.
 *  
 *  The new hyperplane will have the same a_{i} coefficients but a different 
 *  b coefficient, one that satisfies the equation:
 *  \f$ a_{1} * p_{1} + a_{2} * p_{2} + ... + a_{n} * p_{n} = b \f$.
 *  
 *  \sa Hyperplane and Point
 */
Hyperplane 
Hyperplane::parallelThrough(const Point & p) const
{
  double newB = 0.0;
  for (unsigned int i=0; i!=space_dimension(); ++i)
    newB += coefficients_[i] * p[i];

  return Hyperplane(coefficients_.begin(), coefficients_.end(), newB);
}


//! Check if two hyperplanes are parallel.
/*!
 *  \param hyperplane A Hyperplane instance.
 *  \return true if the given hyperplane instance is parallel to the current 
 *          one; false otherwise.
 *  
 *  We should check if the two instances have the same a_{i} coefficients. The 
 *  only problem is that the a_{i} and b coefficients of one instance might 
 *  be scaled:
 *  Let h1 and h2 be two parallel hyperplanes. Scale h2's a_{i}'s and b by a 
 *  constant c. Its slope doesn't change, that is h1 and h2 remain parallel.
 *  
 *  To overcome this problem, we scale each hyperplane's coefficients by the 
 *  opposite hyperplane's a_{1} and expect equal a_{i}'s.
 *  
 *  \sa Hyperplane
 */
bool 
Hyperplane::isParallel(const Hyperplane & hyperplane) const
{
  if (space_dimension() != hyperplane.space_dimension())
    return false;
  for (unsigned int i=0; i!=space_dimension(); ++i) 
    if (coefficients_[i] * hyperplane[0] != hyperplane[i] * coefficients_[0])
      return false;

  return true;
}


//! Reverse the sign of all of the Hyperplane's coefficients.
/*!
 *  Reverse the sign of all the \f$ a_{i} \f$ and b coefficients.
 *  Since we're reversing the sign on both sides of the hyperplane 
 *  equation, the equation stays the same.
 *
 *  Reminder:
 *  A hyperplane on an n-dimensional space can be described by an 
 *  equation of the form:
 *  \f$ a_{1} x_{1} + a_{2} x_{2} + ... + a_{n} x_{n} = b \f$
 *  
 *  \sa Hyperplane
 */
void 
Hyperplane::reverseCoefficientSigns()
{
  for (unsigned int i = 0; i != space_dimension(); ++i)
    coefficients_[i] = -coefficients_[i];
  b_ = -b_;
}


//! Make sure the hyperplane faces away from the given point.
/*!
 *  \param p A Point instance to face away from.
 *  
 *  Reminder:
 *  A hyperplane on an n-dimensional space can be described by an 
 *  equation of the form:
 *  \f$ a_{1} x_{1} + a_{2} x_{2} + ... + a_{n} x_{n} = b \f$
 *
 *  Make sure using the current instance's \f$ a_{i} \f$ coefficients 
 *  as weights in a weighted sum and minimizing "combs" away from p.
 *  
 *  The problem:
 *  If we reverse the sign of all \f$ a_{i} \f$ and b coefficients
 *  of the hyperplane the hyperplane equation doesn't change. The points
 *  that satisfy the equation stay the same. On first glance it seems 
 *  that it doesn't matter which version of the hyperplane equation 
 *  we'll use. 
 *  
 *  But when using the hyperplane's \f$ a_{i} \f$ coefficients as the 
 *  weights of a weighted sum (let's call it WS) their sign matters. 
 *  Minimizing WS graphically means we are "combing" the space with 
 *  the hyperplane i.e. moving the hyperplane in one of the two 
 *  perpendicular to the hyperplane directions. Reversing all 
 *  \f$ a_{i}'s \f$ signs means reversing the direction the hyperplane 
 *  moves ("combs") towards while minimizing.
 *
 *  We'll make sure that minimizing WS will "comb" away from p.
 *  
 *  Possible exceptions:
 *  - May throw a PointIsOnTheHyperplaneException exception if the given 
 *    point satisfies the hyperplane equation. (can't face away from it)
 *  
 *  \sa Hyperplane
 */
void 
Hyperplane::faceAwayFrom(const Point & p)
{
  double ws = arma::as_scalar(this->toRowVec() * p.toVec());
  if (ws == this->b())
    throw PointIsOnTheHyperplaneException();
  // else
  if (ws < this->b())
    this->reverseCoefficientSigns();
}


}  // namespace pareto_approximator


/* @} */
