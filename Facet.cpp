/*! \file Facet.cpp
 *  \brief The implementation of the Facet class.
 *  \author Christos Nitsas
 *  \date 2012
 *  
 *  Won't `include` Facet.h. In fact Facet.h will `include` 
 *  Facet.cpp because it describes a class template (which doesn't allow 
 *  us to split declaration from definition).
 */


#include <assert.h>
#include <iterator>
#include <algorithm>
#include <armadillo>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
namespace pareto_approximator {


//! A Facet constructor.
/*!
 *  \param firstVertex An iterator to the first of the facet vertices.
 *  \param lastVertex An iterator to the past-the-end element in the 
 *                    container of facet vertices.
 *  \param preferPositiveNormalVector While computing the facet's normal 
 *                    vector prefer the all-positive one (if it exists).
 *  
 *  Vertices cannot be null PointAndSolution<S> instances and the 
 *  contained points cannot be null Point instances. 
 *  
 *  Initializes:
 *  - Facet<S>::vertices_ to the sequence of vertices pointed to by 
 *    firstVertex and lastVertex. The range used is 
 *    [firstVertex, lastVertex).
 *  - Calculates the hyperplane passing through the facet's vertices 
 *    and uses its normal vector as the facet's normal vector 
 *    (Facet<S>::normal_). For each set of n vertices there are two 
 *    different n-hyperplanes passing through them with opposite normal 
 *    vectors. This constructor will prefer the all-positive normal 
 *    vector (if one exists) if preferPositiveNormalVector is set to 
 *    true; otherwise it will choose one depending on the order of the 
 *    vertices.
 *  - Facet<S>::approximationErrorUpperBound_ to the distance between 
 *    the Facet and its Lower Distal Point (LDP). Calculates both the 
 *    LDP and the distance.
 *  
 *  Possible exceptions:
 *  - May throw a DifferentDimensionsException exception if the given 
 *    vertices have different dimensions.
 *  - May throw a NullObjectException if some of the given vertices or 
 *    some of the points contained in them are null instances (null 
 *    PointAndSolution<S> and null Point instances respectively).
 *  
 *  \sa Facet
 */
template <class S> 
Facet<S>::Facet(
      typename std::vector< PointAndSolution<S> >::const_iterator firstVertex, 
      typename std::vector< PointAndSolution<S> >::const_iterator lastVertex, 
      bool preferPositiveNormalVector)
{
  spaceDimension_ = firstVertex->dimension();

  // Only accept simplicial facets for now. 
  // - A facet is a simplicial facet if it consists of exactly d vertices, 
  //   where d is the dimension of the space that the facet lives in.
  assert(std::distance(firstVertex, lastVertex) == spaceDimension_);

  // Make sure that all the given vertices are valid.
  // - They all have the correct dimension. (i.e. the dimension of the space 
  //   that the facet lives in)
  // - They (and the points they contain) are not null instances.
  ConstVertexIterator cvi;
  for (cvi = firstVertex; cvi != lastVertex; ++cvi) {
    if (cvi->isNull() or cvi->point.isNull())
      throw exception_classes::NullObjectException();
    if (cvi->dimension() != spaceDimension_)
      throw exception_classes::DifferentDimensionsException();
  }

  // First fill-in Facet<S>::vertices_.
  vertices_.assign(firstVertex, lastVertex);

  // Compute and set the facet's normal vector (Facet<S>::normal_).
  computeAndSetFacetNormal(preferPositiveNormalVector);

  // Compute and set the facet's offset.
  b_ = arma::dot( arma::vec(normal_), vertices_[0].point.toVec() );

  // Compute and set the facet's localApproximationErrorUpperBound_ and 
  // isBoundaryFacet_ attributes.
  computeAndSetLocalApproximationErrorUpperBoundAndIsBoundaryFacet();
}


//! A Facet constructor.
/*!
 *  \param firstVertex An iterator to the first of the facet vertices.
 *  \param lastVertex An iterator to the past-the-end element in the 
 *                    container of facet vertices.
 *  \param firstElemOfFacetNormal An iterator to the first element of 
 *                                the facet's normal.
 *  \param lastElemOfFacetNormal An iterator to the past-the-end 
 *                               element in the container of facet 
 *                               vertices.
 *  
 *  Vertices cannot be null PointAndSolution<S> instances and the 
 *  contained points cannot be null Point instances. 
 *  
 *  Initializes:
 *  - Facet<S>::vertices_ to the sequence of vertices pointed to by 
 *    firstVertex and lastVertex. The range used is 
 *    [firstVertex, lastVertex).
 *  - Facet<S>::normal_ to the sequence of elements pointed to by 
 *    firstElemOfFacetNormal and lastElemOfFacetNormal. The range 
 *    used is [firstElemOfFacetNormal, lastElemOfFacetNormal).
 *  - Facet<S>::approximationErrorUpperBound_ to the distance between 
 *    the Facet and its Lower Distal Point (LDP). Calculates both the 
 *    LDP and the distance.
 *  
 *  We do not verify that the given normal vector is indeed correct, 
 *  i.e. it agrees with the given vertices. The responsibility lies with 
 *  the caller.
 *  
 *  Possible exceptions:
 *  - May throw a DifferentDimensionsException exception if the given 
 *    vertices have different dimensions.
 *  - May throw a NullObjectException if some of the given vertices or 
 *    some of the points contained in them are null instances (null 
 *    PointAndSolution<S> and null Point instances respectively).
 *  
 *  \sa Facet
 */
template <class S> 
Facet<S>::Facet(
      typename std::vector< PointAndSolution<S> >::const_iterator firstVertex, 
      typename std::vector< PointAndSolution<S> >::const_iterator lastVertex, 
      std::vector<double>::const_iterator firstElemOfFacetNormal, 
      std::vector<double>::const_iterator lastElemOfFacetNormal)
{
  spaceDimension_ = std::distance(firstElemOfFacetNormal, 
                                  lastElemOfFacetNormal);

  // Only accept simplicial facets for now. 
  // - A facet is a simplicial facet if it consists of exactly d vertices, 
  //   where d is the dimension of the space that the facet lives in.
  assert(std::distance(firstVertex, lastVertex) == spaceDimension_);

  // Make sure that all the given vertices are valid.
  // - They all have the correct dimension. (i.e. the dimension of the space 
  //   that the facet lives in)
  // - They (and the points they contain) are not null instances.
  ConstVertexIterator cvi;
  for (cvi = firstVertex; cvi != lastVertex; ++cvi) {
    if (cvi->isNull() or cvi->point.isNull())
      throw exception_classes::NullObjectException();
    if (cvi->dimension() != spaceDimension_)
      throw exception_classes::DifferentDimensionsException();
  }

  // First fill-in Facet<S>::vertices_ and Facet<S>::normal_.
  vertices_.assign(firstVertex, lastVertex);
  normal_.assign(firstElemOfFacetNormal, lastElemOfFacetNormal);

  // Compute and set the facet's offset.
  b_ = arma::dot( arma::vec(normal_), vertices_[0].point.toVec() );

  // Compute and set the facet's localApproximationErrorUpperBound_ and 
  // isBoundaryFacet_ attributes.
  computeAndSetLocalApproximationErrorUpperBoundAndIsBoundaryFacet();
}


//! A simple (and empty) destructor.
template <class S> 
Facet<S>::~Facet() { }


//! Return the dimension of the space that the facet lives in.
/*!
 *  \sa Facet
 */
template <class S> 
unsigned int 
Facet<S>::spaceDimension() const
{
  return spaceDimension_;
}


//! Is the facet a boundary facet?
/*!
 *  \return true if the facet is a boundary facet; false otherwise.
 *  
 *  We call a facet a boundary facet if it does not have a Lower Distal 
 *  Point. (The hyperplanes h_{i} associated with its vertices v_{i} 
 *  do not intersect in a unique point. Check the documentation of Facet 
 *  for more info.)
 *  
 *  \sa Facet and Facet<S>::computeLowerDistalPoint()
 */
template <class S> 
bool 
Facet<S>::isBoundaryFacet() const
{
  return isBoundaryFacet_;
}


//! Get an upper bound to the current facet's approximation error.
/*! 
 *  \return the localApproximationErrorUpperBound attribute
 *  
 *  We will use the distance from the facet to it's Lower Distal Point 
 *  as an upper bound to the local approximation error.
 *  
 *  Check the documentation for Facet for a description of what a Lower 
 *  Distal Point is.
 *  
 *  Possible exceptions:
 *  - May throw a BoundaryFacetException if the facet is a 
 *    boundary facet.
 *  
 *  \sa Facet and Facet<S>::computeLowerDistalPoint()
 */
template <class S> 
double 
Facet<S>::getLocalApproximationErrorUpperBound() const
{
  if (isBoundaryFacet())
    throw exception_classes::BoundaryFacetException();
  // else

  return localApproximationErrorUpperBound_;
}


//! Return iterator to the beginning of the vector of facet vertices.
/*! 
 *  \return An iterator pointing to the first vertex in the vector 
 *          of vertices.
 *  
 *  \sa Facet
 */
template <class S> 
typename Facet<S>::ConstVertexIterator 
Facet<S>::beginVertex() const
{
  return vertices_.begin();
}


//! Return iterator to the end of the vector of facet vertices.
/*! 
 *  \return An iterator pointing just after the last vertex in the 
 *          vector of vertices.
 *  
 *  \sa Facet
 */
template <class S> 
typename Facet<S>::ConstVertexIterator 
Facet<S>::endVertex() const
{
  return vertices_.end();
}


//! Return iterator to the beginning of the facet's normal vector.
/*! 
 *  \return An iterator pointing to the first element in the facet's 
 *          normal vector.
 *  
 *  \sa Facet
 */
template <class S> 
typename Facet<S>::ConstFacetNormalIterator 
Facet<S>::beginFacetNormal() const
{
  return normal_.begin();
}


//! Return iterator to the end of the facet's normal vector.
/*! 
 *  \return An iterator pointing just after the last element in the 
 *          facet's normal vector.
 *  
 *  \sa Facet
 */
template <class S> 
typename Facet<S>::ConstFacetNormalIterator 
Facet<S>::endFacetNormal() const
{
  return normal_.end();
}


//! Compute the b coefficient of the hyperplane on which the facet lies.
/*! 
 *  \return The b coefficient of the hyperplane on which the facet lies 
 *          (the one with the same normal vector as the facet).
 *  
 *  We describe a hyperplane using the following linear equation:
 *  \f$ a_{1} \cdot x_{1} + a_{2} \cdot x_{2} + ... + a_{n} \cdot x_{n} 
 *      = b \f$, where \f$ \mathbf{a} = [a_{1} a_{2} ... a_{n}]^{T} \f$ 
 *  is the hyperplane's normal vector and b is the coefficient this 
 *  function will compute.
 *  
 *  b can be found by taking the dot product of the hyperplane's (or 
 *  facet's) normal vector with any point on the hyperplane (or facet).
 *  
 *  \sa Facet
 */
template <class S> 
double 
Facet<S>::b() const
{
  return b_;
}


/*!
 * \brief Compute the mean of all the weight vectors of the facet's 
 *        vertices.
 *
 *  \return A weight vector W of size this->spaceDimension(). Each 
 *          element W_{j} is the mean of all w_{ij}'s, where w_{i} is 
 *          the weight vector inside the facet's i'th vertex's 
 *          PointAndSolution object (PointAndSolution::weightsUsed).
 *  
 *  \sa Facet
 */
template <class S> 
std::vector<double> 
Facet<S>::computeMeanVertexWeights() const
{
  ConstVertexIterator cvi;
  std::vector<double> meanWeights(spaceDimension(), 0.0);
  for (unsigned int i = 0; i != spaceDimension(); ++i) {
    for (cvi = beginVertex(); cvi != endVertex(); ++cvi)
      meanWeights[i] += cvi->weightsUsed[i];
    meanWeights[i] /= spaceDimension();
  }

  return meanWeights;
}


//! Compute the facet's Lower Distal Point (LDP).
/*! 
 *  \return The facet's Lower Distal Point (Point instance) if one 
 *          exists, a null Point instance (i.e. a Point instance whose 
 *          Point::isNull() method returns true) otherwise. 
 *  
 *  Solves (for x) the system of k equations of the form:
 *  \f$ w_{i1} * x_{1} + ... + w_{ik} * x_{k} = w_{i} \dot v_{i} \f$, 
 *  where w_{i} is the weight vector associated with the i'th vertex of 
 *  the facet (the normal of the associated lower-bound hyperplane) and 
 *  v_{i} is the vector of the i'th vectex's coordinates. 
 *  
 *  (\f$ w_{i} \dot v_{i} = b_{i} \f$ is the associated hyperplane's offset 
 *  from the origin)
 *  
 *  The solution, if one (and only one) exists, will be the LDP's coordinates.
 *  If a unique solution does not exist we return a null Point instance 
 *  (i.e. a Point instance whose Point::isNull() method returns true).
 *  
 *  What is the Lower Distal Point (LDP)? Recall that hyperplanes through 
 *  the Pareto Set points with normal vectors equal to the weight vector 
 *  yielding that point are lower bounds of the Pareto Set. 
 *  Let's call the facet's vertices v_{i} and the hyperplane assosiated 
 *  with each vertex h_{i}.
 *  The LDP is the point where the current facet's h_{i}'s intersect, 
 *  provided these N hyperplanes intersect in a unique point.
 *  Intuitively, the LDP is the most distant possible point we might 
 *  generate using the current facet's normal as weights.
 *  
 *  LDP Example:
 *  The LDP in 3 dimensions is the top of the pyramid whose base is the 
 *  current facet. Each of the pyramid's sides lies on a hyperplane 
 *  h_{i}, where h_{i} is the hyperplane used to find the facet's vertex 
 *  v_{i}. (We can recreate h_{i} using v_{i}'s weight vector.)
 *  
 *  The h_{i}'s might not intersect in a unique point. In that case, 
 *  this method returns a null Point instance (i.e. a Point instance 
 *  Point::isNull() method returns true) and the current facet is 
 *  treated as a boundary facet.
 *  
 *  \sa Facet
 */
template <class S> 
Point 
Facet<S>::computeLowerDistalPoint() const
{
  // open a stream to /dev/null (will redirect error messages there)
  // - will redirect error messages there
  std::ofstream f("/dev/null");
  // redirect armadillo error messages to /dev/null
  // - e.g. when arma::solve() finds no solutions
  // - we check arma::solve()'s return value for errors, 
  //   no need for error messages
  arma::set_stream_err2(f);

  arma::mat W;
  arma::vec b;

  ConstVertexIterator cvi;
  // fill in matrix W and vector b
  for (cvi = beginVertex(); cvi != endVertex(); ++cvi) {
    // make sure the weightsUsed field of the current vertex is not empty
    assert(cvi->weightsUsed.size() == spaceDimension());

    arma::rowvec wi(cvi->weightsUsed);
    // fill in row i of the weight (hyperplane-normal) matrix
    W.insert_rows(W.n_rows, wi);
    // fill in element i of the hyperplane-offsets vector
    b.insert_rows(b.n_rows, wi * cvi->point.toVec());
  }
  
  arma::vec x;
  bool hasSolution = arma::solve(x, W, b);
  if (hasSolution) 
    // unique solution
    // - return it as a Point instance
    return Point(x.begin(), x.end());
  else 
    // either no solution or an infinite number of solutions 
    // - return a null Point instance 
    return Point();
}


/*!
 *  \brief Compute a Point instance's distance from the hyperplane 
 *         that the facet lies on.
 *
 *  \param p A Point instance. (should be strictly positive)
 *  \return The distance from p to the hyperplane on which the facet 
 *          lies.
 *  
 *  There are different possible distance metrics we could use (e.g. 
 *  ratio distance, Euclidean distance, additive distance etc.). We use 
 *  the additive distance metric for now.
 *  
 *  Possible exceptions:
 *  - May throw a DifferentDimensionsException exception if the given 
 *    point and the hyperplane belong in spaces of different dimensions.
 *  - May throw an InfiniteRatioDistanceException exception if the given 
 *    point's coordinate vector is perpendicular to the facet's 
 *    normal vector. Multiplying the point by a constant moves it in 
 *    a direction parallel to the facet's supporting hyperplane.
 *  - May throw a NotPositivePointException (or 
 *    NotStrictlyPositivePointException if we are using the multiplicative 
 *    error measure) exception if the given point is not positive (not 
 *    strictly positive, respectively).
 *  - May throw a NullObjectException exception if the given Point 
 *    instance is a null Point instance.
 *  
 *  \sa Point and Facet
 */
template <class S> 
double 
Facet<S>::distance(const Point & p) const
{
  return additiveDistance(p);
}


/*!
 *  \brief Compute the Euclidean distance from the given point to the 
 *         hyperpane on which the facet lies.
 * 
 *  \param p A Point instance.
 *  \return The Euclidean distance from p to the hyperplane on which 
 *          the facet lies. (supporting hyperplane)
 *  
 *  The formula for the Euclidean distance between a d-dimensional point 
 *  p and a d-dimensional facet F with a supporting hyperplane H which 
 *  has normal \f$\mathbf{n}\f$ and is described by the equation 
 *  \f$ \mathbf{n} \dot \mathbf{x} = c \f$ is:
 *  \f$ ED(p, F) = \left|
 *      \frac{ \mathbf{n} \dot \mathbf{p} - c }{ ||\mathbf{n}|| } 
 *      \right| \f$
 *  
 *  \sa Facet and Point
 */
template <class S> 
double 
Facet<S>::euclideanDistance(const Point & p) const
{
  if (p.isNull())
    throw exception_classes::NullObjectException();
  if (spaceDimension() != p.dimension())
    throw exception_classes::DifferentDimensionsException();

  double normOfNormalVector = arma::norm(arma::vec(normal_), 2);

  return std::abs( (arma::dot(arma::vec(normal_), p.toVec()) - b()) / 
                   normOfNormalVector );
}


/*!
 *  \brief Compute a point's ratio distance from the hyperplane that 
 *         the facet lies on.
 *  
 *  \param p A Point instance. (stricty positive)
 *  \return The point's ratio distance from the hyperplane on which the 
 *          facet lies.
 *  
 *  The ratio distance from a point p to a hyperplane H is defined as:
 *  \f$ RD(p, H) = \min_{q \in H} RD(p, q) \f$, where q is a point on H.
 *  The ratio distance from a point p to a point q is defined as:
 *  \f$ RD(p, q) = \max\{ \max_{i}\{(q_{i} - p_{i}) / p_{i}\}, 0.0 \} \f$.
 *  
 *  Intuitively it is the minimum value of \f$ \epsilon \ge 0 \f$ such 
 *  that some point on H \f$\epsilon\f$ -dominates (\f$\epsilon\f$ -covers) 
 *  p in the multiplicative sense.
 *  
 *  In order for the ratio distance to make sense point p must be 
 *  strictly positive, i.e. \f$ p_{i} > 0.0 \f$ must hold for all i.
 *  
 *  Possible exceptions:
 *  - May throw a DifferentDimensionsException exception if the given point 
 *    and the hyperplane belong in spaces of different dimensions.
 *  - May throw an InfiniteRatioDistanceException exception if the given 
 *    point's coordinate vector is perpendicular to the facet's 
 *    normal vector. Multiplying the point by a constant moves it in 
 *    a direction parallel to the facet's supporting hyperplane.
 *  - May throw a NotStrictlyPositivePointException exception if the 
 *    given point is not strictly positive.
 *  - May throw a NullObjectException exception if the given Point 
 *    instance is a null Point instance.
 *  
 *  \sa Point and Facet
 */
template <class S> 
double 
Facet<S>::ratioDistance(const Point & p) const
{
  if (p.isNull())
    throw exception_classes::NullObjectException();
  if (spaceDimension() != p.dimension())
    throw exception_classes::DifferentDimensionsException();
  if (not p.isStrictlyPositive()) 
    throw exception_classes::NotStrictlyPositivePointException();
  // else

  assert(spaceDimension() > 0);

  double dotProduct  = 0.0;
  double facetOffset = b();      // the facet's offset from the origin
  for (unsigned int i=0; i!=spaceDimension(); ++i) {
    dotProduct  += normal_[i] * p[i];
  }

  double result;
  if (dotProduct == facetOffset)
    // the point is on the facet
    // it's okay even if dotProduct == 0.0
    result = 0.0;
  else if (dotProduct == 0.0)
    // multiplying the point by a constant moves it in a direction 
    // parallel to the hyperplane
    throw exception_classes::InfiniteRatioDistanceException();
  else
    result = std::max( (facetOffset - dotProduct) / dotProduct, 0.0 );

  return result;
}


/*!
 *  \brief Compute a point's additive distance from the hyperplane 
 *         that the facet lies on.
 *  
 *  \param p A Point instance.
 *  \return The point's additive distance from the hyperplane that the 
 *          facet lies on. (i.e. the minimum value of \f$\epsilon\f$ 
 *          such that the hyperplane dominates the point in the additive 
 *          sense)
 *  
 *  The additive distance from a point p to a hyperplane H is defined as:
 *  \f$ AD(p, H) = \min_{q \in H} AD(p, q) \f$, where q is a point on H.
 *  The additive distance from a point p to a point q is defined as:
 *  \f$ AD(p, q) = \max\{ \max_{i}\{(q_{i} - p_{i})\}, 0.0 \} \f$.
 *  
 *  Intuitively it is the minimum value of \f$ \epsilon \ge 0 \f$ such 
 *  that some point on H \f$\epsilon\f$ -dominates (\f$\epsilon\f$ -covers) 
 *  p in the additive sense. To calculate it we take advantage of the fact 
 *  that the point (p + \f$\epsilon\f$) will be lying on H.
 *  
 *  Possible exceptions:
 *  - May throw a NullObjectException exception if the given Point 
 *    instance is a null Point instance.
 *  - May throw a DifferentDimensionsException exception if the given 
 *    point and the hyperplane belong in spaces of different dimensions.
 *  - May throw a NotPositivePointException exception if the given point 
 *    is not positive. (i.e. some coordinate is less than 0.0)
 *
 *  \sa Point and Facet
 */
template <class S> 
double 
Facet<S>::additiveDistance(const Point & p) const
{
  if (p.isNull())
    throw exception_classes::NullObjectException();
  if (spaceDimension() != p.dimension())
    throw exception_classes::DifferentDimensionsException();
  if (not p.isPositive())
    throw exception_classes::NotPositivePointException();
  // else

  assert(spaceDimension() > 0);

  double sumOfFacetNormal = 0.0;
  double dotProduct       = 0.0;
  for (unsigned int i=0; i!=spaceDimension(); ++i) {
    sumOfFacetNormal += normal_[i];
    dotProduct       += normal_[i] * p[i];
  }

  // - To calculate the result we take advantage of the fact that point 
  //   (p + \f$\epsilon\f$) will be lying on H (the hyperplane).
  // - sumOfFacetNormal should not be 0.0 - it can only be 0.0 if 
  //   the facet's normal vector is all zero (not a valid facet)
  //   (we assume the facet has an all-positive normal vector; 
  //   if not, there is no point in calling this method)
  assert(sumOfFacetNormal != 0.0);
  return std::max( (b() - dotProduct) / sumOfFacetNormal, 0.0 );
}


//! Check if the Facet approximately dominates the given point.
/*!
 *  \param p A Point instance. (must be positive if we are using the 
 *           additive error measure; strictly positive if we are using 
 *           the multiplicative)
 *  \param eps The approximation factor.
 *  \return true if some point on the facet's supporting hyperplane 
 *          approximately dominates the given point; false otherwise
 *  
 *  There are two different definitions of approximate dominance 
 *  (\f$\epsilon\f$ -dominance) we could use:
 *  - Additive \f$\epsilon\f$ -dominance. Where a point q is 
 *    \f$\epsilon\f$ -dominated by a point p if: 
 *    \f$ p_{i} \le q_{i} + \epsilon \f$ for all i.
 *  - Multiplicative \f$\epsilon\f$ -dominance. Where a point q is 
 *    \f$\epsilon\f$ -dominated by a point p if: 
 *    \f$ p_{i} \le (1 + \epsilon) q_{i} \f$ for all i.
 *  
 *  This method checks if some point on the facet's supporting 
 *  hyperplane H (the hyperplane on which the facet lies, that has the 
 *  same normal vector as the facet) approximately dominates the given 
 *  point.
 *  
 *  Currently using the additive error measure.
 *  
 *  \sa Facet, Point, Point::dominates(), Facet::dominatesAdditive() 
 *      and Facet::dominatesMultiplicative()
 */
template <class S> 
bool 
Facet<S>::dominates(const Point & p, double eps) const
{
  return dominatesAdditive(p, eps);
}


/*! 
 *  \brief Check if the Facet approximately dominates (in the additive 
 *         sense) the given point.
 *  
 *  \param p A Point instance. (must be positive - i.e. all coordinates 
 *           greater than or equal to 0.0)
 *  \param eps The approximation factor.
 *  \return true if some point on the facet's supporting hyperplane 
 *          approximately dominates the given point in the additive 
 *          sense; false otherwise
 *  
 *  Possible exceptions:
 *  - May throw a NullObjectException exception if the given Point 
 *    instance is a null Point instance.
 *  - May throw a DifferentDimensionsException exception if the given 
 *    point and the hyperplane belong in spaces of different dimensions.
 *  - May throw a NegativeApproximationRatioException exception if the 
 *    given approximation ratio/factor/threshold is less than 0.0.
 *  - May throw a NotPositivePointException exception if the given point 
 *    is not positive. (i.e. some coordinate is less than 0.0)
 *  
 *  \sa Facet, Point, Point::dominatesAdditive() and 
 *      Facet::dominates()
 */
template <class S> 
bool 
Facet<S>::dominatesAdditive(const Point & p, double eps) const
{
  if (p.isNull())
    throw exception_classes::NullObjectException();
  if (spaceDimension() != p.dimension())
    throw exception_classes::DifferentDimensionsException();
  if (eps < 0.0)
    throw exception_classes::NegativeApproximationRatioException();
  if (not p.isPositive())
    throw exception_classes::NotPositivePointException();
  // else

  if (additiveDistance(p) <= eps)
    return true;
  else
    return false;
}


/*! 
 *  \brief Check if the Facet approximately dominates (in the 
 *         multiplicative sense) the given point.
 *  
 *  \param p A Point instance. (must be strictly positive - i.e. all 
 *           coordinates greater than 0.0)
 *  \param eps The approximation factor.
 *  \return true if some point on the facet's supporting hyperplane 
 *          approximately dominates the given point in the multiplicative 
 *          sense; false otherwise
 *  
 *  Possible exceptions:
 *  - May throw a NullObjectException exception if the given Point 
 *    instance is a null Point instance.
 *  - May throw a DifferentDimensionsException exception if the given 
 *    point and the hyperplane belong in spaces of different dimensions.
 *  - May throw a NegativeApproximationRatioException exception if the 
 *    given approximation ratio/factor/threshold is less than 0.0.
 *  - May throw a NotStrictlyPositivePointException exception if the 
 *    given point is not strictly positive. (i.e. some coordinate is less 
 *    than or equal to 0.0)
 *  
 *  \sa Facet, Point, Point::dominatesMultiplicative() and 
 *      Facet::dominates()
 */
template <class S> 
bool 
Facet<S>::dominatesMultiplicative(const Point & p, double eps) const
{
  if (p.isNull())
    throw exception_classes::NullObjectException();
  if (spaceDimension() != p.dimension())
    throw exception_classes::DifferentDimensionsException();
  if (not p.isStrictlyPositive()) 
    throw exception_classes::NotStrictlyPositivePointException();
  if (eps < 0.0)
    throw exception_classes::NegativeApproximationRatioException();
  // else

  if (ratioDistance(p) <= eps)
    return true;
  else
    return false;
}


//! Check if every element of the facet's normal vector is non-positive.
/*!
 *  \return true if every element of the facet's normal vector 
 *          (Facet<S>::normal_) is non-positive.
 *  
 *  Each element must be non-positive.
 *  
 *  \sa Facet
 */
template <class S> 
bool 
Facet<S>::hasAllNormalVectorElementsNonPositive() const
{
  for (unsigned int i = 0; i != spaceDimension(); ++i)
    if (normal_[i] > 0.0)
      return false;

  return true;
}


//! Check if every element of the facet's normal vector is non-negative.
/*!
 *  \return true if every element of the facet's normal vector 
 *          (Facet<S>::normal_) is non-negative.
 *  
 *  Each element must be non-negative.
 *  
 *  \sa Facet
 */
template <class S> 
bool 
Facet<S>::hasAllNormalVectorElementsNonNegative() const
{
  for (unsigned int i = 0; i != spaceDimension(); ++i)
    if (normal_[i] < 0.0)
      return false;

  return true;
}


//! Normalizes the facet's normal vector.
/*!
 *  Normalizes the facet's normal vector so that its magnitude 
 *  (i.e. length or L2-norm) becomes 1. 
 *  
 *  First computes "l2Norm", which is the current L2-norm of the 
 *  normal vector. Then divides each normal vector element by "l2Norm".
 *  
 *  \sa Facet
 */
template <class S> 
void 
Facet<S>::normalizeNormalVector()
{
  // get the facet's normal vector inside an armadillo vec
  arma::vec normalVec(normal_);

  // compute its L2-norm
  double l2Norm = arma::norm(normalVec, 2);

  // divide each normal vector element by "l2Norm"
  for (unsigned int i = 0; i != spaceDimension(); ++i)
    normal_[i] /= l2Norm;
}


//! Get a copy of the facet's normal vector.
/*!
 *  \return A copy of the facet's normal vector.
 *
 *  \sa Facet
 */
template <class S> 
std::vector<double> 
Facet<S>::getNormalVector() const
{
  return normal_;
}


//! Compute (and set) the facet's normal vector using the facet's vertices.
/*!
 *  \param preferPositiveNormalVector Should we prefer the all-positive 
 *                                    normal vector (if it exists)?
 *  
 *  Calculates the hyperplane passing through the facet's vertices 
 *  and uses its normal vector as the facet's normal vector. 
 *  
 *  For each set of n vertices there are two different n-hyperplanes passing 
 *  through them with opposite normal vectors. This method will prefer the 
 *  all-positive normal vector (if one exists) if preferPositiveNormalVector 
 *  is set to true; otherwise it will choose one depending on the order of 
 *  the facet vertices.
 *  
 *  \sa Facet
 */
template <class S> 
void 
Facet<S>::computeAndSetFacetNormal(bool preferPositiveNormalVector) 
{
  // fill a matrix will each point's coordinates
  arma::mat M;
  ConstVertexIterator vi;
  for (vi = vertices_.begin(); vi != vertices_.end(); ++vi)
    M.insert_rows(M.n_rows, vi->point.toRowVec());
  // add a column of ones at the end (will make the following easier)
  M.insert_cols(M.n_cols, arma::ones<arma::vec>(spaceDimension()));
  
  // fill in the normal vector's elements
  for (unsigned int i = 0; i != spaceDimension(); ++i) {
    M.swap_cols(i, M.n_cols - 1);
    normal_.push_back(arma::det(M.cols(0, M.n_cols - 2)));
    M.swap_cols(i, M.n_cols - 1);
  }

  if (preferPositiveNormalVector && hasAllNormalVectorElementsNonPositive())
    reverseNormalVectorSign();
}


/*! \brief Compute (and set) the facet's isBoundaryFacet_ and 
 *         localApproximationErrorUpperBound_ attributes.
 *  
 *  Computes the facet's local approximation error upper bound (i.e. 
 *  distance from the facet's Lower Distal Point if (a unique) one 
 *  exists) and sets the facet's localApproximationErrorUpperBound_ 
 *  and isBoundaryFacet_ attributes accordingly.
 *  
 *  We have only created this function in order to erase duplicate 
 *  code from inside the constructors.
 *  
 *  \sa Facet
 */
template <class S> 
void 
Facet<S>::computeAndSetLocalApproximationErrorUpperBoundAndIsBoundaryFacet()
{
  // - First find the facet's Lower Distal Point (LDP).
  Point lowerDistalPoint = computeLowerDistalPoint();

  // - If an LDP exists use it to compute the facet's local approximation 
  //   error, else mark the facet as a boundary facet.
  if (not lowerDistalPoint.isNull()) {
    isBoundaryFacet_ = false;
    if (lowerDistalPoint.isStrictlyPositive()) 
      localApproximationErrorUpperBound_ = euclideanDistance(lowerDistalPoint);
    else {
      // The LDP is not strictly positive.
      // - mark the facet as a boundary facet
      isBoundaryFacet_ = true;
      // localApproximationErrorUpperBound_ is not valid now:
      localApproximationErrorUpperBound_ = -1.0;
    }
  }
  else {
    isBoundaryFacet_ = true;
    // localApproximationErrorUpperBound_ is not valid now:
    localApproximationErrorUpperBound_ = -2.0;
  }
}


//! Reverse the sign of all elements of the facet's normal vector.
/*!
 *  Reverse the sign of all the elements of the facet's normal vector 
 *  (Facet<S>::normal_).
 *  
 *  \sa Facet
 */
template <class S> 
void 
Facet<S>::reverseNormalVectorSign()
{
  for (unsigned int i = 0; i != spaceDimension(); ++i)
    normal_[i] = -normal_[i];
}


}  // namespace pareto_approximator


/* @} */
