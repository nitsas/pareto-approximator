/*! \file Facet.cpp
 *  \brief The implementation of the Facet class.
 *  \author Christos Nitsas
 *  \date 2012
 *  
 *  Won't `#include` "Facet.h". In fact "Facet.h" will `#include` 
 *  "Facet.cpp" because it describes a class template (which doesn't allow 
 *  us to split declaration from definition).
 */


#include <assert.h>
#include <iterator>
#include <algorithm>
#include <armadillo>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Facet's default constructor. (empty)
template <class S> 
Facet<S>::Facet() { }


//! A Facet constructor.
/*!
 *  \param firstVertex An iterator to the first of the facet vertices.
 *  \param lastVertex An iterator to the past-the-end element in the 
 *                    container of facet vertices.
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
 *    vector, if one exists.
 *  - Facet<S>::approximationErrorUpperBound_ to the ratio distance 
 *    between the Facet and its Lower Distal Point (LDP). Calculates 
 *    both the LDP and the ratio distance.
 *  
 *  \sa Facet
 */
template <class S> 
Facet<S>::Facet(typename std::vector< PointAndSolution<S> >::const_iterator firstVertex, 
                typename std::vector< PointAndSolution<S> >::const_iterator lastVertex)
{
  // Make sure that all the given vertices are of the correct dimension.
  // (the dimension of the space that the facet lives in)
  unsigned int thisSpaceDimension = std::distance(firstVertex, lastVertex);
  ConstVertexIterator cvi;
  for (cvi = firstVertex; cvi != lastVertex; ++cvi)
    if (cvi->point.dimension() != thisSpaceDimension)
      throw DifferentDimensionsException();

  // First fill-in Facet<S>::vertices_.
  vertices_.assign(firstVertex, lastVertex);

  // Then compute the facet's normal vector (Facet<S>::normal_).
  computeAndSetFacetNormal();

  // Now compute Facet<S>::localApproximationErrorUpperBound.
  // - First find the facet's Lower Distal Point (LDP).
  Point* lowerDistalPoint = computeLowerDistalPoint();

  // - If an LDP exists use it to compute the facet's local approximation 
  //   error, else mark the facet as a boundary facet.
  if (lowerDistalPoint != NULL) {
    isBoundaryFacet_ = false;
    localApproximationErrorUpperBound_ = ratioDistance(*lowerDistalPoint);
  }
  else {
    isBoundaryFacet_ = true;
    // localApproximationErrorUpperBound_ is not valid now:
    localApproximationErrorUpperBound_ = -1.0;
  }
}


//! Facet constructor.
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
 *  Initializes:
 *  - Facet<S>::vertices_ to the sequence of vertices pointed to by 
 *    firstVertex and lastVertex. The range used is 
 *    [firstVertex, lastVertex).
 *  - Facet<S>::normal_ to the sequence of elements pointed to by 
 *    firstElemOfFacetNormal and lastElemOfFacetNormal. The range 
 *    used is [firstElemOfFacetNormal, lastElemOfFacetNormal).
 *  - Facet<S>::approximationErrorUpperBound_ to the ratio distance between 
 *    the Facet and its Lower Distal Point (LDP). Calculates both the 
 *    LDP and the ratio distance.
 *  
 *  \sa Facet
 */
template <class S> 
Facet<S>::Facet(typename std::vector< PointAndSolution<S> >::const_iterator firstVertex, 
                typename std::vector< PointAndSolution<S> >::const_iterator lastVertex, 
                std::vector<double>::const_iterator firstElemOfFacetNormal, 
                std::vector<double>::const_iterator lastElemOfFacetNormal)
{
  assert(std::distance(firstVertex, lastVertex) == 
         std::distance(firstElemOfFacetNormal, lastElemOfFacetNormal));

  // Make sure that all the given vertices are of the correct dimension.
  // (the dimension of the space that the facet lives in)
  unsigned int thisSpaceDimension = std::distance(firstVertex, lastVertex);
  ConstVertexIterator cvi;
  for (cvi = firstVertex; cvi != lastVertex; ++cvi)
    if (cvi->point.dimension() != thisSpaceDimension)
      throw DifferentDimensionsException();

  // First fill-in Facet<S>::vertices_ and Facet<S>::normal_.
  vertices_.assign(firstVertex, lastVertex);
  normal_.assign(firstElemOfFacetNormal, lastElemOfFacetNormal);

  // Now compute Facet<S>::localApproximationErrorUpperBound.
  // - First find the facet's Lower Distal Point (LDP).
  Point* lowerDistalPoint = computeLowerDistalPoint();

  // - If an LDP exists use it to compute the facet's local approximation 
  //   error, else mark the facet as a boundary facet.
  if (lowerDistalPoint != NULL) {
    isBoundaryFacet_ = false;
    localApproximationErrorUpperBound_ = ratioDistance(*lowerDistalPoint);
  }
  else {
    isBoundaryFacet_ = true;
    // localApproximationErrorUpperBound_ is not valid now:
    localApproximationErrorUpperBound_ = -1.0;
  }
}


//! A simple (and empty) destructor.
template <class S> 
Facet<S>::~Facet() { }


//! Return the dimension of the space the facet lives in.
/*!
 *  \sa Facet
 */
template <class S> 
unsigned int 
Facet<S>::spaceDimension() const
{
  // The space's dimension will be the same as the facet-normal-vector's
  // dimension. 
  return normal_.size();
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
 *  We will use the ratio distance from the facet to it's Lower Distal 
 *  Point as an upper bound to the local approximation error.
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
    throw BoundaryFacetException();
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


//! Compute the mean of all the weight vectors of the facet's vertices.
/*!
 *  \return A weight vector W of size this->spaceDimension(). Each 
 *          element W_{j} is the mean of all w_{ij}'s, where w_{i} is 
 *          the weight vector inside the i'th facet vertex's 
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
 *  \return A pointer to the facet's Lower Distal Point (Point instance) 
 *          if one exists, NULL otherwise.
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
 *  If a unique solution does not exist we return NULL.
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
 *  this method returns NULL and the current facet is treated as a 
 *  boundary facet.
 *  
 *  \sa Facet
 */
template <class S> 
Point* 
Facet<S>::computeLowerDistalPoint() const
{
  // make a stream to redirect errors to
  std::ofstream f("my_log.txt");
  // redirect armadillo error messages to that stream
  // (e.g. when arma::solve() finds no solutions)
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
  Point* result = NULL;
  bool hasSolution = arma::solve(x, W, b);
  if (hasSolution) 
    // unique solution
    result = new Point(x.begin(), x.end());
  else 
    // either no solution or an infinite number of solutions
    result = NULL;
  
  return result;
}


//! Compute a point's ratio distance from the hyperplane the facet lies on.
/*! 
 *  \param p A Point instance. (with non-negative coordinates)
 *  \return The point's ratio distance from the hyperplane on which the 
 *          facet lies.
 *  
 *  The ratio distance from a point p to a hyperplane H is defined as:
 *  \f$ RD(p, H) = \min_{q \in H} RD(p, q) \f$, where q is a point on H.
 *  The ratio distance from a point p to a point q is defined as:
 *  \f$ RD(p, q) = \max\{ \max_{i}\{(q_{i} - p_{i}) / p_{i}\}, 0.0 \} \f$.
 *  
 *  Intuitively it is the minimum value of \f$ \epsilon \ge 0 \f$ such 
 *  that some point on H \f$ \epsilon -covers p \f$.
 *  
 *  In order for the ratio distance to make sense point p must have 
 *  non-negative coordinates.
 *  
 *  Possible exceptions:
 *  - May throw a DifferentDimensionsException exception if the given point 
 *    and the hyperplane belong in spaces of different dimensions.
 *  
 *  \sa Point, Facet and Facet<S>::computeFacetsLowerDistalPoint()
 */
template <class S> 
double 
Facet<S>::ratioDistance(const Point& p) const
{
  assert(spaceDimension() > 0);

  if (spaceDimension() != p.dimension())
    throw DifferentDimensionsException();
  // else
  Point aPointOnTheFacet = vertices_[0].point;
  double dotProduct  = 0.0;
  double facetOffset = 0.0;      // the facet offset from the origin
  for (unsigned int i=0; i!=spaceDimension(); ++i) {
    dotProduct  += normal_[i] * p[i];
    facetOffset += normal_[i] * aPointOnTheFacet[i];
  }

  return std::max( (facetOffset - dotProduct) / dotProduct, 0.0 );
}


//! Check if every element of the facet's normal vector is non-positive.
/*!
 *  \return true if every element of the facet's normal vector 
 *               (Facet<S>::normal_) is non-positive.
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


//! Compute the facet's normal vector.
/*!
 *  Calculates the hyperplane passing through the facet's vertices 
 *  and uses its normal vector as the facet's normal vector. For each 
 *  set of n vertices there are two different n-hyperplanes passing through 
 *  them, with opposite normal vectors. This method will prefer the 
 *  all-positive normal vector, if one exists.
 *  
 *  \sa Facet
 */
template <class S> 
void 
Facet<S>::computeAndSetFacetNormal() 
{
  unsigned int thisSpaceDimension = vertices_.size();

  // fill a matrix will each point's coordinates
  arma::mat M;
  ConstVertexIterator vi;
  for (vi = vertices_.begin(); vi != vertices_.end(); ++vi)
    M.insert_rows(M.n_rows, vi->point.toRowVec());
  // add a column of ones at the end (will make the following easier)
  M.insert_cols(M.n_cols, arma::ones<arma::vec>(thisSpaceDimension));
  
  // fill in the normal vector's elements
  for (unsigned int i = 0; i != thisSpaceDimension; ++i) {
    M.swap_cols(i, M.n_cols - 1);
    normal_.push_back(arma::det(M.cols(0, M.n_cols - 2)));
    M.swap_cols(i, M.n_cols - 1);
  }

  if (hasAllNormalVectorElementsNonPositive())
    reverseNormalVectorSign();
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
