/*! \file Facet.h
 *  \brief The definition of the Facet class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef FACET_H
#define FACET_H

#include <vector>

#include "Point.h"
#include "PointAndSolution.h"
#include "DifferentDimensionsException.h"
#include "NullObjectException.h"
#include "NotStrictlyPositivePointException.h"
#include "BoundaryFacetException.h"
#include "InfiniteRatioDistanceException.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
namespace pareto_approximator {


/*!
 *  \brief A template for a class representing a facet of the current upper 
 *         approximation.
 *  
 *  The template parameter is carried on from the BaseProblem class template.
 *  It is the type of the problem solution.
 *  
 *  In order to represent a facet we will use the vertices of the facet
 *  (Facet::vertices_), the facet normal (Facet::normal_) and the dimension 
 *  of the space that the facet lives in (Facet::spaceDimension_). (The Facet 
 *  class contains two more private variables, Facet::isBoundaryFacet_ and 
 *  Facet::localApproximationErrorUpperBound_ but they are only used to 
 *  store secondary/derived info about the facet.)
 *  
 *  The facet normal is a vector perpendicular to the facet's surface. 
 *  The normal is simply the direction that the facet is facing.
 *  
 *  Each vertex will be a PointAndSolution<S> instance containing:
 *  - the point in objective space (a Point instance)
 *  - the corresponding solution (an S instance)
 *  - the weights we used to find the aforementioned point 
 *    (a std::vector<double>)
 *  
 *  Each Facet instance will also include the distance from the facet to 
 *  its Lower Distal Point (Facet::localApproximationErrorUpperBound). The 
 *  distance from the facet to its LDP is an upper bound to the approximation 
 *  error for the facet.
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
 *  \sa BaseProblem, PointAndSolution and Point
 */
template <class S>
class Facet
{
  public:
    //! The type of a facet vertex.
    typedef PointAndSolution<S> Vertex;

    //! The type of a set (std::vector) of vertices.
    typedef std::vector< PointAndSolution<S> > VerticesVector;

    //! The type of an element of the facet's normal vector.
    typedef double NormalVectorElement;

    //! The type of the facet's normal vector.
    typedef std::vector<double> FacetNormalVector;

    //! Random access iterator to the facet's vertices.
    typedef typename std::vector< PointAndSolution<S> >::iterator 
                     VertexIterator;
    
    //! Constant random access iterator to the facet's vertices.
    typedef typename std::vector< PointAndSolution<S> >::const_iterator 
                     ConstVertexIterator;

    //! Random access iterator to the elements of the facet's normal.
    typedef std::vector<double>::iterator FacetNormalIterator;

    //! Constant random access iterator to the elements of the facet's normal.
    typedef std::vector<double>::const_iterator ConstFacetNormalIterator;

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
     *  - Facet<S>::localApproximationErrorUpperBound_ to the distance 
     *    between the Facet and its Lower Distal Point (LDP). Calculates 
     *    both the LDP and the distance.
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
    Facet(
      typename std::vector< PointAndSolution<S> >::const_iterator firstVertex, 
      typename std::vector< PointAndSolution<S> >::const_iterator lastVertex, 
      bool preferPositiveNormalVector=true);

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
     *  - Facet<S>::localApproximationErrorUpperBound_ to the distance 
     *    between the Facet and its Lower Distal Point (LDP). Calculates 
     *    both the LDP and the distance.
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
    Facet(
      typename std::vector< PointAndSolution<S> >::const_iterator firstVertex, 
      typename std::vector< PointAndSolution<S> >::const_iterator lastVertex, 
      std::vector<double>::const_iterator firstElemOfFacetNormal, 
      std::vector<double>::const_iterator lastElemOfFacetNormal);

    //! A simple (and empty) destructor.
    ~Facet();

    //! Return the dimension of the space the facet lives in.
    /*!
     *  \sa Facet
     */
    unsigned int spaceDimension() const;

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
    bool isBoundaryFacet() const;

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
    double getLocalApproximationErrorUpperBound() const;

    //! Return iterator to the beginning of the vector of facet vertices.
    /*! 
     *  \return An iterator pointing to the first vertex in the vector 
     *          of vertices.
     *  
     *  \sa Facet
     */
    ConstVertexIterator beginVertex() const;

    //! Return iterator to the end of the vector of facet vertices.
    /*! 
     *  \return An iterator pointing just after the last vertex in the 
     *          vector of vertices.
     *  
     *  \sa Facet
     */
    ConstVertexIterator endVertex() const;

    //! Return iterator to the beginning of the facet's normal vector.
    /*! 
     *  \return An iterator pointing to the first element in the facet's 
     *          normal vector.
     *  
     *  \sa Facet
     */
    ConstFacetNormalIterator beginFacetNormal() const;

    //! Return iterator to the end of the facet's normal vector.
    /*! 
     *  \return An iterator pointing just after the last element in the 
     *          facet's normal vector.
     *  
     *  \sa Facet
     */
    ConstFacetNormalIterator endFacetNormal() const;

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
    double b() const;

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
    std::vector<double> computeMeanVertexWeights() const;

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
    Point computeLowerDistalPoint() const;
    
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
    double distance(const Point & p) const;

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
     *  Possible exceptions:
     *  - May throw a DifferentDimensionsException exception if the given point 
     *    and the hyperplane belong in spaces of different dimensions.
     *  - May throw a NullObjectException exception if the given Point 
     *    instance is a null Point instance.
     *  
     *  \sa Facet and Point
     */
    double euclideanDistance(const Point & p) const;

    /*!
     *  \brief Compute a point's ratio distance from the hyperplane that 
     *         the facet lies on.
     *  
     *  \param p A Point instance. (should be strictly positive)
     *  \return The point's ratio distance from the hyperplane on which the 
     *          facet lies.
     *  
     *  The ratio distance from a point p to a hyperplane H is defined as:
     *  \f$ RD(p, H) = \min_{q \in H} RD(p, q) \f$, where q is a point on H.
     *  The ratio distance from a point p to a point q is defined as:
     *  \f$ RD(p, q) = \max\{ \max_{i}\{(q_{i} - p_{i}) / p_{i}\}, 0.0 \} \f$.
     *  
     *  Intuitively it is the minimum value of \f$ \epsilon \ge 0 \f$ such 
     *  that some point on H \f$\epsilon\f$ -dominates (\f$\epsilon\f$ 
     *  -covers) p in the multiplicative sense.
     *  
     *  In order for the ratio distance to make sense point p must be 
     *  strictly positive, i.e. \f$ p_{i} > 0.0 \f$ must hold for all i.
     *  
     *  Possible exceptions:
     *  - May throw a DifferentDimensionsException exception if the given 
     *    point and the hyperplane belong in spaces of different dimensions.
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
    double ratioDistance(const Point & p) const;

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
     *  that some point on H \f$\epsilon\f$ -dominates (\f$\epsilon\f$ 
     *  -covers) p in the additive sense.
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
    double additiveDistance(const Point & p) const;

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
    bool dominates(const Point & p, double eps=0.0) const;

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
    bool dominatesAdditive(const Point & p, double eps=0.0) const;

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
    bool dominatesMultiplicative(const Point & p, double eps=0.0) const;

    //! Check if every element of the facet's normal vector is non-positive.
    /*!
     *  \return true if every element of the facet's normal vector 
     *          (Facet<S>::normal_) is non-positive.
     *  
     *  Each element must be non-positive.
     *  
     *  \sa Facet
     */
    bool hasAllNormalVectorElementsNonPositive() const;

    //! Check if every element of the facet's normal vector is non-negative.
    /*!
     *  \return true if every element of the facet's normal vector 
     *          (Facet<S>::normal_) is non-negative.
     *  
     *  Each element must be non-negative.
     *  
     *  \sa Facet
     */
    bool hasAllNormalVectorElementsNonNegative() const;

    //! Normalizes the facet's normal vector.
    /*!
     *  Normalizes the facet's normal vector so that its magnitude 
     *  (i.e. length or L2-norm) becomes 1. 
     *  
     *  First computes "l2Norm", which is the current L2-norm of the 
     *  normal vector. Then divides each normal vector element with "l2Norm".
     *  
     *  \sa Facet
     */
    void normalizeNormalVector();

    //! Get a copy of the facet's normal vector.
    /*!
     *  \return A copy of the facet's normal vector.
     *
     *  \sa Facet
     */
    std::vector<double> getNormalVector() const;

    //! \brief Check if the facet is coplanar with the given point.
    //!
    //! \param p A point with the same dimensions as the facet.
    //! \return true if the facet and the point are coplanar; false otherwise.
    //!
    //! A point and a facet are coplanar if the point is on the facet's 
    //! supporting hyperplane.
    //!
    //! Possible exceptions:
    //! - May throw a NullObjectException exception if the given Point 
    //!   instance is a null Point instance.
    //! - May throw a DifferentDimensionsException exception if the given 
    //!   point and the hyperplane belong in spaces of different dimensions.
    //! 
    //! \sa Facet and Point
    //!
    bool isCoplanarWith(const Point & p) const;

  private: 

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
    void computeAndSetFacetNormal(bool preferPositiveNormalVector);

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
    void computeAndSetLocalApproximationErrorUpperBoundAndIsBoundaryFacet();

    //! Reverse the sign of all elements of the facet's normal vector.
    /*!
     *  Reverse the sign of all the elements of the facet's normal vector 
     *  (Facet<S>::normal_).
     *  
     *  \sa Facet
     */
    void reverseNormalVectorSign();

    //! The dimension of the space that the facet lives in.
    /*!
     *  \sa Facet
     */
    unsigned int spaceDimension_;

    //! The vertices of the facet.
    /*! 
     *  Each vertex is stored as a PointAndSolution<S> instance containing:
     *  - the point in objective space (a Point instance)
     *  - the corresponding solution (an S instance)
     *  - the weights we used to find the aforementioned point 
     *    (a std::vector<double>)
     *  
     *  \sa Facet
     */
    std::vector< PointAndSolution<S> > vertices_;

    //! The facet normal.
    /*! 
     *  The facet normal is a vector perpendicular to the facet's surface. 
     *  The normal is simply the direction that the facet is facing.
     *  
     *  \sa Facet
     */
    std::vector<double> normal_;

    //! The offset of the hyperplane on which the facet lies..
    /*!
     *  For point on the underlying hyperplane, the dot product of its 
     *  coordinates and the hyperplane/facet normal vector will be equal 
     *  to b_.
     */
    double b_;

    //! An upper bound to the current facet's approximation error.
    /*! 
     *  The distance from the facet to it's Lower Distal Point.
     *  It is an upper bound to the local approximation error.
     *  
     *  Check the documentation for Facet for a description of what a Lower 
     *  Distal Point is.
     *  
     *  \sa Facet and Facet<S>::computeLowerDistalPoint()
     */
    double localApproximationErrorUpperBound_;

    //! Is the facet a boundary facet?
    /*!
     *  We call a facet a boundary facet if it does not have a Lower Distal 
     *  Point. (The hyperplanes h_{i} associated with its vertices v_{i} 
     *  do not intersect in a unique point. Check the documentation of Facet 
     *  for more info.)
     *  
     *  \sa Facet and Facet<S>::computeLowerDistalPoint()
     */
    bool isBoundaryFacet_;
};


}  // namespace pareto_approximator


/* @} */


// We have got to #include the implementation here because we are 
// describing a class template, not a simple class.
#include "Facet.cpp"


#endif  // FACET_H
