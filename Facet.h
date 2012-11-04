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
#include "BoundaryFacetException.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


/*!
 *  \brief A template for a class representing a facet of the current upper 
 *         approximation.
 *  
 *  The template parameter is carried on from the BaseProblem class template.
 *  It is the type of the problem solution.
 *  
 *  In order to represent a facet we will use the vertices of the facet
 *  (Facet::vertices_) plus the facet normal (Facet::normal_). (The Facet 
 *  class contains two more private variables, Facet::isBoundaryFacet_ and 
 *  Facet::localApproximationErrorUpperBound_ but they are only used to 
 *  store secondary info about the facet.)
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
 *  Each Facet instance will also include the ratio distance from the 
 *  facet to it's Lower Distal Point (Facet::localApproximationErrorUpperBound).
 *  The ratio distance from the facet to its LDP is an upper bound to the 
 *  approximation error for the facet.
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
 *  \sa BaseProblem, PointAndSolution, Point and Hyperplane
 */
template <class S>
class Facet
{
  public:
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

    //! Facet's default constructor. (empty)
    Facet();

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
    Facet(typename std::vector< PointAndSolution<S> >::const_iterator firstVertex, 
          typename std::vector< PointAndSolution<S> >::const_iterator lastVertex);

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
     *  Initializes:
     *  - Facet<S>::vertices_ to the sequence of vertices pointed to by 
     *    firstVertex and lastVertex. The range used is 
     *    [firstVertex, lastVertex).
     *  - Facet<S>::normal_ to the sequence of elements pointed to by 
     *    firstElemOfFacetNormal and lastElemOfFacetNormal. The range 
     *    used is [firstElemOfFacetNormal, lastElemOfFacetNormal).
     *  - Facet<S>::approximationErrorUpperBound_ to the ratio distance 
     *    between the Facet and its Lower Distal Point (LDP). Calculates 
     *    both the LDP and the ratio distance.
     *  
     *  \sa Facet
     */
    Facet(typename std::vector< PointAndSolution<S> >::const_iterator firstVertex, 
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

    //! Compute the mean of all the weight vectors of the facet's vertices.
    /*!
     *  \return A weight vector W of size this->spaceDimension(). Each 
     *          element W_{j} is the mean of all w_{ij}'s, where w_{i} is 
     *          the weight vector inside the i'th facet vertex's 
     *          PointAndSolution object (PointAndSolution::weightsUsed).
     *  
     *  \sa Facet
     */
    std::vector<double> computeMeanVertexWeights() const;

    //! Compute the facet's Lower Distal Point (LDP).
    /*! 
     *  \return A pointer to the facet's Lower Distal Point (Point instance) 
     *          if one exists, NULL otherwise.
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
    Point* computeLowerDistalPoint() const;
    
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
    double ratioDistance(const Point& p) const;

    //! Check if every element of the facet's normal vector is non-positive.
    /*!
     *  \return true if every element of the facet's normal vector 
     *               (Facet<S>::normal_) is non-positive.
     *  
     *  Each element must be non-positive.
     *  
     *  \sa Facet
     */
    bool hasAllNormalVectorElementsNonPositive() const;

  private: 
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
    void computeAndSetFacetNormal();

    //! Reverse the sign of all elements of the facet's normal vector.
    /*!
     *  Reverse the sign of all the elements of the facet's normal vector 
     *  (Facet<S>::normal_).
     *  
     *  \sa Facet
     */
    void reverseNormalVectorSign();

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

    //! An upper bound to the current facet's approximation error.
    /*! 
     *  The ratio distance from the facet to it's Lower Distal Point.
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
