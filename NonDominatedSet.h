/*! \file NonDominatedSet.h
 *  \brief The definition of the NonDominatedSet<T> class template.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef NON_DOMINATED_SET_H
#define NON_DOMINATED_SET_H

#include <set>

#include "Point.h"
#include "PointAndSolution.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
namespace pareto_approximator {


//! A container that only keeps non-dominated points. (plus solutions)
/*!
 *  The NonDominatedSet<T> is a type of container/filter. It stores unique 
 *  T instances. T must represent some kind of a point. It must have at 
 *  least the following methods:
 *  - bool dominates(const T & s, double eps): 
 *    Returns true if the current instance eps-covers s; false otherwise.
 *    We say that a point p \f$ \epsilon \f$-covers another point q 
 *    (\f$\epsilon \ge 0 \f$) iff \f$ p_{i} \le (1 + \epsilon) q_{i} \f$, 
 *    for all i. Both p and q must be of the same dimension.
 *    (see Point or PointAndSolution for more info and examples)
 *  - bool operator==(const T & s):
 *    Returns true if the current instance is equal to s; false otherwise.
 *  - bool operator<(const T & s):
 *    Returns true if the current instance is less-than s; false otherwise.
 *    (only needed for the sort method)
 *  
 *  The key difference from a std::set<T> container is the fact that 
 *  NonDominatedSet<T> acts as a filter that only keeps those T instances 
 *  that form a Pareto frontier - that is, a frontier of Pareto-efficient
 *  points (points not dominated by each other or by any point we have 
 *  inserted so far).
 *
 *  Another difference is that we have not exposed all methods of the 
 *  std::set that lies underneath. The NonDominatedSet class is meant for 
 *  a very specific purpose and not as a general container. However, if 
 *  some new functionality is required it shouldn't be hard to add new 
 *  methods that delegate to the underlying set's methods.
 *  
 *  \sa PointAndSolution
 */
template <class T> 
class NonDominatedSet
{
  public:
    //! Bidirectional iterator to the set's contents.
    typedef typename std::set<T>::iterator iterator;

    //! Constant bidirectional iterator to the set's contents.
    typedef typename std::set<T>::const_iterator const_iterator;

    //! Unsigned integral type (usually same as size_t).
    typedef typename std::set<T>::size_type size_type;

    //! Default constructor. Makes an empty set.
    NonDominatedSet();

    //! Iteration constructor. 
    /*!
     *  \param first Input iterator to the initial position in a sequence 
     *               of T instances. 
     *  \param last Input iterator to the past-the-end position in a 
     *              sequence of T instances.
     *  
     *  Reminder: T instances represent some kind of point. (e.g. the Point 
     *            or PointAndSolution<T> classes)
     *  
     *  Iterates between first and last, filtering and inserting the 
     *  elements of the sequence into the container object. The range used 
     *  is [first, last), which includes all the elements between first and 
     *  last, including the element pointed by first but not the element 
     *  pointed by last. 
     *  
     *  Elements will be filtered during the insertion process so that 
     *  only non-dominated T instances are kept.
     *
     *  The function template type can be any type of input iterator 
     *  that points to T instances.
     *
     *  \sa NonDominatedSet
     */
    template <class InputIterator>
    NonDominatedSet(InputIterator first, InputIterator last);

    //! Destructor. (all the contained elements' destructors will be called) 
    ~NonDominatedSet();

    //! Return iterator to beginning.
    /*!
     *  \return An iterator referring to the first T element in the 
     *          container.
     *  
     *  \sa NonDominatedSet
     */
    iterator begin();

    //! Return const_iterator to beginning.
    /*!
     *  \return A const_iterator referring to the first T element in the 
     *          container.
     *  
     *  \sa NonDominatedSet
     */
    const_iterator begin() const;

    //! Return iterator to end.
    /*!
     *  \return An iterator pointing to the past-the-end T element in the 
     *          container.
     *  
     *  \sa NonDominatedSet
     */
    iterator end();

    //! Return const_iterator to end.
    /*!
     *  \return A const_iterator pointing to the past-the-end T element 
     *          in the container.
     *  
     *  \sa NonDominatedSet
     */
    const_iterator end() const;

    //! Test whether container is empty.
    /*!
     *  \return true if the container is empty; false otherwise.
     *  
     *  \sa NonDominatedSet
     */
    bool empty() const;

    //! Returns the number of elements in the container.
    /*!
     *  \return The number of elements that form the set's content.
     *  
     *  Member type size_type is (usually) an unsigned integral type.
     *
     *  \sa NonDominatedSet
     */
    size_type size() const;

    //! Insert element. 
    /*!
     *  \param t The T instance to insert.
     *  \return true if the element was actually inserted (was not 
     *          dominated); false otherwise.
     *
     *  Reminder: T instances represent some kind of point. (e.g. the Point 
     *            or PointAndSolution<T> classes)
     *  
     *  The container is extended by inserting a single new element iff the 
     *  new element is not already dominated.
     *  
     *  The element will be inserted if and only if it is not dominated by 
     *  any other element in the set. Any elements dominated by the newly 
     *  inserted element will be deleted.
     *  
     *  \sa NonDominatedSet
     */
    bool insert(const T & t);

    //! Insert elements.
    /*!
     *  \param first Input iterator to the initial position in a sequence 
     *               of T instances.
     *  \param last Input iterator to the past-the-end position in a 
     *              sequence of T instances.
     *  \return true if at least one element was actually inserted (not 
     *          dominated); false otherwise.
     *  
     *  Reminder: T instances represent some kind of point. (e.g. the Point 
     *            or PointAndSolution<T> classes)
     *  
     *  Iterates between first and last, filtering and inserting the 
     *  elements of the sequence into the container object. The range used 
     *  is [first, last), which includes all the elements between first and 
     *  last, including the element pointed by first but not the element 
     *  pointed by last. 
     *  
     *  Elements will be filtered during the insertion process so that 
     *  only non-dominated T instances are kept.
     *
     *  The function template type can be any type of input iterator 
     *  that points to T instances.
     *
     *  \sa NonDominatedSet
     */
    template <class InputIterator> 
    bool insert(InputIterator first, InputIterator last);

    //! Check if some element in the set dominates the given instance.
    /*!
     *  \param t A T instance.
     *  \return true if some element in the set dominates t; false otherwise.
     *  
     *  For every element (T instance) in the set, check if it dominates t.
     *  Return true if some does; false otherwise.
     *  
     *  \sa NonDominatedSet
     */
    bool dominates(const T & t) const;

    //! Clear content. 
    /*!
     *  All the elements' destructors are called and then they are removed 
     *  from the container, leaving it with a size of 0.
     *  
     *  \sa NonDominatedSet
     */
    void clear();

    //! Get iterator to element.
    /*!
     *  \param t A T instance.
     *  \return An iterator to the T element which is equal to t if one 
     *          exists; an iterator to the past-the-end element otherwise.
     *  
     *  Searches the container for a T instance that is equal to the given 
     *  instance (t) and returns an iterator to it if it exists, otherwise 
     *  it returns an iterator to NonDominatedSet::end().
     *  
     *  \sa NonDominatedSet and NonDominatedSet::end()
     */
    iterator find(const T & t) const;

  private:
    std::set<T> contents_;
};


}  // namespace pareto_approximator


/* @} */


// We've got to #include the implementation here because we are describing 
// a class template, not a simple class.
#include "NonDominatedSet.cpp"


#endif  // NON_DOMINATED_SET_H
