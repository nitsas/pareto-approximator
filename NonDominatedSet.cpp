/*! \file NonDominatedSet.cpp
 *  \brief The definition of the NonDominatedSet<T> class template.
 *  \author Christos Nitsas
 *  \date 2012
 *  
 *  Won't `#include` "NonDominatedSet.h". In fact "NonDominatedSet.h" will 
 *  `#include` "NonDominatedSet.cpp" because it describes a class template 
 *  (which doesn't allow us to split declaration from definition).
 */


#include <assert.h>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Default constructor. Makes an empty set.
template <class T> 
NonDominatedSet<T>::NonDominatedSet() : contents_() { }


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
template <class T> 
template <class InputIterator> 
NonDominatedSet<T>::NonDominatedSet(InputIterator first, InputIterator last)
{
  contents_.insert(first, last);
}


//! Destructor. (all the contained elements' destructors will be called) 
template <class T> 
NonDominatedSet<T>::~NonDominatedSet() { }


//! Return iterator to beginning.
/*!
 *  \return An iterator referring to the first T element in the 
 *          container.
 *  
 *  \sa NonDominatedSet
 */
template <class T> 
typename NonDominatedSet<T>::iterator 
NonDominatedSet<T>::begin()
{
  return contents_.begin();
}


//! Return const_iterator to beginning.
/*!
 *  \return A const_iterator referring to the first T element in the 
 *          container.
 *  
 *  \sa NonDominatedSet
 */
template <class T> 
typename NonDominatedSet<T>::const_iterator 
NonDominatedSet<T>::begin() const
{
  return contents_.begin();
}


//! Return iterator to end.
/*!
 *  \return An iterator pointing to the past-the-end T element in the 
 *          container.
 *  
 *  \sa NonDominatedSet
 */
template <class T> 
typename NonDominatedSet<T>::iterator 
NonDominatedSet<T>::end()
{
  return contents_.end();
}


//! Return const_iterator to end.
/*!
 *  \return A const_iterator pointing to the past-the-end T element 
 *          in the container.
 *  
 *  \sa NonDominatedSet
 */
template <class T> 
typename NonDominatedSet<T>::const_iterator 
NonDominatedSet<T>::end() const
{
  return contents_.end();
}


//! Test whether container is empty.
/*!
 *  \return true if the container is empty; false otherwise.
 *  
 *  \sa NonDominatedSet
 */
template <class T> 
bool 
NonDominatedSet<T>::empty() const
{
  return contents_.empty();
}


//! Returns the number of elements in the container.
/*!
 *  \return The number of elements that form the set's content.
 *  
 *  Member type size_type is (usually) an unsigned integral type.
 *
 *  \sa NonDominatedSet
 */
template <class T> 
typename NonDominatedSet<T>::size_type 
NonDominatedSet<T>::size() const
{
  return contents_.size();
}


//! Insert element.
/*!
 *  \param t The T instance to insert.
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
template <class T> 
void 
NonDominatedSet<T>::insert(const T & t)
{
  bool tIsDominated = this->dominates(t);
  if (!tIsDominated) {
    // first remove any elements dominated by t
    for (iterator it = contents_.begin(); it != contents_.end(); )
      if (t.dominates(*it)) 
        // increment "it" before removing its previous value from the set
        contents_.erase(it++);
      else
        // just increment "it"
        it++;
    // now insert the new element
    contents_.insert(t);
  }
}


//! Insert elements.
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
template <class T> 
template <class InputIterator> 
void 
NonDominatedSet<T>::insert(InputIterator first, InputIterator last)
{
  for ( ; first != last; ++first)
    this->insert(*first);
}


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
template <class T> 
bool 
NonDominatedSet<T>::dominates(const T & t) const
{
  iterator it;
  for (it = contents_.begin(); it != contents_.end(); ++it) 
    if (it->dominates(t))
      return true;

  return false;
}


//! Clear content. 
/*!
 *  All the elements' destructors are called and then they are removed 
 *  from the container, leaving it with a size of 0.
 *  
 *  \sa NonDominatedSet
 */
template <class T> 
void 
NonDominatedSet<T>::clear()
{
  contents_.clear();
}


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
template <class T> 
typename NonDominatedSet<T>::iterator 
NonDominatedSet<T>::find(const T & t) const
{
  iterator it;
  for (it = contents_.begin(); it != contents_.end(); ++it)
    if (*it == t)
      break;

  return it;
}


}  // namespace pareto_approximator


/* @} */
