/*! \file BoundaryFacetException.h
 *  \brief The declaration and definition of the BoundaryFacetException 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef BOUNDARY_FACET_EXCEPTION_H
#define BOUNDARY_FACET_EXCEPTION_H

#include <exception>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the Pareto set approximation algorithms.
 *  @{
 */


//! The namespace containing everything needed for the Pareto set approximation algorithms.
namespace pareto_approximator {


//! The namespace containing all the exception classes.
namespace exception_classes {


/*! 
 *  \brief Exception thrown by Facet::getLocalApproximationErrorUpperBound().
 *
 *  An exception thrown when:
 *  - The localApproximationErrorUpperBound_ attribute of a Facet instance 
 *    was requested (through the getter) but the facet is a boundary facet 
 *    (localApproximationErrorUpperBound_ not defined).
 */
class BoundaryFacetException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw()
    {
      return "The facet is a boundary facet.";
    }
};


}  // namespace exception_classes


}  // namespace pareto_approximator


/*! @} */


#endif  // BOUNDARY_FACET_EXCEPTION_H
