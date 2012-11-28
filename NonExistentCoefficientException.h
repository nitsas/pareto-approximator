/*! \file NonExistentCoefficientException.h
 *  \brief The declaration and definition of the NonExistentCoefficient 
 *         exception class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef NON_EXISTENT_COEFFICIENT_EXCEPTION_H
#define NON_EXISTENT_COEFFICIENT_EXCEPTION_H

#include <exception>


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Exception thrown by Hyperplane::operator[]().
/*! 
 *  An exception thrown when the requested Hyperplane coefficient does not 
 *  exist. (out of bounds)
 */
class NonExistentCoefficientException : public std::exception
{
  public:
    //! Return a simple char* message.
    const char* what() const throw()
    {
      return "The requested coefficient does not exist. (out of bounds)";
    }
};


}  // namespace pareto_approximator


/* @} */


#endif  // NON_EXISTENT_COEFFICIENT_EXCEPTION_H
