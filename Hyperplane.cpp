/*! \file Hyperplane.cpp
 *  \brief The declaration of the Hyperplane class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <assert.h>
#include <string>
#include <sstream>
#include <vector>

#include "Hyperplane.h"


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
Hyperplane::Hyperplane(std::vector<int>::iterator first, 
                       std::vector<int>::iterator last, int b) : b_(b)
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
Hyperplane::Hyperplane(std::vector<double>::iterator first, 
                       std::vector<double>::iterator last, double b) : b_(b)
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
Hyperplane::Hyperplane(int* first, int* last, int b) : b_(b)
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
Hyperplane::Hyperplane(double* first, double* last, double b) : b_(b)
{
  assert(coefficients_.size() == 0);
  coefficients_.assign(first, last);
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


}  // namespace pareto_approximator
