/*! \file Hyperplane.h
 *  \brief The declaration of the Hyperplane class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef HYPERPLANE_H
#define HYPERPLANE_H

#include <iostream>
#include <string>
#include <vector>

#include "Point.h"
#include "NonExistentCoefficientException.h"


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! A class representing a hyperplane on an n-dimensional space.
/*!
 *  A hyperplane on an n-dimensional space can be described by an equation 
 *  of the form:
 *  \f$ a_{1} x_{1} + a_{2} x_{2} + ... + a_{n} x_{n} = b \f$
 *  
 *  We use the a_{i} coefficients and b to represent the hyperplane.
 *  
 *  \sa Hyperplane() and ~Hyperplane()
 */
class Hyperplane
{
  public:
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
    Hyperplane(double a1, double a2, double b);

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
    Hyperplane(double a1, double a2, double a3, double b);

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
    Hyperplane(std::vector<int>::iterator first, 
               std::vector<int>::iterator last, int b);

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
    Hyperplane(std::vector<double>::iterator first, 
               std::vector<double>::iterator last, double b);

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
    Hyperplane(int* first, int* last, int b);

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
    Hyperplane(double* first, double* last, double b);

    //! A simple (and empty) destructor.
    ~Hyperplane();

    //! Access operator for the hyperplane's coefficients (except b).
    /*!
     *  \param pos The position (in the coefficients vector) to access.
     *  \return The hyperplane's a_{pos+1} coefficient. Remember coefficients 
     *          are labeled a_{1}, a_{2}, ..., a_{n}.
     *  
     *  - May throw a NonExistentCoefficientException if the requested 
     *    coefficient does not exist (pos is out of bounds).
     *  
     *  \sa Hyperplane
     */
    double operator[](unsigned int pos) const;

    /*!
     *  \brief Return the Hyperplane instance's b coefficient. (the hyperplane 
     *         equation's right hand side)
     */
    double b() const;

    //! Return the dimension of the space the hyperplane lives in.
    unsigned int space_dimension() const;

  private:
    //! A std::vector of all the a_{i} coefficients.
    std::vector<double> coefficients_;
    //! The right hand side of the hyperplane equation.
    double b_;
};


}  // namespace pareto_approximator


#endif  // HYPERPLANE_H
