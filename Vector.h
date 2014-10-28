/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: A simple, parallel Vector class with the utility needed to
 *              support the basis generation methods of this library.  A
 *              distributed Vector has its rows distributed across processors.
 *
 *****************************************************************************/

#ifndef included_Vector_h
#define included_Vector_h

#include "Utilities.h"

namespace CAROM {

/**
 * A simple vector class in which the dimensions may be distributed across
 * multiple processes.  This class supports only the basic operations that are
 * needed by the SVD library.
 */
class Vector
{
   public:
      /**
       * @brief Constructor creating a Vector with uninitialized values.
       *
       * @pre dim > 0
       *
       * @param[in] dim When undistributed, the total dimension of the Vector.
       *                When distributed, the part of the total dimension of
       *                the Vector on this processor.
       * @param[in] distributed If true the dimensions of the Vector are spread
       *                        over all processors.
       */
      Vector(
         int dim,
         bool distributed);

      /**
       * @brief Constructor in which the Vector is given its initial values.
       *
       * @pre vec != 0
       * @pre dim > 0
       *
       * @param[in] vec The initial values of the Vector.
       * @param[in] dim When undistributed, the total dimension of the Vector.
       *                When distributed, the part of the total dimension of
       *                the Vector on this processor.
       * @param[in] distributed If true the dimensions of the Vector are spread
       *                        over all processors.
       */
      Vector(
         const double* vec,
         int dim,
         bool distributed);

      /**
       * @brief Copy constructor.
       *
       * @param[in] other The Vector to copy.
       */
      Vector(
         const Vector& other);

      /**
       * @brief Destructor.
       */
      ~Vector();

      /**
       * @brief Assignment operator.
       *
       * @param[in] rhs The Vector to assign to this.
       *
       * @return This after rhs has been assigned to it.
       */
      Vector&
      operator = (
         const Vector& rhs);

      /**
       * @brief Returns true if the Vector is distributed.
       *
       * @return True if the Vector is distributed.
       */
      bool
      distributed() const
      {
         return d_distributed;
      }

      /**
       * @brief Returns the dimension of the Vector on this processor.
       *
       * @return The part of the Vector's dimension on this processor.
       */
      int
      dim() const
      {
         return d_dim;
      }

      /**
       * @brief Inner product, reference form.
       *
       * For distributed Vectors this is a parallel operation.
       *
       * @pre d_dim == other.d_dim
       * @pre d_distributed == other.d_distributed
       *
       * @param[in] other The Vector to form the inner product with this.
       *
       * @return The inner product of this and other.
       */
      double
      inner_product(
         const Vector& other) const;

      /**
       * @brief Inner product, pointer version.
       *
       * For distributed Vectors this is a parallel operation.
       *
       * @pre other != 0
       * @pre d_dim == other->d_dim
       * @pre d_distributed == other.d_distributed
       *
       * @param[in] other The Vector to form the inner product with this.
       *
       * @return The inner product of this and other.
       */
      double
      inner_product(
         const Vector* const other) const
      {
         CAROM_ASSERT(other != 0);
         return inner_product(*other);
      }

      /**
       * @brief Form the norm of this.
       *
       * For a distributed Vector this is a parallel operation.
       *
       * @return The norm of this.
       */
      double
      norm() const;

      /**
       * @brief Normalizes the Vector and returns its norm.
       *
       * For a distributed Vector this is a parallel operation.
       *
       * @return The norm of this.
       */
      double
      normalize();

      /**
       * @brief Adds other and this and returns the result, reference version.
       *
       * @pre d_distributed == other.d_distributed
       * @pre d_dim == other.d_dim
       *
       * @param[in] other The other summand.
       *
       * @return this + other
       */
      Vector*
      plus(
         const Vector& other) const;

      /**
       * @brief Adds other and this and returns the result, pointer version.
       *
       * @pre other != 0
       * @pre d_distributed == other->d_distributed
       * @pre d_dim == other->d_dim
       *
       * @param[in] other The other summand.
       *
       * @return this + other
       */
      Vector*
      plus(
         const Vector* const other) const
      {
         CAROM_ASSERT(other != 0);
         return plus(*other);
      }

      /**
       * @brief Subtracts other and this and returns the result, reference
       * version.
       *
       * @pre d_distributed == other.d_distributed
       * @pre d_dim == other.d_dim
       *
       * @param[in] other The other subtrahand.
       *
       * @return this - other
       */
      Vector*
      minus(
         const Vector& other) const;

      /**
       * @brief Subtracts other and this and returns the result, pointer
       * version.
       *
       * @pre other != 0
       * @pre d_distributed == other->d_distributed
       * @pre d_dim == other->d_dim
       *
       * @param[in] other The other subtrahand.
       *
       * @return this - other
       */
      Vector*
      minus(
         const Vector* const other) const
      {
         CAROM_ASSERT(other != 0);
         return minus(*other);
      }

      /**
       * @brief Const Vector member access.
       *
       * @pre (0 <= i) && (i < d_dim)
       *
       * @param[in] i The component of the Vector on this processor requested.
       *
       * @return The requested component of the Vector on this processor.
       */
      const double&
      item(
         int i) const
      {
         CAROM_ASSERT((0 <= i) && (i < d_dim));
         return d_vec[i];
      }

      /**
       * @brief Non-const Vector member access.
       *
       * Allows constructs of the form vec[i] = val;
       *
       * @pre (0 <= i) && (i < d_dim)
       *
       * @param[in] i The component of the Vector on this processor requested.
       *
       * @return The requested component of the Vector on this processor.
       */
      double&
      item(
         int i)
      {
         CAROM_ASSERT((0 <= i) && (i < d_dim));
         return d_vec[i];
      }
      
   private:
      /**
       * @brief Default constructor is not implemented.
       */
      Vector();

      /**
       * @brief The storage for the Vector's values on this processor.
       */
      double* d_vec;

      /**
       * @brief The part of the Vector's dimension on this processor.
       */
      int d_dim;

      /**
       * @brief If true, the Vector's dimensions are distributed over all
       * processors.
       *
       * Each processor does not need to hold the same number of dimensions.
       */
      bool d_distributed;
};

}

#endif
