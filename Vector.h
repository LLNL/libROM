/******************************************************************************
 *
 * Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory
 * Written by William Arrighi wjarrighi@llnl.gov
 * CODE-686965
 * All rights reserved.
 *
 * This file is part of libROM.
 * For details, see https://computation.llnl.gov/librom
 * Please also read README_BSD_NOTICE.
 *
 * Redistribution and use in source and binary forms, with or without
 * modifications, are permitted provided that the following conditions are met:
 *
 *    o Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the disclaimer below.
 *    o Redistribution in binary form must reproduce the above copyright
 *      notice, this list of conditions and the disclaimer (as noted below) in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *    o Neither the name of the LLNS/LLNL nor the names of its contributors may
 *      be used to endorse or promote products derived from this software
 *      without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
 * LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR SUCH DAMAGE.
 *
 *****************************************************************************/

// Description: A simple, parallel Vector class with the utility needed to
//              support the basis generation methods of this library.  A
//              distributed Vector has its rows distributed across processors.

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
       * @brief Sets the length of the vector and reallocates storage if
       * needed.
       *
       * @param[in] num_rows New number of rows
       * @param[in] num_cols New number of cols
       */
      void
      setSize(
         int dim)
      {
         if (dim > d_alloc_size) {
            if (d_vec) {
               delete [] d_vec;
            }
            d_vec = new double [dim];
            d_alloc_size = dim;
         }
         d_dim = dim;
      }

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
       * @pre dim() == other.dim()
       * @pre distributed() == other.distributed()
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
       * @pre dim() == other->dim()
       * @pre distributed() == other->distributed()
       *
       * @param[in] other The Vector to form the inner product with this.
       *
       * @return The inner product of this and other.
       */
      double
      inner_product(
         const Vector* other) const
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
       * @pre distributed() == other.distributed()
       * @pre dim() == other.dim()
       *
       * @param[in] other The other summand.
       *
       * @return this + other
       */
      Vector*
      plus(
         const Vector& other) const
      {
         Vector* result = 0;
         plus(other, result);
         return result;
      }

      /**
       * @brief Adds other and this and returns the result, pointer version.
       *
       * @pre other != 0
       * @pre distributed() == other->distributed()
       * @pre dim() == other->dim()
       *
       * @param[in] other The other summand.
       *
       * @return this + other
       */
      Vector*
      plus(
         const Vector* other) const
      {
         CAROM_ASSERT(other != 0);
         return plus(*other);
      }

      /**
       * @brief Adds other and this and fills result with the answer.
       *
       * @pre result == 0 || result->distributed() == distributed()
       * @pre distributed() == other.distributed()
       * @pre dim() == other.dim()
       *
       * @param[in] other The other summand.
       * @param[out] result this + other
       */
      void
      plus(
         const Vector& other,
         Vector*& result) const;

      /**
       * @brief Subtracts other and this and returns the result, reference
       * version.
       *
       * @pre distributed() == other.distributed()
       * @pre dim() == other.dim()
       *
       * @param[in] other The other subtrahand.
       *
       * @return this - other
       */
      Vector*
      minus(
         const Vector& other) const
      {
         Vector* result = 0;
         minus(other, result);
         return result;
      }

      /**
       * @brief Subtracts other and this and returns the result, pointer
       * version.
       *
       * @pre other != 0
       * @pre distributed() == other->distributed()
       * @pre dim() == other->dim()
       *
       * @param[in] other The other subtrahand.
       *
       * @return this - other
       */
      Vector*
      minus(
         const Vector* other) const
      {
         CAROM_ASSERT(other != 0);
         return minus(*other);
      }

      /**
       * @brief Subtracts other and this and fills result with the answer.
       *
       * @pre result == 0 || result->distributed() == distributed()
       * @pre distributed() == other.distributed()
       * @pre dim() == other.dim()
       *
       * @param[in] other The other subtrahend.
       * @param[out] result this - other
       */
      void
      minus(
         const Vector& other,
         Vector*& result) const;

      /**
       * @brief Const Vector member access.
       *
       * @pre (0 <= i) && (i < dim())
       *
       * @param[in] i The component of the Vector on this processor requested.
       *
       * @return The requested component of the Vector on this processor.
       */
      const double&
      item(
         int i) const
      {
         CAROM_ASSERT((0 <= i) && (i < dim()));
         return d_vec[i];
      }

      /**
       * @brief Non-const Vector member access.
       *
       * Allows constructs of the form vec[i] = val;
       *
       * @pre (0 <= i) && (i < dim())
       *
       * @param[in] i The component of the Vector on this processor requested.
       *
       * @return The requested component of the Vector on this processor.
       */
      double&
      item(
         int i)
      {
         CAROM_ASSERT((0 <= i) && (i < dim()));
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
       * @brief The currently allocated size.
       *
       * d_dim <= d_alloc_size
       */
      int d_alloc_size;

      /**
       * @brief If true, the Vector's dimensions are distributed over all
       * processors.
       *
       * Each processor does not need to hold the same number of dimensions.
       */
      bool d_distributed;

      /**
       * @brief The number of processors being run on.
       */
      int d_num_procs;
};

}

#endif
