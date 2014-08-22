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

// A simple vector class in which the dimensions may be distributed evenly
// across multiple processes.  Supports only the basic inner product and norm
// related operations that we need in the CAROM project.
class Vector
{
   public:
      // Constructor in which vector manages storage.
      Vector(
         int dim,
         bool distributed,
         int rank,
         int num_procs);

      // Constructor in which vector is given its values.
      Vector(
         const double* vec,
         int dim,
         bool distributed,
         int rank,
         int num_procs);

      // Destructor.
      ~Vector();

      // Returns true if vector is distributed.
      bool
      distributed() const
      {
         return d_distributed;
      }

      // Returns the dimension of the vector.
      int
      dim() const
      {
         return d_dim;
      }

      // Dot product.
      double
      dot(const Vector& other) const;

      // Dot product.
      double
      dot(const Vector* const other) const
      {
         return dot(*other);
      }

      // Takes the norm of this.
      double
      norm() const;

      // Normalizes the vector and returns its norm.
      double
      normalize();

      // Adds other and this and returns the result.
      Vector*
      add(const Vector& other) const;

      // Adds other and this and returns the result.
      Vector*
      add(const Vector* const other) const
      {
         return add(*other);
      }

      // Subtracts other and this and returns the result.
      Vector*
      subtract(const Vector& other) const;

      // Subtracts other and this and returns the result.
      Vector*
      subtract(const Vector* const other) const
      {
         return subtract(*other);
      }

      // Const vector member access.
      const double&
      item(
         int i) const
      {
         CAROM_ASSERT((0 <= i) && (i < d_dim));
         return d_vec[i];
      }

      // Non-const vector member access.
      double&
      item(
         int i)
      {
         CAROM_ASSERT((0 <= i) && (i < d_dim));
         return d_vec[i];
      }
      
   private:
      // Default constructor is not implemented.
      Vector();

      // Copy constructor is not implemented.
      Vector(
         const Vector& other);

      // Assignment operator is not implemented.
      Vector&
      operator = (
         const Vector& rhs);

      // The storage for the vector.
      double* d_vec;

      // The dimension of the vector.
      int d_dim;

      // If true, the vector's dimensions are distributed over all processors.
      // Each processor holds the same number of dimensions.
      bool d_distributed;

      // The MPI rank of the process owning this object.
      int d_rank;

      // The number of MPI processes.
      int d_num_procs;
};

}

#endif
