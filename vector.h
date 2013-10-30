#ifndef included_vector_h
#define included_vector_h

#include "CAROM_config.h"

#include <assert.h>

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

      // Constructor in which vector is given storage managed by some other
      // entity.
      Vector(
         double* vec,
         int dim,
         bool distributed,
         int rank,
         int num_procs);

      // Destructor.
      ~Vector();

      // Dot product.
      double
      dot(const Vector& other) const;

      // Takes the norm of this.
      double
      norm() const;

      // Normalizes the vector and returns its norm.
      double
      normalize();

      // Const vector member access.
      const double&
      item(const int i) const
      {
#ifdef DEBUG_CHECK_ASSERTIONS
         assert((0 <= i) && (i < d_dim));
#endif
         return d_vec[i];
      }

      // Non-const vector member access.
      double&
      item(const int i)
      {
#ifdef DEBUG_CHECK_ASSERTIONS
         assert((0 <= i) && (i < d_dim));
#endif
         return d_vec[i];
      }
      
   private:
      friend class Matrix;

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
      int d_distributed;

      // The MPI rank of the process owning this object.
      int d_rank;

      // The number of MPI processes.
      int d_num_procs;

      // True if the object manages its storage.
      bool d_manages_storage;
};

}

#endif
