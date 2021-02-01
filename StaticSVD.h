/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A class implementing interface of SVD for the static SVD
//              algorithm.

#ifndef included_StaticSVD_h
#define included_StaticSVD_h

#include "SVD.h"
#include "Options.h"
#include "scalapack_wrapper.h"

#include <limits>
#include <memory>
#include <vector>

namespace CAROM {

/**
 * StaticSVD implements the interface of class SVD for the static SVD
 * algorithm.  This algorithm is not scalable and is intended primarily as a
 * sanity check of the incremental svd algorithm.
 */
class StaticSVD : public SVD
{
   public:
      /**
       * Destructor.
       */
      ~StaticSVD();

      /**
       * @brief Collect the new sample, u_in at supplied time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The new sample.
       * @param[in] time The simulation time of the new sample.
       * @param[in] add_without_increase If true, then addLinearlyDependent will be invoked
       *
       * @return True if the sampling was successful.
       */
      virtual
      bool
      takeSample(
         double* u_in,
         double time,
         bool add_without_increase = false);

      /**
       * @brief Returns the basis vectors for the current time interval.
       *
       * @post thisIntervalBasisCurrent()
       *
       * @return The basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getSpatialBasis();

      /**
       * @brief Returns the temporal basis vectors for the current time interval.
       *
       * @post thisIntervalBasisCurrent()
       *
       * @return The temporal basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getTemporalBasis();

      /**
       * @brief Returns the singular values for the current time interval.
       *
       * @post thisIntervalBasisCurrent()
       *
       * @return The singular values for the current time interval.
       */
      virtual
      const Vector*
      getSingularValues();

      /**
       * @brief Returns the snapshot matrix for the current time interval.
       *
       * @return The snapshot matrix for the current time interval.
       */
      virtual
      const Matrix*
      getSnapshotMatrix();

   protected:

      /**
        * @brief Constructor.
        *
        * @param[in] options The struct containing the options for this SVD
        *                    implementation.
        */
       StaticSVD(
          Options options
          );

      /**
       * @brief Gathers samples from all other processors to form complete
       * sample of system and computes the SVD.
       */
      virtual void
      computeSVD();

      /**
       * @brief Tells if the basis vectors for this time interval are up to
       * date.
       *
       * @return True if the basis vectors for this time interval are up to
       * date.
       */
      bool
      thisIntervalBasisCurrent()
      {
         return d_this_interval_basis_current;
      }

      /**
       * @brief Current samples of the system.
       */
      std::unique_ptr<SLPK_Matrix> d_samples;

      /**
       * @brief Factorization manager object used to compute the SVD
       */
      std::unique_ptr<SVDManager> d_factorizer;

      /**
       * @brief Flag to indicate if the basis vectors for the current time
       * interval are up to date.
       */
      bool d_this_interval_basis_current;

      /**
       * @brief The rank of the process this object belongs to.
       */
      int d_rank;

      /**
       * @brief The number of processors being run on.
       */
      int d_num_procs;

      /**
       * @brief The starting row (0-based) of the matrix that each process owns.
       * Stored to avoid an MPI operation to get this operation every time we
       * scatter a sample.
       */
      std::vector<int> d_istarts;

      /**
       * @brief The number of rows that each process owns. Stored to avoid
       * an MPI operation to get this operation every time we scatter a sample.
       */
      std::vector<int> d_dims;

      /**
       * @brief The total dimension of the system (row dimension)
       */
      int d_total_dim;

      /**
       * @brief The number of processor rows and processor columns in the grid.
       */
      int d_nprow;
      int d_npcol;

      /**
       * @brief The block size used internally for computing the SVD.
       */
      int d_blocksize;

      /**
       * @brief Get the system's total row dimension and where my rows sit in
       * the matrix.
       */
      void get_global_info();

      /**
       * @brief The max number of basis vectors to return.
       */
      int d_max_basis_dimension;

      /**
       * @brief The tolerance for singular values below which to drop vectors
       */
      double d_singular_value_tol;

      void delete_samples();
      void delete_factorizer();

      void broadcast_sample(const double* u_in);

    private:

       friend class BasisGenerator;

       /**
        * @brief Unimplemented default constructor.
        */
       StaticSVD();

       /**
        * @brief Unimplemented copy constructor.
        */
       StaticSVD(
          const StaticSVD& other);

       /**
        * @brief Unimplemented assignment operator.
        */
       StaticSVD&
       operator = (
          const StaticSVD& rhs);
};

}

#endif
