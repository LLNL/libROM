/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
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
#include "scalapack_wrapper.h"

/* Use C++11 built-in shared pointers if available; else fallback to Boost. */
#if __cplusplus >= 201103L
#include <memory>
using std::shared_ptr;
#else
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;
#endif

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
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre sample_per_time_interval > 0
       *
       * @param[in] dim The dimension of the system distributed to this
       *                processor.
       * @param[in] samples_per_time_interval The maximium number of samples
       *                                      collected in a time interval.
       * @param[in] debug_algorithm If true results of the algorithm will be
       *                            printed to facilitate debugging.
       */
      StaticSVD(
         int dim,
         int samples_per_time_interval,
         bool debug_algorithm = false);

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
      const Matrix*
      getSingularValues();

   private:
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

      /**
       * @brief Gathers samples from all other processors to form complete
       * sample of system and computes the SVD.
       */
      void
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
      shared_ptr<SLPK_Matrix> d_samples;

      /**
       * @brief Factorization manager object used to compute the SVD
       */
      shared_ptr<SVDManager> d_factorizer;

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
       * @brief The starting row (0-based) of the matrix that I own.
       */
      int d_istart;

      /**
       * @brief The total dimension of the system (row dimension)
       */
      int d_total_dim;

      /**
       * @brief The number of processor rows and processor columns in the grid.
       */
      int d_nprow;
      int d_npcol;
      int d_blocksize;

      /**
       * @brief Get the system's total row dimension and where my rows sit in
       * the matrix.
       */
      void get_total_dim(int*, int*);

      /**
       * @brief MPI message tag.
       */
      static const int COMMUNICATE_A;

      void delete_samples();
      void delete_factorizer();

      void broadcast_sample(const double* u_in);
};

}

#endif
