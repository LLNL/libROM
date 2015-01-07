/******************************************************************************
 *
 * This file is part of the CAROM distribution.  For full copyright
 * information, see COPYRIGHT.
 *
 * Copyright:   (c) 2013-2014 Lawrence Livermore National Security, LLC
 * Description: A class implementing interface of SVD for the static SVD
 *              algorithm.
 *
 *****************************************************************************/

#ifndef included_StaticSVD_h
#define included_StaticSVD_h

#include "SVD.h"

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
       * @param[in] sample_per_time_interval The maximium number of samples
       *                                     collected in a time interval.
       * @param[in] debug_rom If true results of svd will be printed to
       *                      facilitate debugging.
       */
      StaticSVD(
         int dim,
         int samples_per_time_interval,
         bool debug_rom = false);

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
       */
      virtual
      void
      takeSample(
         const double* u_in,
         double time);

      /**
       * @brief Returns the basis vectors for the current time interval.
       *
       * @post thisIntervalBasisCurrent()
       *
       * @return The basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getBasis();

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
       * @brief Preforms the actual SVD via lapack.
       *
       * Given a matrix, A, computes the 3 components of the singular value
       * decomposition.
       *
       * @pre A != 0
       * @pre total_dim > 0
       *
       * @param[in] A The globalized system whose SVD will be computed.
       * @param[in] total_dim The total dimension of the system that has been
       *                      distributed of multiple processors.
       */
      void
      svd(
         double* A,
         int total_dim);

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
      std::vector<double*> d_samples;

      /**
       * @brief The globalized matrix U.
       *
       * U is large and each process owns all of U.
       */
      Matrix* d_U;

      /**
       * @brief The globalized matrix S.
       *
       * S is small and each process owns all of S.
       */
      Matrix* d_S;

      /**
       * @brief The globalized matrix L.
       *
       * L is small and each process owns all of L.
       */
      Matrix* d_V;

      /**
       * @brief Flag to indicate if the basis vectors for the current time
       * interval are up to date.
       */
      bool d_this_interval_basis_current;

      /**
       * @brief MPI message tag.
       */
      static const int COMMUNICATE_A;
};

}

#endif
