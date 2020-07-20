/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete implementation of the incremental SVD algorithm
//              using Matthew Brand's "fast update" method.

#ifndef included_IncrementalSVDFastUpdate_h
#define included_IncrementalSVDFastUpdate_h

#include "IncrementalSVD.h"

namespace CAROM {

/**
 * IncrementalSVDFastUpdate implements Brand's fast update incremental SVD
 * algorithm by implementing the pure virtual methods of the IncrementalSVD
 * base class.
 */
class IncrementalSVDFastUpdate : public IncrementalSVD
{
   public:
      /**
       * @brief Constructor.
       *
       * @pre dim > 0
       * @pre linearity_tol > 0.0
       * @pre samples_per_time_interval > 0
       *
       * @param[in] dim The dimension of the system on this processor.
       * @param[in] linearity_tol Tolerance to determine whether or not a
       *                          sample is linearly dependent.
       * @param[in] skip_linearly_dependent If true skip linearly dependent
       *                                    samples.
       * @param[in] samples_per_time_interval The number of samples to be
       *                                      collected for each time interval.
       * @param[in] basis_file_name The base part of the name of the file
       *                            containing the basis vectors.  Each process
       *                            will append its process ID to this base
       *                            name.
       * @param[in] save_state If true the state of the SVD will be written to
       *                       disk when the object is deleted.  If there are
       *                       multiple time intervals then the state will not
       *                       be saved as restoring such a state makes no
       *                       sense.
       * @param[in] restore_state If true the state of the SVD will be restored
       *                          when the object is created.
       * @param[in] debug_algorithm If true results of the algorithm will be
       *                            printed to facilitate debugging.
       */
      IncrementalSVDFastUpdate(
         int dim,
         double linearity_tol,
         double singular_value_tol,
         bool skip_linearly_dependent,
         int max_basis_dimension,
         int samples_per_time_interval,
         const std::string& basis_file_name,
         bool save_state = false,
         bool restore_state = false,
         bool updateRightSV = false,
         bool debug_algorithm = false);

      /**
       * @brief Destructor.
       */
      ~IncrementalSVDFastUpdate();

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      IncrementalSVDFastUpdate();

      /**
       * @brief Unimplemented copy constructor.
       */
      IncrementalSVDFastUpdate(
         const IncrementalSVDFastUpdate& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      IncrementalSVDFastUpdate&
      operator = (
         const IncrementalSVDFastUpdate& rhs);

      /**
       * @brief Constructs the first svd.
       *
       * @pre u != 0
       * @pre time >= 0.0
       *
       * @param[in] u The first state.
       * @param[in] time The simulation time for the first state.
       */
      virtual
      void
      buildInitialSVD(
         double* u,
         double time);

      /**
       * @brief Computes the current basis vectors.
       */
      virtual
      void
      computeBasis();

      /**
       * @brief Add a linearly dependent sample to the svd.
       *
       * @pre A != 0
       * @pre sigma != 0
       *
       * @param[in] A The left singular vectors.
       * @param[in] W The right singular vectors.
       * @param[in] sigma The singular values.
       */
      void
      addLinearlyDependentSample(
         const Matrix* A,
         const Matrix* W,
         const Matrix* sigma);

      /**
       * @brief Add a new, unique sample to the svd.
       *
       * @pre j != 0
       * @pre A != 0
       * @pre W != 0
       * @pre sigma != 0
       *
       * @param[in] j The new column of d_U.
       * @param[in] A The left singular vectors.
       * @param[in] W The right singular vectors.
       * @param[in] sigma The singular values.
       */
      void
      addNewSample(
         const Vector* j,
         const Matrix* A,
         const Matrix* W,
         Matrix* sigma);

      /**
       * @brief The matrix U'.
       *
       * U' is not distributed and the entire matrix exists on each processor.
       */
      Matrix* d_Up;

      /**
       * @brief 
       *
       * singular_value_tol is a tolerance value used to remove small singular values
       */
      double d_singular_value_tol;

};

}

#endif
