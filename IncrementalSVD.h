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

// Description: The abstract incremental SVD algorithm defines algorithm
//              interface.

#ifndef included_IncrementalSVD_h
#define included_IncrementalSVD_h

#include "SVD.h"
#include "Database.h"

namespace CAROM {

/**
 * IncrementalSVD is an abstract class defining the internal API of the
 * incremental SVD algorithm.
 */
  class IncrementalSVD : public SVD
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
       * @param[in] debug_algorithm If true results of algorithm will be
       *                            printed to facilitate debugging.
       */
      IncrementalSVD(
         int dim,
         double linearity_tol,
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
      virtual
      ~IncrementalSVD();

      /**
       * @brief Sample new state, u_in, at the given time.
       *
       * @pre u_in != 0
       * @pre time >= 0.0
       *
       * @param[in] u_in The state at the specified time.
       * @param[in] time The simulation time for the state.
       * @param[in] add_without_increase If true, addLinearlyDependent is invoked.  
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
       * @return The basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getBasis();

      /**
       * @brief Returns the temporal basis vectors for the current time interval.
       *
       * @return The temporal basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getTBasis();

      /**
       * @brief Returns the singular values for the current time interval.
       *
       * @return The singular values for the current time interval.
       */
      virtual
      const Matrix*
      getSingularValues();

   protected:
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
         double time) = 0;

      /**
       * @brief Adds the new sampled the state vector, u, to the system.
       *
       * @pre u != 0
       *
       * @param[in] u The new state.
       *
       * @return True if building the incremental svd was successful.
       */
      virtual
      bool
      buildIncrementalSVD(
         double* u, bool add_without_increase = false);

      /**
       * @brief Computes the current basis vectors.
       */
      virtual
      void
      computeBasis() = 0;

      /**
       * @brief Construct the matrix Q whose svd is needed.
       *
       * @pre l != 0
       * @pre l.dim() == numSamples()
       *
       * @param[out] Q The matrix to be constructed. [d_S,l; 0,k]
       * @param[in] l The last column of Q.
       * @param[in] k The lower right element of Q.
       */
      void
      constructQ(
         double*& Q,
         const Vector* l,
         double k);

      /**
       * @brief Given a matrix, A, returns 2 of the 3 components of its
       * singular value decomposition.
       *
       * The right singular vectors are not needed and therefore not returned.
       *
       * @pre A != 0
       *
       * @param[in] A The matrix whose svd is needed.
       * @param[out] U The left singular vectors of A.
       * @param[out] S The singular values of A.
       * @param[out] V The right singular vectors of A.
       *
       * @return True if the svd succeeded.
       */
      bool
      svd(
         double* A,
         Matrix*& U,
         Matrix*& S,
         Matrix*& V);

      /**
       * Add a linearly dependent sample to the svd.
       *
       * @pre A != 0
       * @pre sigma != 0
       *
       * @param[in] A The left singular vectors.
       * @param[in] W The right singular vectors.
       * @param[in] sigma The singular values.
       */
      virtual
      void
      addLinearlyDependentSample(
         const Matrix* A,
         const Matrix* W,
         const Matrix* sigma) = 0;

      /**
       * @brief Add a new, unique sample to the svd.
       *
       * @pre j != 0
       * @pre A != 0
       * @pre sigma != 0
       *
       * @param[in] j The new column of d_U.
       * @param[in] A The left singular vectors.
       * @param[in] sigma The singular values.
       */
      virtual
      void
      addNewSample(
         const Vector* j,
         const Matrix* A,
         const Matrix* W,
         Matrix* sigma) = 0;

      /**
       * @brief The number of samples stored.
       *
       * @return The number of samples stored.
       */
      int
      numSamples()
      {
         return d_num_samples;
      }

      /**
       * @brief Computes and returns the orthogonality of m.
       *
       * @pre m != 0
       *
       * @param[in] m The matrix to check.
       *
       * @return The orthogonality of m.
       */
      double
      checkOrthogonality(
        const Matrix* m);

      /**
       * @brief Reorthogonalizes m.
       *
       * @pre m != 0
       *
       * @param[in,out] m The matrix to reorthogonalize.
       */
      void
      reOrthogonalize(
         Matrix* m);

      /**
       * @brief Tolerance to determine whether or not a sample is linearly
       * dependent.
       */
      double d_linearity_tol;

      /**
       * @brief If true, skip linearly dependent samples.
       */
      bool d_skip_linearly_dependent;

      /**
       * @brief the maximum basis dimension
       */
      int d_max_basis_dimension;

      /**
       * @brief Total number of processors.
       */
      int d_size;

      /**
       * @brief The rank of the processor owning this object.
       */
      int d_rank;

      /**
       * @brief Dimension of the system on each processor.
       */
      std::vector<int> d_proc_dims;

      /**
       * @brief The total dimension of the system.
       */
      long int d_total_dim;

      /**
       * @brief If true the state of the SVD will be written to disk when the
       * object is deleted.
       *
       * If there are multiple time intervals then the state will not be saved
       * as restoring such a state makes no sense.
       */
      bool d_save_state;

      /**
       * @brief If true the right singular vectors will be updated 
       */
      bool d_updateRightSV;

      /**
       * @brief Pointer to the database that will hold saved state data if the
       * state is to be saved.
       */
      Database* d_state_database;

      /**
       * @brief Name of file to which state is save to or restored from.
       */
      std::string d_state_file_name;

      /**
       * @brief MPI message tag.
       */
      static const int COMMUNICATE_U;

   private:
      /**
       * @brief Unimplemented default constructor.
       */
      IncrementalSVD();

      /**
       * @brief Unimplemented copy constructor.
       */
      IncrementalSVD(
         const IncrementalSVD& other);

      /**
       * @brief Unimplemented assignment operator.
       */
      IncrementalSVD&
      operator = (
         const IncrementalSVD& rhs);
};

}

#endif
