/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The abstract incremental SVD algorithm defines algorithm
//              interface.

#ifndef included_IncrementalSVD_h
#define included_IncrementalSVD_h

#include "SVD.h"
#include "Database.h"

namespace CAROM {

struct IncrementalSVDOptions : virtual public SVDOptions
{
  /**
   * @brief Constructor.
   *
   * @pre dim > 0
   * @pre max_basis_dimension > 0
   * @pre max_basis_dimension <= dim
   * @pre linearity_tol > 0.0
   * @pre initial_dt > 0.0
   * @pre samples_per_time_interval > 0
   * @pre sampling_tol > 0.0
   * @pre max_time_between_samples > 0.0
   * @pre min_sampling_time_step_scale >= 0.0
   * @pre sampling_time_step_scale >= 0.0
   * @pre max_sampling_time_step_scale >= 0.0
   * @pre min_sampling_time_step_scale <= max_sampling_time_step_scale
   *
   * @param[in] dim The dimension of the system on this processor.
   * @param[in] linearity_tol Tolerance to determine whether or not a
   *                          sample is linearly dependent.
   * @param[in] skip_linearly_dependent If true skip linearly dependent
   *                                    samples.
   * @param[in] fast_update If true use the fast update algorithm.
   * @param[in] maximum basis dimension
   * @param[in] initial_dt Initial simulation time step.
   * @param[in] samples_per_time_interval The maximum number of samples in
   *                                      each time interval.
   * @param[in] sampling_tol Sampling control tolerance.  Limits error in
   *                         projection of solution into reduced order
   *                         space.
   * @param[in] max_time_between_samples Upper bound on time between
   *                                     samples.
   * @param[in] save_state If true the state of the SVD will be written to
   *                       disk when the object is deleted.  If there are
   *                       multiple time intervals then the state will not
   *                       be saved as restoring such a state makes no
   *                       sense.
   * @param[in] restore_state If true the state of the SVD will be restored
   *                          when the object is created.
   * @param[in] update_rightSV Whether to update the right SV or not.
   * @param[in] min_sampling_time_step_scale Minimum overall scale factor
   *                                         to apply to time step.
   * @param[in] sampling_time_step_scale Scale factor to apply to sampling
   *                                     algorithm.
   * @param[in] max_sampling_time_step_scale Maximum overall scale factor
   *                                         to apply to time step.
   * @param[in] debug_algorithm If true results of incremental svd
   *                            algorithm will be printed to facilitate
   *                            debugging.
   * @param[in] singular_value_tol Tolerance to determine whether or to include
   *                               a singular value in the SVD.
   * @param[in] max_time_intervals The maximum number of time intervals.
   * @param[in] write_snapshots Whether to automatically write snapshots matrices
   *                        instead of basis matrices.
   */

   IncrementalSVDOptions() = delete;

   IncrementalSVDOptions(int dim_,
     int samples_per_time_interval_,
     double linearity_tol_,
     int max_basis_dimension_,
     double initial_dt_,
     double sampling_tol_,
     double max_time_between_samples_,
     bool skip_linearly_dependent_ = false,
     bool fast_update_ = false,
     bool save_state_ = false,
     bool restore_state_ = false,
     bool updateRightSV_ = false,
     double min_sampling_time_step_scale_ = 0.1,
     double sampling_time_step_scale_ = 0.8,
     double max_sampling_time_step_scale_ = 5.0,
     bool debug_algorithm_ = false,
     double singular_value_tol_ = 0.0,
     int max_time_intervals_ = -1,
     bool write_snapshots_ = false
   ) : SVDOptions(dim_, samples_per_time_interval_, debug_algorithm_,
     max_time_intervals, write_snapshots_),
   linearity_tol(linearity_tol_),
   max_basis_dimension(max_basis_dimension_),
   initial_dt(initial_dt_),
   sampling_tol(sampling_tol_),
   max_time_between_samples(max_time_between_samples_),
   skip_linearly_dependent(skip_linearly_dependent_),
   fast_update(fast_update_),
   save_state(save_state_),
   restore_state(restore_state_),
   updateRightSV(updateRightSV_),
   min_sampling_time_step_scale(min_sampling_time_step_scale_),
   sampling_time_step_scale(sampling_time_step_scale_),
   max_sampling_time_step_scale(max_sampling_time_step_scale_),
   singular_value_tol(singular_value_tol_) {};

   double linearity_tol;
   int max_basis_dimension;
   double initial_dt;
   double sampling_tol;
   double max_time_between_samples;
   bool skip_linearly_dependent;
   bool fast_update;
   bool save_state;
   bool restore_state;
   bool updateRightSV;
   double min_sampling_time_step_scale;
   double sampling_time_step_scale;
   double max_sampling_time_step_scale;
   double singular_value_tol;
};

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
       * @param[in] options The struct containing the options for this basis
       *                    generator.
       * @param[in] basis_file_name The base part of the name of the file
       *                            containing the basis vectors.  Each process
       *                            will append its process ID to this base
       *                            name.
       */
      IncrementalSVD(
         IncrementalSVDOptions options,
         const std::string& basis_file_name);

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
      getSpatialBasis();

      /**
       * @brief Returns the temporal basis vectors for the current time interval.
       *
       * @return The temporal basis vectors for the current time interval.
       */
      virtual
      const Matrix*
      getTemporalBasis();

      /**
       * @brief Returns the singular values for the current time interval.
       *
       * @return The singular values for the current time interval.
       */
      virtual
      const Matrix*
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
