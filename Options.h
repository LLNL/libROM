/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A class defining the parameters to the BasisGenerator class and
// SVD algorithm.

#ifndef included_Options_h
#define included_Options_h

namespace CAROM {

/**
 * A class defining the parameters to the BasisGenerator class and SVD algorithm.
 */
class Options
{
  /**
   * @brief Constructor.
   *
   * @pre dim > 0
   * @pre samples_per_time_interval > 0
   * @pre max_basis_dimension > 0
   * @pre singular_value_tol >= 0
   *
   * @param[in] dim The dimension of the system on this processor.
   * @param[in] samples_per_time_interval The maximum number of samples in
   *                                      each time interval.
   * @param[in] max_time_intervals The maximum number of time intervals.
   * @param[in] update_right_SV Whether to update the right SV or not.
   * @param[in] write_snapshots Whether to automatically write snapshots matrices
   *                        instead of basis matrices.
   * @param[in] max_basis_dimension (typemax(int)) The maximum number of
   *                                vectors returned in the basis.
   * @param[in] singular_value_tol Tolerance to determine whether to include
   *                               a singular value in the SVD.
   * @param[in] debug_algorithm If true results of static Options algorithm
   *                            will be printed to facilitate debugging.
   *
   * Static SVD
   *
   * If both max_basis_dimension and singular_value_tol would result in
   * truncating the basis, the dimension of the returned basis will be the
   * *minimum* of the number of vectors that is computed from each.
   *
   *
   * Incremental SVD
   *
   * @pre linearity_tol > 0.0
   * @pre initial_dt > 0.0
   * @pre sampling_tol > 0.0
   * @pre max_time_between_samples > 0.0
   * @pre min_sampling_time_step_scale >= 0.0
   * @pre sampling_time_step_scale >= 0.0
   * @pre max_sampling_time_step_scale >= 0.0
   * @pre min_sampling_time_step_scale <= max_sampling_time_step_scale
   *
   * @param[in] linearity_tol Tolerance to determine whether or not a
   *                          sample is linearly dependent.
   * @param[in] initial_dt Initial simulation time step.
   * @param[in] sampling_tol Sampling control tolerance.  Limits error in
   *                         projection of solution into reduced order
   *                         space.
   * @param[in] max_time_between_samples Upper bound on time between
   *                                     samples.
   * @param[in] fast_update If true use the fast update algorithm.
   * @param[in] skip_linearly_dependent If true skip linearly dependent
   *                                    samples.
   * @param[in] save_state If true the state of the SVD will be written to
   *                       disk when the object is deleted.  If there are
   *                       multiple time intervals then the state will not
   *                       be saved as restoring such a state makes no
   *                       sense.
   * @param[in] restore_state If true the state of the SVD will be restored
   *                          when the object is created.
   * @param[in] min_sampling_time_step_scale Minimum overall scale factor
   *                                         to apply to time step.
   * @param[in] sampling_time_step_scale Scale factor to apply to sampling
   *                                     algorithm.
   * @param[in] max_sampling_time_step_scale Maximum overall scale factor
   *                                         to apply to time step.
   */
public:

   Options(int dim_,
     int samples_per_time_interval_,
     int max_time_intervals_ = -1,
     bool update_right_SV_ = false,
     bool write_snapshots_ = false
   ): dim(dim_),
   max_basis_dimension(dim_),
   samples_per_time_interval(samples_per_time_interval_),
   max_time_intervals(max_time_intervals_),
   update_right_SV(update_right_SV_),
   write_snapshots(write_snapshots_) {};

   Options setMaxBasisDimension(
      int max_basis_dimension_
   )
   {
       max_basis_dimension = max_basis_dimension_;
       return *this;
   }

   Options setSingularValueTol(
       double singular_value_tol_
   )
   {
       singular_value_tol = singular_value_tol_;
       return *this;
   }

   Options setDebugMode(
        bool debug_algorithm_
      )
   {
      debug_algorithm = debug_algorithm_;
      return *this;
   }

   Options setIncrementalSVD(
        double linearity_tol_,
        double initial_dt_,
        double sampling_tol_,
        double max_time_between_samples_,
        bool fast_update_ = false,
        bool skip_linearly_dependent_ = false
      )
   {
      linearity_tol = linearity_tol_;
      initial_dt = initial_dt_;
      sampling_tol = sampling_tol_;
      max_time_between_samples = max_time_between_samples_;
      fast_update = fast_update_;
      skip_linearly_dependent = skip_linearly_dependent_;
      return *this;
   }

   Options setStateIO(
        bool save_state_,
        bool restore_state_
      )
   {
      save_state = save_state_;
      restore_state = restore_state_;
      return *this;
   }

   Options setSamplingTimeStepScale(
        double min_sampling_time_step_scale_,
        double sampling_time_step_scale_,
        double max_sampling_time_step_scale_
      )
   {
      min_sampling_time_step_scale = min_sampling_time_step_scale_;
      sampling_time_step_scale = sampling_time_step_scale_;
      max_sampling_time_step_scale = max_sampling_time_step_scale_;
      return *this;
   }

   int dim = -1;
   int samples_per_time_interval = -1;
   int max_time_intervals = -1;
   bool update_right_SV = false;
   bool write_snapshots = false;

   int max_basis_dimension = -1;
   double singular_value_tol = 0;
   bool debug_algorithm = false;

   // Incremental SVD
   double linearity_tol = -1;
   double initial_dt = -1;
   double sampling_tol = -1;
   double max_time_between_samples = -1;
   bool fast_update = false;
   bool skip_linearly_dependent = false;

   bool save_state = false;
   bool restore_state= false;

   double min_sampling_time_step_scale = 0.1;
   double sampling_time_step_scale = 0.8;
   double max_sampling_time_step_scale = 5.0;
};

}

#endif
