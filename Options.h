/******************************************************************************
 *
 * Copyright (c) 2013-2019, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: An abstract class defining the interface to the generic Options
//              algorithm.

#ifndef included_Options_h
#define included_Options_h

#include <limits>

namespace CAROM {

class Options
{
  /**
   * @brief Constructor.
   *
   * @pre dim > 0
   * @pre samples_per_time_interval > 0
   * @pre max_basis_dimension > 0
   * @pre max_basis_dimension <= dim
   *
   * @param[in] dim The dimension of the system on this processor.
   * @param[in] samples_per_time_interval The maximum number of samples in
   *                                      each time interval.
   * @param[in] maximum basis dimension (typemax(int)) The maximum number of
   *                                vectors returned in the basis.
   * @param[in] update_rightSV Whether to update the right SV or not.
   * @param[in] max_time_intervals The maximum number of time intervals.
   * @param[in] debug_algorithm If true results of static Options algorithm
   *                            will be printed to facilitate debugging.
   * @param[in] write_snapshots Whether to automatically write snapshots matrices
   *                        instead of basis matrices.
   *
   * Static SVD
   * @pre sigma_tolerance >= 0
   *
   * @param[in] sigma_tolerance This tolerance is based on the ratio of
   *                            singular values to the largest singular
   *                            value. If sigma[i] / sigma[0] < sigma_tolerance,
   *                            the associated vector is dropped from the
   *                            basis.
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
   * @param[in] linearity_tol Tolerance to determine whether or not a
   *                          sample is linearly dependent.
   * @param[in] skip_linearly_dependent If true skip linearly dependent
   *                                    samples.
   * @param[in] fast_update If true use the fast update algorithm.
   * @param[in] initial_dt Initial simulation time step.
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
   */
public:

   Options(int dim_,
     int samples_per_time_interval_,
     int max_basis_dimension_ = 0,
     double singular_value_tol_ = 0.0,
     bool update_right_SV_ = false,
     int max_time_intervals_ = -1,
     bool write_snapshots_ = false
   ): dim(dim_),
   samples_per_time_interval(samples_per_time_interval_),
   singular_value_tol(singular_value_tol_),
   update_right_SV(update_right_SV_),
   max_time_intervals(max_time_intervals_),
   write_snapshots(write_snapshots_) {
     if (max_basis_dimension_ > 0) {
       max_basis_dimension = max_basis_dimension_;
     }
   };

   Options setDebugMode(
        bool debug_algorithm_
      )
   {
      debug_algorithm = debug_algorithm_;
      return *this;
   }

   Options setIncrementalSVD(
        double linearity_tol_,
        int max_basis_dimension_,
        double initial_dt_,
        double sampling_tol_,
        double max_time_between_samples_,
        bool skip_linearly_dependent_ = false,
        bool fast_update_ = false
      )
   {
      linearity_tol = linearity_tol_;
      max_basis_dimension = max_basis_dimension_;
      initial_dt = initial_dt_;
      sampling_tol = sampling_tol_;
      max_time_between_samples = max_time_between_samples_;
      skip_linearly_dependent = skip_linearly_dependent_;
      fast_update = fast_update_;
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
   int max_basis_dimension = std::numeric_limits<int>::max();
   double singular_value_tol = 0;
   bool update_right_SV = false;
   int max_time_intervals = -1;
   bool write_snapshots = false;
   bool debug_algorithm = false;

   // Incremental SVD
   double linearity_tol = -1;
   double initial_dt = -1;
   double sampling_tol = -1;
   double max_time_between_samples = -1;
   bool skip_linearly_dependent = false;
   bool fast_update = false;
   bool save_state = false;
   bool restore_state= false;
   double min_sampling_time_step_scale = 0.1;
   double sampling_time_step_scale = 0.8;
   double max_sampling_time_step_scale = 5.0;
};

}

#endif
