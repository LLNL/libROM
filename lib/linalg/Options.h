/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
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
 * Class Options defines the parameters to the BasisGenerator class and SVD
 *       algorithm. StaticSVD, RandomizedSVD, and IncrementalSVD parameters
 *       can be set.
 */
class Options
{
public:
    /**
     * @brief Constructor.
     *
     * @pre dim_ > 0
     * @pre samples_per_time_interval_ > 0
     * @pre max_time_intervals == -1 || max_time_intervals > 0
     *
     * @param[in] dim_ The dimension of the system on this processor.
     * @param[in] samples_per_time_interval_ The maximum number of samples in
     *                                      each time interval.
     * @param[in] max_time_intervals_ The maximum number of time intervals.
     * @param[in] update_right_SV_ Whether to update the right SV or not.
     * @param[in] write_snapshots_ Whether to automatically write snapshots matrices
     *                        instead of basis matrices.
     *
     */
    Options(int dim_,
            int samples_per_time_interval_,
            int max_time_intervals_ = -1,
            bool update_right_SV_ = false,
            bool write_snapshots_ = false
           ): dim(dim_),
        max_basis_dimension(samples_per_time_interval_),
        samples_per_time_interval(samples_per_time_interval_),
        max_time_intervals(max_time_intervals_),
        update_right_SV(update_right_SV_),
        write_snapshots(write_snapshots_) {};

    /**
     * @brief Sets the maximum basis dimension of the SVD algorithm.
     *
     * @pre max_basis_dimension_ > 0
     *
     * @param[in] max_basis_dimension_ The maximum basis dimension of the
     *                                 basis.
     */
    Options setMaxBasisDimension(
        int max_basis_dimension_
    )
    {
        max_basis_dimension = max_basis_dimension_;
        return *this;
    }

    /**
     * @brief Sets the singular value tolerance of the SVD algorithm.
     *
     * @pre singular_value_tol_ >= 0
     *
     * @param[in] singular_value_tol_ The singular value tolerance of the SVD
     *                                algorithm. If both max_basis_dimension and
     *                                singular_value_tol would result in
     *                                truncating the basis, the dimension of the
     *                                returned basis will be the minimum of the
     *                                number of vectors that is computed from each.
     */
    Options setSingularValueTol(
        double singular_value_tol_
    )
    {
        singular_value_tol = singular_value_tol_;
        return *this;
    }

    /**
     * @brief Sets whether debugging is turned on.
     *
     * @param[in] debug_algorithm_ Whether to turn on debugging.
     */
    Options setDebugMode(
        bool debug_algorithm_
    )
    {
        debug_algorithm = debug_algorithm_;
        return *this;
    }

    /**
     * @brief Sets the parameters of the randomized SVD algorithm.
     *
     * @param[in] randomized_ Whether to use randomization.
     * @param[in] randomized_subspace_dim_ The dimension of the randomized
                                           subspace
     * @param[in] random_seed_ The random seed used to initialize the algorithm.
     */
    Options setRandomizedSVD(
        bool randomized_,
        int randomized_subspace_dim_ = -1,
        int random_seed_ = 1
    )
    {
        randomized = randomized_;
        randomized_subspace_dim = randomized_subspace_dim_;
        random_seed = random_seed_;
        return *this;
    }

    /**
     * @brief Sets the essential parameters of the incremental SVD algorithm.
     *
     * @pre linearity_tol_ > 0.0
     * @pre initial_dt_ > 0.0
     * @pre sampling_tol_ > 0.0
     *
     * @param[in] linearity_tol_ Tolerance to determine whether or not a
     *                          sample is linearly dependent.
     * @param[in] initial_dt_ Initial simulation time step.
     * @param[in] sampling_tol_ Sampling control tolerance.  Limits error in
     *                         projection of solution into reduced order
     *                         space.
     * @param[in] max_time_between_samples_ Upper bound on time between
     *                                     samples.
     * @param[in] fast_update_ If true use the fast update algorithm.
     * @param[in] skip_linearly_dependent_ If true skip linearly dependent
     *                                    samples.
     */
    Options setIncrementalSVD(
        double linearity_tol_,
        double initial_dt_,
        double sampling_tol_,
        double max_time_between_samples_,
        bool fast_update_ = false,
        bool fast_update_brand_ = false,
        bool skip_linearly_dependent_ = false
    )
    {
        linearity_tol = linearity_tol_;
        initial_dt = initial_dt_;
        sampling_tol = sampling_tol_;
        max_time_between_samples = max_time_between_samples_;
        fast_update = fast_update_;
        fast_update_brand = fast_update_brand_;
        skip_linearly_dependent = skip_linearly_dependent_;
        return *this;
    }

    /**
     * @brief Sets the state IO parameters of the incremental SVD algorithm.
     *
     * @param[in] save_state_ If true the state of the SVD will be written to
     *                       disk when the object is deleted.  If there are
     *                       multiple time intervals then the state will not
     *                       be saved as restoring such a state makes no
     *                       sense.
     * @param[in] restore_state_ If true the state of the SVD will be restored
     *                          when the object is created.
     */
    Options setStateIO(
        bool save_state_,
        bool restore_state_
    )
    {
        save_state = save_state_;
        restore_state = restore_state_;
        return *this;
    }

    /**
     * @brief Sets the sampling time step scaling of the incremental SVD
     *        algorithm.
     *
     * @pre max_time_between_samples_ > 0.0
     * @pre min_sampling_time_step_scale_ >= 0.0
     * @pre sampling_time_step_scale_ >= 0.0
     * @pre max_sampling_time_step_scale_ >= 0.0
     * @pre min_sampling_time_step_scale_ <= max_sampling_time_step_scale
     *
     * @param[in] min_sampling_time_step_scale_ Minimum overall scale factor
     *                                         to apply to time step.
     * @param[in] sampling_time_step_scale_ Scale factor to apply to sampling
     *                                     algorithm.
     * @param[in] max_sampling_time_step_scale_ Maximum overall scale factor
     *                                         to apply to time step.
     */
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

    /**
     * @brief The dimension of the system on this processor.
     */
    int dim = -1;

    /**
     * @brief The number of samples per time interval.
     */
    int samples_per_time_interval = -1;

    /**
     * @brief The maximum number of time intervals.
     */
    int max_time_intervals = -1;

    /**
     * @brief Whether to update the right singular values.
     */
    bool update_right_SV = false;

    /**
     * @brief Whether to write snapshots or bases.
     */
    bool write_snapshots = false;

    /**
     * @brief The maximum dimension of the basis.
     */
    int max_basis_dimension = -1;

    /**
     * @brief The singular value tolerance used in the SVD algorithm.
     */
    double singular_value_tol = 0;

    /**
     * @brief Whether debugging is turned on (any randomness is turned off).
     */
    bool debug_algorithm = false;

    // Randomized SVD

    /**
     * @brief Whether to use the randomized SVD algorithm.
     */
    bool randomized = false;

    /**
     * @brief The dimension of the randomized subspace used in the randomized SVD
     *        algogrithm
     */
    int randomized_subspace_dim = -1;

    /**
     * @brief The random seed used to initialize the algorithm.
     */
    int random_seed = 1;

    // Incremental SVD

    /**
     * @brief The tolerance of the incremental SVD algorithm to determine whether
     *        or not a sample is linearly dependent.
     */
    double linearity_tol = -1;

    /**
     * @brief The initial simulation time step of the
     *        incremental SVD algorithm.
     */
    double initial_dt = -1;

    /**
     * @brief The sampling control tolerance of the
     *        incremental SVD algorithm. Limits error in projection of
     *        solution into reduced order space.
     */
    double sampling_tol = -1;

    /**
     * @brief The upper bound on time between samples of the
     *        incremental SVD algorithm.
     */
    double max_time_between_samples = -1;

    /**
     * @brief If true use the fast update version of the
     *        incremental SVD algorithm.
     */
    bool fast_update = false;

    /**
     * @brief If true use the exact Brand's algorithm for the
     *        incremental SVD.
     */
    bool fast_update_brand = false;

    /**
     * @brief If true skip linearly dependent samples of the
     *        incremental SVD algorithm.
     */
    bool skip_linearly_dependent = false;

    /**
     * @brief If true the state of the incremental SVD will be written to
     *        disk when the object is deleted.  If there are
     *        multiple time intervals then the state will not
     *        be saved as restoring such a state makes no
     *        sense.
     */
    bool save_state = false;

    /**
     * @brief If true the state of the incremental SVD will be restored when the
     *        object is created.
     */
    bool restore_state = false;

    /**
     * @brief The minimum overall scale factor to apply to time step of the
     *        incremental SVD algorithm.
     */
    double min_sampling_time_step_scale = 0.1;

    /**
     * @brief The scaling factor to apply to sampling algorithm of the
     *        incremental SVD algorithm.
     */
    double sampling_time_step_scale = 0.8;

    /**
     * @brief The maximum overall scale factor to apply to time step of the
     *        incremental SVD algorithm.
     */
    double max_sampling_time_step_scale = 5.0;
};

}

#endif
