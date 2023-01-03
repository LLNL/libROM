/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete implementation of the Choi et. al's greedy algorithm
//              using random sampling or latin-hypercube sampling.

#ifndef included_GreedyRandomSampler_h
#define included_GreedyRandomSampler_h

#include "GreedySampler.h"

namespace CAROM {

/**
 * Class GreedyRandomSampler implements a variant of
 *               Choi et. al's greedy algorithm using random sampling or
 *               latin-hypercube sampling.
 */
class GreedyRandomSampler : public GreedySampler
{
public:
    /**
     * @brief Constructor.
     *
     * @param[in] param_space_min A CAROM::Vector representing the minimum
                                  of the parameter space domain.
     * @param[in] param_space_max A CAROM::Vector representing the maximum
                                  of the parameter space domain.
     * @param[in] num_parameter_points The maximum number of parameter points
                                       to sample.
     * @param[in] check_local_rom Compute local ROM error indicator each iteration.
     * @param[in] relative_error_tolerance The relative error tolerance value
                                           for which to end the algorithm.
     * @param[in] alpha A alpha constant to increase greedy algorithm by each
                        iteration.
     * @param[in] max_clamp A scalar factor representing the maximum amount
                            the error indicator tolerance can change per iteration.
     * @param[in] subset_size The size of the random subset.
     * @param[in] convergence_subset_size The size of the convergence subset.
     * @param[in] use_latin_hypercube Whether to use latin-hypercube sampling
                                      instead of random sampling.
     * @param[in] output_log_path The path to the output log file. If not used,
                                  outputs to stdout.
     * @param[in] warm_start_file_name The path to the HDF5 file of a previous
                                       database to use as a warm start.
     * @param[in] use_centroid Whether to use the centroid heuristic when
                               determining the first parameter point to sample.
     * @param[in] random_seed A random seed.
     * @param[in] debug_algorithm Whether to turn off all randomness for
     *                            debugging purposes.
     */
    GreedyRandomSampler(
        Vector param_space_min,
        Vector param_space_max,
        int num_parameter_points,
        bool check_local_rom,
        double relative_error_tolerance,
        double alpha,
        double max_clamp,
        int subset_size,
        int convergence_subset_size,
        bool use_latin_hypercube,
        std::string output_log_path = "",
        std::string warm_start_file_name = "",
        bool use_centroid = true,
        int random_seed = 1,
        bool debug_algorithm = false);

    /**
     * @brief Constructor.
     *
     * @param[in] param_space_min A double representing the minimum
                                  of the parameter space domain.
     * @param[in] param_space_max A double representing the maximum
                                  of the parameter space domain.
     * @param[in] num_parameter_points The maximum number of parameter points
                                       to sample.
     * @param[in] check_local_rom Compute local ROM error indicator each iteration.
     * @param[in] relative_error_tolerance The relative error tolerance value
                                           for which to end the algorithm.
     * @param[in] alpha A alpha constant to increase greedy algorithm by each
                        iteration.
     * @param[in] max_clamp A scalar factor representing the maximum amount
                            the error indicator tolerance can change per iteration.
     * @param[in] subset_size The size of the random subset.
     * @param[in] convergence_subset_size The size of the convergence subset.
     * @param[in] use_latin_hypercube Whether to use latin-hypercube sampling
                                      instead of random sampling.
     * @param[in] output_log_path The path to the output log file. If not used,
                                  outputs to stdout.
     * @param[in] warm_start_file_name The path to the HDF5 file of a previous
                                       database to use as a warm start.
     * @param[in] use_centroid Whether to use the centroid heuristic when
                               determining the first parameter point to sample.
     * @param[in] random_seed A random seed.
     * @param[in] debug_algorithm Whether to turn off all randomness for
     *                            debugging purposes.
     */
    GreedyRandomSampler(
        double param_space_min,
        double param_space_max,
        int num_parameter_points,
        bool check_local_rom,
        double relative_error_tolerance,
        double alpha,
        double max_clamp,
        int subset_size,
        int convergence_subset_size,
        bool use_latin_hypercube,
        std::string output_log_path = "",
        std::string warm_start_file_name = "",
        bool use_centroid = true,
        int random_seed = 1,
        bool debug_algorithm = false);

    /**
     * @brief Constructor.
     *
     * @param[in] base_file_name The base part of the file of the
                                 database to load when restarting from a save.
     * @param[in] output_log_path The path to the output log file. If not used,
                                  outputs to stdout.
     */
    GreedyRandomSampler(
        std::string base_file_name,
        std::string output_log_path = "");

    /**
     * @brief Save the object state to a file.
     *
     * @param[in] base_file_name The base part of the file to save the
     *                           database to.
     */
    void
    save(std::string base_file_name);

    /**
     * @brief Destructor.
     */
    ~GreedyRandomSampler();

protected:

    /**
     * @brief Load the object state from a file.
     *
     * @param[in] base_file_name The base part of the file to load the
     *                           database from.
     */
    void load(std::string base_file_name);

    /**
     * @brief Construct the list of parameter point candidates to sample.
     */
    void constructParameterPoints();

    /**
     * @brief Get the next parameter point to sample after a convergence failure.
     */
    void getNextParameterPointAfterConvergenceFailure();

    /**
     * @brief Use latin hypercube sampling instead of fixed uniform.
     */
    bool d_use_latin_hypercube;
};

}

#endif
