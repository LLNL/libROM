/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete implementation of the Choi et. al's greedy algorithm
//              using pre-defined parameter points.

#ifndef included_GreedyParameterPointPreDefinedSampler_h
#define included_GreedyParameterPointPreDefinedSampler_h

#include "GreedyParameterPointSampler.h"

namespace CAROM {

/**
 * Class GreedyParameterPointRandomSampler implements a variant of
 *               Choi et. al's greedy algorithm using pre-defined
 *               parameter points.
 */
class GreedyParameterPointPreDefinedSampler : public GreedyParameterPointSampler
{
public:
    /**
     * @brief Constructor.
     *
     * @param[in] parameter_points A vector of CAROM::Vectors containing
                                   the different parameter points.
     * @param[in] check_local_rom Compute local ROM error indicator each iteration.
     * @param[in] relative_error_tolerance The relative error tolerance value
                                           for which to end the algorithm.
     * @param[in] alpha A alpha constant to increase greedy algorithm by each
                        iteration.
     * @param[in] max_clamp A scalar factor representing the maximum amount
                            the error indicator tolerance can change per iteration.
     * @param[in] subset_size The size of the random subset.
     * @param[in] convergence_subset_size The size of the convergence subset.
     * @param[in] use_centroid Whether to use the centroid heuristic when
                               determining the first parameter point to sample.
     * @param[in] random_seed A random seed.
     * @param[in] debug_algorithm Whether to turn off all randomness for
     *                            debugging purposes.
     */
    GreedyParameterPointPreDefinedSampler(
        std::vector<Vector> parameter_points,
        bool check_local_rom,
        double relative_error_tolerance,
        double alpha,
        double max_clamp,
        int subset_size,
        int convergence_subset_size,
        std::string output_log_path = "",
        std::string warm_start_file_name = "",
        bool use_centroid = true,
        int random_seed = 1,
        bool debug_algorithm = false);

    /**
     * @brief Constructor.
     *
     * @param[in] parameter_points A vector of doubles containing
                                   the different parameter points.
     * @param[in] check_local_rom Compute local ROM error indicator each iteration.
     * @param[in] relative_error_tolerance The relative error tolerance value
                                           for which to end the algorithm.
     * @param[in] alpha A alpha constant to increase greedy algorithm by each
                        iteration.
     * @param[in] max_clamp A scalar factor representing the maximum amount
                            the error indicator tolerance can change per iteration.
     * @param[in] subset_size The size of the random subset.
     * @param[in] convergence_subset_size The size of the convergence subset.
     * @param[in] use_centroid Whether to use the centroid heuristic when
                               determining the first parameter point to sample.
     * @param[in] random_seed A random seed.
     * @param[in] debug_algorithm Whether to turn off all randomness for
     *                            debugging purposes.
     */
    GreedyParameterPointPreDefinedSampler(
        std::vector<double> parameter_points,
        bool check_local_rom,
        double relative_error_tolerance,
        double alpha,
        double max_clamp,
        int subset_size,
        int convergence_subset_size,
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
    GreedyParameterPointPreDefinedSampler(
        std::string base_file_name,
        std::string output_log_path = "");

    /**
     * @brief Destructor.
     */
    ~GreedyParameterPointPreDefinedSampler();

protected:

    /**
     * @brief Construct the list of parameter point candidates to sample.
     */
    void constructParameterPoints();

    /**
     * @brief Get the next parameter point to sample after a convergence failure.
     */
    void getNextParameterPointAfterConvergenceFailure();
};

}

#endif
