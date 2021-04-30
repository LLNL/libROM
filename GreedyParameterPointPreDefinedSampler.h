/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: This class greedily selects parameter points
//              for the construction of a ROM database. The implemented
//              greedy algorithm is obtained from algorithm 2 of Choi et. al's
//              paper "Gradient-based constrained optimization using a database
//              of linear reduced-order models": https://arxiv.org/abs/1506.07849
//
//              The greedy algorithm workflow is as follows:
//              1. Construct the GreedyParameterPointPreDefinedSampler by giving it a
//                 domain of parameter points.
//              2. Request a parameter point to sample.
//              3. Request a parameter point to compute an error residual for.
//              4. Request the nearest ROM to the parameter point requiring a
//                 residual.
//              5. Give the computed residual to the GreedyParameterPointPreDefinedSampler.
//              6. Repeat steps 4 and 5 until the GreedyParameterPointPreDefinedSampler
//                 no longer requires more residuals to be computed.
//              7. Repeat steps 2 to 6 until the GreedyParameterPointPreDefinedSampler
//                 no longer requires more parameter points to be sampled.
//              8. The ROM database is now complete, meeting the error tolerance
//                 for all parameter points within the domain.

#ifndef included_GreedyParameterPointPreDefinedSampler_h
#define included_GreedyParameterPointPreDefinedSampler_h

#include "GreedyParameterPointSampler.h"

namespace CAROM {

/**
 * Class GreedyParameterPointPreDefinedSampler is a class defining the interface of a
 *       greedy algorithm that given a domain of parameter points, iteratively
 *       returns the next best parameter point to sample in order to create
 *       a ROM database efficiently.
 */
class GreedyParameterPointPreDefinedSampler : public GreedyParameterPointSampler
{
public:
    /**
     * @brief Constructor.
     *
     * @param[in] parameter_points A vector of CAROM::Vectors containing
                                   the different parameter points.
     * @param[in] check_local_rom Compute local ROM residual each iteration.
     * @param[in] tolerance A tolerance value for which to end the algorithm.
     * @param[in] saturation A saturation constant.
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
        double saturation,
        int subset_size,
        int convergence_subset_size,
        std::string output_log_path = "",
        std::string warm_start_file_name = "",
        bool use_centroid = true,
        int random_seed = 1,
        bool debug_algorithm = false);

    GreedyParameterPointPreDefinedSampler(
        std::vector<double> parameter_points,
        bool check_local_rom,
        double relative_error_tolerance,
        double saturation,
        int subset_size,
        int convergence_subset_size,
        std::string output_log_path = "",
        std::string warm_start_file_name = "",
        bool use_centroid = true,
        int random_seed = 1,
        bool debug_algorithm = false);

    GreedyParameterPointPreDefinedSampler(
        std::string base_file_name,
        std::string output_log_path = "");

    /**
     * @brief Destructor.
     */
    ~GreedyParameterPointPreDefinedSampler();

protected:

    void constructParameterPoints();

    void getNextParameterPointAfterConvergenceFailure();

    /**
     * @brief Use latin hypercube sampling instead of fixed uniform.
     */
    bool d_use_latin_hypercube;
};

}

#endif
