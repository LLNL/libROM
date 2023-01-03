/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: This class greedily selects parameter points
//              for the construction of a ROM database.

#include "GreedyCustomSampler.h"
#include <cmath>
#include <algorithm>
#include <limits.h>
#include <fstream>

namespace CAROM {

GreedyCustomSampler::GreedyCustomSampler(
    std::vector<double> parameter_points,
    bool check_local_rom,
    double relative_error_tolerance,
    double alpha,
    double max_clamp,
    int subset_size,
    int convergence_subset_size,
    std::string output_log_path,
    std::string warm_start_file_name,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm) :
    GreedySampler(
        parameter_points,
        check_local_rom,
        relative_error_tolerance,
        alpha,
        max_clamp,
        subset_size,
        convergence_subset_size,
        output_log_path,
        warm_start_file_name,
        use_centroid,
        random_seed,
        debug_algorithm
    )
{
    printSamplingType("pre-defined");
    constructParameterPoints();
    initializeParameterPoints();

    if (warm_start_file_name != "")
    {
        addDatabaseFromFile(warm_start_file_name);
    }
}

GreedyCustomSampler::GreedyCustomSampler(
    std::vector<Vector> parameter_points,
    bool check_local_rom,
    double relative_error_tolerance,
    double alpha,
    double max_clamp,
    int subset_size,
    int convergence_subset_size,
    std::string output_log_path,
    std::string warm_start_file_name,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm) :
    GreedySampler(
        parameter_points,
        check_local_rom,
        relative_error_tolerance,
        alpha,
        max_clamp,
        subset_size,
        convergence_subset_size,
        output_log_path,
        warm_start_file_name,
        use_centroid,
        random_seed,
        debug_algorithm
    )
{
    printSamplingType("pre-defined");
    constructParameterPoints();
    initializeParameterPoints();

    if (warm_start_file_name != "")
    {
        addDatabaseFromFile(warm_start_file_name);
    }
}

GreedyCustomSampler::GreedyCustomSampler(
    std::string base_file_name,
    std::string output_log_path) :
    GreedySampler(
        base_file_name,
        output_log_path
    ) {}

void
GreedyCustomSampler::constructParameterPoints()
{
    d_min_param_point = d_parameter_points[0];
    d_max_param_point = d_parameter_points[0];

    for (int i = 0; i < d_parameter_points.size(); i++)
    {
        for (int j = 0; j < d_parameter_points[i].dim(); j++)
        {
            d_min_param_point.item(j) = std::min(d_min_param_point.item(j),
                                                 d_parameter_points[i].item(j));
            d_max_param_point.item(j) = std::max(d_max_param_point.item(j),
                                                 d_parameter_points[i].item(j));
        }
    }

    checkParameterPointInput();
}

void
GreedyCustomSampler::getNextParameterPointAfterConvergenceFailure()
{
    d_next_point_to_sample = getNearestNonSampledPoint(
                                 d_convergence_points[d_counter]);
}

GreedyCustomSampler::~GreedyCustomSampler()
{
}

}
