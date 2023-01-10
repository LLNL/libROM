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

#include "GreedyRandomSampler.h"
#include "utils/HDFDatabase.h"
#include "mpi.h"
#include <cmath>
#include <algorithm>
#include <limits.h>
#include <fstream>

namespace CAROM {

GreedyRandomSampler::GreedyRandomSampler(
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
    std::string output_log_path,
    std::string warm_start_file_name,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm) :
    GreedySampler(
        param_space_min,
        param_space_max,
        num_parameter_points,
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
    d_use_latin_hypercube = use_latin_hypercube;

    if (d_use_latin_hypercube)
    {
        printSamplingType("latin-hypercube");
    }
    else
    {
        printSamplingType("random");
    }

    constructParameterPoints();
    initializeParameterPoints();

    if (warm_start_file_name != "")
    {
        addDatabaseFromFile(warm_start_file_name);
    }
}

GreedyRandomSampler::GreedyRandomSampler(
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
    std::string output_log_path,
    std::string warm_start_file_name,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm) :
    GreedySampler(
        param_space_min,
        param_space_max,
        num_parameter_points,
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
    d_use_latin_hypercube = use_latin_hypercube;

    if (d_use_latin_hypercube)
    {
        printSamplingType("latin-hypercube");
    }
    else
    {
        printSamplingType("random");
    }

    constructParameterPoints();
    initializeParameterPoints();

    if (warm_start_file_name != "")
    {
        addDatabaseFromFile(warm_start_file_name);
    }
}

GreedyRandomSampler::GreedyRandomSampler(
    std::string base_file_name,
    std::string output_log_path) :
    GreedySampler(
        base_file_name,
        output_log_path
    ) {}

void
GreedyRandomSampler::constructParameterPoints()
{
    checkParameterPointInput();

    Vector vec(d_min_param_point.dim(), false);

    if (d_use_latin_hypercube)
    {
        std::vector<double> frequencies;
        std::vector<std::vector<double>> point_coordinates;
        std::vector<int> point_indices;
        for (int i = 0; i < d_min_param_point.dim(); i++)
        {
            frequencies.push_back(std::abs(d_max_param_point.item(i) -
                                           d_min_param_point.item(i)) / d_num_parameter_points);
        }

        for (int i = 0; i < d_num_parameter_points; i++) {
            point_coordinates.push_back(std::vector<double>());
            for (int j = 0; j < d_min_param_point.dim(); j++)
            {
                double low = d_min_param_point.item(j) + i * frequencies[j];
                double high = d_min_param_point.item(j) + ((i + 1) * frequencies[j]);
                std::uniform_real_distribution<double> unif(low, high);
                point_coordinates[i].push_back(unif(rng));
            }
        }

        for (int i = 0; i < std::pow(d_num_parameter_points, d_min_param_point.dim());
                i++)
        {
            point_indices.push_back(i);
        }

        std::shuffle(point_indices.begin(), point_indices.end(), rng);

        for (int i = 0; i < d_num_parameter_points; i++)
        {
            int point_index = point_indices[i];
            for (int j = 0; j < d_min_param_point.dim(); j++)
            {
                vec.item(j) = point_coordinates[point_index % d_num_parameter_points][j];
                point_index /= d_num_parameter_points;
            }
            d_parameter_points.push_back(vec);
        }
    }
    else
    {
        d_parameter_points = generateRandomPoints(d_num_parameter_points);
    }
}

void
GreedyRandomSampler::load(std::string base_file_name)
{
    GreedySampler::load(base_file_name);

    char tmp[100];
    std::string full_file_name = base_file_name;
    HDFDatabase database;
    database.open(full_file_name, "r");

    if (!d_procedure_completed)
    {
        int bool_int_temp;
        sprintf(tmp, "use_latin_hypercube");
        database.getInteger(tmp, bool_int_temp);
        d_use_latin_hypercube = (bool) bool_int_temp;
    }

    database.close();
}

void
GreedyRandomSampler::save(std::string base_file_name)
{
    GreedySampler::save(base_file_name);

    if (d_rank == 0)
    {
        char tmp[100];
        std::string full_file_name = base_file_name;
        HDFDatabase database;
        database.open(full_file_name, "wr");

        if (!d_procedure_completed)
        {
            sprintf(tmp, "use_latin_hypercube");
            database.putInteger(tmp, d_use_latin_hypercube);
        }

        database.close();
    }
}

void
GreedyRandomSampler::getNextParameterPointAfterConvergenceFailure()
{
    d_parameter_points.push_back(d_convergence_points[d_counter]);
    d_parameter_point_errors.push_back(d_max_error);
    d_parameter_point_local_rom.push_back(getNearestROMIndexToParameterPoint(
            d_next_point_to_sample, true));
    d_parameter_point_random_indices.push_back(
        d_parameter_point_random_indices.size());
    d_next_point_to_sample = d_parameter_points.size() - 1;
}

GreedyRandomSampler::~GreedyRandomSampler()
{
}

}
