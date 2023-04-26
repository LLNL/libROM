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

#include "GreedySampler.h"
#include "utils/HDFDatabase.h"
#include "mpi.h"
#include <cmath>
#include <algorithm>
#include <limits.h>
#include <fstream>

namespace CAROM {

struct GreedyErrorIndicatorPoint
createGreedyErrorIndicatorPoint(Vector* point, Vector* localROM)
{
    struct GreedyErrorIndicatorPoint result;
    result.point = std::shared_ptr<Vector>(point);
    result.localROM = std::shared_ptr<Vector>(localROM);
    return result;
}

struct GreedyErrorIndicatorPoint
createGreedyErrorIndicatorPoint(Vector* point,
                                std::shared_ptr<Vector>& localROM)
{
    struct GreedyErrorIndicatorPoint result;
    result.point = std::shared_ptr<Vector>(point);
    result.localROM = localROM;
    return result;
}

Vector
getNearestPoint(std::vector<Vector> param_points, Vector point)
{

    int closest_point_index = getNearestPointIndex(param_points, point);
    return param_points[closest_point_index];
}

int
getNearestPointIndex(std::vector<Vector> param_points, Vector point)
{

    double closest_dist_to_points = INT_MAX;
    int closest_point_index = -1;

    for (int i = 0; i < param_points.size(); i++)
    {
        Vector diff;
        point.minus(param_points[i], diff);
        double dist = diff.norm();
        if (dist < closest_dist_to_points)
        {
            closest_dist_to_points = dist;
            closest_point_index = i;
        }
    }

    return closest_point_index;
}

double
getNearestPoint(std::vector<double> param_points, double point)
{

    int closest_point_index = getNearestPointIndex(param_points, point);
    return param_points[closest_point_index];
}

int
getNearestPointIndex(std::vector<double> param_points, double point)
{

    double closest_dist_to_points = INT_MAX;
    int closest_point_index = -1;

    for (int i = 0; i < param_points.size(); i++)
    {
        double dist = std::abs(point - param_points[i]);
        if (dist < closest_dist_to_points)
        {
            closest_dist_to_points = dist;
            closest_point_index = i;
        }
    }

    return closest_point_index;
}

GreedySampler::GreedySampler(
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
    bool debug_algorithm)
{

    d_num_parameter_points = parameter_points.size();
    CAROM_VERIFY(d_num_parameter_points >= 1);
    d_parameter_points = parameter_points;

    constructObject(check_local_rom, relative_error_tolerance, alpha, max_clamp,
                    subset_size, convergence_subset_size, output_log_path, use_centroid,
                    random_seed, debug_algorithm);
}

GreedySampler::GreedySampler(
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
    bool debug_algorithm)
{

    d_num_parameter_points = parameter_points.size();

    std::vector<Vector> parameter_points_vec;
    for (int i = 0; i < parameter_points.size(); i++)
    {
        Vector vec(1, false);
        vec.item(0) = parameter_points[i];
        parameter_points_vec.push_back(vec);
    }
    CAROM_VERIFY(d_num_parameter_points >= 1);
    d_parameter_points = parameter_points_vec;

    constructObject(check_local_rom, relative_error_tolerance, alpha, max_clamp,
                    subset_size, convergence_subset_size, output_log_path, use_centroid,
                    random_seed, debug_algorithm);
}

GreedySampler::GreedySampler(
    Vector param_space_min,
    Vector param_space_max,
    int num_parameter_points,
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
    bool debug_algorithm)
{

    d_min_param_point = param_space_min;
    d_max_param_point = param_space_max;
    d_num_parameter_points = num_parameter_points;

    constructObject(check_local_rom, relative_error_tolerance, alpha, max_clamp,
                    subset_size,
                    convergence_subset_size, output_log_path, use_centroid, random_seed,
                    debug_algorithm);
}

GreedySampler::GreedySampler(
    double param_space_min,
    double param_space_max,
    int num_parameter_points,
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
    bool debug_algorithm)
{

    Vector param_space_min_vec(1, false);
    param_space_min_vec.item(0) = param_space_min;
    Vector param_space_max_vec(1, false);
    param_space_max_vec.item(0) = param_space_max;

    d_min_param_point = param_space_min_vec;
    d_max_param_point = param_space_max_vec;
    d_num_parameter_points = num_parameter_points;

    constructObject(check_local_rom, relative_error_tolerance, alpha, max_clamp,
                    subset_size,
                    convergence_subset_size, output_log_path, use_centroid, random_seed,
                    debug_algorithm);
}

GreedySampler::GreedySampler(
    std::string base_file_name,
    std::string output_log_path)
{
    CAROM_ASSERT(!base_file_name.empty());

    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init) {
        MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    }
    else {
        d_rank = 0;
    }

    d_output_log_path = output_log_path;

    load(base_file_name);
}

void
GreedySampler::addDatabaseFromFile(
    std::string const& warm_start_file_name)
{
    char tmp[100];
    std::string full_file_name = warm_start_file_name;
    HDFDatabase database;
    database.open(full_file_name, "r");

    sprintf(tmp, "num_parameter_sampled_indices");
    int num_parameter_sampled_indices;
    database.getInteger(tmp, num_parameter_sampled_indices);
    if (num_parameter_sampled_indices > 0)
    {
        int temp_parameter_sampled_indices[num_parameter_sampled_indices];
        sprintf(tmp, "parameter_sampled_indices");
        database.getIntegerArray(tmp, &temp_parameter_sampled_indices[0],
                                 num_parameter_sampled_indices);
        for (int i = 0; i < num_parameter_sampled_indices; i++)
        {
            std::string vec_path = warm_start_file_name + "_" + std::to_string(i);
            Vector point;
            point.local_read(vec_path, 0);
            for (int j = 0; j < d_parameter_points.size(); j++)
            {
                Vector diff;
                point.minus(d_parameter_points[j], diff);
                double dist = diff.norm();
                if (dist > 1e-12)
                {
                    d_parameter_points.push_back(point);
                    d_parameter_sampled_indices.insert(d_parameter_points.size() - 1);
                }
            }
        }
    }
}

void
GreedySampler::load(
    std::string base_file_name)
{
    char tmp[100];
    std::string full_file_name = base_file_name;
    HDFDatabase database;
    database.open(full_file_name, "r");

    sprintf(tmp, "num_parameter_points");
    int num_parameter_points;
    database.getInteger(tmp, num_parameter_points);
    for (int i = 0; i < num_parameter_points; i++)
    {
        std::string vec_path = base_file_name + "_" + std::to_string(i);
        Vector point;
        point.local_read(vec_path, 0);
        d_parameter_points.push_back(point);
    }

    std::string vec_path = base_file_name + "_min_point";
    d_min_param_point.local_read(vec_path, 0);

    vec_path = base_file_name + "_max_point";
    d_max_param_point.local_read(vec_path, 0);

    sprintf(tmp, "num_parameter_sampled_indices");
    int num_parameter_sampled_indices;
    database.getInteger(tmp, num_parameter_sampled_indices);
    if (num_parameter_sampled_indices > 0)
    {
        int temp_parameter_sampled_indices[num_parameter_sampled_indices];
        sprintf(tmp, "parameter_sampled_indices");
        database.getIntegerArray(tmp, &temp_parameter_sampled_indices[0],
                                 num_parameter_sampled_indices);
        for (int i = 0; i < num_parameter_sampled_indices; i++)
        {
            d_parameter_sampled_indices.insert(temp_parameter_sampled_indices[i]);
        }
    }

    int bool_int_temp;
    sprintf(tmp, "procedure_completed");
    database.getInteger(tmp, bool_int_temp);
    d_procedure_completed = (bool) bool_int_temp;

    if (!d_procedure_completed)
    {
        sprintf(tmp, "max_error");
        database.getDouble(tmp, d_max_error);
        sprintf(tmp, "curr_relative_error");
        database.getDouble(tmp, d_curr_relative_error);
        sprintf(tmp, "error_indicator_tol");
        database.getDouble(tmp, d_error_indicator_tol);
        sprintf(tmp, "relative_error_tol");
        database.getDouble(tmp, d_relative_error_tol);
        sprintf(tmp, "alpha");
        database.getDouble(tmp, d_alpha);
        sprintf(tmp, "max_clamp");
        database.getDouble(tmp, d_max_clamp);
        sprintf(tmp, "max_num_parameter_points");
        database.getInteger(tmp, d_num_parameter_points);
        sprintf(tmp, "subset_size");
        database.getInteger(tmp, d_subset_size);
        sprintf(tmp, "convergence_subset_size");
        database.getInteger(tmp, d_convergence_subset_size);
        sprintf(tmp, "next_point_to_sample");
        database.getInteger(tmp, d_next_point_to_sample);
        sprintf(tmp, "next_point_requiring_error_indicator");
        database.getInteger(tmp, d_next_point_requiring_error_indicator);
        sprintf(tmp, "check_local_rom");
        database.getInteger(tmp, bool_int_temp);
        d_check_local_rom = (bool) bool_int_temp;
        sprintf(tmp, "use_centroid");
        database.getInteger(tmp, bool_int_temp);
        d_use_centroid = (bool) bool_int_temp;
        sprintf(tmp, "iteration_started");
        database.getInteger(tmp, bool_int_temp);
        d_iteration_started = (bool) bool_int_temp;
        sprintf(tmp, "convergence_started");
        database.getInteger(tmp, bool_int_temp);
        d_convergence_started = (bool) bool_int_temp;
        sprintf(tmp, "next_parameter_point_computed");
        database.getInteger(tmp, bool_int_temp);
        d_next_parameter_point_computed = (bool) bool_int_temp;
        sprintf(tmp, "point_requiring_error_indicator_computed");
        database.getInteger(tmp, bool_int_temp);
        d_point_requiring_error_indicator_computed = (bool) bool_int_temp;
        sprintf(tmp, "subset_created");
        database.getInteger(tmp, bool_int_temp);
        d_subset_created = (bool) bool_int_temp;
        sprintf(tmp, "debug_algorithm");
        database.getInteger(tmp, bool_int_temp);
        d_debug_algorithm = (bool) bool_int_temp;
        sprintf(tmp, "counter");
        database.getInteger(tmp, d_counter);
        sprintf(tmp, "subset_counter");
        database.getInteger(tmp, d_subset_counter);
        sprintf(tmp, "random_seed");
        database.getInteger(tmp, d_random_seed);

        sprintf(tmp, "parameter_point_random_indices");
        d_parameter_point_random_indices.resize(d_parameter_points.size());
        database.getIntegerArray(tmp, &d_parameter_point_random_indices[0],
                                 d_parameter_points.size());

        sprintf(tmp, "parameter_point_errors");
        d_parameter_point_errors.resize(d_parameter_points.size());
        database.getDoubleArray(tmp, &d_parameter_point_errors[0],
                                d_parameter_points.size());

        sprintf(tmp, "parameter_point_local_rom");
        d_parameter_point_local_rom.resize(d_parameter_points.size());
        database.getIntegerArray(tmp, &d_parameter_point_local_rom[0],
                                 d_parameter_points.size());

        for (int i = 0; i < d_convergence_subset_size; i++)
        {
            std::string vec_path = base_file_name + "_conv_" + std::to_string(i);
            Vector point;
            point.local_read(vec_path, 0);
            d_convergence_points.push_back(point);
        }
    }
    database.close();

    rng.seed(d_random_seed + d_parameter_sampled_indices.size());
}

void
GreedySampler::checkParameterPointInput()
{
    CAROM_VERIFY(d_min_param_point.dim() == d_max_param_point.dim());
    CAROM_VERIFY(d_num_parameter_points >= 1);

    bool isGreater = false;
    for (int i = 0; i < d_min_param_point.dim(); i++)
    {
        if (d_num_parameter_points == 1)
        {
            if (d_max_param_point.item(i) >= d_min_param_point.item(i))
            {
                isGreater = true;
                break;
            }
        }
        else
        {
            if (d_max_param_point.item(i) > d_min_param_point.item(i))
            {
                isGreater = true;
                break;
            }
        }
    }
    CAROM_VERIFY(isGreater);

    if (d_rank == 0)
    {
        std::string str;
        str += "Total number of sample points: " + std::to_string(
                   d_num_parameter_points) + "\n";
        str += "Parameter space minimum: [ ";
        for (int i = 0 ; i < d_min_param_point.dim(); i++)
        {
            str += std::to_string(d_min_param_point.item(i)) + " ";
        }
        str += "]\n";
        str += "Parameter space maximum: [ ";
        for (int i = 0 ; i < d_max_param_point.dim(); i++)
        {
            str += std::to_string(d_max_param_point.item(i)) + " ";
        }
        str += "]\n";
        agnosticPrint(str);
    }
}

void
GreedySampler::constructObject(
    bool check_local_rom,
    double relative_error_tolerance,
    double alpha,
    double max_clamp,
    int subset_size,
    int convergence_subset_size,
    std::string output_log_path,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm)
{
    CAROM_VERIFY(relative_error_tolerance > 0.0);
    CAROM_VERIFY(alpha >= 1.0);
    CAROM_VERIFY(max_clamp >= 1.0);
    CAROM_VERIFY(subset_size > 0);
    CAROM_VERIFY(convergence_subset_size > 0);
    CAROM_VERIFY(subset_size < convergence_subset_size);
    CAROM_VERIFY(random_seed > 0);

    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init) {
        MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    }
    else {
        d_rank = 0;
    }

    d_check_local_rom = check_local_rom;
    d_error_indicator_tol = 0.0;
    d_relative_error_tol = relative_error_tolerance;
    d_alpha = alpha;
    d_max_clamp = max_clamp;
    d_subset_size = subset_size;
    d_convergence_subset_size = convergence_subset_size;
    d_output_log_path = output_log_path;
    d_use_centroid = use_centroid;
    d_max_error = 0;
    d_next_point_to_sample = -1;
    d_next_point_requiring_error_indicator = -1;
    d_next_parameter_point_computed = false;
    d_point_requiring_error_indicator_computed = false;
    d_iteration_started = false;
    d_convergence_started = false;
    d_subset_created = false;
    d_procedure_completed = false;
    d_subset_counter = 0;
    d_counter = -1;
    d_debug_algorithm = debug_algorithm;
    d_random_seed = random_seed;

    rng.seed(d_random_seed);

    if (d_rank == 0)
    {
        std::string str;
        str += "Relative error tolerance: " + std::to_string(d_relative_error_tol) +
               "\n";
        str += "Alpha constant: " + std::to_string(d_alpha) + "\n";
        str += "Max clamp constant: " + std::to_string(d_max_clamp) + "\n";
        str += "Iteration subset size: " + std::to_string(d_subset_size) + "\n";
        str += "Convergence subset size: " + std::to_string(d_convergence_subset_size) +
               "\n";
        agnosticPrint(str);
    }
}

void
GreedySampler::initializeParameterPoints()
{
    CAROM_VERIFY(d_parameter_points.size() > 0);
    CAROM_VERIFY(d_subset_size <= d_parameter_points.size());
    CAROM_VERIFY(d_convergence_subset_size <= d_parameter_points.size());

    for (int i = 0; i < d_parameter_points.size() - 1; i++) {
        CAROM_VERIFY(d_parameter_points[i].dim() == d_parameter_points[i + 1].dim());
    }

    for (int i = 0 ; i < d_parameter_points.size(); i++) {
        d_parameter_point_errors.push_back(INT_MAX);
        d_parameter_point_local_rom.push_back(-1);
        d_parameter_point_random_indices.push_back(i);
    }

    generateConvergenceSubset();
}

std::shared_ptr<Vector>
GreedySampler::getNextParameterPoint()
{
    if (isComplete())
    {
        return std::shared_ptr<Vector>(nullptr);
    }
    if (d_iteration_started)
    {
        return std::shared_ptr<Vector>(nullptr);
    }
    if (d_next_parameter_point_computed)
    {
        Vector* result = new Vector(d_parameter_points[d_next_point_to_sample]);
        return std::shared_ptr<Vector>(result);
    }

    if (d_parameter_sampled_indices.size() == 0)
    {
        d_next_point_to_sample = getCenterPoint(d_parameter_points, d_use_centroid);
    }

    d_max_error = 0;

    int curr_point_to_sample = d_next_point_to_sample;

    d_convergence_started = false;
    d_subset_created = false;
    d_subset_counter = 0;
    d_counter = -1;

    auto search = d_parameter_sampled_indices.find(curr_point_to_sample);
    CAROM_VERIFY(search == d_parameter_sampled_indices.end());

    d_parameter_sampled_indices.insert(curr_point_to_sample);

    if (d_rank == 0)
    {
        std::string str;
        str += "\nPoint sampled at [ ";
        for (int i = 0 ; i < d_parameter_points[curr_point_to_sample].dim(); i++)
        {
            str += std::to_string(d_parameter_points[curr_point_to_sample].item(i)) + " ";
        }
        str += "]\n";
        agnosticPrint(str);
    }

    d_next_parameter_point_computed = true;

    Vector* result = new Vector(d_parameter_points[curr_point_to_sample]);
    return std::shared_ptr<Vector>(result);
}

struct GreedyErrorIndicatorPoint
GreedySampler::getNextPointRequiringRelativeError()
{
    if (isComplete())
    {
        return createGreedyErrorIndicatorPoint(nullptr, nullptr);
    }
    if (!d_next_parameter_point_computed)
    {
        return createGreedyErrorIndicatorPoint(nullptr, nullptr);
    }

    Vector* result1 = new Vector(d_parameter_points[d_next_point_to_sample]);
    Vector* result2 = NULL;

    if (d_parameter_sampled_indices.size() == 1)
    {
        result2 = new Vector(d_parameter_points[getNearestROMIndexToParameterPoint(
                d_next_point_to_sample, false)]);
    }
    else
    {
        result2 = new Vector(d_parameter_points[getNearestROMIndexToParameterPoint(
                d_next_point_to_sample, true)]);
    }

    return createGreedyErrorIndicatorPoint(result1, result2);
}

struct GreedyErrorIndicatorPoint
GreedySampler::getNextPointRequiringErrorIndicator()
{
    if (isComplete())
    {
        return createGreedyErrorIndicatorPoint(nullptr, nullptr);
    }
    if (!d_iteration_started)
    {
        return createGreedyErrorIndicatorPoint(nullptr, nullptr);
    }

    if (d_convergence_started)
    {
        return getNextConvergencePointRequiringErrorIndicator();
    }
    else
    {
        return getNextSubsetPointRequiringErrorIndicator();
    }
}

struct GreedyErrorIndicatorPoint
GreedySampler::getNextSubsetPointRequiringErrorIndicator()
{
    if (d_point_requiring_error_indicator_computed)
    {
        Vector* result1 = new Vector(
            d_parameter_points[d_next_point_requiring_error_indicator]);
        Vector* result2 = new Vector(
            d_parameter_points[getNearestROMIndexToParameterPoint(
                                   d_next_point_requiring_error_indicator, false)]);
        return createGreedyErrorIndicatorPoint(result1, result2);
    }
    if (d_subset_counter == d_subset_size)
    {
        return createGreedyErrorIndicatorPoint(nullptr, nullptr);
    }

    if(!d_subset_created)
    {
        // generate random shuffle
        if (!d_debug_algorithm)
        {
            std::shuffle(d_parameter_point_random_indices.begin(),
                         d_parameter_point_random_indices.end(), rng);
        }
        d_subset_created = true;
    }

    d_next_point_requiring_error_indicator = -1;

    while (d_counter < (int) d_parameter_points.size() - 1)
    {
        d_counter++;
        if (d_subset_counter == d_subset_size)
        {
            break;
        }
        auto search = d_parameter_sampled_indices.find(
                          d_parameter_point_random_indices[d_counter]);
        if (search == d_parameter_sampled_indices.end())
        {
            d_subset_counter++;
            double curr_error =
                d_parameter_point_errors[d_parameter_point_random_indices[d_counter]];
            if (curr_error > d_max_error)
            {
                // if we have already computed this error indicator at the same local rom, the error indicator will not improve
                // no need to calculate the error indicator again
                if (d_parameter_point_local_rom[d_parameter_point_random_indices[d_counter]] ==
                        getNearestROMIndexToParameterPoint(d_parameter_point_random_indices[d_counter],
                                false))
                {
                    d_max_error = curr_error;
                    d_next_point_to_sample = d_parameter_point_random_indices[d_counter];
                    if (d_rank == 0)
                    {
                        std::string str;
                        str += "Error indicator at [ ";
                        for (int i = 0 ;
                                i < d_parameter_points[d_parameter_point_random_indices[d_counter]].dim(); i++)
                        {
                            str += std::to_string(
                                       d_parameter_points[d_parameter_point_random_indices[d_counter]].item(i)) + " ";
                        }
                        str += "] skipped.\n";
                        str += "Error indicator " + std::to_string(curr_error) +
                               " already computed at the same local ROM.\n";
                        agnosticPrint(str);
                    }
                }
                else
                {
                    d_next_point_requiring_error_indicator =
                        d_parameter_point_random_indices[d_counter];
                    d_point_requiring_error_indicator_computed = true;
                    Vector* result1 = new Vector(
                        d_parameter_points[d_next_point_requiring_error_indicator]);
                    Vector* result2 = new Vector(
                        d_parameter_points[getNearestROMIndexToParameterPoint(
                                               d_next_point_requiring_error_indicator, false)]);
                    return createGreedyErrorIndicatorPoint(result1, result2);
                }
            }
            else
            {
                if (d_rank == 0)
                {
                    std::string str;
                    str += "Error indicator at [ ";
                    for (int i = 0 ;
                            i < d_parameter_points[d_parameter_point_random_indices[d_counter]].dim(); i++)
                    {
                        str += std::to_string(
                                   d_parameter_points[d_parameter_point_random_indices[d_counter]].item(i)) + " ";
                    }
                    str += "] skipped.\n";
                    str += "Error indicator " + std::to_string(curr_error) +
                           " is less than current max error " + std::to_string(d_max_error) + "\n";
                    agnosticPrint(str);
                }
            }
        }
    }
    if (d_next_point_requiring_error_indicator == -1)
    {
        if (d_rank == 0)
        {
            std::string str;
            str += "Ran out of points to calculate error indicator for in this iteration.\n";
            agnosticPrint(str);
        }
        if (d_max_error < d_error_indicator_tol)
        {
            startConvergence();
        }
        else
        {
            printErrorIndicatorToleranceNotMet();
            d_iteration_started = false;
        }
    }

    return createGreedyErrorIndicatorPoint(nullptr, nullptr);
}

struct GreedyErrorIndicatorPoint
GreedySampler::getNextConvergencePointRequiringErrorIndicator()
{
    if (d_point_requiring_error_indicator_computed)
    {
        Vector* result1 = new Vector(
            d_convergence_points[d_next_point_requiring_error_indicator]);
        std::shared_ptr<Vector> result2 = getNearestROM(
                                              d_convergence_points[d_next_point_requiring_error_indicator]);
        return createGreedyErrorIndicatorPoint(result1, result2);
    }
    if (d_counter == d_convergence_subset_size)
    {
        return createGreedyErrorIndicatorPoint(nullptr, nullptr);
    }

    d_next_point_requiring_error_indicator = -1;

    //get next point requiring error indicator
    if (d_counter < (int) d_convergence_points.size())
    {
        d_next_point_requiring_error_indicator = d_counter;
        d_point_requiring_error_indicator_computed = true;
        Vector* result1 = new Vector(
            d_convergence_points[d_next_point_requiring_error_indicator]);
        std::shared_ptr<Vector> result2 = getNearestROM(
                                              d_convergence_points[d_next_point_requiring_error_indicator]);
        return createGreedyErrorIndicatorPoint(result1, result2);
    }

    if (d_next_point_requiring_error_indicator == -1)
    {
        printConvergenceAchieved();
        d_procedure_completed = true;
    }

    return createGreedyErrorIndicatorPoint(nullptr, nullptr);
}

void
GreedySampler::printSamplingType(std::string sampling_type)
{
    if (d_rank == 0)
    {
        std::string str;
        str += "Sampling type: " + sampling_type + "\n";
        agnosticPrint(str);
    }
}

void
GreedySampler::printConvergenceAchieved()
{
    if (d_rank == 0)
    {
        std::string str;
        str += "Convergence achieved.\n";
        agnosticPrint(str);

        str = "\nSampled Parameter Points\n";
        std::vector<std::pair<double, int>> first_dim_of_sampled_points;
        for (auto itr = d_parameter_sampled_indices.begin();
                itr != d_parameter_sampled_indices.end(); ++itr) {
            first_dim_of_sampled_points.push_back(std::make_pair(
                    d_parameter_points[*itr].item(0), *itr));
        }
        sort(first_dim_of_sampled_points.begin(), first_dim_of_sampled_points.end());
        for (int i = 0; i < first_dim_of_sampled_points.size(); i++)
        {
            str += "[ ";
            for (int j = 0 ;
                    j < d_parameter_points[first_dim_of_sampled_points[i].second].dim(); j++)
            {
                str += std::to_string(
                           d_parameter_points[first_dim_of_sampled_points[i].second].item(j)) + " ";
            }
            str += "]\n";
        }
        agnosticPrint(str);
    }
}

void
GreedySampler::setPointRelativeError(double error)
{
    CAROM_VERIFY(error >= 0);
    CAROM_VERIFY(d_next_parameter_point_computed);

    if (!std::isfinite(error))
    {
        error = INT_MAX;
    }

    if (d_rank == 0)
    {
        if (d_output_log_path == "")
        {
            std::string str;
            str += "Relative error computed at [ ";
            for (int i = 0 ; i < d_parameter_points[d_next_point_to_sample].dim(); i++)
            {
                str += std::to_string(d_parameter_points[d_next_point_to_sample].item(i)) + " ";
            }
            str += "]\n";
            str += "Relative error: " + std::to_string(error) + "\n";
            agnosticPrint(str);
        }
    }

    d_curr_relative_error = error;
    d_next_parameter_point_computed = false;
    d_iteration_started = true;

    double old_error_indicator_tol = d_error_indicator_tol;

    if (d_parameter_sampled_indices.size() > 1)
    {
        double max1 = d_alpha * d_error_indicator_tol;
        double min1 = d_max_clamp * d_parameter_point_errors[d_next_point_to_sample];
        double min2 = d_relative_error_tol *
                      d_parameter_point_errors[d_next_point_to_sample] / d_curr_relative_error;

        if (d_curr_relative_error <= d_relative_error_tol)
        {
            if (d_rank == 0)
            {
                std::string str;
                str += "The relative error was smaller than the relative error tolerance. The error indicator tolerance must increase.\n";
                str += "Alpha constant * relative error tolerance: " + std::to_string(
                           max1) + "\n";
                str += "The minimum value the error indicator tolerance can take is: " +
                       std::to_string(max1) + "\n";
                str += "Max clamp constant * current max error indicator: " + std::to_string(
                           min1) + "\n";
                str += "Relative error tolerance * current max error indicator / current relative error: "
                       + std::to_string(min2) + "\n";
                str += "The maximum value the error indicator tolerance can take is the minimum of the previous two values: "
                       + std::to_string(std::min(min1, min2)) + "\n";
                agnosticPrint(str);
            }
            d_error_indicator_tol = std::max(max1, std::min(min1, min2));
        }
        else
        {
            double min1 = d_relative_error_tol *
                          d_parameter_point_errors[d_next_point_to_sample] / d_curr_relative_error;

            if (d_rank == 0)
            {
                std::string str;
                str += "The relative error was larger than the relative error tolerance. The error indicator tolerance must decrease.\n";
                str += "Current error indicator tolerance: " + std::to_string(
                           d_error_indicator_tol) + "\n";
                str += "Relative error tolerance * current max error indicator / current relative error: "
                       + std::to_string(min1) + "\n";
                str += "The minimum value the error indicator tolerance can take is: " +
                       std::to_string(std::min(d_error_indicator_tol, min1)) + "\n";
                agnosticPrint(str);
            }
            d_error_indicator_tol = std::min(d_error_indicator_tol, min1);
        }
        if (d_rank == 0)
        {
            if (old_error_indicator_tol != d_error_indicator_tol)
            {
                std::string str;
                str += "Error indicator tolerance was adaptively changed from " +
                       std::to_string(old_error_indicator_tol) + " to " + std::to_string(
                           d_error_indicator_tol) + "\n";
                agnosticPrint(str);
            }
        }
    }

    d_parameter_point_errors[d_next_point_to_sample] = 0;
    d_parameter_point_local_rom[d_next_point_to_sample] = d_next_point_to_sample;

    if (d_parameter_sampled_indices.size() > 1
            && d_curr_relative_error <= d_relative_error_tol)
    {
        startConvergence();
    }
    else
    {
        if (d_check_local_rom || d_parameter_sampled_indices.size() == 1)
        {
            d_next_point_requiring_error_indicator = d_next_point_to_sample;
            d_point_requiring_error_indicator_computed = true;
        }
        else
        {
            // Precompute next error indicator point
            // This will allow us to figure out if the greedy algorithm has terminated
            // early without needing an extra call to the error indicator function.
            getNextPointRequiringErrorIndicator();
        }
    }
}

void
GreedySampler::setPointErrorIndicator(double error, int vec_size)
{
    CAROM_VERIFY(error >= 0);
    CAROM_VERIFY(d_point_requiring_error_indicator_computed);

    double proc_errors = pow(error, 2);
    int total_vec_size = vec_size;
    CAROM_VERIFY(MPI_Allreduce(MPI_IN_PLACE,
                               &proc_errors,
                               1,
                               MPI_DOUBLE,
                               MPI_SUM,
                               MPI_COMM_WORLD) == MPI_SUCCESS);
    CAROM_VERIFY(MPI_Allreduce(MPI_IN_PLACE,
                               &total_vec_size,
                               1,
                               MPI_INT,
                               MPI_SUM,
                               MPI_COMM_WORLD) == MPI_SUCCESS);
    proc_errors = sqrt(proc_errors);
    proc_errors /= total_vec_size;

    if (!std::isfinite(proc_errors))
    {
        proc_errors = INT_MAX;
    }

    if (d_convergence_started)
    {
        setConvergenceErrorIndicator(proc_errors);
    }
    else
    {
        setSubsetErrorIndicator(proc_errors);
    }
}

void
GreedySampler::printErrorIndicator(Vector errorIndicatorPoint,
                                   double proc_errors)
{
    if (d_rank == 0)
    {
        std::string str;
        str += "Error indicator computed at [ ";
        for (int i = 0 ; i < errorIndicatorPoint.dim(); i++)
        {
            str += std::to_string(errorIndicatorPoint.item(i)) + " ";
        }
        str += "]\n";
        str += "Error indicator: " + std::to_string(proc_errors) + "\n";
        agnosticPrint(str);
    }
}

void
GreedySampler::agnosticPrint(std::string str)
{
    if (d_output_log_path == "")
    {
        std::cout << str;
    }
    else
    {
        std::ofstream database_history;
        database_history.open(d_output_log_path, std::ios::app);
        database_history << str;
        database_history.close();
    }
}

void
GreedySampler::printErrorIndicatorToleranceNotMet()
{
    if (d_rank == 0)
    {
        std::string str;
        str += "Error indicator tolerance " + std::to_string(d_error_indicator_tol) +
               " not met.\n";
        agnosticPrint(str);
    }
}

void
GreedySampler::setSubsetErrorIndicator(double proc_errors)
{
    if (d_check_local_rom || d_parameter_sampled_indices.size() == 1)
    {
        auto search = d_parameter_sampled_indices.find(
                          d_next_point_requiring_error_indicator);
        if (search != d_parameter_sampled_indices.end())
        {
            d_parameter_point_errors[d_next_point_requiring_error_indicator] = proc_errors;
            d_parameter_point_local_rom[d_next_point_requiring_error_indicator] =
                d_next_point_requiring_error_indicator;

            double old_error_indicator_tol = d_error_indicator_tol;
            if (d_parameter_sampled_indices.size() == 1)
            {
                d_error_indicator_tol = std::max(d_error_indicator_tol, proc_errors);
            }
            if (d_rank == 0)
            {
                std::string str;
                str += "Local ROM error indicator computed at [ ";
                for (int i = 0 ;
                        i < d_parameter_points[d_next_point_requiring_error_indicator].dim(); i++)
                {
                    str += std::to_string(
                               d_parameter_points[d_next_point_requiring_error_indicator].item(i)) + " ";
                }
                str += "]\n";
                str += "Local ROM error indicator (tolerance unchecked): " + std::to_string(
                           proc_errors) + "\n";
                if (old_error_indicator_tol != d_error_indicator_tol)
                {
                    str += "Error indicator at the local ROM was higher than the previous tolerance.\n";
                    str += "The error indicator tolerance should always be at least the error indicator at the local ROM.\n";
                    str += "Error indicator tolerance was adaptively changed from " +
                           std::to_string(old_error_indicator_tol) + " to " + std::to_string(
                               d_error_indicator_tol) + "\n";
                }
                agnosticPrint(str);
            }

            d_point_requiring_error_indicator_computed = false;

            // Precompute next error indicator point
            // This will allow us to figure out if the greedy algorithm has terminated
            // early without needing an extra call to the error indicator function.
            getNextPointRequiringErrorIndicator();

            return;
        }
    }

    if (proc_errors <
            d_parameter_point_errors[d_parameter_point_random_indices[d_counter]])
    {
        d_parameter_point_errors[d_parameter_point_random_indices[d_counter]] =
            proc_errors;
        d_parameter_point_local_rom[d_parameter_point_random_indices[d_counter]] =
            getNearestROMIndexToParameterPoint(d_parameter_point_random_indices[d_counter],
                                               false);
    }

    printErrorIndicator(
        d_parameter_points[d_parameter_point_random_indices[d_counter]], proc_errors);

    d_point_requiring_error_indicator_computed = false;

    if (proc_errors > d_max_error)
    {
        d_max_error = proc_errors;
        d_next_point_to_sample = d_parameter_point_random_indices[d_counter];
    }

    if (d_subset_counter == d_subset_size
            || d_counter == (int) d_parameter_points.size() - 1)
    {
        if (d_max_error < d_error_indicator_tol)
        {
            startConvergence();
        }
        else
        {
            d_iteration_started = false;
            printErrorIndicatorToleranceNotMet();
        }
    }

    // Precompute next error indicator point
    // This will allow us to figure out if the greedy algorithm has terminated
    // early without needing an extra call to the error indicator function.
    getNextPointRequiringErrorIndicator();

    return;
}

void
GreedySampler::setConvergenceErrorIndicator(double proc_errors)
{
    printErrorIndicator(d_convergence_points[d_counter], proc_errors);

    d_point_requiring_error_indicator_computed = false;

    if (proc_errors > d_max_error)
    {
        d_max_error = proc_errors;
    }

    if (proc_errors >= d_error_indicator_tol)
    {
        d_iteration_started = false;
        printErrorIndicatorToleranceNotMet();
        getNextParameterPointAfterConvergenceFailure();
        generateConvergenceSubset();
    }
    else
    {
        d_counter++;
    }

    if (d_counter == d_convergence_subset_size)
    {
        printConvergenceAchieved();
        d_procedure_completed = true;
    }

    // Precompute next error indicator point
    // This will allow us to figure out if the greedy algorithm has terminated
    // early without needing an extra call to the error indicator function.
    getNextPointRequiringErrorIndicator();
}

void
GreedySampler::generateConvergenceSubset()
{
    d_convergence_points.clear();
    d_convergence_points = generateRandomPoints(d_convergence_subset_size);
}

void
GreedySampler::startConvergence()
{
    d_convergence_started = true;
    d_max_error = 0;
    d_counter = 0;
    d_subset_counter = 0;

    if (d_rank == 0)
    {
        std::string str;
        str += "Error indicator tolerance " + std::to_string(d_error_indicator_tol) +
               " met. Computing convergence.\n";
        agnosticPrint(str);
    }

    // Precompute next error indicator point
    // This will allow us to figure out if the greedy algorithm has terminated
    // early without needing an extra call to the error indicator function.
    getNextPointRequiringErrorIndicator();
}

std::vector<Vector>
GreedySampler::generateRandomPoints(int num_points)
{
    std::vector<Vector> random_points;

    std::vector<std::uniform_real_distribution<double>> unif;
    for (int i = 0; i < d_min_param_point.dim(); i++)
    {
        unif.push_back(std::uniform_real_distribution<double>(d_min_param_point.item(i),
                       d_max_param_point.item(i)));
    }

    for (int i = 0; i < num_points; i++)
    {
        Vector point(d_min_param_point.dim(), false);
        for (int j = 0; j < point.dim(); j++)
        {
            point.item(j) = unif[j](rng);
        }
        random_points.push_back(point);
    }
    return random_points;
}

std::shared_ptr<Vector>
GreedySampler::getNearestROM(Vector point)
{

    CAROM_VERIFY(point.dim() == d_parameter_points[0].dim());

    double closest_dist_to_points = INT_MAX;
    int closest_point_index = -1;

    for (auto itr = d_parameter_sampled_indices.begin();
            itr != d_parameter_sampled_indices.end(); ++itr) {
        Vector diff;
        point.minus(d_parameter_points[*itr], diff);
        double dist = diff.norm();
        if (dist < closest_dist_to_points)
        {
            closest_dist_to_points = dist;
            closest_point_index = *itr;
        }
    }

    if (closest_point_index == -1)
    {
        return std::shared_ptr<Vector>(nullptr);
    }

    Vector* result = new Vector(d_parameter_points[closest_point_index]);
    return std::shared_ptr<Vector>(result);
}

int
GreedySampler::getNearestNonSampledPoint(Vector point)
{
    double closest_dist_to_points = INT_MAX;
    int closest_point_index = -1;

    for (int i = 0; i < d_parameter_points.size(); i++)
    {
        auto search = d_parameter_sampled_indices.find(i);
        if (search == d_parameter_sampled_indices.end())
        {
            Vector diff;
            point.minus(d_parameter_points[i], diff);
            double dist = diff.norm();
            if (dist < closest_dist_to_points)
            {
                closest_dist_to_points = dist;
                closest_point_index = i;
            }
        }
    }

    return closest_point_index;
}

int
GreedySampler::getNearestROMIndexToParameterPoint(int index, bool ignore_self)
{

    CAROM_VERIFY(index >= 0 && index < d_parameter_points.size());

    double closest_dist_to_points = INT_MAX;
    int closest_point_index = -1;

    for (auto itr = d_parameter_sampled_indices.begin();
            itr != d_parameter_sampled_indices.end(); ++itr) {
        if (index != *itr)
        {
            Vector diff;
            d_parameter_points[index].minus(d_parameter_points[*itr], diff);
            double dist = diff.norm();
            if (dist < closest_dist_to_points ||
                    (dist == closest_dist_to_points && *itr == d_parameter_point_local_rom[index]))
            {
                closest_dist_to_points = dist;
                closest_point_index = *itr;
            }
        }
        else if (!ignore_self)
        {
            closest_dist_to_points = 0;
            closest_point_index = *itr;
            break;
        }
    }

    return closest_point_index;
}

std::vector<Vector>
GreedySampler::getParameterPointDomain()
{
    return d_parameter_points;
}

std::vector<Vector>
GreedySampler::getSampledParameterPoints()
{
    std::vector<Vector> sampled_points;
    for (auto itr = d_parameter_sampled_indices.begin();
            itr != d_parameter_sampled_indices.end(); ++itr) {
        sampled_points.push_back(d_parameter_points[*itr]);
    }
    return sampled_points;
}

void
GreedySampler::save(std::string base_file_name)
{
    CAROM_ASSERT(!base_file_name.empty());

    if (d_rank == 0)
    {
        char tmp[100];
        std::string full_file_name = base_file_name;
        HDFDatabase database;
        database.create(full_file_name);

        sprintf(tmp, "num_parameter_points");
        database.putInteger(tmp, d_parameter_points.size());
        for (int i = 0; i < d_parameter_points.size(); i++)
        {
            std::string vec_path = base_file_name + "_" + std::to_string(i);
            d_parameter_points[i].write(vec_path);
        }
        for (int i = 0; i < d_convergence_points.size(); i++)
        {
            std::string vec_path = base_file_name + "_conv_" + std::to_string(i);
            d_convergence_points[i].write(vec_path);
        }

        std::string vec_path = base_file_name + "_min_point";
        d_min_param_point.write(vec_path);

        vec_path = base_file_name + "_max_point";
        d_max_param_point.write(vec_path);

        sprintf(tmp, "num_parameter_sampled_indices");
        database.putInteger(tmp, d_parameter_sampled_indices.size());
        if (d_parameter_sampled_indices.size() > 0)
        {
            sprintf(tmp, "parameter_sampled_indices");
            std::vector<int> d_parameter_sampled_indices_vec(
                d_parameter_sampled_indices.begin(), d_parameter_sampled_indices.end());
            database.putIntegerArray(tmp, &d_parameter_sampled_indices_vec[0],
                                     d_parameter_sampled_indices.size());
        }

        sprintf(tmp, "procedure_completed");
        database.putInteger(tmp, d_procedure_completed);

        if (!d_procedure_completed)
        {
            sprintf(tmp, "max_error");
            database.putDouble(tmp, d_max_error);
            sprintf(tmp, "curr_relative_error");
            database.putDouble(tmp, d_curr_relative_error);
            sprintf(tmp, "error_indicator_tol");
            database.putDouble(tmp, d_error_indicator_tol);
            sprintf(tmp, "relative_error_tol");
            database.putDouble(tmp, d_relative_error_tol);
            sprintf(tmp, "alpha");
            database.putDouble(tmp, d_alpha);
            sprintf(tmp, "max_clamp");
            database.putDouble(tmp, d_max_clamp);
            sprintf(tmp, "max_num_parameter_points");
            database.putInteger(tmp, d_num_parameter_points);
            sprintf(tmp, "subset_size");
            database.putInteger(tmp, d_subset_size);
            sprintf(tmp, "convergence_subset_size");
            database.putInteger(tmp, d_convergence_subset_size);
            sprintf(tmp, "next_point_to_sample");
            database.putInteger(tmp, d_next_point_to_sample);
            sprintf(tmp, "next_point_requiring_error_indicator");
            database.putInteger(tmp, d_next_point_requiring_error_indicator);
            sprintf(tmp, "check_local_rom");
            database.putInteger(tmp, d_check_local_rom);
            sprintf(tmp, "use_centroid");
            database.putInteger(tmp, d_use_centroid);
            sprintf(tmp, "iteration_started");
            database.putInteger(tmp, d_iteration_started);
            sprintf(tmp, "convergence_started");
            database.putInteger(tmp, d_convergence_started);
            sprintf(tmp, "next_parameter_point_computed");
            database.putInteger(tmp, d_next_parameter_point_computed);
            sprintf(tmp, "point_requiring_error_indicator_computed");
            database.putInteger(tmp, d_point_requiring_error_indicator_computed);
            sprintf(tmp, "subset_created");
            database.putInteger(tmp, d_subset_created);
            sprintf(tmp, "debug_algorithm");
            database.putInteger(tmp, d_debug_algorithm);
            sprintf(tmp, "counter");
            database.putInteger(tmp, d_counter);
            sprintf(tmp, "subset_counter");
            database.putInteger(tmp, d_subset_counter);
            sprintf(tmp, "random_seed");
            database.putInteger(tmp, d_random_seed);

            sprintf(tmp, "parameter_point_random_indices");
            database.putIntegerArray(tmp, &d_parameter_point_random_indices[0],
                                     d_parameter_point_random_indices.size());
            sprintf(tmp, "parameter_point_errors");
            database.putDoubleArray(tmp, &d_parameter_point_errors[0],
                                    d_parameter_point_errors.size());
            sprintf(tmp, "parameter_point_local_rom");
            database.putIntegerArray(tmp, &d_parameter_point_local_rom[0],
                                     d_parameter_point_local_rom.size());
        }
        database.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

bool
GreedySampler::isComplete()
{
    if (!d_procedure_completed
            && d_parameter_sampled_indices.size() == d_num_parameter_points)
    {
        if (d_rank == 0)
        {
            std::string str;
            str += "Maximum number of parameter points reached. Stopping now.\n";
            agnosticPrint(str);
        }
        d_procedure_completed = true;
    }
    return d_procedure_completed;
}

GreedySampler::~GreedySampler()
{
}

}
