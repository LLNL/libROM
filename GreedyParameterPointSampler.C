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
//              for the construction of a ROM database.

#include "GreedyParameterPointSampler.h"
#include "HDFDatabase.h"
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
createGreedyErrorIndicatorPoint(Vector* point, std::shared_ptr<Vector>& localROM)
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

GreedyParameterPointSampler::GreedyParameterPointSampler(
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
                    subset_size, convergence_subset_size, output_log_path, use_centroid, random_seed, debug_algorithm);
}

GreedyParameterPointSampler::GreedyParameterPointSampler(
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
                    subset_size, convergence_subset_size, output_log_path, use_centroid, random_seed, debug_algorithm);
}

GreedyParameterPointSampler::GreedyParameterPointSampler(
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

    constructObject(check_local_rom, relative_error_tolerance, alpha, max_clamp, subset_size,
                    convergence_subset_size, output_log_path, use_centroid, random_seed, debug_algorithm);
}

GreedyParameterPointSampler::GreedyParameterPointSampler(
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

    constructObject(check_local_rom, relative_error_tolerance, alpha, max_clamp, subset_size,
                    convergence_subset_size, output_log_path, use_centroid, random_seed, debug_algorithm);
}

GreedyParameterPointSampler::GreedyParameterPointSampler(
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
GreedyParameterPointSampler::addDatabaseFromFile(
    std::string const& warm_start_file_name)
{
    char tmp[100];
    sprintf(tmp, ".%06d", d_rank);
    std::string full_file_name = warm_start_file_name + tmp;
    HDFDatabase database;
    database.open(full_file_name);

    sprintf(tmp, "num_parameter_sampled_indices");
    int num_parameter_sampled_indices;
    database.getInteger(tmp, num_parameter_sampled_indices);
    if (num_parameter_sampled_indices > 0)
    {
        int temp_parameter_sampled_indices[num_parameter_sampled_indices];
        sprintf(tmp, "parameter_sampled_indices");
        database.getIntegerArray(tmp, &temp_parameter_sampled_indices[0], num_parameter_sampled_indices);
        for (int i = 0; i < num_parameter_sampled_indices; i++)
        {
            std::string vec_path = warm_start_file_name + "_" + std::to_string(i);
            Vector point;
            point.read(vec_path);
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
GreedyParameterPointSampler::load(
    std::string base_file_name)
{
    char tmp[100];
    sprintf(tmp, ".%06d", d_rank);
    std::string full_file_name = base_file_name + tmp;
    HDFDatabase database;
    database.open(full_file_name);

    sprintf(tmp, "num_parameter_points");
    int num_parameter_points;
    database.getInteger(tmp, num_parameter_points);
    for (int i = 0; i < num_parameter_points; i++)
    {
        std::string vec_path = base_file_name + "_" + std::to_string(i);
        Vector point;
        point.read(vec_path);
        d_parameter_points.push_back(point);
    }

    std::string vec_path = base_file_name + "_min_point";
    d_min_param_point.read(vec_path);

    vec_path = base_file_name + "_max_point";
    d_max_param_point.read(vec_path);

    sprintf(tmp, "num_parameter_sampled_indices");
    int num_parameter_sampled_indices;
    database.getInteger(tmp, num_parameter_sampled_indices);
    if (num_parameter_sampled_indices > 0)
    {
        int temp_parameter_sampled_indices[num_parameter_sampled_indices];
        sprintf(tmp, "parameter_sampled_indices");
        database.getIntegerArray(tmp, &temp_parameter_sampled_indices[0], num_parameter_sampled_indices);
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
        database.getIntegerArray(tmp, &d_parameter_point_random_indices[0], d_parameter_points.size());

        sprintf(tmp, "parameter_point_errors");
        d_parameter_point_errors.resize(d_parameter_points.size());
        database.getDoubleArray(tmp, &d_parameter_point_errors[0], d_parameter_points.size());

        sprintf(tmp, "parameter_point_local_rom");
        d_parameter_point_local_rom.resize(d_parameter_points.size());
        database.getIntegerArray(tmp, &d_parameter_point_local_rom[0], d_parameter_points.size());

        for (int i = 0; i < d_convergence_subset_size; i++)
        {
            std::string vec_path = base_file_name + "_conv_" + std::to_string(i);
            Vector point;
            point.read(vec_path);
            d_convergence_points.push_back(point);
        }
    }
    database.close();

    rng.seed(d_random_seed + d_parameter_sampled_indices.size());
}

void
GreedyParameterPointSampler::constructParameterPoints()
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
        if (d_output_log_path == "")
        {
            std::cout << "Total number of sample points: " << d_num_parameter_points << std::endl;
            std::cout << "Parameter space minimum: [ ";
            for (int i = 0 ; i < d_min_param_point.dim(); i++)
            {
                std::cout << d_min_param_point.item(i) << " ";
            }
            std::cout << "]" << std::endl;
            std::cout << "Parameter space maximum: [ ";
            for (int i = 0 ; i < d_max_param_point.dim(); i++)
            {
                std::cout << d_max_param_point.item(i) << " ";
            }
            std::cout << "]" << std::endl;
        }
        else
        {
            std::ofstream database_history;
            database_history.open(d_output_log_path, std::ios::app);
            database_history << "Total number of sample points: " << d_num_parameter_points << std::endl;
            database_history << "Parameter space minimum: [ ";
            for (int i = 0 ; i < d_min_param_point.dim(); i++)
            {
                database_history << d_min_param_point.item(i) << " ";
            }
            database_history << "]" << std::endl;
            database_history << "Parameter space maximum: [ ";
            for (int i = 0 ; i < d_max_param_point.dim(); i++)
            {
                database_history << d_max_param_point.item(i) << " ";
            }
            database_history << "]" << std::endl;
            database_history.close();
        }
    }
}

void
GreedyParameterPointSampler::constructObject(
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
    CAROM_VERIFY(max_clamp >= 0.0 || max_clamp == -1);
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
        if (d_output_log_path == "")
        {
            std::cout << "Greedy relative error tolerance: " << d_relative_error_tol << std::endl;
            std::cout << "Greedy alpha constant: " << d_alpha << std::endl;
            std::cout << "Greedy max clamp constant: " << d_max_clamp << std::endl;
            std::cout << "Greedy iteration subset size: " << d_subset_size << std::endl;
            std::cout << "Greedy convergence subset size: " << d_convergence_subset_size << std::endl;
        }
        else
        {
            std::ofstream database_history;
            database_history.open(d_output_log_path, std::ios::app);
            database_history << "Greedy relative error tolerance: " << d_relative_error_tol << std::endl;
            database_history << "Greedy alpha constant: " << d_alpha << std::endl;
            database_history << "Greedy max clamp constant: " << d_max_clamp << std::endl;
            database_history << "Greedy iteration subset size: " << d_subset_size << std::endl;
            database_history << "Greedy convergence subset size: " << d_convergence_subset_size << std::endl;
            database_history.close();
        }
    }
}

void
GreedyParameterPointSampler::initializeParameterPoints()
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
GreedyParameterPointSampler::getNextParameterPoint()
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
        // get center-most point
        // obtain the centroid and find the point closest to centroid
        if (d_use_centroid)
        {
            Vector centroid(d_parameter_points[0]);
            for (int i = 1; i < d_parameter_points.size(); i++) {
                centroid += d_parameter_points[i];
            }
            centroid.mult(1.0 / d_parameter_points.size(), centroid);

            double min_dist_to_centroid;
            for (int i = 0; i < d_parameter_points.size(); i++) {
                Vector diff;
                centroid.minus(d_parameter_points[i], diff);
                double dist_to_centroid = diff.norm();
                if (i == 0 || dist_to_centroid < min_dist_to_centroid)
                {
                    min_dist_to_centroid = dist_to_centroid;
                    d_next_point_to_sample = i;
                }
            }
        }

        // otherwise, do a brute force search to obtain
        // the point closest to all other points
        else
        {
            std::vector<double> dist_to_points(d_parameter_points.size(), 0);
            for (int i = 0; i < d_parameter_points.size(); i++) {
                for (int j = i + 1; j < d_parameter_points.size(); j++) {
                    Vector diff;
                    d_parameter_points[i].minus(d_parameter_points[j], diff);
                    double dist = diff.norm();
                    dist_to_points[i] += dist;
                    dist_to_points[j] += dist;
                }
            }

            double closest_dist_to_points = INT_MAX;
            for (int i = 0; i < dist_to_points.size(); i++) {
                if (dist_to_points[i] < closest_dist_to_points)
                {
                    closest_dist_to_points = dist_to_points[i];
                    d_next_point_to_sample = i;
                }
            }
        }
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
        if (d_output_log_path == "")
        {
            std::cout << "Point sampled at [ ";
            for (int i = 0 ; i < d_parameter_points[curr_point_to_sample].dim(); i++)
            {
                std::cout << d_parameter_points[curr_point_to_sample].item(i) << " ";
            }
            std::cout << "]" << std::endl;
        }
        else
        {
            std::ofstream database_history;
            database_history.open(d_output_log_path, std::ios::app);
            database_history << "Point sampled at [ ";
            for (int i = 0 ; i < d_parameter_points[curr_point_to_sample].dim(); i++)
            {
                database_history << d_parameter_points[curr_point_to_sample].item(i) << " ";
            }
            database_history << "]" << std::endl;
            database_history.close();
        }
    }

    d_next_parameter_point_computed = true;

    Vector* result = new Vector(d_parameter_points[curr_point_to_sample]);
    return std::shared_ptr<Vector>(result);
}

struct GreedyErrorIndicatorPoint
GreedyParameterPointSampler::getNextPointRequiringRelativeError()
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
        result2 = new Vector(d_parameter_points[getNearestROMIndex(d_next_point_to_sample, false)]);
    }
    else
    {
        result2 = new Vector(d_parameter_points[getNearestROMIndex(d_next_point_to_sample, true)]);
    }

    return createGreedyErrorIndicatorPoint(result1, result2);
}

struct GreedyErrorIndicatorPoint
GreedyParameterPointSampler::getNextPointRequiringErrorIndicator()
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
GreedyParameterPointSampler::getNextSubsetPointRequiringErrorIndicator()
{
    if (d_point_requiring_error_indicator_computed)
    {
        Vector* result1 = new Vector(d_parameter_points[d_next_point_requiring_error_indicator]);
        Vector* result2 = new Vector(d_parameter_points[getNearestROMIndex(d_next_point_requiring_error_indicator, false)]);
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
            std::shuffle(d_parameter_point_random_indices.begin(), d_parameter_point_random_indices.end(), rng);
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
        auto search = d_parameter_sampled_indices.find(d_parameter_point_random_indices[d_counter]);
        if (search == d_parameter_sampled_indices.end())
        {
            d_subset_counter++;
            double curr_error = d_parameter_point_errors[d_parameter_point_random_indices[d_counter]];
            if (curr_error > d_max_error)
            {
                // if we have already computed this error indicator at the same local rom, the error indicator will not improve
                // no need to calculate the error indicator again
                if (d_parameter_point_local_rom[d_parameter_point_random_indices[d_counter]] == getNearestROMIndex(d_parameter_point_random_indices[d_counter], false))
                {
                    d_max_error = curr_error;
                    d_next_point_to_sample = d_parameter_point_random_indices[d_counter];
                    if (d_rank == 0)
                    {
                        if (d_output_log_path == "")
                        {
                            std::cout << "Error indicator at [ ";
                            for (int i = 0 ; i < d_parameter_points[d_parameter_point_random_indices[d_counter]].dim(); i++)
                            {
                                std::cout << d_parameter_points[d_parameter_point_random_indices[d_counter]].item(i) << " ";
                            }
                            std::cout << "] skipped." << std::endl;
                            std::cout << "Error indicator " << curr_error << " already computed at the same local ROM." << std::endl;
                        }
                        else
                        {
                            std::ofstream database_history;
                            database_history.open(d_output_log_path, std::ios::app);
                            database_history << "Error indicator at [ ";
                            for (int i = 0 ; i < d_parameter_points[d_parameter_point_random_indices[d_counter]].dim(); i++)
                            {
                                database_history << d_parameter_points[d_parameter_point_random_indices[d_counter]].item(i) << " ";
                            }
                            database_history << "] skipped." << std::endl;
                            database_history << "Error indicator " << curr_error << " already computed at the same local ROM." << std::endl;
                            database_history.close();
                        }
                    }
                }
                else
                {
                    d_next_point_requiring_error_indicator = d_parameter_point_random_indices[d_counter];
                    d_point_requiring_error_indicator_computed = true;
                    Vector* result1 = new Vector(d_parameter_points[d_next_point_requiring_error_indicator]);
                    Vector* result2 = new Vector(d_parameter_points[getNearestROMIndex(d_next_point_requiring_error_indicator, false)]);
                    return createGreedyErrorIndicatorPoint(result1, result2);
                }
            }
            else
            {
                if (d_rank == 0)
                {
                    if (d_output_log_path == "")
                    {
                        std::cout << "Error indicator at [ ";
                        for (int i = 0 ; i < d_parameter_points[d_parameter_point_random_indices[d_counter]].dim(); i++)
                        {
                            std::cout << d_parameter_points[d_parameter_point_random_indices[d_counter]].item(i) << " ";
                        }
                        std::cout << "] skipped." << std::endl;
                        std::cout << "Error indicator " << curr_error
                                  << " is less than current max error " << d_max_error << std::endl;
                    }
                    else
                    {
                        std::ofstream database_history;
                        database_history.open(d_output_log_path, std::ios::app);
                        database_history << "Error indicator at [ ";
                        for (int i = 0 ; i < d_parameter_points[d_parameter_point_random_indices[d_counter]].dim(); i++)
                        {
                            database_history << d_parameter_points[d_parameter_point_random_indices[d_counter]].item(i) << " ";
                        }
                        database_history << "] skipped." << std::endl;
                        database_history << "Error indicator " << curr_error
                                         << " is less than current max error " << d_max_error << std::endl;
                        database_history.close();
                    }
                }
            }
        }
    }
    if (d_next_point_requiring_error_indicator == -1)
    {
        if (d_rank == 0)
        {
            if (d_output_log_path == "")
            {
                std::cout << "Ran out of points to calculate error indicator for in this iteration." << std::endl;
            }
            else
            {
                std::ofstream database_history;
                database_history.open(d_output_log_path, std::ios::app);
                database_history << "Ran out of points to calculate error indicator for in this iteration." << std::endl;
                database_history.close();
            }
        }
        if (d_max_error < d_error_indicator_tol)
        {
            startConvergence();
        }
        else
        {
            d_iteration_started = false;
        }
    }

    return createGreedyErrorIndicatorPoint(nullptr, nullptr);
}

struct GreedyErrorIndicatorPoint
GreedyParameterPointSampler::getNextConvergencePointRequiringErrorIndicator()
{
    if (d_point_requiring_error_indicator_computed)
    {
        Vector* result1 = new Vector(d_convergence_points[d_next_point_requiring_error_indicator]);
        std::shared_ptr<Vector> result2 = getNearestROM(d_convergence_points[d_next_point_requiring_error_indicator]);
        return createGreedyErrorIndicatorPoint(result1, result2);
    }
    if (d_counter == d_convergence_subset_size)
    {
        return createGreedyErrorIndicatorPoint(nullptr, nullptr);
    }

    d_next_point_requiring_error_indicator = -1;

    //get next point requiring error indicator
    while (d_counter < (int) d_convergence_points.size())
    {
        d_next_point_requiring_error_indicator = d_counter;
        d_point_requiring_error_indicator_computed = true;
        Vector* result1 = new Vector(d_convergence_points[d_next_point_requiring_error_indicator]);
        std::shared_ptr<Vector> result2 = getNearestROM(d_convergence_points[d_next_point_requiring_error_indicator]);
        return createGreedyErrorIndicatorPoint(result1, result2);
    }

    if (d_next_point_requiring_error_indicator == -1)
    {
        if (d_rank == 0)
        {
            if (d_output_log_path == "")
            {
                std::cout << "Convergence achieved.";
            }
            else
            {
                std::ofstream database_history;
                database_history.open(d_output_log_path, std::ios::app);
                database_history << "Convergence achieved." << std::endl;
                database_history.close();
            }
        }
        d_procedure_completed = true;
    }

    return createGreedyErrorIndicatorPoint(nullptr, nullptr);
}

void
GreedyParameterPointSampler::setPointRelativeError(double error)
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
            std::cout << "Relative error computed at [ ";
            for (int i = 0 ; i < d_parameter_points[d_next_point_to_sample].dim(); i++)
            {
                std::cout << d_parameter_points[d_next_point_to_sample].item(i) << " ";
            }
            std::cout << "]" << std::endl;
            std::cout << "Relative error: " << error << std::endl;
        }
        else
        {
            std::ofstream database_history;
            database_history.open(d_output_log_path, std::ios::app);
            database_history << "Relative error computed at [ ";
            for (int i = 0 ; i < d_parameter_points[d_next_point_to_sample].dim(); i++)
            {
                database_history << d_parameter_points[d_next_point_to_sample].item(i) << " ";
            }
            database_history << "]" << std::endl;
            database_history << "Relative error: " << error << std::endl;
            database_history.close();
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
        double min2 = d_relative_error_tol * d_parameter_point_errors[d_next_point_to_sample] / d_curr_relative_error;

        if (d_curr_relative_error <= d_relative_error_tol)
        {
            if (d_rank == 0)
            {
                if (d_output_log_path == "")
                {
                    std::cout << "The relative error was smaller than the Greedy relative error tolerance. The error indicator tolerance must increase." << std::endl;
                    std::cout << "Greedy alpha constant * Greedy relative error tolerance: " << max1 << std::endl;
                    std::cout << "The minimum value the greedy error indicator tolerance must increase to is: " << max1 << std::endl;
                    std::cout << "Greedy max clamp constant * Current max error indicator: " << min1 << std::endl;
                    std::cout << "Greedy relative error tolerance * Current max error indicator / Current relative error: " << min2 << std::endl;
                    std::cout << "The maximum value the greedy error indicator tolerance is allowed to change is the minimum of the previous two values: " << std::min(min1, min2) << std::endl;
                }
                else
                {
                    std::ofstream database_history;
                    database_history.open(d_output_log_path, std::ios::app);
                    database_history << "The relative error was smaller than the Greedy relative error tolerance. The error indicator tolerance must increase." << std::endl;
                    database_history << "Greedy alpha constant * Greedy relative error tolerance: " << max1 << std::endl;
                    database_history << "The minimum value the greedy error indicator tolerance must increase to is: " << max1 << std::endl;
                    database_history << "Greedy max clamp constant * Current max error indicator: " << min1 << std::endl;
                    database_history << "Greedy relative error tolerance * Current max error indicator / Current relative error: " << min2 << std::endl;
                    database_history << "The maximum value the greedy error indicator tolerance is allowed to change is the minimum of the previous two values: " << std::min(min1, min2) << std::endl;
                    database_history.close();
                }
            }
            d_error_indicator_tol = std::max(d_alpha * d_error_indicator_tol, std::min(d_max_clamp * d_parameter_point_errors[d_next_point_to_sample], d_relative_error_tol * d_parameter_point_errors[d_next_point_to_sample] / d_curr_relative_error));
        }
        else
        {
            double min1 = d_relative_error_tol * d_parameter_point_errors[d_next_point_to_sample] / d_curr_relative_error;

            if (d_rank == 0)
            {
                if (d_output_log_path == "")
                {
                    std::cout << "The relative error was larger than the Greedy relative error tolerance. The error indicator tolerance must decrease." << std::endl;
                    std::cout << "Greedy alpha constant * Greedy relative error tolerance: " << d_error_indicator_tol<< std::endl;
                    std::cout << "Greedy relative error tolerance * Current max error indicator / Current relative error: " << min1 << std::endl;
                    std::cout << "The minimum value the greedy error indicator tolerance is allowed to change is: " << min1 << std::endl;
                }
                else
                {
                    std::ofstream database_history;
                    database_history.open(d_output_log_path, std::ios::app);
                    database_history << "The relative error was larger than the Greedy relative error tolerance. The error indicator tolerance must decrease." << std::endl;
                    database_history << "Greedy alpha constant * Greedy relative error tolerance: " << d_error_indicator_tol<< std::endl;
                    database_history << "Greedy relative error tolerance * Current max error indicator / Current relative error: " << min1 << std::endl;
                    database_history << "The minimum value the greedy error indicator tolerance is allowed to change is: " << min1 << std::endl;
                    database_history.close();
                }
            }
            d_error_indicator_tol = std::min(d_error_indicator_tol, d_relative_error_tol * d_parameter_point_errors[d_next_point_to_sample] / d_curr_relative_error);
        }
        if (d_rank == 0)
        {
            if (d_output_log_path == "")
            {
                if (old_error_indicator_tol != d_error_indicator_tol)
                {
                    std::cout << "Error indicator tolerance was adaptively changed from " << old_error_indicator_tol << " to " << d_error_indicator_tol << std::endl;
                }
            }
            else
            {
                std::ofstream database_history;
                database_history.open(d_output_log_path, std::ios::app);
                if (old_error_indicator_tol != d_error_indicator_tol)
                {
                    database_history << "Error indicator tolerance was adaptively changed from " << old_error_indicator_tol << " to " << d_error_indicator_tol << std::endl;
                }
                database_history.close();
            }
        }
    }

    d_parameter_point_errors[d_next_point_to_sample] = 0;
    d_parameter_point_local_rom[d_next_point_to_sample] = d_next_point_to_sample;

    if (d_parameter_sampled_indices.size() > 1 && d_curr_relative_error <= d_relative_error_tol)
    {
        if (d_rank == 0)
        {
            if (d_output_log_path == "")
            {
                std::cout << "The relative error was smaller than the Greedy relative error tolerance. Computing convergence." << std::endl;
            }
            else
            {
                std::ofstream database_history;
                database_history.open(d_output_log_path, std::ios::app);
                database_history << "The relative error was smaller than the Greedy relative error tolerance. Computing convergence." << std::endl;
                database_history.close();
            }
        }
        startConvergence();
    }
    else
    {
        if (d_check_local_rom)
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
GreedyParameterPointSampler::setPointErrorIndicator(double error, int vec_size)
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
GreedyParameterPointSampler::printErrorIndicator(Vector errorIndicatorPoint, double proc_errors)
{
    if (d_rank == 0)
    {
        if (d_output_log_path == "")
        {
            std::cout << "Error indicator computed at [ ";
            for (int i = 0 ; i < errorIndicatorPoint.dim(); i++)
            {
                std::cout << errorIndicatorPoint.item(i) << " ";
            }
            std::cout << "]" << std::endl;
            std::cout << "Error indicator: " << proc_errors << std::endl;
        }
        else
        {
            std::ofstream database_history;
            database_history.open(d_output_log_path, std::ios::app);
            database_history << "Error indicator computed at [ ";
            for (int i = 0 ; i < errorIndicatorPoint.dim(); i++)
            {
                database_history << errorIndicatorPoint.item(i) << " ";
            }
            database_history << "]" << std::endl;
            database_history << "Error indicator: " << proc_errors << std::endl;
            database_history.close();
        }
    }
}

void
GreedyParameterPointSampler::printErrorIndicatorToleranceNotMet()
{
    if (d_rank == 0)
    {
        if (d_output_log_path == "")
        {
            std::cout << "Error indicator tolerance " << d_error_indicator_tol << " not met." << std::endl;
        }
        else
        {
            std::ofstream database_history;
            database_history.open(d_output_log_path, std::ios::app);
            database_history << "Error indicator tolerance " << d_error_indicator_tol << " not met." << std::endl;
            database_history.close();
        }
    }
}

void
GreedyParameterPointSampler::setSubsetErrorIndicator(double proc_errors)
{
    if (d_check_local_rom)
    {
        auto search = d_parameter_sampled_indices.find(d_next_point_requiring_error_indicator);
        if (search != d_parameter_sampled_indices.end())
        {
            d_parameter_point_errors[d_next_point_requiring_error_indicator] = proc_errors;
            d_parameter_point_local_rom[d_next_point_requiring_error_indicator] = d_next_point_requiring_error_indicator;

            if (d_parameter_sampled_indices.size() == 1)
            {
                double old_error_indicator_tol = d_error_indicator_tol;
                d_error_indicator_tol = std::max(d_error_indicator_tol, proc_errors);

                if (d_rank == 0)
                {
                    if (d_output_log_path == "")
                    {
                        std::cout << "Local ROM Error indicator computed at [ ";
                        for (int i = 0 ; i < d_parameter_points[d_next_point_requiring_error_indicator].dim(); i++)
                        {
                            std::cout << d_parameter_points[d_next_point_requiring_error_indicator].item(i) << " ";
                        }
                        std::cout << "]" << std::endl;
                        std::cout << "Local ROM Error indicator (tolerance unchecked): " << proc_errors << std::endl;
                        if (old_error_indicator_tol != d_error_indicator_tol)
                        {
                            std::cout << "Error indicator at the local ROM was higher than the previous tolerance." << std::endl;
                            std::cout << "The error indicator tolerance should always be at least the error indicator at the local ROM since this error indicator is a lower bound." << std::endl;
                            std::cout << "Error indicator tolerance was adaptively changed from " << old_error_indicator_tol << " to " << d_error_indicator_tol << std::endl;
                        }
                    }
                    else
                    {
                        std::ofstream database_history;
                        database_history.open(d_output_log_path, std::ios::app);
                        database_history << "Local ROM Error indicator computed at [ ";
                        for (int i = 0 ; i < d_parameter_points[d_next_point_requiring_error_indicator].dim(); i++)
                        {
                            database_history << d_parameter_points[d_next_point_requiring_error_indicator].item(i) << " ";
                        }
                        database_history << "]" << std::endl;
                        database_history << "Local ROM Error indicator (tolerance unchecked): " << proc_errors << std::endl;
                        if (old_error_indicator_tol != d_error_indicator_tol)
                        {
                            database_history << "Error indicator at the local ROM was higher than the previous tolerance." << std::endl;
                            database_history << "The error indicator tolerance should always be at least the error indicator at the local ROM since this error indicator is a lower bound." << std::endl;
                            database_history << "Error indicator tolerance was adaptively changed from " << old_error_indicator_tol << " to " << d_error_indicator_tol << std::endl;
                        }
                        database_history.close();
                    }
                }
            }

            d_point_requiring_error_indicator_computed = false;

            // Precompute next error indicator point
            // This will allow us to figure out if the greedy algorithm has terminated
            // early without needing an extra call to the error indicator function.
            getNextPointRequiringErrorIndicator();

            return;
        }
    }

    if (proc_errors < d_parameter_point_errors[d_parameter_point_random_indices[d_counter]])
    {
        d_parameter_point_errors[d_parameter_point_random_indices[d_counter]] = proc_errors;
        d_parameter_point_local_rom[d_parameter_point_random_indices[d_counter]] = getNearestROMIndex(d_parameter_point_random_indices[d_counter], false);
    }

    printErrorIndicator(d_parameter_points[d_parameter_point_random_indices[d_counter]], proc_errors);

    d_point_requiring_error_indicator_computed = false;

    if (proc_errors > d_max_error)
    {
        d_max_error = proc_errors;
        d_next_point_to_sample = d_parameter_point_random_indices[d_counter];
    }

    if (d_subset_counter == d_subset_size || d_counter == (int) d_parameter_points.size() - 1)
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
GreedyParameterPointSampler::setConvergenceErrorIndicator(double proc_errors)
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
        d_procedure_completed = true;
    }

    // Precompute next error indicator point
    // This will allow us to figure out if the greedy algorithm has terminated
    // early without needing an extra call to the error indicator function.
    getNextPointRequiringErrorIndicator();
}

void
GreedyParameterPointSampler::generateConvergenceSubset()
{
    d_convergence_points.clear();
    d_convergence_points = generateRandomPoints(d_convergence_subset_size);
}

void
GreedyParameterPointSampler::startConvergence()
{
    d_convergence_started = true;
    d_max_error = 0;
    d_counter = 0;
    d_subset_counter = 0;

    if (d_rank == 0)
    {
        if (d_output_log_path == "")
        {
            std::cout << "Error indicator tolerance " << d_error_indicator_tol << " met. Computing convergence." << std::endl;
        }
        else
        {
            std::ofstream database_history;
            database_history.open(d_output_log_path, std::ios::app);
            database_history << "Error indicator tolerance " << d_error_indicator_tol << " met. Computing convergence." << std::endl;
            database_history.close();
        }
    }

    // Precompute next error indicator point
    // This will allow us to figure out if the greedy algorithm has terminated
    // early without needing an extra call to the error indicator function.
    getNextPointRequiringErrorIndicator();
}

std::vector<Vector>
GreedyParameterPointSampler::generateRandomPoints(int num_points)
{
    std::vector<Vector> random_points;

    std::vector<std::uniform_real_distribution<double>> unif;
    for (int i = 0; i < d_min_param_point.dim(); i++)
    {
        unif.push_back(std::uniform_real_distribution<double>(d_min_param_point.item(i), d_max_param_point.item(i)));
    }

    for (int i = 0; i < d_convergence_subset_size; i++)
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
GreedyParameterPointSampler::getNearestROM(Vector point)
{

    CAROM_VERIFY(point.dim() == d_parameter_points[0].dim());

    double closest_dist_to_points = INT_MAX;
    int closest_point_index = -1;

    for (auto itr = d_parameter_sampled_indices.begin(); itr != d_parameter_sampled_indices.end(); ++itr) {
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
GreedyParameterPointSampler::getNearestNonSampledPoint(Vector point)
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
GreedyParameterPointSampler::getNearestROMIndex(int index, bool ignore_self)
{

    CAROM_VERIFY(index >= 0 && index < d_parameter_points.size());

    double closest_dist_to_points = INT_MAX;
    int closest_point_index = -1;

    for (auto itr = d_parameter_sampled_indices.begin(); itr != d_parameter_sampled_indices.end(); ++itr) {
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
GreedyParameterPointSampler::getParameterPointDomain()
{
    return d_parameter_points;
}

std::vector<Vector>
GreedyParameterPointSampler::getSampledParameterPoints()
{
    std::vector<Vector> sampled_points;
    for (auto itr = d_parameter_sampled_indices.begin(); itr != d_parameter_sampled_indices.end(); ++itr) {
        sampled_points.push_back(d_parameter_points[*itr]);
    }
    return sampled_points;
}

void
GreedyParameterPointSampler::save(std::string base_file_name)
{
    CAROM_ASSERT(!base_file_name.empty());

    char tmp[100];
    sprintf(tmp, ".%06d", d_rank);
    std::string full_file_name = base_file_name + tmp;
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
        std::vector<int> d_parameter_sampled_indices_vec(d_parameter_sampled_indices.begin(), d_parameter_sampled_indices.end());
        database.putIntegerArray(tmp, &d_parameter_sampled_indices_vec[0], d_parameter_sampled_indices.size());

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
        database.putIntegerArray(tmp, &d_parameter_point_random_indices[0], d_parameter_point_random_indices.size());
        sprintf(tmp, "parameter_point_errors");
        database.putDoubleArray(tmp, &d_parameter_point_errors[0], d_parameter_point_errors.size());
        sprintf(tmp, "parameter_point_local_rom");
        database.putIntegerArray(tmp, &d_parameter_point_local_rom[0], d_parameter_point_local_rom.size());
    }
    database.close();
}

bool
GreedyParameterPointSampler::isComplete()
{
    if (d_parameter_sampled_indices.size() == d_num_parameter_points)
    {
        d_procedure_completed = true;
    }
    return d_procedure_completed;
}

GreedyParameterPointSampler::~GreedyParameterPointSampler()
{
}

}
