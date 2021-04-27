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

#include "GreedyParameterPointSelector.h"
#include "Vector.h"
#include "HDFDatabase.h"
#include "mpi.h"
#include <cmath>
#include <limits.h>
#include <fstream>

namespace CAROM {

struct GreedyResidualPoint
createGreedyResidualPoint(Vector* point, Vector* localROM)
{
    struct GreedyResidualPoint result;
    result.point = std::shared_ptr<Vector>(point);
    result.localROM = std::shared_ptr<Vector>(localROM);
    return result;
}

struct GreedyResidualPoint
createGreedyResidualPoint(Vector* point, std::shared_ptr<Vector>& localROM)
{
    struct GreedyResidualPoint result;
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

GreedyParameterPointSelector::GreedyParameterPointSelector(
    std::vector<Vector> parameter_points,
    bool check_local_rom,
    double tolerance,
    double saturation,
    int subset_size,
    int convergence_subset_size,
    std::string output_log_path,
    std::string warm_start_file_name,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm)
{
    constructObject(check_local_rom, tolerance, saturation,
                    subset_size, convergence_subset_size, output_log_path, use_centroid, random_seed, debug_algorithm);

    initializeParameterPoints(parameter_points);

    if (warm_start_file_name != "")
    {
        addDatabaseFromFile(warm_start_file_name);
    }
}

GreedyParameterPointSelector::GreedyParameterPointSelector(
    std::vector<double> parameter_points,
    bool check_local_rom,
    double tolerance,
    double saturation,
    int subset_size,
    int convergence_subset_size,
    std::string output_log_path,
    std::string warm_start_file_name,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm)
{
    std::vector<Vector> parameter_points_vec;
    for (int i = 0; i < parameter_points.size(); i++) {
        Vector vec(1, false);
        vec.item(0) = parameter_points[i];
        parameter_points_vec.push_back(vec);
    }

    constructObject(check_local_rom, tolerance, saturation,
                    subset_size, convergence_subset_size, output_log_path, use_centroid, random_seed, debug_algorithm);

    initializeParameterPoints(parameter_points_vec);

    if (warm_start_file_name != "")
    {
        addDatabaseFromFile(warm_start_file_name);
    }
}

GreedyParameterPointSelector::GreedyParameterPointSelector(
    double param_space_min,
    double param_space_max,
    int param_space_size,
    bool check_local_rom,
    double tolerance,
    double saturation,
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

    constructObject(check_local_rom, tolerance, saturation, subset_size,
                    convergence_subset_size, output_log_path, use_centroid, random_seed, debug_algorithm);

    std::vector<Vector> parameter_points_vec =
        constructParameterPoints(param_space_min_vec, param_space_max_vec, param_space_size);

    initializeParameterPoints(parameter_points_vec);

    if (warm_start_file_name != "")
    {
        addDatabaseFromFile(warm_start_file_name);
    }
}

GreedyParameterPointSelector::GreedyParameterPointSelector(
    Vector param_space_min,
    Vector param_space_max,
    int param_space_size,
    bool check_local_rom,
    double tolerance,
    double saturation,
    int subset_size,
    int convergence_subset_size,
    std::string output_log_path,
    std::string warm_start_file_name,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm)
{
    constructObject(check_local_rom, tolerance, saturation, subset_size,
                    convergence_subset_size, output_log_path, use_centroid, random_seed, debug_algorithm);

    std::vector<Vector> parameter_points_vec =
        constructParameterPoints(param_space_min, param_space_max, param_space_size);

    initializeParameterPoints(parameter_points_vec);

    if (warm_start_file_name != "")
    {
        addDatabaseFromFile(warm_start_file_name);
    }
}

GreedyParameterPointSelector::GreedyParameterPointSelector(
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
GreedyParameterPointSelector::addDatabaseFromFile(
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
GreedyParameterPointSelector::load(
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
        sprintf(tmp, "tol");
        database.getDouble(tmp, d_tol);
        sprintf(tmp, "sat");
        database.getDouble(tmp, d_sat);
        sprintf(tmp, "subset_size");
        database.getInteger(tmp, d_subset_size);
        sprintf(tmp, "convergence_subset_size");
        database.getInteger(tmp, d_convergence_subset_size);
        sprintf(tmp, "next_point_to_sample");
        database.getInteger(tmp, d_next_point_to_sample);
        sprintf(tmp, "next_point_requiring_residual");
        database.getInteger(tmp, d_next_point_requiring_residual);
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
        sprintf(tmp, "point_requiring_residual_computed");
        database.getInteger(tmp, bool_int_temp);
        d_point_requiring_residual_computed = (bool) bool_int_temp;
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
}

std::vector<Vector>
GreedyParameterPointSelector::constructParameterPoints(
    Vector param_space_min,
    Vector param_space_max,
    int param_space_size)
{
    CAROM_VERIFY(param_space_min.dim() == param_space_max.dim());
    CAROM_VERIFY(param_space_size >= 1);

    bool isGreater = false;
    for (int i = 0; i < param_space_min.dim(); i++)
    {
        if (param_space_size == 1)
        {
            if (param_space_max.item(i) >= param_space_min.item(i))
            {
                isGreater = true;
                break;
            }
        }
        else
        {
            if (param_space_max.item(i) > param_space_min.item(i))
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
            std::cout << "Total number of sample points: " << param_space_size << std::endl;
            std::cout << "Parameter space minimum: [ ";
            for (int i = 0 ; i < param_space_min.dim(); i++)
            {
                std::cout << param_space_min.item(i) << " ";
            }
            std::cout << "]" << std::endl;
            std::cout << "Parameter space maximum: [ ";
            for (int i = 0 ; i < param_space_max.dim(); i++)
            {
                std::cout << param_space_max.item(i) << " ";
            }
            std::cout << "]" << std::endl;
        }
        else
        {
            std::ofstream database_history;
            database_history.open(d_output_log_path, std::ios::app);
            database_history << "Total number of sample points: " << param_space_size << std::endl;
            database_history << "Parameter space minimum: [ ";
            for (int i = 0 ; i < param_space_min.dim(); i++)
            {
                database_history << param_space_min.item(i) << " ";
            }
            database_history << "]" << std::endl;
            database_history << "Parameter space maximum: [ ";
            for (int i = 0 ; i < param_space_max.dim(); i++)
            {
                database_history << param_space_max.item(i) << " ";
            }
            database_history << "]" << std::endl;
            database_history.close();
        }
    }

    std::vector<Vector> parameter_points_vec;
    std::vector<double> frequencies;
    Vector vec(param_space_min.dim(), false);
    for (int i = 0; i < param_space_min.dim(); i++)
    {
        frequencies.push_back(std::abs(param_space_max.item(i) - param_space_min.item(i)) / (param_space_size - 1));
    }
    for (int i = 0; i < param_space_size; i++) {
        for (int j = 0; j < param_space_min.dim(); j++)
        {
            vec.item(j) = param_space_min.item(j) + i * frequencies[j];
        }
        parameter_points_vec.push_back(vec);
    }

    return parameter_points_vec;
}

void
GreedyParameterPointSelector::constructObject(
    bool check_local_rom,
    double tolerance,
    double saturation,
    int subset_size,
    int convergence_subset_size,
    std::string output_log_path,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm)
{
    CAROM_VERIFY(tolerance > 0.0);
    CAROM_VERIFY(saturation > 0.0);
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
    d_tol = tolerance;
    d_sat = saturation;
    d_subset_size = subset_size;
    d_convergence_subset_size = convergence_subset_size;
    d_output_log_path = output_log_path;
    d_use_centroid = use_centroid;
    d_max_error = 0;
    d_next_point_to_sample = -1;
    d_next_point_requiring_residual = -1;
    d_point_requiring_residual_computed = false;
    d_iteration_started = false;
    d_convergence_started = false;
    d_subset_created = false;
    d_procedure_completed = false;
    d_subset_counter = 0;
    d_counter = -1;
    d_debug_algorithm = debug_algorithm;

    rng.seed(random_seed);

    if (d_rank == 0)
    {
        if (d_output_log_path == "")
        {
            std::cout << "Greedy tolerance: " << d_tol << std::endl;
            std::cout << "Greedy saturation constant: " << d_sat << std::endl;
            std::cout << "Greedy iteration subset size: " << d_subset_size << std::endl;
            std::cout << "Greedy convergence subset size: " << d_convergence_subset_size << std::endl;
        }
        else
        {
            std::ofstream database_history;
            database_history.open(d_output_log_path, std::ios::app);
            database_history << "Greedy tolerance: " << d_tol << std::endl;
            database_history << "Greedy saturation constant: " << d_sat << std::endl;
            database_history << "Greedy iteration subset size: " << d_subset_size << std::endl;
            database_history << "Greedy convergence subset size: " << d_convergence_subset_size << std::endl;
            database_history.close();
        }
    }
}

void
GreedyParameterPointSelector::initializeParameterPoints(
    std::vector<Vector> parameter_points)
{
    CAROM_VERIFY(parameter_points.size() > 0);
    CAROM_VERIFY(d_subset_size <= parameter_points.size());
    CAROM_VERIFY(d_convergence_subset_size <= parameter_points.size());

    d_parameter_points = parameter_points;

    for (int i = 0; i < d_parameter_points.size() - 1; i++) {
        CAROM_VERIFY(d_parameter_points[i].dim() == d_parameter_points[i + 1].dim());
    }

    for (int i = 0 ; i < d_parameter_points.size(); i++) {
        d_parameter_point_errors.push_back(INT_MAX);
        d_parameter_point_local_rom.push_back(-1);
        d_parameter_point_random_indices.push_back(i);
    }

    d_min_param_point = d_parameter_points.front();
    d_max_param_point = d_parameter_points.back();

    for (int i = 0; i < d_parameter_points.size(); i++)
    {
        for (int j = 0; j < d_parameter_points[i].dim(); j++)
        {
            d_min_param_point.item(j) = std::min(d_min_param_point.item(j), d_parameter_points[i].item(j));
            d_max_param_point.item(j) = std::max(d_max_param_point.item(j), d_parameter_points[i].item(j));
        }
    }

    generateConvergenceSubset();
}

std::shared_ptr<Vector>
GreedyParameterPointSelector::getNextParameterPoint()
{
    if (isComplete())
    {
        return std::shared_ptr<Vector>(nullptr);
    }
    if (d_iteration_started)
    {
        return std::shared_ptr<Vector>(nullptr);
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

    d_iteration_started = true;
    d_max_error = 0;

    int curr_point_to_sample = d_next_point_to_sample;
    d_next_point_to_sample = -1;

    d_parameter_point_errors[curr_point_to_sample] = 0;
    d_parameter_point_local_rom[curr_point_to_sample] = curr_point_to_sample;

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

    if (d_check_local_rom)
    {
        d_next_point_requiring_residual = curr_point_to_sample;
        d_point_requiring_residual_computed = true;
    }
    else
    {
        // Precompute next residual point
        // This will allow us to figure out if the greedy algorithm has terminated
        // early without needing an extra call to the residual function.
        getNextPointRequiringResidual();
    }

    Vector* result = new Vector(d_parameter_points[curr_point_to_sample]);
    return std::shared_ptr<Vector>(result);
}

struct GreedyResidualPoint
GreedyParameterPointSelector::getNextPointRequiringResidual()
{
    if (isComplete())
    {
        return createGreedyResidualPoint(nullptr, nullptr);
    }
    if (!d_iteration_started)
    {
        return createGreedyResidualPoint(nullptr, nullptr);
    }

    if (d_convergence_started)
    {
        return getNextConvergencePointRequiringResidual();
    }
    else
    {
        return getNextSubsetPointRequiringResidual();
    }
}

struct GreedyResidualPoint
GreedyParameterPointSelector::getNextSubsetPointRequiringResidual()
{
    if (d_point_requiring_residual_computed)
    {
        Vector* result1 = new Vector(d_parameter_points[d_next_point_requiring_residual]);
        Vector* result2 = new Vector(d_parameter_points[getNearestROMIndex(d_next_point_requiring_residual)]);
        return createGreedyResidualPoint(result1, result2);
    }
    if (d_subset_counter == d_subset_size)
    {
        return createGreedyResidualPoint(nullptr, nullptr);
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

    d_next_point_requiring_residual = -1;

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
            double error_with_sat_factor = d_sat * d_parameter_point_errors[d_parameter_point_random_indices[d_counter]];
            if (error_with_sat_factor > d_max_error)
            {
                // if we have already computed this residual at the same local rom, the residual will not improve
                // no need to calculate the residual again
                if (d_parameter_point_local_rom[d_parameter_point_random_indices[d_counter]] == getNearestROMIndex(d_parameter_point_random_indices[d_counter]))
                {
                    d_max_error = error_with_sat_factor;
                    d_next_point_to_sample = d_parameter_point_random_indices[d_counter];
                    if (d_rank == 0)
                    {
                        if (d_output_log_path == "")
                        {
                            std::cout << "Residual at [ ";
                            for (int i = 0 ; i < d_parameter_points[d_parameter_point_random_indices[d_counter]].dim(); i++)
                            {
                                std::cout << d_parameter_points[d_parameter_point_random_indices[d_counter]].item(i) << " ";
                            }
                            std::cout << "] skipped." << std::endl;
                            std::cout << "Residual multiplied by saturation factor " << error_with_sat_factor << " already computed at the same local ROM." << std::endl;
                        }
                        else
                        {
                            std::ofstream database_history;
                            database_history.open(d_output_log_path, std::ios::app);
                            database_history << "Residual at [ ";
                            for (int i = 0 ; i < d_parameter_points[d_parameter_point_random_indices[d_counter]].dim(); i++)
                            {
                                database_history << d_parameter_points[d_parameter_point_random_indices[d_counter]].item(i) << " ";
                            }
                            database_history << "] skipped." << std::endl;
                            database_history << "Residual multiplied by saturation factor " << error_with_sat_factor << " already computed at the same local ROM." << std::endl;
                            database_history.close();
                        }
                    }
                }
                else
                {
                    d_next_point_requiring_residual = d_parameter_point_random_indices[d_counter];
                    d_point_requiring_residual_computed = true;
                    Vector* result1 = new Vector(d_parameter_points[d_next_point_requiring_residual]);
                    Vector* result2 = new Vector(d_parameter_points[getNearestROMIndex(d_next_point_requiring_residual)]);
                    return createGreedyResidualPoint(result1, result2);
                }
            }
            else
            {
                if (d_rank == 0)
                {
                    if (d_output_log_path == "")
                    {
                        std::cout << "Residual at [ ";
                        for (int i = 0 ; i < d_parameter_points[d_parameter_point_random_indices[d_counter]].dim(); i++)
                        {
                            std::cout << d_parameter_points[d_parameter_point_random_indices[d_counter]].item(i) << " ";
                        }
                        std::cout << "] skipped." << std::endl;
                        std::cout << "Residual multiplied by saturation factor " << error_with_sat_factor
                                  << " is less than current max error " << d_max_error << std::endl;
                    }
                    else
                    {
                        std::ofstream database_history;
                        database_history.open(d_output_log_path, std::ios::app);
                        database_history << "Residual at [ ";
                        for (int i = 0 ; i < d_parameter_points[d_parameter_point_random_indices[d_counter]].dim(); i++)
                        {
                            database_history << d_parameter_points[d_parameter_point_random_indices[d_counter]].item(i) << " ";
                        }
                        database_history << "] skipped." << std::endl;
                        database_history << "Residual multiplied by saturation factor " << error_with_sat_factor
                                         << " is less than current max error " << d_max_error << std::endl;
                        database_history.close();
                    }
                }
            }
        }
    }
    if (d_next_point_requiring_residual == -1)
    {
        if (d_rank == 0)
        {
            if (d_output_log_path == "")
            {
                std::cout << "Ran out of points to calculate residual for in this iteration." << std::endl;
            }
            else
            {
                std::ofstream database_history;
                database_history.open(d_output_log_path, std::ios::app);
                database_history << "Ran out of points to calculate residual for in this iteration." << std::endl;
                database_history.close();
            }
        }
        if (d_max_error < d_tol)
        {
            startConvergence();
        }
        else
        {
            d_iteration_started = false;
        }
    }

    return createGreedyResidualPoint(nullptr, nullptr);
}

struct GreedyResidualPoint
GreedyParameterPointSelector::getNextConvergencePointRequiringResidual()
{
    if (d_point_requiring_residual_computed)
    {
        Vector* result1 = new Vector(d_convergence_points[d_next_point_requiring_residual]);
        std::shared_ptr<Vector> result2 = getNearestROM(d_convergence_points[d_next_point_requiring_residual]);
        return createGreedyResidualPoint(result1, result2);
    }
    if (d_counter == d_convergence_subset_size)
    {
        return createGreedyResidualPoint(nullptr, nullptr);
    }

    d_next_point_requiring_residual = -1;

    //get next point requiring residual
    while (d_counter < (int) d_convergence_points.size())
    {
        d_next_point_requiring_residual = d_counter;
        d_point_requiring_residual_computed = true;
        Vector* result1 = new Vector(d_convergence_points[d_next_point_requiring_residual]);
        std::shared_ptr<Vector> result2 = getNearestROM(d_convergence_points[d_next_point_requiring_residual]);
        return createGreedyResidualPoint(result1, result2);
    }

    if (d_next_point_requiring_residual == -1)
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

    return createGreedyResidualPoint(nullptr, nullptr);
}

void
GreedyParameterPointSelector::setPointResidual(double error, int vec_size)
{
    CAROM_VERIFY(error >= 0);
    CAROM_VERIFY(d_point_requiring_residual_computed);

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
        setConvergenceResidual(proc_errors);
    }
    else
    {
        setSubsetResidual(proc_errors);
    }
}

void
GreedyParameterPointSelector::printResidual(Vector residualPoint, double proc_errors)
{
    if (d_rank == 0)
    {
        if (d_output_log_path == "")
        {
            std::cout << "Residual computed at [ ";
            for (int i = 0 ; i < residualPoint.dim(); i++)
            {
                std::cout << residualPoint.item(i) << " ";
            }
            std::cout << "]" << std::endl;
            std::cout << "Residual: " << proc_errors << std::endl;
        }
        else
        {
            std::ofstream database_history;
            database_history.open(d_output_log_path, std::ios::app);
            database_history << "Residual computed at [ ";
            for (int i = 0 ; i < residualPoint.dim(); i++)
            {
                database_history << residualPoint.item(i) << " ";
            }
            database_history << "]" << std::endl;
            database_history << "Residual: " << proc_errors << std::endl;
            database_history.close();
        }
    }
}

void
GreedyParameterPointSelector::printToleranceNotMet()
{
    if (d_rank == 0)
    {
        if (d_output_log_path == "")
        {
            std::cout << "Tolerance " << d_tol << " not met." << std::endl;
        }
        else
        {
            std::ofstream database_history;
            database_history.open(d_output_log_path, std::ios::app);
            database_history << "Tolerance " << d_tol << " not met." << std::endl;
            database_history.close();
        }
    }
}

void
GreedyParameterPointSelector::setSubsetResidual(double proc_errors)
{
    if (d_check_local_rom)
    {
        auto search = d_parameter_sampled_indices.find(d_next_point_requiring_residual);
        if (search != d_parameter_sampled_indices.end())
        {
            d_parameter_point_errors[d_next_point_requiring_residual] = proc_errors;
            d_parameter_point_local_rom[d_next_point_requiring_residual] = d_next_point_requiring_residual;

            if (d_rank == 0)
            {
                if (d_output_log_path == "")
                {
                    std::cout << "Local ROM Residual computed at [ ";
                    for (int i = 0 ; i < d_parameter_points[d_next_point_requiring_residual].dim(); i++)
                    {
                        std::cout << d_parameter_points[d_next_point_requiring_residual].item(i) << " ";
                    }
                    std::cout << "]" << std::endl;
                    std::cout << "Local ROM Residual (tolerance unchecked): " << proc_errors << std::endl;
                }
                else
                {
                    std::ofstream database_history;
                    database_history.open(d_output_log_path, std::ios::app);
                    database_history << "Local ROM Residual computed at [ ";
                    for (int i = 0 ; i < d_parameter_points[d_next_point_requiring_residual].dim(); i++)
                    {
                        database_history << d_parameter_points[d_next_point_requiring_residual].item(i) << " ";
                    }
                    database_history << "]" << std::endl;
                    database_history << "Local ROM Residual (tolerance unchecked): " << proc_errors << std::endl;
                    database_history.close();
                }
            }

            d_point_requiring_residual_computed = false;

            // Precompute next residual point
            // This will allow us to figure out if the greedy algorithm has terminated
            // early without needing an extra call to the residual function.
            getNextPointRequiringResidual();

            return;
        }
    }

    if (proc_errors < d_parameter_point_errors[d_parameter_point_random_indices[d_counter]])
    {
        d_parameter_point_errors[d_parameter_point_random_indices[d_counter]] = proc_errors;
        d_parameter_point_local_rom[d_parameter_point_random_indices[d_counter]] = getNearestROMIndex(d_parameter_point_random_indices[d_counter]);
    }

    printResidual(d_parameter_points[d_parameter_point_random_indices[d_counter]], proc_errors);

    d_point_requiring_residual_computed = false;

    if (proc_errors > d_max_error)
    {
        d_max_error = proc_errors;
        d_next_point_to_sample = d_parameter_point_random_indices[d_counter];
    }

    if (d_subset_counter == d_subset_size || d_counter == (int) d_parameter_points.size() - 1)
    {
        if (d_max_error < d_tol)
        {
            startConvergence();
        }
        else
        {
            d_iteration_started = false;
            printToleranceNotMet();
        }
    }

    // Precompute next residual point
    // This will allow us to figure out if the greedy algorithm has terminated
    // early without needing an extra call to the residual function.
    getNextPointRequiringResidual();

    return;
}

void
GreedyParameterPointSelector::setConvergenceResidual(double proc_errors)
{
    printResidual(d_convergence_points[d_counter], proc_errors);

    d_point_requiring_residual_computed = false;

    if (proc_errors >= d_tol)
    {
        d_iteration_started = false;
        double curr_max_error = 0.0;
        for (int i = 0; i < d_parameter_points.size(); i++)
        {
            auto search = d_parameter_sampled_indices.find(i);
            if (search == d_parameter_sampled_indices.end())
            {
                if (d_parameter_point_errors[i] > curr_max_error)
                {
                    curr_max_error = d_parameter_point_errors[i];
                    d_next_point_to_sample = i;
                }
            }
        }
        printToleranceNotMet();
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

    // Precompute next residual point
    // This will allow us to figure out if the greedy algorithm has terminated
    // early without needing an extra call to the residual function.
    getNextPointRequiringResidual();

    return;
}

void
GreedyParameterPointSelector::generateConvergenceSubset()
{
    d_convergence_points.clear();

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
        d_convergence_points.push_back(point);
    }
}

void
GreedyParameterPointSelector::startConvergence()
{
    d_convergence_started = true;
    d_max_error = 0;
    d_counter = 0;
    d_subset_counter = 0;

    if (d_rank == 0)
    {
        if (d_output_log_path == "")
        {
            std::cout << "Tolerance " << d_tol << " met. Computing convergence." << std::endl;
        }
        else
        {
            std::ofstream database_history;
            database_history.open(d_output_log_path, std::ios::app);
            database_history << "Tolerance " << d_tol << " met. Computing convergence." << std::endl;
            database_history.close();
        }
    }

    // Precompute next residual point
    // This will allow us to figure out if the greedy algorithm has terminated
    // early without needing an extra call to the residual function.
    getNextPointRequiringResidual();
}

std::shared_ptr<Vector>
GreedyParameterPointSelector::getNearestROM(Vector point)
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
GreedyParameterPointSelector::getNearestROMIndex(int index)
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
        else
        {
            closest_dist_to_points = 0;
            closest_point_index = *itr;
            break;
        }
    }

    return closest_point_index;
}

std::vector<Vector>
GreedyParameterPointSelector::getParameterPointDomain()
{
    return d_parameter_points;
}

std::vector<Vector>
GreedyParameterPointSelector::getSampledParameterPoints()
{
    std::vector<Vector> sampled_points;
    for (auto itr = d_parameter_sampled_indices.begin(); itr != d_parameter_sampled_indices.end(); ++itr) {
        sampled_points.push_back(d_parameter_points[*itr]);
    }
    return sampled_points;
}

void
GreedyParameterPointSelector::save(std::string base_file_name)
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
        sprintf(tmp, "tol");
        database.putDouble(tmp, d_tol);
        sprintf(tmp, "sat");
        database.putDouble(tmp, d_sat);
        sprintf(tmp, "subset_size");
        database.putInteger(tmp, d_subset_size);
        sprintf(tmp, "convergence_subset_size");
        database.putInteger(tmp, d_convergence_subset_size);
        sprintf(tmp, "next_point_to_sample");
        database.putInteger(tmp, d_next_point_to_sample);
        sprintf(tmp, "next_point_requiring_residual");
        database.putInteger(tmp, d_next_point_requiring_residual);
        sprintf(tmp, "check_local_rom");
        database.putInteger(tmp, d_check_local_rom);
        sprintf(tmp, "use_centroid");
        database.putInteger(tmp, d_use_centroid);
        sprintf(tmp, "iteration_started");
        database.putInteger(tmp, d_iteration_started);
        sprintf(tmp, "convergence_started");
        database.putInteger(tmp, d_convergence_started);
        sprintf(tmp, "point_requiring_residual_computed");
        database.putInteger(tmp, d_point_requiring_residual_computed);
        sprintf(tmp, "subset_created");
        database.putInteger(tmp, d_subset_created);
        sprintf(tmp, "debug_algorithm");
        database.putInteger(tmp, d_debug_algorithm);
        sprintf(tmp, "counter");
        database.putInteger(tmp, d_counter);
        sprintf(tmp, "subset_counter");
        database.putInteger(tmp, d_subset_counter);

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
GreedyParameterPointSelector::isComplete()
{
    if (d_parameter_sampled_indices.size() == d_parameter_points.size())
    {
        d_procedure_completed = true;
    }
    return d_procedure_completed;
}

GreedyParameterPointSelector::~GreedyParameterPointSelector()
{
}

}
