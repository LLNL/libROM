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
#include "mpi.h"
#include <cmath>
#include <limits.h>
#include <fstream>

namespace CAROM {

GreedyParameterPointSelector::GreedyParameterPointSelector(
    std::vector<Vector> parameter_points,
    double tolerance,
    double saturation,
    int subset_size,
    int convergence_subset_size,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm)
{

    constructObject(parameter_points, tolerance, saturation,
        subset_size, convergence_subset_size, use_centroid, random_seed, debug_algorithm);
}

GreedyParameterPointSelector::GreedyParameterPointSelector(
    std::vector<double> parameter_points,
    double tolerance,
    double saturation,
    int subset_size,
    int convergence_subset_size,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm)
{
    //convert parameter_points from double to Vector
    CAROM_VERIFY(parameter_points.size() > 0);

    std::vector<Vector> parameter_points_vec;
    for (int i = 0; i < parameter_points.size(); i++) {
        Vector vec(1, false);
        vec.item(0) = parameter_points[i];
        parameter_points_vec.push_back(vec);
    }

    constructObject(parameter_points_vec, tolerance, saturation,
        subset_size, convergence_subset_size, use_centroid, random_seed, debug_algorithm);
}

GreedyParameterPointSelector::GreedyParameterPointSelector(
    double param_space_min,
    double param_space_max,
    double param_space_size,
    double tolerance,
    double saturation,
    int subset_size,
    int convergence_subset_size,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm)
{
    //convert parameter_points from double to Vector
    CAROM_VERIFY(param_space_max > param_space_min);
    CAROM_VERIFY(param_space_size > 0);

    std::vector<Vector> parameter_points_vec;
    Vector vec(1, false);
    if (param_space_size == 1)
    {
        vec.item(0) = (param_space_min + param_space_max) / 2;
        parameter_points_vec.push_back(vec);
    }
    else if (param_space_size == 2)
    {
        double frequency = abs(param_space_max - param_space_min) / (param_space_size + 1);
        for (int i = 0; i < param_space_size; i++) {
            vec.item(0) = param_space_min + (i + 1) * frequency;
            parameter_points_vec.push_back(vec);
        }
    }
    else
    {
        double frequency = abs(param_space_max - param_space_min) / param_space_size;
        for (int i = 0; i < param_space_size; i++) {
            vec.item(0) = param_space_min + i * frequency;
            parameter_points_vec.push_back(vec);
        }
    }

    constructObject(parameter_points_vec, tolerance, saturation,
        subset_size, convergence_subset_size, use_centroid, random_seed, debug_algorithm);
}

void
GreedyParameterPointSelector::constructObject(
    std::vector<Vector> parameter_points,
    double tolerance,
    double saturation,
    int subset_size,
    int convergence_subset_size,
    bool use_centroid,
    int random_seed,
    bool debug_algorithm)
{
    CAROM_VERIFY(parameter_points.size() > 0);
    CAROM_VERIFY(tolerance > 0.0);
    CAROM_VERIFY(saturation > 0.0);
    CAROM_VERIFY(subset_size > 0 && subset_size <= parameter_points.size());
    CAROM_VERIFY(convergence_subset_size > 0 && convergence_subset_size <= parameter_points.size());
    CAROM_VERIFY(random_seed > 0);

    d_parameter_points = parameter_points;
    d_tol = tolerance;
    d_sat = saturation;
    d_subset_size = subset_size;
    d_convergence_subset_size = convergence_subset_size;
    d_use_centroid = use_centroid;
    d_max_error = 0;
    d_next_point_to_sample = -1;
    d_next_point_requiring_residual = -1;
    d_point_requiring_residual_computed = false;
    d_iteration_started = false;
    d_subset_created = false;
    d_procedure_completed = false;
    d_subset_counter = 0;
    d_counter = -1;
    d_debug_algorithm = debug_algorithm;

    for (int i = 0; i < d_parameter_points.size() - 1; i++) {
        CAROM_VERIFY(d_parameter_points[i].dim() == d_parameter_points[i + 1].dim());
    }

    for (int i = 0 ; i < d_parameter_points.size(); i++) {
        d_parameter_point_errors.push_back(INT_MAX);
        d_parameter_point_random_indices.push_back(i);
    }
    rng.seed(random_seed);
}

int
GreedyParameterPointSelector::getNextParameterPoint()
{
    if (d_parameter_sampled_indices.size() == d_parameter_points.size())
    {
        d_procedure_completed = true;
        return -1;
    }
    if (d_iteration_started || d_procedure_completed)
    {
        return -1;
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

            double min_dist_to_centroid = INT_MAX;
            for (int i = 0; i < d_parameter_points.size(); i++) {
                Vector diff;
                centroid.minus(d_parameter_points[i], diff);
                double dist_to_centroid = diff.norm();
                if (dist_to_centroid < min_dist_to_centroid)
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

    d_subset_created = false;
    d_subset_counter = 0;
    d_counter = -1;

    auto search = d_parameter_sampled_indices.find(curr_point_to_sample);
    CAROM_VERIFY(search == d_parameter_sampled_indices.end());

    d_parameter_sampled_indices.insert(curr_point_to_sample);

    return curr_point_to_sample;
}

int
GreedyParameterPointSelector::getNextPointRequiringResidual()
{
    if (d_parameter_sampled_indices.size() == d_parameter_points.size())
    {
        d_procedure_completed = true;
        return -1;
    }
    if (d_subset_counter == d_subset_size)
    {
        return -1;
    }
    if (!d_iteration_started || d_procedure_completed)
    {
        return -1;
    }
    if (d_point_requiring_residual_computed)
    {
        return -1;
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

    //get next point requiring residual
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
            if (d_sat * d_parameter_point_errors[d_parameter_point_random_indices[d_counter]] > d_max_error)
            {
                d_next_point_requiring_residual = d_parameter_point_random_indices[d_counter];
                d_point_requiring_residual_computed = true;
                return d_next_point_requiring_residual;
            }
        }
    }

    if (d_next_point_requiring_residual == -1)
    {
        d_iteration_started = false;
    }

    return d_next_point_requiring_residual;
}

void
GreedyParameterPointSelector::setPointResidual(double error, int rank, int num_procs)
{
    CAROM_VERIFY(error >= 0);
    CAROM_VERIFY(d_point_requiring_residual_computed);

    double *proc_errors = new double[num_procs];
    proc_errors[rank] = error;
    CAROM_VERIFY(MPI_Allgather(MPI_IN_PLACE,
            1,
            MPI_DOUBLE,
            proc_errors,
            1,
            MPI_DOUBLE,
            MPI_COMM_WORLD) == MPI_SUCCESS);

    double total_error = 0;
    for (int i = 0; i < num_procs; i++) {
        total_error += pow(proc_errors[i], 2);
    }
    total_error = sqrt(total_error);

    d_parameter_point_errors[d_parameter_point_random_indices[d_counter]] = total_error;
    if (total_error > d_max_error)
    {
        d_max_error = total_error;
        d_next_point_to_sample = d_parameter_point_random_indices[d_counter];
    }

    d_point_requiring_residual_computed = false;

    if (d_subset_counter == d_subset_size || d_counter == (int) d_parameter_points.size() - 1)
    {
        d_iteration_started = false;
        if (d_max_error < d_tol)
        {
            if (!d_debug_algorithm)
            {
              std::shuffle(d_parameter_point_random_indices.begin(), d_parameter_point_random_indices.end(), rng);
            }

            d_procedure_completed = true;
            for (int i = 0; i < d_convergence_subset_size; i++)
            {
                if (d_parameter_point_errors[d_parameter_point_random_indices[i]] >= d_tol)
                {
                    d_procedure_completed = false;
                }
            }
        }
    }
}

int
GreedyParameterPointSelector::getNearestROM(int index)
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
            if (dist < closest_dist_to_points)
            {
                closest_dist_to_points = dist;
                closest_point_index = *itr;
            }
        }
    }

    return closest_point_index;
}

std::vector<Vector>
GreedyParameterPointSelector::getParameterPointDomain()
{
    return d_parameter_points;
}

void
GreedyParameterPointSelector::writeSampledPoints(std::string const& path)
{
    std::ofstream out_path;
    out_path.open(path);
    out_path << d_parameter_sampled_indices.size() << std::endl;
    out_path.close();

    int counter = 0;
    for (auto itr = d_parameter_sampled_indices.begin(); itr != d_parameter_sampled_indices.end(); ++itr) {
        std::string vec_path = path + "_" + std::to_string(counter);
        d_parameter_points[*itr].write(vec_path);
        counter++;
    }
}

GreedyParameterPointSelector::~GreedyParameterPointSelector()
{
}

int getNearestPoint(std::vector<Vector> paramPoints, Vector point)
{

    double closest_dist_to_points = INT_MAX;
    int closest_point_index = -1;

    for (int i = 0; i < paramPoints.size(); i++)
    {
        Vector diff;
        point.minus(paramPoints[i], diff);
        double dist = diff.norm();
        if (dist < closest_dist_to_points)
        {
            closest_dist_to_points = dist;
            closest_point_index = i;
        }
    }

    return closest_point_index;
}

int getNearestPoint(std::vector<double> paramPoints, double point)
{

    double closest_dist_to_points = INT_MAX;
    int closest_point_index = -1;

    for (int i = 0; i < paramPoints.size(); i++)
    {
        double dist = std::abs(point - paramPoints[i]);
        if (dist < closest_dist_to_points)
        {
            closest_dist_to_points = dist;
            closest_point_index = i;
        }
    }

    return closest_point_index;
}

std::vector<Vector>
readSampledPoints(std::string const& path)
{

    std::ifstream in_path;
    in_path.open(path);
    int num_sampled_points = 0;
    in_path >> num_sampled_points;
    in_path.close();

    std::vector<Vector> sampled_points;
    for (int i = 0; i < num_sampled_points; i++)
    {
        std::string vec_path = path + "_" + std::to_string(i);
        Vector point;
        point.read(vec_path);
        sampled_points.push_back(point);
    }

    return sampled_points;
}

}
