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

#include "ParameterPointGreedySelector.h"
#include "Vector.h"
#include <random>
#include <limits.h>

namespace CAROM {

ParameterPointGreedySelector::ParameterPointGreedySelector(
    std::vector<Vector> parameter_points,
    double tolerance,
    double saturation,
    int subset_size,
    int convergence_subset_size,
    int random_seed) :
   d_parameter_points(parameter_points),
   d_tol(tolerance),
   d_sat(saturation),
   d_subset_size(subset_size),
   d_convergence_subset_size(convergence_subset_size),
   d_max_error(0.0),
   d_iteration_completed(false),
   d_procedure_completed(false)
{
    CAROM_VERIFY(parameter_points.size() > 0);
    CAROM_VERIFY(tolerance > 0.0);
    CAROM_VERIFY(saturation > 0.0);
    CAROM_VERIFY(subset_size > 0 && subset_size <= parameter_points.size());
    CAROM_VERIFY(convergence_subset_size > 0 && convergence_subset_size <= parameter_points.size());
    CAROM_VERIFY(random_seed > 0);

    for (int i = 0 ; i < d_parameter_points.size(); i++) {
        d_parameter_point_errors.push_back(INT_MAX);
        d_parameter_point_random_indices.push_back(i);
    }
    rng.seed(random_seed);
}

ParameterPointGreedySelector::ParameterPointGreedySelector(
    std::vector<double> parameter_points,
    double tolerance,
    double saturation,
    int subset_size,
    int convergence_subset_size,
    int random_seed)
{
    //convert parameter_points from double to Vector

    std::vector<Vector> parameter_points_vec;

    ParameterPointGreedySelector(parameter_points_vec, tolerance, saturation,
        subset_size, convergence_subset_size, random_seed);
}

int
ParameterPointGreedySelector::computeNextSampleParameterPoint()
{
    if (!d_iteration_completed)
    {

        // generate random shuffle
        std::shuffle(d_parameter_point_random_indices.begin(), d_parameter_point_random_indices.end(), rng);

        if (d_rom_database.size() == 0)
        {
            //get center-most point
            d_iteration_completed = true;
        }
        else
        {
            //get error of next point
            for (int i = 0; i < d_subset_size; i++) {
                if (d_tol * d_parameter_point_errors[d_parameter_point_random_indices[i]] > d_max_error)
                {
                    // compute error here
                    double error;
                    d_parameter_point_errors[d_parameter_point_random_indices[i]] = error;
                    if (error > d_max_error)
                    {
                        d_max_error = error;
                        d_next_point = d_parameter_point_random_indices[i];
                    }
                }
            }
        }

        d_iteration_completed = true;
        if (d_max_error < d_tol)
        {
            std::shuffle(d_parameter_point_random_indices.begin(), d_parameter_point_random_indices.end(), rng);
            for (int i = 0; i < d_convergence_subset_size; i++)
            {
                if (d_tol * d_parameter_point_errors[d_parameter_point_random_indices[i]] > d_max_error)
                {
                    return d_next_point;
                }
            }
            d_next_point = -1;
            d_procedure_completed = true;
        }

        return d_next_point;
    }

    return d_next_point;
}

void
ParameterPointGreedySelector::addROMToDatabase(BasisGenerator* rom)
{
    CAROM_VERIFY(d_iteration_completed);
    CAROM_VERIFY(!d_procedure_completed);
    if (d_iteration_completed)
    {
        d_rom_database.push_back(rom);
    }
    d_iteration_completed = false;
}

ParameterPointGreedySelector::~ParameterPointGreedySelector()
{
}

}
