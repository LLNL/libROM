/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the AdaptiveDMD algorithm.

#include "AdaptiveDMD.h"
#include "manifold_interp/VectorInterpolator.h"

#include "linalg/Matrix.h"
#include "linalg/Vector.h"

namespace CAROM {

AdaptiveDMD::AdaptiveDMD(int dim, double desired_dt, std::string rbf, std::string interp_method, double epsilon) : DMD(dim, desired_dt)
{
    CAROM_VERIFY(desired_dt > 0.0);
    CAROM_VERIFY(rbf == "G" || rbf == "IQ" || rbf == "MQ" || rbf == "IMQ");
    CAROM_VERIFY(interp_method == "LS" || interp_method == "IDW" || interp_method == "LP");
    d_dt = desired_dt;
    d_interp_method = interp_method;
    d_rbf = rbf;
    d_epsilon = epsilon;
}

void AdaptiveDMD::takeSample(double* u_in, double t)
{
    CAROM_VERIFY(u_in != 0);
    CAROM_VERIFY(t >= 0.0);
    Vector* sample = new Vector(u_in, d_dim, true);
    if (d_snapshots.empty())
    {
        d_t_offset = t;
        t = 0.0;
    }
    else
    {
        t -= d_t_offset;
    }

    // If we have sampled another snapshot at the same timestep, replace
    // the previous sample with the new one.
    if (!d_sampled_times.empty() && d_sampled_times.back()->item(0) == t)
    {
        Vector* last_snapshot = d_snapshots.back();
        delete last_snapshot;
        d_snapshots.pop_back();
        d_snapshots.push_back(sample);
    }
    else
    {
        d_snapshots.push_back(sample);
        Vector* sampled_time = new Vector(&t, 1, false);
        d_sampled_times.push_back(sampled_time);
    }
}

void AdaptiveDMD::train(double energy_fraction)
{
    const Matrix* f_snapshots = getInterpolatedSnapshots();
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(energy_fraction > 0 && energy_fraction <= 1);
    d_energy_fraction = energy_fraction;
    constructDMD(f_snapshots, d_rank, d_num_procs);

    delete f_snapshots;
}

void AdaptiveDMD::train(int k)
{
    const Matrix* f_snapshots = getInterpolatedSnapshots();
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(k > 0 && k <= f_snapshots->numColumns() - 1);
    d_energy_fraction = -1.0;
    d_k = k;
    constructDMD(f_snapshots, d_rank, d_num_procs);

    delete f_snapshots;
}

Vector* AdaptiveDMD::predict(double t)
{
    const std::pair<Vector*, Vector*> d_projected_init_pair(d_projected_init_real, d_projected_init_imaginary);
    return predict(d_projected_init_pair, t);
}

Vector* AdaptiveDMD::predict(const std::pair<Vector*, Vector*> init, double t)
{
    t -= d_t_offset;
    return DMD::predict(init, t);
}

void AdaptiveDMD::interpolateSnapshots()
{
    CAROM_VERIFY(d_sampled_times.back()->item(0) > d_dt);
    CAROM_VERIFY(d_interp_snapshots.size() == 0);

    // Find the nearest dt that evenly divides the snapshots.
    int num_time_steps = std::round(d_sampled_times.back()->item(0) / d_dt);
    double new_dt = d_sampled_times.back()->item(0) / num_time_steps;
    if (new_dt != d_dt)
    {
        d_dt = new_dt;
        std::cout << "Setting desired dt to " << d_dt << " to ensure a constant dt given the final sampled time." << std::endl;
    }

    // Solve the linear system if required.
    Matrix* f_T = NULL;
    if (d_interp_method == "LS")
    {
        if (d_epsilon <= 0.0) d_epsilon = 0.5 / d_dt;
        f_T = solveLinearSystem(d_sampled_times, d_snapshots, d_interp_method, d_rbf, d_epsilon);
        std::cout << "Epsilon auto-corrected by the linear solve to " << d_epsilon << std::endl;
    }

    // Create interpolated snapshots using d_dt as the desired dt.
    for (int i = 0; i <= num_time_steps; i++)
    {
        double curr_time = i * d_dt;
        std::cout << "Creating new interpolated sample at: " << curr_time << std::endl;
        CAROM::Vector* point = new Vector(&curr_time, 1, false);

        // Obtain distances from database points to new point
        std::vector<double> rbf = obtainRBFToTrainingPoints(d_sampled_times, d_interp_method, d_rbf, d_epsilon, point);

        // Obtain the interpolated snapshot.
        CAROM::Vector* curr_interpolated_snapshot = obtainInterpolatedVector(d_snapshots, f_T, d_interp_method, rbf);
        d_interp_snapshots.push_back(curr_interpolated_snapshot);

        delete point;
    }
}

double AdaptiveDMD::getTruedt()
{
    return d_dt;
}

const Matrix* AdaptiveDMD::getInterpolatedSnapshots()
{
    if (d_interp_snapshots.size() == 0) interpolateSnapshots();
    return createSnapshotMatrix(d_interp_snapshots);
}
}
