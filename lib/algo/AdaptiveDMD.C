/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
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
#include <algorithm>

namespace CAROM {

AdaptiveDMD::AdaptiveDMD(int dim, double desired_dt, std::string rbf, std::string interp_method,
                         double closest_rbf_val) : DMD(dim, desired_dt)
{
    CAROM_VERIFY(rbf == "G" || rbf == "IQ" || rbf == "IMQ");
    CAROM_VERIFY(interp_method == "LS" || interp_method == "IDW" || interp_method == "LP");
    CAROM_VERIFY(closest_rbf_val >= 0.0 && closest_rbf_val <= 1.0);
    d_dt = desired_dt;
    d_interp_method = interp_method;
    d_rbf = rbf;
    d_closest_rbf_val = closest_rbf_val;
}

void AdaptiveDMD::takeSample(double* u_in, double t)
{
    CAROM_VERIFY(u_in != 0);
    CAROM_VERIFY(t >= 0.0);
    Vector* sample = new Vector(u_in, d_dim, true);

    double orig_t = t;
    if (d_snapshots.empty())
    {
        d_t_offset = t;
        t = 0.0;
    }
    else
    {
        t -= d_t_offset;
    }

    // Erase any snapshots taken at the same or later time
    while (!d_sampled_times.empty() && d_sampled_times.back()->item(0) >= t)
    {
        if (d_rank == 0) std::cout << "Removing existing snapshot at time: " << d_t_offset + d_sampled_times.back()->item(0) << std::endl;
        Vector* last_snapshot = d_snapshots.back();
        delete last_snapshot;
        d_snapshots.pop_back();
        d_sampled_times.pop_back();
    }

    if (d_snapshots.empty())
    {
        d_t_offset = orig_t;
        t = 0.0;
    }
    else
    {
        CAROM_VERIFY(d_sampled_times.back()->item(0) < t);
    }

    d_snapshots.push_back(sample);
    Vector* sampled_time = new Vector(&t, 1, false);
    d_sampled_times.push_back(sampled_time);
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

void AdaptiveDMD::interpolateSnapshots()
{
    CAROM_VERIFY(d_interp_snapshots.size() == 0);
    CAROM_VERIFY(d_snapshots.size() == d_sampled_times.size());
    CAROM_VERIFY(d_sampled_times.size() > 1);

    if (d_rank == 0) std::cout << "Number of snapshots is: " << d_snapshots.size() << std::endl;

    bool automate_dt = false;
    if (d_dt <= 0.0)
    {
        automate_dt = true;
        std::vector<double> d_sampled_dts;
        for (int i = 1; i < d_sampled_times.size(); i++)
        {
            d_sampled_dts.push_back(d_sampled_times[i]->item(0) - d_sampled_times[i - 1]->item(0));
        }

        auto m = d_sampled_dts.begin() + d_sampled_dts.size() / 2;
        std::nth_element(d_sampled_dts.begin(), m, d_sampled_dts.end());
        if (d_rank == 0) std::cout << "Setting desired dt to the median dt between the samples: " << d_sampled_dts[d_sampled_dts.size() / 2] << std::endl;
        d_dt = d_sampled_dts[d_sampled_dts.size() / 2];
    }
    CAROM_VERIFY(d_sampled_times.back()->item(0) > d_dt);

    // Find the nearest dt that evenly divides the snapshots.
    int num_time_steps = std::round(d_sampled_times.back()->item(0) / d_dt);
    if (automate_dt && num_time_steps < d_sampled_times.size())
    {
        num_time_steps = d_sampled_times.size();
        if (d_rank == 0) std::cout << "There will be less interpolated snapshots than FOM snapshots. dt will be decreased." << std::endl;
    }
    double new_dt = d_sampled_times.back()->item(0) / num_time_steps;
    if (new_dt != d_dt)
    {
        d_dt = new_dt;
        if (d_rank == 0) std::cout << "Setting desired dt to " << d_dt << " to ensure a constant dt given the final sampled time." << std::endl;
    }

    // Solve the linear system if required.
    Matrix* f_T = NULL;
    double epsilon = convertClosestRBFToEpsilon(d_sampled_times, d_rbf, d_closest_rbf_val);
    if (d_interp_method == "LS")
    {
        f_T = solveLinearSystem(d_sampled_times, d_snapshots, d_interp_method, d_rbf, epsilon);
    }

    // Create interpolated snapshots using d_dt as the desired dt.
    for (int i = 0; i <= num_time_steps; i++)
    {
        double curr_time = i * d_dt;
        if (d_rank == 0) std::cout << "Creating new interpolated sample at: " << d_t_offset + curr_time << std::endl;
        CAROM::Vector* point = new Vector(&curr_time, 1, false);

        // Obtain distances from database points to new point
        std::vector<double> rbf = obtainRBFToTrainingPoints(d_sampled_times, d_interp_method, d_rbf, epsilon, point);

        // Obtain the interpolated snapshot.
        CAROM::Vector* curr_interpolated_snapshot = obtainInterpolatedVector(d_snapshots, f_T, d_interp_method, rbf);
        d_interp_snapshots.push_back(curr_interpolated_snapshot);

        delete point;
    }

    if (d_rank == 0) std::cout << "Number of interpolated snapshots is: " << d_interp_snapshots.size() << std::endl;
}

double AdaptiveDMD::getTrueDt() const
{
    return d_dt;
}

const Matrix* AdaptiveDMD::getInterpolatedSnapshots()
{
    if (d_interp_snapshots.size() == 0) interpolateSnapshots();
    return createSnapshotMatrix(d_interp_snapshots);
}

}
