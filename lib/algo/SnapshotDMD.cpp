/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of the AdaptiveDMD algorithm.

#include "manifold_interp/SnapshotInterpolator.h"
#include "SnapshotDMD.h"
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "utils/CSVDatabase.h"
#include <algorithm>

namespace CAROM {

SnapshotDMD::~SnapshotDMD()
{
}

void SnapshotDMD::train(int k, const Matrix* W0, double linearity_tol)
{
    CAROM_VERIFY(d_snapshots.size() > 0);

    if(k >= d_snapshots.size())
    {
        interpolateToNSnapshots(k + 1);
    }

    const Matrix* f_snapshots = getSnapshotMatrix();
    CAROM_VERIFY(f_snapshots->numColumns() > 1);
    CAROM_VERIFY(k > 0 && k <= f_snapshots->numColumns() - 1);
    d_energy_fraction = -1.0;
    d_k = k;
    constructDMD(f_snapshots, d_rank, d_num_procs, W0, linearity_tol);

    delete f_snapshots;
}

void SnapshotDMD::train(double energy_fraction, const Matrix* W0,
                        double linearity_tol)
{
    DMD::train(energy_fraction,W0,linearity_tol);
}

void SnapshotDMD::interpolateToNSnapshots(int n)
{
    SnapshotInterpolator* interp = new SnapshotInterpolator();
    std::vector<Vector*> new_snapshots;
    std::vector<Vector*> new_times;

    interp->interpolate(d_sampled_times,d_snapshots,n,new_times,new_snapshots);
    d_snapshots = std::move(new_snapshots);
    d_sampled_times = std::move(new_times);
    d_dt = d_sampled_times[2]->getData()[0]-d_sampled_times[1]->getData()[0];

}

void SnapshotDMD::interpolateToNSnapshots(int n, int run, int window)
{
    SnapshotInterpolator* interp = new SnapshotInterpolator();
    std::vector<Vector*> new_snapshots;
    std::vector<Vector*> new_times;
    CSVDatabase* csv_db(new CSVDatabase);
    int dim = d_snapshots[0]->dim();
    std::vector<double> print_times;
    std::string debugPath = "snapshot_interpolation";

    if(!d_rank)
    {
        std::cout << "Interpolating DMD snapshots with dimension " << dim << std::endl;
        std::cout << "Input times: ";
        for(int i = 0; i < d_snapshots.size(); ++i)
        {
            csv_db->putDoubleArray(debugPath+"/snapshot_" + std::to_string(
                                       i) + "pre_interp_run" +
                                   std::to_string(run) + "_win" + std::to_string(window) +
                                   ".csv",d_snapshots[i]->getData(),dim);
            print_times.push_back(d_sampled_times[i]->getData()[0]);
            std::cout << print_times[i] << ", ";
        }
        std::cout << std::endl;
        csv_db->putDoubleArray("/times_pre_interp_run_" + std::to_string(
                                   run) + "_win" + std::to_string(window) + ".csv", &print_times[0],
                               d_snapshots.size());
        std::cout << "Window " << window << " start with " << d_snapshots.size() <<
                  " snapshots" <<  std::endl;

        std::cout << "Window " << window << " spans: [" << print_times[0] << "," <<
                  print_times[d_snapshots.size()-1] << " with dt = " << d_dt << std::endl;
    }

    interp->interpolate(d_sampled_times,d_snapshots,n,new_times,new_snapshots);
    d_snapshots = std::move(new_snapshots);
    d_sampled_times = std::move(new_times);

    d_dt = d_sampled_times[2]->getData()[0]-d_sampled_times[1]->getData()[0];


    if(!d_rank)
    {
        std::cout << "output times: ";
        for(int i = 0; i < d_snapshots.size(); ++i)
        {
            csv_db->putDoubleArray(debugPath+"/snapshot_" + std::to_string(
                                       i) + "post_interp_run" +
                                   std::to_string(run) + "_win" + std::to_string(window) +
                                   ".csv",d_snapshots[i]->getData(),dim);
            print_times.push_back(d_sampled_times[i]->getData()[0]);
            std::cout << d_sampled_times[i]->getData()[0] << ", ";
        }
        std::cout << std::endl;
        csv_db->putDoubleArray("/times_post_interp_run_" + std::to_string(
                                   run) + "_win" + std::to_string(window) + ".csv", &print_times[0],
                               d_snapshots.size());
        std::cout << "Window " << window << " ends with " << d_snapshots.size() <<
                  " snapshots" << std::endl;
        std::cout << "Window " << window << " spans: [" <<
                  d_sampled_times[0]->getData()[0] << "," <<
                  d_sampled_times[d_snapshots.size()-1]->getData()[0] << " with dt = " << d_dt <<
                  std::endl;
    }


}

}
