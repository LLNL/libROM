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

#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "linalg/scalapack_wrapper.h"
#include "mpi.h"

/* Use C++11 built-in shared pointers if available; else fallback to Boost. */
#if __cplusplus >= 201103L
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#endif

/* Use automatically detected Fortran name-mangling scheme */
#define zgetrf CAROM_FC_GLOBAL(zgetrf, ZGETRF)
#define zgetri CAROM_FC_GLOBAL(zgetri, ZGETRI)

extern "C" {
    // LU decomposition of a general matrix.
    void zgetrf(int*, int*, double*, int*, int*, int*);

    // Generate inverse of a matrix given its LU decomposition.
    void zgetri(int*, double*, int*, int*, double*, int*, int*);
}

using namespace std;

namespace CAROM {

AdaptiveDMD::AdaptiveDMD(int dim) : DMD(dim)
{
}

void AdaptiveDMD::takeSample(double* u_in, double t)
{
    CAROM_VERIFY(u_in != 0);
    Vector sample(u_in, d_dim, true);
    if (d_snapshots.empty())
    {
        CAROM_VERIFY(t == 0.0);
    }

    // If we have sampled another snapshot at the same timestep, replace
    // the previous sample with the new one.
    if (!d_sampled_times.empty() && d_sampled_times.back() == t)
    {
        d_snapshots.pop_back();
        d_snapshots.push_back(sample);
    }
    else
    {
        d_snapshots.push_back(sample);
        d_sampled_times.push_back(t);
    }
}

double
AdaptiveDMD::interpolateSampledTime(double n)
{
    for (int i = 0; i < d_sampled_times.size(); i++)
    {
        if (n == d_sampled_times[i])
        {
            return (double) i;
        }
        else if (n < d_sampled_times[i])
        {
            double range = d_sampled_times[i] - d_sampled_times[i - 1];
            double offset = (n - d_sampled_times[i - 1]) / range;
            return (double) ((i - 1) + offset);
        }
    }

    // We are unable to handle this case and will error out.
    std::cout << "The inputted time was greater than the final sampled time. Interpolation can not occur. Aborting." << std::endl;
    abort();
}

Vector*
AdaptiveDMD::predict(const std::pair<Vector*, Vector*> init,
                     double n)
{
    CAROM_VERIFY(d_trained);
    CAROM_VERIFY(n >= 0.0);
    CAROM_VERIFY(n <= d_sampled_times.back());
    double interpolated_sampled_time = interpolateSampledTime(n);
    DMD::predict(init, interpolated_sampled_time);
}

}
