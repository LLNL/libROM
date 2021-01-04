#include "../BasisGenerator.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring> // for memcpy
#include <random>

#include "mpi.h"

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int totaldim = nprocs * 12;
    const double entry = 1 / sqrt(static_cast<double>(3 * nprocs));

    double* sample1 = new double[5] {0.5377, 1.8339, -2.2588, 0.8622, 0.3188};
    double* sample2 = new double[5] {-1.3077, -0.4336, 0.3426, 3.5784, 2.7694};
    double* sample3 = new double[5] {-1.3499, 3.0349, 0.7254, -0.0631, 0.7147};

    { /* Wrap the sampler in its own scope; ScaLAPACK will call MPI_Comm_finalize
         in the destructor after MPI_Finalize otherwise. */

    // Construct an SVDSampler to send our matrix. I take V = I for simplicity,
    // so the matrix A that we factor is just the columns scaled by the sigmas.
    CAROM::Options randomized_svd_options = CAROM::Options(5, 3, 1);
    randomized_svd_options.setRandomizedSVD(true);
    CAROM::BasisGenerator sampler(randomized_svd_options, false);
    sampler.takeSample(sample1, 0, 0);
    sampler.takeSample(sample2, 0, 0);
    sampler.takeSample(sample3, 0, 0);

    auto distU = sampler.getSpatialBasis();

    }

    MPI_Finalize();

    return 0;
}
