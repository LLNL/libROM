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

    // Columns are constructed to be orthogonal with unit norm.
    std::vector<std::vector<double>> columns;
    for (int i = 0; i < 4; ++i) {
        columns.emplace_back(12u);
        for (int j = i; j < 12; j += 4)
            columns[i][j] = entry;
    }

    std::vector<double> sigmas(4);

    if (rank == 0) {
        std::normal_distribution<double> dist;
        auto generator = std::default_random_engine(std::random_device()());
        for (auto& x: sigmas) {
            x = dist(generator);
            x = x < 0 ? -x : x;
        }
        std::sort(sigmas.begin(), sigmas.end(), std::greater<double>());
    }
    MPI_Bcast(sigmas.data(), 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Generated singular values: ");
        for (int i = 0; i < 4; ++i)
            printf("%8.4E  ", sigmas[i]);
        printf("\n");
    }

    { /* Wrap the sampler in its own scope; ScaLAPACK will call MPI_Comm_finalize
         in the destructor after MPI_Finalize otherwise. */

    // Construct an SVDSampler to send our matrix. I take V = I for simplicity,
    // so the matrix A that we factor is just the columns scaled by the sigmas.
    CAROM::Options static_svd_options = CAROM::Options(12, 4, 4).setDebugMode(true);
    CAROM::BasisGenerator sampler(static_svd_options, false);

    for (unsigned j = 0; j < 4; ++j) {
        std::vector<double> similar(columns[j]);
        for (unsigned i = 0; i < 12; ++i)
            similar[i] *= sigmas[j];
        sampler.takeSample(similar.data(), 0);
    }

    auto distU = sampler.getSpatialBasis();
    if (rank == 0) {
        CAROM::Matrix U(totaldim, 4, false);
        memcpy(&U.item(0, 0), &distU->item(0, 0), 48*sizeof(double));
        for (int src = 1; src < nprocs; ++src)
            MPI_Recv(&U.item(12*src, 0), 48, MPI_DOUBLE, src, MPI_ANY_TAG,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf("Actual U                                        Computed U\n");
        printf("========                                        ==========\n");
        for (int i = 0; i < totaldim; ++i) {
            for (int j = 0; j < 4; ++j)
                printf("%8.4E ", columns[j][i%12]);
            printf("    ");
            for (int j = 0; j < 4; ++j)
                printf("%8.4E ", U.item(i, j));
            printf("\n");
        }
    } else {
        MPI_Send(&distU->item(0, 0), 48, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    }

    MPI_Finalize();

    return 0;
}
