/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Interface to the DMD algorithm.

#include "DMD.h"

#include "Matrix.h"
#include "Vector.h"
#include "scalapack_wrapper.h"
#include "mpi.h"

/* Use C++11 built-in shared pointers if available; else fallback to Boost. */
#if __cplusplus >= 201103L
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#endif

using namespace std;

namespace CAROM {

    DMD::DMD(const Matrix* f_snapshots,
             double energy_fraction,
             int rank,
             int num_procs)
    {
        CAROM_VERIFY(f_snapshots->numColumns() > 1);
        CAROM_VERIFY(energy_fraction > 0 && energy_fraction <= 1);

        Matrix* f_snapshots_minus = new Matrix(f_snapshots->numRows(),
            f_snapshots->numColumns() - 1, f_snapshots->distributed());
        Matrix* f_snapshots_plus = new Matrix(f_snapshots->numRows(),
            f_snapshots->numColumns() - 1, f_snapshots->distributed());

        for (int i = 0; i < f_snapshots->numRows(); i++)
        {
            f_snapshots_minus->item(i, 0) = f_snapshots->item(i, 0);
            for (int j = 1; j < f_snapshots->numColumns() - 1; j++)
            {
                f_snapshots_minus->item(i, j) = f_snapshots->item(i, j);
                f_snapshots_plus->item(i, j - 1) = f_snapshots->item(i, j);
            }
            f_snapshots_minus->item(i, f_snapshots->numColumns() - 2) =
                f_snapshots->item(i, f_snapshots->numColumns() - 1);
        }

        int *row_offset = new int[num_procs + 1];
        row_offset[num_procs] = f_snapshots_minus->numDistributedRows();
        row_offset[rank] = f_snapshots_minus->numRows();

        CAROM_VERIFY(MPI_Allgather(MPI_IN_PLACE,
                                   1,
                                   MPI_INT,
                                   row_offset,
                                   1,
                                   MPI_INT,
                                   MPI_COMM_WORLD) == MPI_SUCCESS);
        for (int i = num_procs - 1; i >= 0; i--) {
            row_offset[i] = row_offset[i + 1] - row_offset[i];
        }

        CAROM_VERIFY(row_offset[0] == 0);

        int d_blocksize = row_offset[num_procs] / num_procs;
        if (row_offset[num_procs] % num_procs != 0) d_blocksize += 1;

        SLPK_Matrix svd_input;
        initialize_matrix(&svd_input, f_snapshots_minus->numColumns(),
                          f_snapshots_minus->numDistributedRows(),
                          1, num_procs, d_blocksize, d_blocksize);
        scatter_block(&svd_input, 1, 1,
                      f_snapshots_minus->getData(),
                      f_snapshots_minus->numColumns(),
                      f_snapshots_minus->numDistributedRows(),0);

        std::unique_ptr<SVDManager> d_factorizer;

        // This block does the actual ScaLAPACK call to do the factorization.
        svd_init(d_factorizer.get(), &svd_input);

        d_factorizer->dov = 1;

        factorize(d_factorizer.get());
        free_matrix_data(&svd_input);

        delete f_snapshots_minus;
        delete f_snapshots_plus;
    }

}
