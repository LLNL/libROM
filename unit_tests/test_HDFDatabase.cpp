/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

#ifdef CAROM_HAS_GTEST

#include<gtest/gtest.h>
#include "linalg/Matrix.h"
#include "linalg/Vector.h"
#include "utils/HDFDatabase.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring> // for memcpy
#include <random>
#include "mpi.h"
#include "utils/mpi_utils.h"
#include <chrono>
using namespace std::chrono;

/**
 * Simple smoke test to make sure Google Test is properly linked
 */
TEST(GoogleTestFramework, GoogleTestFrameworkFound) {
    SUCCEED();
}

TEST(HDF5, Test_parallel_writing)
{
    int nproc, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    hid_t plist_id;
    hid_t file_id;
    herr_t errf = 0;
    
    /*
     * Set up file access property list with parallel I/O access
     */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    /*
     * OPTIONAL: It is generally recommended to set collective
     *           metadata reads/writes on FAPL to perform metadata reads
     *           collectively, which usually allows datasets
     *           to perform better at scale, although it is not
     *           strictly necessary.
     */
    H5Pset_all_coll_metadata_ops(plist_id, true);
    H5Pset_coll_metadata_write(plist_id, true);
    
    /*
     * Create a new file collectively and release property list identifier.
     */
    file_id = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    /*
     * Create the dataspace for the dataset.
     */
    const int dim_rank = 1;
    const int nelements = 10000000;
    hsize_t dim_global[dim_rank] = { static_cast<hsize_t>(nelements) };
    hid_t filespace = H5Screate_simple(dim_rank, dim_global, NULL);

    /*
     * Create the dataset with default properties and close filespace.
     */
    hid_t dset_id =
        H5Dcreate(file_id, "test-int", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    const int nelem_local = CAROM::split_dimension(nelements, MPI_COMM_WORLD);
    std::vector<int> offsets;
    int dummy = CAROM::get_global_offsets(nelem_local, offsets, MPI_COMM_WORLD);

    /* hyperslab selection parameters */
    hsize_t count[dim_rank];            
    hsize_t offset[dim_rank];
    count[0]  = nelem_local;
    // count[1]  = dimsf[1];
    offset[0] = offsets[rank];
    // offset[1] = 0;
    hid_t memspace  = H5Screate_simple(dim_rank, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    /*
     * Initialize data buffer
     */
    int data[nelem_local];
    for (int d = 0; d < nelem_local; d++)
        data[d] = d + offsets[rank];

    hid_t space = H5Screate_simple(1, dim_global, 0);
    CAROM_VERIFY(space >= 0);

    /*
     * Create property list for collective dataset write.
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    /*
     * Write the data collectively.
     */
    MPI_Barrier(MPI_COMM_WORLD);
    const auto start = steady_clock::now();
    errf = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
    CAROM_VERIFY(errf >= 0);
    MPI_Barrier(MPI_COMM_WORLD);
    const auto stop = steady_clock::now();
    const auto duration = duration_cast<milliseconds>(stop-start);
    printf("rank: %d, duration: %dms\n", rank, duration);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Fclose(file_id);
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}

#else // #ifndef CAROM_HAS_GTEST

int main()
{
    std::cout << "libROM was compiled without Google Test support, so unit "
              << "tests have been disabled. To enable unit tests, compile "
              << "libROM with Google Test support." << std::endl;
}

#endif // #endif CAROM_HAS_GTEST
