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
#include "linalg/BasisGenerator.h"
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

    const int dim_rank = 2;
    const int nrow = 1000, ncol = 400;
    const int nrow_local = CAROM::split_dimension(nrow, MPI_COMM_WORLD);
    std::vector<int> offsets;
    int dummy = CAROM::get_global_offsets(nrow_local, offsets, MPI_COMM_WORLD);

    /*
     * Initialize data buffer
     */
    int data[nrow_local * ncol];
    for (int d = 0; d < nrow_local * ncol; d++)
        data[d] = d + offsets[rank] * ncol;

    MPI_Barrier(MPI_COMM_WORLD);
    const auto start = steady_clock::now();

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
    hsize_t dim_global[dim_rank] = { static_cast<hsize_t>(nrow), static_cast<hsize_t>(ncol) };
    hid_t filespace = H5Screate_simple(dim_rank, dim_global, NULL);

    /*
     * Create the dataset with default properties and close filespace.
     */
    hid_t dset_id =
        H5Dcreate(file_id, "test-int", H5T_NATIVE_INT, filespace, H5P_DEFAULT,
                  H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    /* hyperslab selection parameters */
    hsize_t count[dim_rank];
    hsize_t offset[dim_rank];
    count[0]  = nrow_local;
    count[1]  = ncol;
    offset[0] = offsets[rank];
    offset[1] = 0;
    hid_t memspace  = H5Screate_simple(dim_rank, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    /*
     * Create property list for collective dataset write.
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    /*
     * Write the data collectively.
     */
    errf = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
    CAROM_VERIFY(errf >= 0);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Fclose(file_id);

    MPI_Barrier(MPI_COMM_WORLD);
    const auto stop = steady_clock::now();
    const auto duration = duration_cast<milliseconds>(stop-start);
    printf("rank: %d, duration: %dms\n", rank, duration);
}

TEST(HDF5, Test_parallel_reading)
{
    int nproc, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int dim_rank = 2;
    const int nrow = 1000, ncol = 400;
    const int nrow_local = CAROM::split_dimension(nrow, MPI_COMM_WORLD);
    std::vector<int> offsets;
    int dummy = CAROM::get_global_offsets(nrow_local, offsets, MPI_COMM_WORLD);

    /*
     * Initialize data buffer
     */
    int data[nrow_local * ncol];

    MPI_Barrier(MPI_COMM_WORLD);
    const auto start = steady_clock::now();

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
     * Open a new file collectively and release property list identifier.
     */
    file_id = H5Fopen("test.h5", H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    /*
     * Open the dataset with default properties.
     */
    hid_t dset_id = H5Dopen(file_id, "test-int", H5P_DEFAULT);

    /*
     * Get filespace and read the global size.
     */
    hid_t filespace = H5Dget_space(dset_id);
    CAROM_VERIFY(filespace >= 0);
    int ndims = H5Sget_simple_extent_ndims(filespace);
    CAROM_VERIFY(ndims == dim_rank);
    hsize_t dim_global[dim_rank];
    errf = H5Sget_simple_extent_dims(filespace, dim_global, NULL);
    CAROM_VERIFY(errf >= 0);
    CAROM_VERIFY(dim_global[0] == nrow);
    CAROM_VERIFY(dim_global[1] == ncol);

    H5Sclose(filespace);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    /* hyperslab selection parameters */
    hsize_t count[dim_rank];
    hsize_t offset[dim_rank];
    count[0]  = nrow_local;
    count[1]  = ncol;
    offset[0] = offsets[rank];
    offset[1] = 0;
    hid_t memspace  = H5Screate_simple(dim_rank, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    /*
     * Create property list for collective dataset write.
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    /*
     * Read the data collectively.
     */
    errf = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
    CAROM_VERIFY(errf >= 0);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Fclose(file_id);

    MPI_Barrier(MPI_COMM_WORLD);
    const auto stop = steady_clock::now();
    const auto duration = duration_cast<milliseconds>(stop-start);
    printf("rank: %d, duration: %dms\n", rank, duration);

    for (int d = 0; d < nrow_local * ncol; d++)
        EXPECT_TRUE(data[d] == d + offsets[rank] * ncol);
}

TEST(HDF5, Test_selective_parallel_writing)
{
    int nproc, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Only the first 2 processes will participate in I/O. */
    const int ioproc = std::min(2, nproc);
    int iorank[ioproc];
    for (int r = 0; r < ioproc; r++)
        iorank[r] = r;
    MPI_Group world_grp, io_grp;
    MPI_Comm_group(MPI_COMM_WORLD, &world_grp);
    MPI_Group_incl(world_grp, ioproc, iorank, &io_grp);
    MPI_Comm io_comm;
    MPI_Comm_create(MPI_COMM_WORLD, io_grp, &io_comm);

    const int dim_rank = 2;
    const int nrow = 10, ncol = 4;
    int nrow_local = 0;
    if (rank < ioproc)
        nrow_local = CAROM::split_dimension(nrow, io_comm);
    printf("rank %d, I/O row size: %d\n", rank, nrow_local);
    std::vector<int> offsets;
    int dummy = CAROM::get_global_offsets(nrow_local, offsets, MPI_COMM_WORLD);
    for (int k = offsets.size(); k < nproc; k++)
        offsets.push_back(offsets.back());

    /*
     * Initialize data buffer
     */
    int data[nrow_local * ncol];
    for (int d = 0; d < nrow_local * ncol; d++)
        data[d] = d + offsets[rank] * ncol;

    MPI_Barrier(MPI_COMM_WORLD);
    const auto start = steady_clock::now();

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
    hsize_t dim_global[dim_rank] = { static_cast<hsize_t>(nrow), static_cast<hsize_t>(ncol) };
    hid_t filespace = H5Screate_simple(dim_rank, dim_global, NULL);

    /*
     * Create the dataset with default properties and close filespace.
     */
    hid_t dset_id =
        H5Dcreate(file_id, "test-int", H5T_NATIVE_INT, filespace, H5P_DEFAULT,
                  H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    /* hyperslab selection parameters */
    hsize_t count[dim_rank];
    hsize_t offset[dim_rank];
    count[0]  = nrow_local;
    count[1]  = ncol;
    offset[0] = offsets[rank];
    offset[1] = 0;
    hid_t memspace = H5Screate_simple(dim_rank, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dset_id);
    if (rank < ioproc)
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    else
        H5Sselect_none(filespace);

    /*
     * Create property list for collective dataset write.
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    /*
     * Write the data collectively.
     */
    errf = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
    CAROM_VERIFY(errf >= 0);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Fclose(file_id);

    MPI_Barrier(MPI_COMM_WORLD);
    const auto stop = steady_clock::now();
    const auto duration = duration_cast<milliseconds>(stop-start);
    printf("rank: %d, duration: %dms\n", rank, duration);
}

TEST(HDF5, Test_selective_parallel_reading)
{
    int nproc, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Only the first 2 processes will participate in I/O. */
    const int ioproc = std::min(2, nproc);
    int iorank[ioproc];
    for (int r = 0; r < ioproc; r++)
        iorank[r] = r;
    MPI_Group world_grp, io_grp;
    MPI_Comm_group(MPI_COMM_WORLD, &world_grp);
    MPI_Group_incl(world_grp, ioproc, iorank, &io_grp);
    MPI_Comm io_comm;
    MPI_Comm_create(MPI_COMM_WORLD, io_grp, &io_comm);

    const int dim_rank = 2;
    const int nrow = 10, ncol = 4;
    int nrow_local = 0;
    if (rank < ioproc)
        nrow_local = CAROM::split_dimension(nrow, io_comm);
    printf("rank %d, I/O row size: %d\n", rank, nrow_local);
    std::vector<int> offsets;
    int dummy = CAROM::get_global_offsets(nrow_local, offsets, MPI_COMM_WORLD);

    /*
     * Initialize data buffer
     */
    int data[nrow_local * ncol];

    MPI_Barrier(MPI_COMM_WORLD);
    const auto start = steady_clock::now();

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
     * Open a new file collectively and release property list identifier.
     */
    file_id = H5Fopen("test.h5", H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    /*
     * Open the dataset with default properties.
     */
    hid_t dset_id = H5Dopen(file_id, "test-int", H5P_DEFAULT);

    /*
     * Get filespace and read the global size.
     */
    hid_t filespace = H5Dget_space(dset_id);
    CAROM_VERIFY(filespace >= 0);
    int ndims = H5Sget_simple_extent_ndims(filespace);
    CAROM_VERIFY(ndims == dim_rank);
    hsize_t dim_global[dim_rank];
    errf = H5Sget_simple_extent_dims(filespace, dim_global, NULL);
    CAROM_VERIFY(errf >= 0);
    CAROM_VERIFY(dim_global[0] == nrow);
    CAROM_VERIFY(dim_global[1] == ncol);

    H5Sclose(filespace);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    /* hyperslab selection parameters */
    hsize_t count[dim_rank];
    hsize_t offset[dim_rank];
    count[0]  = nrow_local;
    count[1]  = ncol;
    offset[0] = offsets[rank];
    offset[1] = 0;
    hid_t memspace  = H5Screate_simple(dim_rank, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dset_id);
    if (rank < ioproc)
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    else
        H5Sselect_none(filespace);

    /*
     * Create property list for collective dataset write.
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    /*
     * Read the data collectively.
     */
    errf = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
    CAROM_VERIFY(errf >= 0);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Fclose(file_id);

    MPI_Barrier(MPI_COMM_WORLD);
    const auto stop = steady_clock::now();
    const auto duration = duration_cast<milliseconds>(stop-start);
    printf("rank: %d, duration: %dms\n", rank, duration);

    for (int d = 0; d < nrow_local * ncol; d++)
        EXPECT_TRUE(data[d] == d + offsets[rank] * ncol);
}

TEST(BasisGeneratorIO, HDFDatabase)
{
    // Get the rank of this process, and the number of processors.
    int mpi_init, d_rank, d_num_procs;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        MPI_Init(nullptr, nullptr);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &d_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &d_num_procs);

    const int nrow = 100, ncol = 10;
    int nrow_local = CAROM::split_dimension(nrow, MPI_COMM_WORLD);
    std::vector<int> row_offset(d_num_procs + 1);
    const int dummy = CAROM::get_global_offsets(nrow_local, row_offset,
                                                MPI_COMM_WORLD);
    EXPECT_EQ(nrow, dummy);

    std::default_random_engine generator;
    generator.seed(
        1234); // fix the seed to keep the same result for different nproc.
    std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);
    std::normal_distribution<double> normal_distribution(0.0, 1.0);

    // distribute from a global matrix to keep the same system for different nproc.
    CAROM::Matrix snapshots(nrow, ncol, false);
    for (int i = 0; i < nrow; i++)
        for (int j = 0; j < ncol; j++)
            snapshots(i, j) = normal_distribution(generator);
    snapshots.distribute(nrow_local);

    CAROM::Options svd_options = CAROM::Options(nrow_local, ncol, 1);
    svd_options.setMaxBasisDimension(nrow);
    svd_options.setRandomizedSVD(false);
    CAROM::BasisGenerator sampler(svd_options, false, "test_basis");
    CAROM::Vector sample(nrow_local, true);
    for (int s = 0; s < ncol; s++)
    {
        for (int d = 0; d < nrow_local; d++)
            sample(d) = snapshots(d, s);

        sampler.takeSample(sample.getData(), 0, 0);
    }
    sampler.endSamples();
    sampler.writeSnapshot();

    const CAROM::Matrix *spatial_basis = sampler.getSpatialBasis();
    const CAROM::Matrix *snapshot = sampler.getSnapshotMatrix();

    CAROM::BasisReader basis_reader("test_basis");
    const CAROM::Matrix *spatial_basis1 = basis_reader.getSpatialBasis();

    EXPECT_EQ(spatial_basis->numRows(), spatial_basis1->numRows());
    EXPECT_EQ(spatial_basis->numColumns(), spatial_basis1->numColumns());
    for (int i = 0; i < spatial_basis->numRows(); i++)
        for (int j = 0; j < spatial_basis->numColumns(); j++)
            EXPECT_NEAR((*spatial_basis)(i, j), (*spatial_basis1)(i, j), 1.0e-15);

    CAROM::BasisReader snapshot_reader("test_basis_snapshot");
    const CAROM::Matrix *snapshot1 = snapshot_reader.getSnapshotMatrix();

    EXPECT_EQ(snapshot->numRows(), snapshot1->numRows());
    EXPECT_EQ(snapshot->numColumns(), snapshot1->numColumns());
    for (int i = 0; i < snapshot->numRows(); i++)
        for (int j = 0; j < snapshot->numColumns(); j++)
            EXPECT_NEAR((*snapshot)(i, j), (*snapshot1)(i, j), 1.0e-15);
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
