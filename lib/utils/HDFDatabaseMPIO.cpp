/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete database implementation using HDF5.

#include "HDFDatabaseMPIO.h"
#include "utils/mpi_utils.h"
#include <iostream>

namespace CAROM {

HDFDatabaseMPIO::HDFDatabaseMPIO() :
#if HDF5_IS_PARALLEL
    HDFDatabase(),
    d_rank(-1),
    d_comm(MPI_COMM_NULL)
#else
    HDFDatabase()
#endif
{}

#if HDF5_IS_PARALLEL

bool
HDFDatabaseMPIO::create(
    const std::string& file_name,
    const MPI_Comm comm)
{
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init) {
        CAROM_VERIFY(comm != MPI_COMM_NULL);
        d_comm = comm;
        MPI_Comm_rank(d_comm, &d_rank);
    }
    else
        d_rank = 0;

    /* Create a single file with a name compatible to HDFDatabase. */
    std::string file_name_ext(file_name + ".000000");
    if (d_rank == 0)
        Database::create(file_name_ext);
    CAROM_VERIFY(!file_name.empty());

    hid_t plist_id;
    /*
     * Set up file access property list with parallel I/O access
     */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, d_comm, MPI_INFO_NULL);
    /*
     * OPTIONAL: It is generally recommended to set collective
     *           metadata reads/writes on FAPL to perform metadata reads
     *           collectively, which usually allows datasets
     *           to perform better at scale, although it is not
     *           strictly necessary.
     */
    H5Pset_all_coll_metadata_ops(plist_id, true);
    H5Pset_coll_metadata_write(plist_id, true);

    hid_t file_id = H5Fcreate(file_name_ext.c_str(),
                              H5F_ACC_TRUNC,
                              H5P_DEFAULT,
                              plist_id);
    bool result = file_id >= 0;
    CAROM_VERIFY(result);
    d_is_file = true;
    d_file_id = file_id;
    d_group_id = file_id;

    H5Pclose(plist_id);

    return result;
}

bool
HDFDatabaseMPIO::open(
    const std::string& file_name,
    const std::string& type,
    const MPI_Comm comm)
{
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init) {
        CAROM_VERIFY(comm != MPI_COMM_NULL);
        d_comm = comm;
        MPI_Comm_rank(d_comm, &d_rank);
    }
    else
        d_rank = 0;

    std::string file_name_ext(file_name + ".000000");
    if (d_rank == 0)
        Database::open(file_name_ext, type);
    CAROM_VERIFY(!file_name.empty());
    CAROM_VERIFY(type == "r" || type == "wr");

    hid_t plist_id;
    /*
     * Set up file access property list with parallel I/O access
     */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, d_comm, MPI_INFO_NULL);
    /*
     * OPTIONAL: It is generally recommended to set collective
     *           metadata reads/writes on FAPL to perform metadata reads
     *           collectively, which usually allows datasets
     *           to perform better at scale, although it is not
     *           strictly necessary.
     */
    H5Pset_all_coll_metadata_ops(plist_id, true);
    H5Pset_coll_metadata_write(plist_id, true);

    hid_t file_id;
    if (type == "r")
    {
        file_id = H5Fopen(file_name_ext.c_str(),
                          H5F_ACC_RDONLY,
                          plist_id);
    }
    else if (type == "wr")
    {
        file_id = H5Fopen(file_name_ext.c_str(),
                          H5F_ACC_RDWR,
                          plist_id);
    }
    bool result = file_id >= 0;
    CAROM_VERIFY(result);
    d_is_file = true;
    d_file_id = file_id;
    d_group_id = file_id;

    H5Pclose(plist_id);

    return result;
}

void
HDFDatabaseMPIO::putIntegerArray_parallel(
    const std::string& key,
    const int* const data,
    int nelem_local)
{
    CAROM_VERIFY(!key.empty());
    CAROM_VERIFY(data != nullptr);
    CAROM_VERIFY(nelem_local >= 0);

    /* determine global nelements and offsets */
    std::vector<int> offsets;
    int nelements = CAROM::get_global_offsets(nelem_local, offsets, d_comm);
    CAROM_VERIFY(nelements > 0);

    const int dim_rank = 1;
    hsize_t dim[dim_rank] = { static_cast<hsize_t>(nelements) };
    hid_t filespace = H5Screate_simple(dim_rank, dim, NULL);
    CAROM_VERIFY(filespace >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
    hid_t dataset = H5Dcreate(d_group_id,
                              key.c_str(),
                              H5T_NATIVE_INT,
                              filespace,
                              H5P_DEFAULT,
                              H5P_DEFAULT,
                              H5P_DEFAULT);
#else
    CAROM_ERROR("HDFDatabaseMPIO is not compatible with current HDF5 version!\n");
#endif
    CAROM_VERIFY(dataset >= 0);
    H5Sclose(filespace);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    /* hyperslab selection parameters */
    hsize_t count[dim_rank];
    hsize_t offset[dim_rank];
    count[0]  = nelem_local;
    offset[0] = offsets[d_rank];
    hid_t memspace = H5Screate_simple(dim_rank, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dataset);
    if (nelem_local > 0)
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    else
        H5Sselect_none(filespace);

    /*
     * Create property list for collective dataset write.
     */
    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    herr_t errf = H5Dwrite(dataset,
                           H5T_NATIVE_INT,
                           memspace,
                           filespace,
                           plist_id,
                           data);
    CAROM_VERIFY(errf >= 0);

    // Write attribute so we know what kind of data this is.
    writeAttribute(KEY_INT_ARRAY, dataset);

    errf = H5Sclose(filespace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Sclose(memspace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Pclose(plist_id);
    CAROM_VERIFY(errf >= 0);

    errf = H5Dclose(dataset);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif
}

void
HDFDatabaseMPIO::putDoubleArray_parallel(
    const std::string& key,
    const double* const data,
    int nelem_local)
{
    CAROM_VERIFY(!key.empty());
    CAROM_VERIFY(data != nullptr);
    CAROM_VERIFY(nelem_local >= 0);

    /* determine global nelements and offsets */
    std::vector<int> offsets;
    int nelements = CAROM::get_global_offsets(nelem_local, offsets, d_comm);
    CAROM_VERIFY(nelements > 0);

    const int dim_rank = 1;
    hsize_t dim[dim_rank] = { static_cast<hsize_t>(nelements) };
    hid_t filespace = H5Screate_simple(dim_rank, dim, 0);
    CAROM_VERIFY(filespace >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
    hid_t dataset = H5Dcreate(d_group_id,
                              key.c_str(),
                              H5T_NATIVE_DOUBLE,
                              filespace,
                              H5P_DEFAULT,
                              H5P_DEFAULT,
                              H5P_DEFAULT);
#else
    CAROM_ERROR("HDFDatabaseMPIO is not compatible with current HDF5 version!\n");
#endif
    CAROM_VERIFY(dataset >= 0);
    H5Sclose(filespace);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    /* hyperslab selection parameters */
    hsize_t count[dim_rank];
    hsize_t offset[dim_rank];
    count[0]  = nelem_local;
    offset[0] = offsets[d_rank];
    hid_t memspace = H5Screate_simple(dim_rank, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dataset);
    if (nelem_local > 0)
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    else
        H5Sselect_none(filespace);

    /*
     * Create property list for collective dataset write.
     */
    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    herr_t errf = H5Dwrite(dataset,
                           H5T_NATIVE_DOUBLE,
                           memspace,
                           filespace,
                           plist_id,
                           data);
    CAROM_VERIFY(errf >= 0);

    // Write attribute so we know what kind of data this is.
    writeAttribute(KEY_DOUBLE_ARRAY, dataset);

    errf = H5Sclose(filespace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Sclose(memspace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Pclose(plist_id);
    CAROM_VERIFY(errf >= 0);

    errf = H5Dclose(dataset);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif
}

void
HDFDatabaseMPIO::getIntegerArray_parallel(
    const std::string& key,
    int* data,
    int nelem_local)
{
    CAROM_VERIFY(nelem_local >= 0);
    /* determine global nelements and offsets */
    std::vector<int> offsets;
    int nelements = CAROM::get_global_offsets(nelem_local, offsets, d_comm);
    if (nelements == 0) return;

    CAROM_VERIFY(!key.empty());
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(nelem_local);
#endif

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
    hid_t dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
    CAROM_ERROR("HDFDatabaseMPIO is not compatible with current HDF5 version!\n");
#endif
    CAROM_VERIFY(dset >= 0);

    hid_t filespace = H5Dget_space(dset);
    CAROM_VERIFY(filespace >= 0);

    hsize_t nsel = H5Sget_select_npoints(filespace);
    CAROM_VERIFY(static_cast<int>(nsel) == nelements);

    const int dim_rank = 1;
    H5Sclose(filespace);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    /* hyperslab selection parameters */
    hsize_t count[dim_rank];
    hsize_t offset[dim_rank];
    count[0]  = nelem_local;
    offset[0] = offsets[d_rank];
    hid_t memspace  = H5Screate_simple(dim_rank, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dset);
    if (nelem_local > 0)
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    else
        H5Sselect_none(filespace);

    /*
     * Create property list for collective dataset write.
     */
    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    herr_t errf;
    errf = H5Dread(dset, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
    CAROM_VERIFY(errf >= 0);

    errf = H5Sclose(filespace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Sclose(memspace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Pclose(plist_id);
    CAROM_VERIFY(errf >= 0);

    errf = H5Dclose(dset);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif
}

void
HDFDatabaseMPIO::getDoubleArray_parallel(
    const std::string& key,
    double* data,
    int nelem_local)
{
    CAROM_VERIFY(nelem_local >= 0);
    /* determine global nelements and offsets */
    std::vector<int> offsets;
    int nelements = CAROM::get_global_offsets(nelem_local, offsets, d_comm);
    if (nelements == 0) return;

    CAROM_VERIFY(!key.empty());
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(nelements);
#endif

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
    hid_t dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
    hid_t dset = H5Dopen(d_group_id, key.c_str());
#endif
    CAROM_VERIFY(dset >= 0);

    hid_t filespace = H5Dget_space(dset);
    CAROM_VERIFY(filespace >= 0);

    hsize_t nsel = H5Sget_select_npoints(filespace);
    CAROM_VERIFY(static_cast<int>(nsel) == nelements);

    const int dim_rank = 1;
    H5Sclose(filespace);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    /* hyperslab selection parameters */
    hsize_t count[dim_rank];
    hsize_t offset[dim_rank];
    count[0]  = nelem_local;
    offset[0] = offsets[d_rank];
    hid_t memspace  = H5Screate_simple(dim_rank, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dset);
    if (nelem_local > 0)
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    else
        H5Sselect_none(filespace);

    /*
     * Create property list for collective dataset write.
     */
    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    herr_t errf;
    errf = H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
    CAROM_VERIFY(errf >= 0);

    errf = H5Sclose(filespace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Sclose(memspace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Pclose(plist_id);
    CAROM_VERIFY(errf >= 0);

    errf = H5Dclose(dset);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif
}

void
HDFDatabaseMPIO::getDoubleArray_parallel(
    const std::string& key,
    double* data,
    int nelem_local,
    const std::vector<int>& idx_local)
{
    std::vector<int> idx_tmp = idx_local;
    if (idx_local.size() == 0)
    {
        for (int i = 0; i < nelem_local; i++)
            idx_tmp.push_back(i);
    }

    std::vector<double> alldata(nelem_local);
    getDoubleArray_parallel(key, alldata.data(), nelem_local);
    int k = 0;
    for (int i = 0; i < nelem_local; ++i)
    {
        if (idx_tmp[k] == i)
        {
            data[k++] = alldata[i];
        }
        if (k == idx_tmp.size())
        {
            break;
        }
    }
    CAROM_VERIFY(k == idx_tmp.size());
}

void
HDFDatabaseMPIO::getDoubleArray_parallel(
    const std::string& key,
    double* data,
    int nelem_local,
    int block_offset_global,
    int block_size_global,
    int stride_global)
{
    CAROM_VERIFY(nelem_local >= 0);
    CAROM_VERIFY(CAROM::is_same(stride_global, d_comm));
    CAROM_VERIFY(CAROM::is_same(block_size_global, d_comm));
    CAROM_VERIFY(CAROM::is_same(block_offset_global, d_comm));
    /* We assume stride sets the maximum block size */
    CAROM_VERIFY(block_offset_global + block_size_global <= stride_global);
    /* determine global nelements and offsets */
    hsize_t num_local_blocks = nelem_local / block_size_global;
    std::vector<int> global_offsets;
    int nelements = CAROM::get_global_offsets(num_local_blocks * stride_global,
                    global_offsets, d_comm);

    CAROM_VERIFY(!key.empty());
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(nelements);
#endif

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
    hid_t dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
    hid_t dset = H5Dopen(d_group_id, key.c_str());
#endif
    CAROM_VERIFY(dset >= 0);

    hid_t filespace = H5Dget_space(dset);
    CAROM_VERIFY(filespace >= 0);

    hsize_t nsel = H5Sget_select_npoints(filespace);
    CAROM_VERIFY((nsel == 0) || (nsel == nelements));

    const int dim_rank = 1;
    H5Sclose(filespace);

    herr_t errf;
    if (nsel > 0) {
        /*
         * Each process defines dataset in memory and writes it to the hyperslab
         * in the file.
         */
        /* hyperslab selection parameters */
        hsize_t num_blocks[dim_rank] = {num_local_blocks};
        hsize_t buffer_array_size[dim_rank] = {(hsize_t) nelem_local};
        hsize_t offset[dim_rank] = {(hsize_t) global_offsets[d_rank] + block_offset_global};
        hsize_t strides[dim_rank] = {(hsize_t) stride_global};
        hsize_t block_sizes[1] = {(hsize_t) block_size_global};
        hid_t memspace  = H5Screate_simple(dim_rank, buffer_array_size, NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset);
        if (nelem_local > 0)
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset,
                                strides, num_blocks, block_sizes);
        else
            H5Sselect_none(filespace);

        /*
         * Create property list for collective dataset write.
         */
        hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        errf = H5Dread(dset, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
        CAROM_VERIFY(errf >= 0);

        errf = H5Sclose(memspace);
        CAROM_VERIFY(errf >= 0);

        errf = H5Pclose(plist_id);
        CAROM_VERIFY(errf >= 0);
    }

    errf = H5Sclose(filespace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Dclose(dset);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif
}

void
HDFDatabaseMPIO::writeAttribute(
    int type_key,
    hid_t dataset_id)
{
    CAROM_VERIFY(CAROM::is_same(type_key, d_comm));

    hid_t attr_id = H5Screate(H5S_SCALAR);
    CAROM_VERIFY(attr_id >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
    hid_t attr = H5Acreate(dataset_id,
                           "Type",
                           H5T_NATIVE_INT,
                           attr_id,
                           H5P_DEFAULT,
                           H5P_DEFAULT);
#else
    CAROM_ERROR("HDFDatabaseMPIO is not compatible with current HDF5 version!\n");
#endif
    CAROM_VERIFY(attr >= 0);

    /* only the last process writes the attribute */
    herr_t errf = H5Awrite(attr, H5T_NATIVE_INT, &type_key);
    CAROM_VERIFY(errf >= 0);

    errf = H5Aclose(attr);
    CAROM_VERIFY(errf >= 0);

    errf = H5Sclose(attr_id);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif
}

#endif

}
