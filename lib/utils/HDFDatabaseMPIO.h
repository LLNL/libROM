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

#ifndef included_HDFDatabaseMPI_h
#define included_HDFDatabaseMPI_h

#include "HDFDatabase.h"

namespace CAROM {

/**
 * HDFDatabaseMPIO implements Multi-Path Input/Output (MPIO) of HDF5.
 */
class HDFDatabaseMPIO : public HDFDatabase
{
public:
    /**
     * @brief Default constructor.
     */
    HDFDatabaseMPIO();

    /**
     * @brief Destructor.
     */
    virtual
    ~HDFDatabaseMPIO() {}

    /**
     * @brief Creates a new HDF5 database file for distributed data
     *        with the supplied name.
     *
     * @param[in] file_name Name of HDF5 database file to create.
     *
     * @return True if file create was successful.
     */
    bool
    create_parallel(
        const std::string& file_name,
        const MPI_Comm comm) override;

    /**
     * @brief Creates a new HDF5 database file for non-distributed data
     *        with the supplied name.
     *        For HDFDatabaseMPIO, the behavior is equilvalent to create_parallel.
     *
     * @param[in] file_name Name of HDF5 database file to create.
     *
     * @return True if file create was successful.
     */
    bool
    create(
        const std::string& file_name) override
    {
        return create_parallel(file_name, MPI_COMM_WORLD);
    }

    /**
     * @brief Opens an existing HDF5 database file for distributed data
     *        with the supplied name.
     *
     * @param[in] file_name Name of existing HDF5 database file to open.
     * @param[in] type Read/write type ("r"/"wr")
     *
     * @return True if file open was successful.
     */
    bool
    open_parallel(
        const std::string& file_name,
        const std::string& type,
        const MPI_Comm comm) override;

    /**
     * @brief Opens an existing HDF5 database file for non-distributed data
     *        with the supplied name.
     *        For HDFDatabaseMPIO, the behavior is equilvalent to open_parallel.
     *
     * @param[in] file_name Name of existing HDF5 database file to open.
     * @param[in] type Read/write type ("r"/"wr")
     *
     * @return True if file open was successful.
     */
    bool
    open(
        const std::string& file_name,
        const std::string& type) override
    {
        return open_parallel(file_name, type, MPI_COMM_WORLD);
    }

    /**
     * @brief Writes a distributed array of integers
     *        associated with the supplied key
     *        to the currently open HDF5 database file.
     *
     * @pre !key.empty()
     * @pre data != nullptr
     * @pre nelem_local >= 0
     * @pre nelements > 0
     *
     * @param[in] key The key associated with the array of values to be
     *                written.
     * @param[in] data The array of integer values to be written.
     * @param[in] nelem_local The local number of integers in the array.
     */
    virtual
    void
    putIntegerArray_parallel(
        const std::string& key,
        const int* const data,
        int nelem_local);

    /**
     * @brief Writes a local array of integers only for root rank,
     * with associated supplied key to
     * the currently open HDF5 database file.
     *
     * @pre !key.empty()
     * @pre data != nullptr
     * @pre nelements > 0
     *
     * @param[in] key The key associated with the array of values to be
     *                written.
     * @param[in] data The array of integer values to be written.
     * @param[in] nelements The number of integers in the array.
     */
    void
    putIntegerArray(
        const std::string& key,
        const int* const data,
        int nelements) override
    {
        if (d_rank != 0)
            nelements = 0;
        putIntegerArray_parallel(key, data, nelements);
    }

    /**
     * @brief Writes a distributed array of doubles
     *        associated with the supplied key to
     *        the currently open HDF5 database file.
     *
     * @pre !key.empty()
     * @pre data != nullptr
     * @pre nelem_local >= 0
     * @pre nelements > 0
     *
     * @param[in] key The key associated with the array of values to be
     *                written.
     * @param[in] data The array of double values to be written.
     * @param[in] nelem_local The number of doubles in the array.
     */
    virtual
    void
    putDoubleArray_parallel(
        const std::string& key,
        const double* const data,
        int nelem_local);

    /**
     * @brief Writes an array of doubles in the root rank
     *        associated with the supplied key to
     *        the currently open HDF5 database file.
     *
     * @pre !key.empty()
     * @pre data != nullptr
     * @pre nelements > 0
     *
     * @param[in] key The key associated with the array of values to be
     *                written.
     * @param[in] data The array of double values to be written.
     * @param[in] nelements The number of doubles in the array.
     */
    void
    putDoubleArray(
        const std::string& key,
        const double* const data,
        int nelements) override
    {
        if (d_rank != 0)
            nelements = 0;
        putDoubleArray_parallel(key, data, nelements);
    }

    /**
     * @brief Reads a distributed array of integers
     *        associated with the supplied key
     *        from the currently open HDF5 database file.
     *
     * @pre !key.empty()
     * @pre data != nullptr || nelements == 0
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated array of integer values to be read.
     * @param[in] nelem_local The local number of integers in the array.
     */
    virtual
    void
    getIntegerArray_parallel(
        const std::string& key,
        int* data,
        int nelem_local);

    /**
     * @brief Reads an array of integers associated with the supplied key
     * from the currently open HDF5 database file.
     * All processes share the same non-distributed integer array.
     *
     * @pre !key.empty()
     * @pre data != nullptr || nelements == 0
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated array of integer values to be read.
     * @param[in] nelements The number of integers in the array.
     */
    void
    getIntegerArray(
        const std::string& key,
        int* data,
        int nelements) override
    {
        int read_size = (d_rank == 0) ? nelements : 0;
        getIntegerArray_parallel(key, data, read_size);

        MPI_Bcast(data, nelements, MPI_INT, 0, MPI_COMM_WORLD);
    }

    /**
     * @brief Reads a distributed array of doubles
     * associated with the supplied key
     * from the currently open HDF5 database file.
     *
     * @pre !key.empty()
     * @pre data != nullptr || nelements == 0
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated array of double values to be read.
     * @param[in] nelem_local The number of doubles in the array.
     */
    virtual
    void
    getDoubleArray_parallel(
        const std::string& key,
        double* data,
        int nelem_local);

    /**
     * @brief Reads an array of doubles associated with the supplied key
     * from the currently open HDF5 database file.
     * All processes share the same non-distributed double array.
     *
     * @pre !key.empty()
     * @pre data != nullptr || nelements == 0
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated array of double values to be read.
     * @param[in] nelements The number of doubles in the array.
     */
    void
    getDoubleArray(
        const std::string& key,
        double* data,
        int nelements) override
    {
        int read_size = (d_rank == 0) ? nelements : 0;
        getDoubleArray_parallel(key, data, read_size);

        MPI_Bcast(data, nelements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    /**
     * @brief Reads a distributed sub-array of doubles
     * associated with the supplied key
     * from the currently open HDF5 database file.
     *
     * @pre !key.empty()
     * @pre data != nullptr || nelements == 0
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated sub-array of double values to be read.
     * @param[in] nelem_local The number of doubles in the full array.
     * @param[in] idx_local The set of indices in the sub-array.
     */
    virtual
    void
    getDoubleArray_parallel(
        const std::string& key,
        double* data,
        int nelem_local,
        const std::vector<int>& idx_local);

    /**
     * @brief Reads a sub-array of doubles associated with the supplied key
     * from the currently open HDF5 database file.
     * All processes share the same non-distributed double array.
     *
     * @pre !key.empty()
     * @pre data != nullptr || nelements == 0
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated sub-array of double values to be read.
     * @param[in] nelements The number of doubles in the full array.
     * @param[in] idx The set of indices in the sub-array.
     */
    void
    getDoubleArray(
        const std::string& key,
        double* data,
        int nelements,
        const std::vector<int>& idx) override
    {
        int read_size = (d_rank == 0) ? nelements : 0;
        getDoubleArray_parallel(key, data, read_size, idx);

        MPI_Bcast(data, nelements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    /**
     * @brief Reads a distributed array of doubles
     * associated with the supplied key
     * from the currently open HDF5 database file.
     *
     * @pre !key.empty()
     * @pre data != nullptr || nelements == 0
     * @pre block_offset_global + block_size_global <= stride_global
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated array of double values to be read.
     * @param[in] nelem_local The number of doubles in the array at the local process.
     * @param[in] block_offset_global The initial offset in the global array.
     *                                Typically, this is a global column index of the matrix data.
     * @param[in] block_size_global The total block size to read from the HDF5 dataset.
     *                                Typically, this is a number of columns (in global) of the matrix data.
     * @param[in] stride_global The global stride to read from the HDF5 dataset.
     *                          Typically, this is the total number of columns of the matrix data.
     */
    virtual
    void
    getDoubleArray_parallel(
        const std::string& key,
        double* data,
        int nelem_local,
        int block_offset_global,
        int block_size_global,
        int stride_global);

    /**
     * @brief Reads an array of doubles associated with the supplied key
     * from the currently open HDF5 database file.
     * All processes share the same non-distributed double array.
     *
     * @pre !key.empty()
     * @pre data != nullptr || nelements == 0
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated array of double values to be read.
     * @param[in] nelements The number of doubles in the array.
     * @param[in] offset The initial offset in the array.
     * @param[in] block_size The block size to read from the HDF5 dataset.
     * @param[in] stride The stride to read from the HDF5 dataset.
     */
    void
    getDoubleArray(
        const std::string& key,
        double* data,
        int nelements,
        int offset,
        int block_size,
        int stride) override
    {
        int read_size = (d_rank == 0) ? nelements : 0;
        getDoubleArray_parallel(key, data, read_size, offset, block_size, stride);

        MPI_Bcast(data, nelements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    void
    writeAttribute(
        int type_key,
        hid_t dataset_id) override;

private:
    MPI_Comm d_comm;
    int d_rank;
};

}

#endif
