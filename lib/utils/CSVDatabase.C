/******************************************************************************
 *
 * Copyright (c) 2013-2021, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete database implementation using CSV.

#include "CSVDatabase.h"
#include "Utilities.h"

namespace CAROM {

CSVDatabase::CSVDatabase() :
    d_file_name("")
{
}

CSVDatabase::~CSVDatabase()
{
}

bool
CSVDatabase::create(
    const std::string& file_name)
{
    return true;
}

bool
CSVDatabase::open(
    const std::string& file_name,
    const std::string& type)
{
    return true;
}

bool
CSVDatabase::close()
{
    return true;
}

void
CSVDatabase::putIntegerArray(
    const std::string& key,
    const int* const data,
    int nelements)
{
    CAROM_ASSERT(!key.empty());
    CAROM_ASSERT(data != 0);
    CAROM_ASSERT(nelements > 0);

    std::ofstream d_fs(key.c_str());
    for (int i = 0; i < nelements; ++i)
    {
        d_fs << data[i] << std::endl;
    }
}

void
CSVDatabase::putDoubleArray(
    const std::string& key,
    const double* const data,
    int nelements)
{
    CAROM_ASSERT(!key.empty());
    CAROM_ASSERT(data != 0);
    CAROM_ASSERT(nelements > 0);

    std::ofstream d_fs(key.c_str());
    for (int i = 0; i < nelements; ++i)
    {
        d_fs << data[i] << std::endl;
    }
}

void
CSVDatabase::getIntegerArray(
    const std::string& key,
    int* data,
    int nelements)
{
    CAROM_ASSERT(!key.empty());
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(nelements);
#endif

    std::ifstream d_fs(key.c_str());
    std::vector<int> tmp;
    int data_entry = 0.0;
    while (d_fs >> data_entry)
    {
        tmp.push_back(data_entry);
    }
    data = new int[tmp.Ssze()];
    for (int i = 0; i < tmp.size(); ++i)
    {
        data[i] = tmp[i];
    }
}

void
CSVDatabase::getDoubleArray(
    const std::string& key,
    double* data,
    int nelements)
{
    CAROM_ASSERT(!key.empty());
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(nelements);
#endif

    std::ifstream d_fs(key.c_str());
    std::vector<double> tmp;
    double data_entry = 0.0;
    while (d_fs >> data_entry)
    {
        tmp.push_back(data_entry);
    }
    data = new double[tmp.size()];
    for (int i = 0; i < tmp.size(); ++i)
    {
        data[i] = tmp[i];
    }
}

void
CSVDatabase::getDoubleArray(
    const std::string& key,
    double* data,
    int nelements,
    int offset,
    int block_size,
    int stride)
{
    CAROM_ASSERT(!key.empty());
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(nelements);
#endif

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
    hid_t dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
    hid_t dset = H5Dopen(d_group_id, key.c_str());
#endif
    CAROM_ASSERT(dset >= 0);

    hid_t dspace = H5Dget_space(dset);
    CAROM_ASSERT(dspace >= 0);

    hsize_t nsel = H5Sget_select_npoints(dspace);

    herr_t errf;
    if (nsel > 0) {

        hsize_t num_blocks[1] = {(hsize_t) nelements/block_size};
        hsize_t buffer_array_size[1] = {(hsize_t) nelements};
        hsize_t offsets[1] = {(hsize_t) offset};
        hsize_t strides[1] = {(hsize_t) stride};
        hsize_t block_sizes[1] = {(hsize_t) block_size};
        hid_t nodespace = H5Screate_simple(1,buffer_array_size,NULL);

        // select hyperslab
        H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offsets,
                            strides, num_blocks, block_sizes);

        errf = H5Dread(dset, H5T_NATIVE_DOUBLE, nodespace, dspace, H5P_DEFAULT, data);
        CAROM_ASSERT(errf >= 0);
    }

    errf = H5Sclose(dspace);
    CAROM_ASSERT(errf >= 0);

    errf = H5Dclose(dset);
    CAROM_ASSERT(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif
}

}
