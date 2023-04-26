/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: The concrete database implementation using HDF5.

#include "HDFDatabase.h"
#include "Utilities.h"

namespace CAROM {

const int HDFDatabase::KEY_DOUBLE_ARRAY = 0;
const int HDFDatabase::KEY_INT_ARRAY = 1;

HDFDatabase::HDFDatabase() :
    d_is_file(false),
    d_file_id(-1),
    d_group_id(-1)
{
}

HDFDatabase::~HDFDatabase()
{
    if (d_is_file) {
        close();
    }

    if (d_group_id != -1) {
        herr_t errf = H5Gclose(d_group_id);
        CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
        CAROM_NULL_USE(errf);
#endif
    }
}

bool
HDFDatabase::create(
    const std::string& file_name)
{
    CAROM_VERIFY(!file_name.empty());
    hid_t file_id = H5Fcreate(file_name.c_str(),
                              H5F_ACC_TRUNC,
                              H5P_DEFAULT,
                              H5P_DEFAULT);
    bool result = file_id >= 0;
    CAROM_VERIFY(result);
    d_is_file = true;
    d_file_id = file_id;
    d_group_id = file_id;

    return result;
}

bool
HDFDatabase::open(
    const std::string& file_name,
    const std::string& type)
{
    CAROM_VERIFY(!file_name.empty());
    CAROM_VERIFY(type == "r" || type == "wr");
    hid_t file_id;
    if (type == "r")
    {
        file_id = H5Fopen(file_name.c_str(),
                          H5F_ACC_RDONLY,
                          H5P_DEFAULT);
    }
    else if (type == "wr")
    {
        file_id = H5Fopen(file_name.c_str(),
                          H5F_ACC_RDWR,
                          H5P_DEFAULT);
    }
    bool result = file_id >= 0;
    CAROM_VERIFY(result);
    d_is_file = true;
    d_file_id = file_id;
    d_group_id = file_id;

    return result;
}

bool
HDFDatabase::close()
{
    herr_t errf = 0;
    if (d_is_file) {
        errf = H5Fclose(d_file_id);
        CAROM_VERIFY(errf >= 0);

        if (d_group_id == d_file_id) {
            d_group_id = -1;
        }
        d_file_id = -1;
        d_is_file = false;
    }

    return errf >= 0;
}

void
HDFDatabase::putIntegerArray(
    const std::string& key,
    const int* const data,
    int nelements)
{
    CAROM_VERIFY(!key.empty());
    CAROM_VERIFY(data != nullptr);
    CAROM_VERIFY(nelements > 0);

    hsize_t dim[] = { static_cast<hsize_t>(nelements) };
    hid_t space = H5Screate_simple(1, dim, 0);
    CAROM_VERIFY(space >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
    hid_t dataset = H5Dcreate(d_group_id,
                              key.c_str(),
                              H5T_STD_I32BE,
                              space,
                              H5P_DEFAULT,
                              H5P_DEFAULT,
                              H5P_DEFAULT);
#else
    hid_t dataset = H5Dcreate(d_group_id,
                              key.c_str(),
                              H5T_STD_I32BE,
                              space,
                              H5P_DEFAULT);
#endif
    CAROM_VERIFY(dataset >= 0);

    herr_t errf = H5Dwrite(dataset,
                           H5T_NATIVE_INT,
                           H5S_ALL,
                           H5S_ALL,
                           H5P_DEFAULT,
                           data);
    CAROM_VERIFY(errf >= 0);

    // Write attribute so we know what kind of data this is.
    writeAttribute(KEY_INT_ARRAY, dataset);

    errf = H5Sclose(space);
    CAROM_VERIFY(errf >= 0);

    errf = H5Dclose(dataset);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif
}

void
HDFDatabase::putDoubleArray(
    const std::string& key,
    const double* const data,
    int nelements)
{
    CAROM_VERIFY(!key.empty());
    CAROM_VERIFY(data != nullptr);
    CAROM_VERIFY(nelements > 0);

    hsize_t dim[] = { static_cast<hsize_t>(nelements) };
    hid_t space = H5Screate_simple(1, dim, 0);
    CAROM_VERIFY(space >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
    hid_t dataset = H5Dcreate(d_group_id,
                              key.c_str(),
                              H5T_IEEE_F64BE,
                              space,
                              H5P_DEFAULT,
                              H5P_DEFAULT,
                              H5P_DEFAULT);
#else
    hid_t dataset = H5Dcreate(d_group_id,
                              key.c_str(),
                              H5T_IEEE_F64BE,
                              space,
                              H5P_DEFAULT);
#endif
    CAROM_VERIFY(dataset >= 0);

    herr_t errf = H5Dwrite(dataset,
                           H5T_NATIVE_DOUBLE,
                           H5S_ALL,
                           H5S_ALL,
                           H5P_DEFAULT,
                           data);
    CAROM_VERIFY(errf >= 0);

    // Write attribute so we know what kind of data this is.
    writeAttribute(KEY_DOUBLE_ARRAY, dataset);

    errf = H5Sclose(space);
    CAROM_VERIFY(errf >= 0);

    errf = H5Dclose(dataset);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif
}

void
HDFDatabase::putDoubleVector(
    const std::string& key,
    const std::vector<double>& data,
    int nelements)
{
    putDoubleArray(key, data.data(), nelements);
}

void
HDFDatabase::getIntegerArray(
    const std::string& key,
    int* data,
    int nelements)
{
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

    hid_t dspace = H5Dget_space(dset);
    CAROM_VERIFY(dspace >= 0);

    hsize_t nsel = H5Sget_select_npoints(dspace);
    CAROM_VERIFY(static_cast<int>(nsel) == nelements);

    herr_t errf;
    if (nsel > 0) {
        errf = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        CAROM_VERIFY(errf >= 0);
    }

    errf = H5Sclose(dspace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Dclose(dset);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif
}

int HDFDatabase::getDoubleArraySize(const std::string& key)
{
    CAROM_VERIFY(!key.empty());

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
    hid_t dset = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
    hid_t dset = H5Dopen(d_group_id, key.c_str());
#endif
    CAROM_VERIFY(dset >= 0);

    hid_t dspace = H5Dget_space(dset);
    CAROM_VERIFY(dspace >= 0);

    hsize_t nsel = H5Sget_select_npoints(dspace);

    herr_t errf = H5Sclose(dspace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Dclose(dset);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif

    return static_cast<int>(nsel);
}

void
HDFDatabase::getDoubleArray(
    const std::string& key,
    double* data,
    int nelements)
{
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

    hid_t dspace = H5Dget_space(dset);
    CAROM_VERIFY(dspace >= 0);

    hsize_t nsel = H5Sget_select_npoints(dspace);
    CAROM_VERIFY(static_cast<int>(nsel) == nelements);

    herr_t errf;
    if (nsel > 0) {
        errf = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        CAROM_VERIFY(errf >= 0);
    }

    errf = H5Sclose(dspace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Dclose(dset);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif
}

void
HDFDatabase::getDoubleArray(
    const std::string& key,
    double* data,
    int nelements,
    const std::vector<int>& idx)
{
    if (idx.size() == 0)
    {
        getDoubleArray(key, data, nelements);
    }
    else
    {
        std::vector<double> alldata(nelements);
        getDoubleArray(key, alldata.data(), nelements);
        int k = 0;
        for (int i = 0; i < nelements; ++i)
        {
            if (idx[k] == i)
            {
                data[k++] = alldata[i];
            }
            if (k == idx.size())
            {
                break;
            }
        }
        CAROM_VERIFY(k == idx.size());
    }
}

void
HDFDatabase::getDoubleArray(
    const std::string& key,
    double* data,
    int nelements,
    int offset,
    int block_size,
    int stride)
{
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

    hid_t dspace = H5Dget_space(dset);
    CAROM_VERIFY(dspace >= 0);

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
        CAROM_VERIFY(errf >= 0);
    }

    errf = H5Sclose(dspace);
    CAROM_VERIFY(errf >= 0);

    errf = H5Dclose(dset);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif
}

bool
HDFDatabase::isInteger(
    const std::string& key)
{
    bool is_int = false;

    if (!key.empty()) {
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
        hid_t this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
        hid_t this_set = H5Dopen(d_group_id, key.c_str());
#endif
        if (this_set > 0) {
            int type_key = readAttribute(this_set);
            if (type_key == KEY_INT_ARRAY) {
                is_int = true;
            }
            herr_t errf = H5Dclose(this_set);
            CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
            CAROM_NULL_USE(errf);
#endif
        }
    }

    return is_int;
}

bool
HDFDatabase::isDouble(
    const std::string& key)
{
    bool is_double = false;

    if (!key.empty()) {
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
        hid_t this_set = H5Dopen(d_group_id, key.c_str(), H5P_DEFAULT);
#else
        hid_t this_set = H5Dopen(d_group_id, key.c_str());
#endif
        if (this_set > 0) {
            int type_key = readAttribute(this_set);
            if (type_key == KEY_DOUBLE_ARRAY) {
                is_double = true;
            }
            herr_t errf = H5Dclose(this_set);
            CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
            CAROM_NULL_USE(errf);
#endif
        }
    }

    return is_double;
}

void
HDFDatabase::writeAttribute(
    int type_key,
    hid_t dataset_id)
{
    hid_t attr_id = H5Screate(H5S_SCALAR);
    CAROM_VERIFY(attr_id >= 0);

#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 6))
    hid_t attr = H5Acreate(dataset_id,
                           "Type",
                           H5T_STD_I8BE,
                           attr_id,
                           H5P_DEFAULT,
                           H5P_DEFAULT);
#else
    hid_t attr = H5Acreate(dataset_id,
                           "Type",
                           H5T_STD_I8BE,
                           attr_id,
                           H5P_DEFAULT);
#endif
    CAROM_VERIFY(attr >= 0);

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

int
HDFDatabase::readAttribute(
    hid_t dataset_id)
{
    hid_t attr = H5Aopen_name(dataset_id, "Type");
    CAROM_VERIFY(attr >= 0);

    int type_key;
    herr_t errf = H5Aread(attr, H5T_NATIVE_INT, &type_key);
    CAROM_VERIFY(errf >= 0);

    errf = H5Aclose(attr);
    CAROM_VERIFY(errf >= 0);
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(errf);
#endif

    return type_key;
}

}
