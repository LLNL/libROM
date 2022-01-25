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
#include <vector>

namespace CAROM {

CSVDatabase::CSVDatabase()
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
    d_fs.close();
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
    d_fs.close();
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
    for (int i = 0; i < nelements; ++i)
    {
        d_fs >> data_entry;
        data[i] = data_entry;
    }
    CAROM_ASSERT(d_fs.eof());
    d_fs.close();
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
    double data_entry = 0.0;
    for (int i = 0; i < nelements; ++i)
    {
        d_fs >> data_entry;
        data[i] = data_entry;
    }
    CAROM_ASSERT(d_fs.eof());
    d_fs.close();
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

    std::ifstream d_fs(key.c_str());
    std::string line, data_entry;
    int count = 0;
    while (count < nelements && d_fs >> line)
    {
        std::stringstream d_ss(line);
        while (std::getline(d_ss, data_entry, ','))
        {
            if (offset-- > stride)
            {
                data[count++] = std::stod(data_entry); 
            }
            if (offset == 0)
            {
                offset = stride + block_size;
            }
        }
    }
    d_fs.close();
}

}
