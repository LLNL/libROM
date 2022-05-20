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
#include <complex>

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
    const std::string& file_name,
    const int* const data,
    int nelements)
{
    CAROM_ASSERT(!file_name.empty());
    CAROM_ASSERT(data != 0);
    CAROM_ASSERT(nelements > 0);

    std::ofstream d_fs(file_name.c_str());
    for (int i = 0; i < nelements; ++i)
    {
        d_fs << data[i] << std::endl;
    }
    d_fs.close();
}

void
CSVDatabase::putDoubleArray(
    const std::string& file_name,
    const double* const data,
    int nelements)
{
    CAROM_ASSERT(!file_name.empty());
    CAROM_ASSERT(data != 0);
    CAROM_ASSERT(nelements > 0);

    std::ofstream d_fs(file_name.c_str());
    for (int i = 0; i < nelements; ++i)
    {
        d_fs << data[i] << std::endl;
    }
    d_fs.close();
}

void
CSVDatabase::putDoubleVector(
    const std::string& file_name,
    const std::vector<double> data,
    int nelements)
{
    CAROM_ASSERT(!file_name.empty());
    CAROM_ASSERT(data != 0);
    CAROM_ASSERT(nelements > 0);

    std::ofstream d_fs(file_name.c_str());
    for (int i = 0; i < nelements; ++i)
    {
        d_fs << data[i] << std::endl;
    }
    d_fs.close();
}

void
CSVDatabase::putComplexVector(
    const std::string& file_name,
    const std::vector<std::complex<double>> data,
    int nelements)
{
    CAROM_ASSERT(!file_name.empty());
    CAROM_ASSERT(data != 0);
    CAROM_ASSERT(nelements > 0);

    std::ofstream d_fs(file_name.c_str());
    for (int i = 0; i < nelements; ++i)
    {
        d_fs << std::real(data[i]) << "," << std::imag(data[i]) << std::endl;
    }
    d_fs.close();
}

void
CSVDatabase::putStringVector(
    const std::string& file_name,
    const std::vector<std::string> data,
    int nelements)
{
    CAROM_ASSERT(!file_name.empty());
    CAROM_ASSERT(data != 0);
    CAROM_ASSERT(nelements > 0);

    std::ofstream d_fs(file_name.c_str());
    for (int i = 0; i < nelements; ++i)
    {
        d_fs << data[i] << std::endl;
    }
    d_fs.close();
}

void
CSVDatabase::getIntegerArray(
    const std::string& file_name,
    int* data,
    int nelements)
{
    CAROM_ASSERT(!file_name.empty());
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(nelements);
#endif

    std::ifstream d_fs(file_name.c_str());
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
CSVDatabase::getIntegerVector(
    const std::string& file_name,
    std::vector<int> &data,
    bool append)
{
    CAROM_ASSERT(!file_name.empty());
    if (!append) data.clear();

    std::ifstream d_fs(file_name.c_str());
    int data_entry;
    while (d_fs >> data_entry)
    {
        data.push_back(data_entry);
    }
    d_fs.close();
}

void
CSVDatabase::getDoubleArray(
    const std::string& file_name,
    double* data,
    int nelements)
{
    CAROM_ASSERT(!file_name.empty());
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(nelements);
#endif

    std::ifstream d_fs(file_name.c_str());
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
    const std::string& file_name,
    double* data,
    int nelements,
    std::vector<int> idx)
{
    CAROM_ASSERT(!file_name.empty());
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(nelements);
#endif

    if (idx.size() == 0)
    {
        getDoubleArray(file_name, data, nelements);
    }
    else
    {
        std::ifstream d_fs(file_name.c_str());
        int k = 0;
        double data_entry = 0.0;
        for (int i = 0; i < nelements; ++i)
        {
            d_fs >> data_entry;
            if (idx[k] == i)
            {
                data[k++] = data_entry;
            }
        }
        CAROM_ASSERT(d_fs.eof());
        CAROM_ASSERT(k == idx.size());
        d_fs.close();
    }
}

void
CSVDatabase::getDoubleArray(
    const std::string& file_name,
    double* data,
    int nelements,
    int offset,
    int block_size,
    int stride)
{
    CAROM_ASSERT(!file_name.empty());
#ifndef DEBUG_CHECK_ASSERTIONS
    CAROM_NULL_USE(nelements);
#endif

    std::ifstream d_fs(file_name.c_str());
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

void
CSVDatabase::getDoubleVector(
    const std::string& file_name,
    std::vector<double> &data,
    bool append)
{
    CAROM_ASSERT(!file_name.empty());
    if (!append) data.clear();

    std::ifstream d_fs(file_name.c_str());
    double data_entry;
    while (d_fs >> data_entry)
    {
        data.push_back(data_entry);
    }
    d_fs.close();
}

void
CSVDatabase::getStringVector(
    const std::string& file_name,
    std::vector<std::string> &data,
    bool append)
{
    CAROM_ASSERT(!file_name.empty());
    if (!append) data.clear();

    std::ifstream d_fs(file_name.c_str());
    std::string data_entry;
    while (std::getline(d_fs, data_entry))
    {
        data.push_back(data_entry);
    }
    d_fs.close();
}

int
CSVDatabase::getLineCount(
    const std::string& file_name)
{
    int count = 0;
    CAROM_ASSERT(!file_name.empty());
    std::ifstream d_fs(file_name.c_str());
    std::string data_entry;
    while (std::getline(d_fs, data_entry))
    {
        count += 1;
    }
    return count;
}

}
