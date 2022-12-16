/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Implementation of DMDData class.

#include "DMDData.h"

#ifndef _WIN32
#include <sys/stat.h>  // mkdir
#else
#include <direct.h>    // _mkdir
#define mkdir(dir, mode) _mkdir(dir)
#endif

#include "mpi.h"

#include "utils/Utilities.h"
#include "utils/CSVDatabase.h"
#include "utils/HDFDatabase.h"

namespace CAROM {

DMDData::DMDData(
    const char* output_path,
    DataFormat data_format
)
:   d_output_path { output_path },
    d_data_format { data_format }
{
    switch (data_format)
    {
        case csv:
            d_db = std::make_unique<CSVDatabase>();
            mkdir(output_path, 0777);
            break;
        case hdf5:
            d_db = std::make_unique<HDFDatabase>();
            d_db->create(d_output_path + ".hdf");
            break;
    }
}

DMDData::DMDData(DMDData&& other)
{
    moveData(std::move(other));
}

DMDData& DMDData::operator=(DMDData&& rhs)
{
    moveData(std::move(rhs));
    return *this;
}

DMDData::~DMDData()
{
    if (!d_steps.empty())
    {
        std::string output_prefix;
        std::string output_postfix;
        if (d_data_format == csv)
        {
            output_prefix = d_output_path + "/";
            output_postfix = ".csv";
        }
        d_db->putDoubleVector(output_prefix + "tval" + output_postfix, d_t_vals, 
            d_t_vals.size());
        d_db->putInteger(output_prefix + "numsnap", d_steps.size());
        if (d_data_format == hdf5)
        {
            d_db->putInteger("snap_bound_size", 0);
        }
        d_db->putIntegerArray(output_prefix + "snap_list" + output_postfix, 
            d_steps.data(), d_steps.size());
    }
}

void DMDData::addSnapshot(double t, int step, const double *data, int size)
{
    std::string output_file;
    if (d_data_format == csv)
    {
        output_file = d_output_path + "/step" + std::to_string(step);
        mkdir(output_file.c_str(), 0777);
        output_file = output_file + "/sol.csv";
    }
    else // data_format == hdf5
    {
        output_file = "step" + std::to_string(step) + "sol";
    }
    d_db->putDoubleArray(output_file, data, size);
    d_t_vals.push_back(t);
    d_steps.push_back(step);
}

void DMDData::moveData(DMDData&& other)
{
    d_output_path = std::move(other.d_output_path);
    d_data_format = other.d_data_format;
    d_db = std::move(other.d_db);
    d_t_vals = std::move(other.d_t_vals);
    d_steps = std::move(other.d_steps);
}

}
