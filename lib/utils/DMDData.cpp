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
:   output_path_ { output_path },
    data_format_ { data_format }
{
    switch (data_format)
    {
        case csv:
            db_ = std::make_unique<CSVDatabase>();
            break;
        case hdf5:
            db_ = std::make_unique<HDFDatabase>();
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
    std::string output_prefix;
    std::string output_postfix;
    if (data_format_ == csv)
    {
        output_prefix = output_path_ + "/";
        output_postfix = ".csv";
    }
    db_->putDoubleVector(output_prefix + "tval" + output_postfix, t_vals_, 
        t_vals_.size());
    db_->putInteger(output_prefix + "numsnap", steps_.size());
    if (data_format_ == hdf5)
    {
        db_->putInteger("snap_bound_size", 0);
    }
    db_->putIntegerArray(output_prefix + "snap_list" + output_postfix, 
        steps_.data(), steps_.size());
}

void DMDData::addSnapshot(double t, int step, const double *data, int size)
{
    std::string output_file;
    if (data_format_ == csv)
    {
        output_file = output_path_ + "/step" + std::to_string(step);
        mkdir(output_file.c_str(), 0777);
        output_file = output_file + "/sol.csv";
    }
    else // data_format == hdf5
    {
        output_file = "step" + std::to_string(step) + "sol";
    }
    db_->putDoubleArray(output_file, data, size);
    t_vals_.push_back(t);
    steps_.push_back(step);
}

void DMDData::moveData(DMDData&& other)
{
    output_path_ = std::move(other.output_path_);
    data_format_ = other.data_format_;
    db_ = std::move(other.db_);
    t_vals_ = std::move(other.t_vals_);
    steps_ = std::move(other.steps_);
}

}
