/******************************************************************************
 *
 * Copyright (c) 2013-2022, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: Creates HDF5 data for training a parametric DMD with time
//              windowing.

#ifndef included_DMDData_h
#define included_DMDData_h

#include <memory>
#include <string>
#include <vector>

namespace CAROM {

class Database;

/**
 * Class DMDData implements methods to save training data for a time or distance
 * windowed (possibly parametric) DMD using local_dw_csv, local_tw_csv,
 * parametric_dw_csv, or parametric_tw_csv
 */
class DMDData
{
public:

    enum DMDType {parametric, local};
    enum DataFormat {hdf5, csv};

    /**
     * @brief Constructor.
     *
     * @param[in] output_path  Directory to store output.
     * @param[in] data_format  Output format
     */
    DMDData(const char* output_path,
        DataFormat data_format = csv
    );

    /**
     * @brief Copy constructor (deleted)
     */
    DMDData(const DMDData& other) = delete;

    /**
     * @brief Move constructor
     * 
     * @param other The DMDData object to move
     */
    DMDData(DMDData&& other);

    /**
     * @brief Copy assignment operator (deleted)
     */
    DMDData& operator=(const DMDData& other) = delete;

    /**
     * @brief Move assignment operator
     * 
     * @param rhs The DMDData object to move
     * @return This after rhs has been moved
     */
    DMDData& operator=(DMDData&& rhs);

    /**
     * @brief Destroy the DMDData object and write summary data
     */
    ~DMDData();

    /**
     * @brief Adds a new snapshot to the output
     * 
     * @param t    Analysis time
     * @param step Analysis step
     * @param data Pointer to snapshot data array
     * @param size Size of snapshot data array
     */
    void addSnapshot(double t, int step, const double* data, int size);

private:

    /**
     * @brief Moves data from other to this
     * 
     * @param other DMDData object to move data from
     */
    void moveData(DMDData&& other);

    /**
     * @brief Relative path to store CSV output data
     */
    std::string d_output_path;

    /**
     * @brief Specifies HDF5 or CSV output
     */
    DataFormat d_data_format;

    /**
     * @brief Data writing object (base class object)
     */
    std::unique_ptr<CAROM::Database> d_db;

    /**
     * @brief Output times
     */
    std::vector<double> d_t_vals;

    /**
     * @brief Output steps
     */
    std::vector<int> d_steps;

};

}

#endif
