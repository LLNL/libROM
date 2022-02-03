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

#ifndef included_CSVDatabase_h
#define included_CSVDatabase_h

#include "Database.h"
#include <string>
#include <fstream>
#include <vector>

namespace CAROM {

/**
 * CSVDatabase implements the interface of Database for CSV database files.
 */
class CSVDatabase : public Database
{
public:
    /**
     * @brief Default constructor.
     */
    CSVDatabase();

    /**
     * @brief Destructor.
     */
    virtual
    ~CSVDatabase();

    /**
     * @brief Creates a new CSV database file with the supplied name.
     *
     * @param[in] file_name Name of CSV database file to create.
     *
     * @return True if file create was successful.
     */
    virtual
    bool
    create(
        const std::string& file_name);

    /**
     * @brief Opens an existing CSV database file with the supplied name.
     *
     * @param[in] file_name Name of existing CSV database file to open.
     * @param[in] type Read/write type ("r"/"wr")
     *
     * @return True if file open was successful.
     */
    virtual
    bool
    open(
        const std::string& file_name,
        const std::string& type);

    /**
     * @brief Closes the currently open CSV database file.
     *
     * @return True if the file close was successful.
     */
    virtual
    bool
    close();

    /**
     * @brief Writes an array of integers associated with the supplied key to
     * the currently open CSV database file.
     *
     * @pre !key.empty()
     * @pre data != 0
     * @pre nelements > 0
     *
     * @param[in] key The key associated with the array of values to be
     *                written.
     * @param[in] data The array of integer values to be written.
     * @param[in] nelements The number of integers in the array.
     */
    virtual
    void
    putIntegerArray(
        const std::string& key,
        const int* const data,
        int nelements);

    /**
     * @brief Writes an array of doubles associated with the supplied key to
     * the currently open CSV database file.
     *
     * @pre !key.empty()
     * @pre data != 0
     * @pre nelements > 0
     *
     * @param[in] key The key associated with the array of values to be
     *                written.
     * @param[in] data The array of double values to be written.
     * @param[in] nelements The number of doubles in the array.
     */
    virtual
    void
    putDoubleArray(
        const std::string& key,
        const double* const data,
        int nelements);

    /**
     * @brief Reads an array of integers associated with the supplied key
     * from the currently open CSV database file.
     *
     * @pre !key.empty()
     * @pre data != 0 || nelements == 0
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated array of integer values to be read.
     * @param[in] nelements The number of integers in the array.
     */
    virtual
    void
    getIntegerArray(
        const std::string& key,
        int* data,
        int nelements);

    /**
     * @brief Reads an array of integers associated with the supplied key
     * from the currently open CSV database file.
     *
     * @pre !key.empty()
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated vector of integer values to be read.
     * @param[in] append True if append to the list, otherwise overwite.
     */
    void
    getIntegerArray(
        const std::string& key,
        std::vector<int> &data, 
        bool append);

    /**
     * @brief Reads an array of doubles associated with the supplied key
     * from the currently open CSV database file.
     *
     * @pre !key.empty()
     * @pre data != 0 || nelements == 0
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated array of double values to be read.
     * @param[in] nelements The number of doubles in the array.
     */
    virtual
    void
    getDoubleArray(
        const std::string& key,
        double* data,
        int nelements);

    /**
     * @brief Reads a sub-array of doubles associated with the supplied key
     * from the currently open CSV database file.
     *
     * @pre !key.empty()
     * @pre data != 0 || nelements == 0
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated sub-array of double values to be read.
     * @param[in] nelements The number of doubles in the full array.
     * @param[in] idx The set of indices in the sub-array.
     */
    virtual
    void
    getDoubleArray(
        const std::string& key,
        double* data,
        int nelements, 
        std::vector<int> idx);

    /**
     * @brief Reads an array of doubles associated with the supplied key
     * from the currently open CSV database file.
     *
     * @pre !key.empty()
     * @pre data != 0 || nelements == 0
     *
     * @param[in] key The key associated with the array of values to be
     *                read.
     * @param[out] data The allocated array of double values to be read.
     * @param[in] nelements The number of doubles in the array.
     * @param[in] offset The initial offset in the array.
     * @param[in] block_size The block size to read from the CSV dataset.
     * @param[in] stride The stride to read from the CSV dataset.
     */
    virtual
    void
    getDoubleArray(
        const std::string& key,
        double* data,
        int nelements,
        int offset,
        int block_size,
        int stride);

    /**
     * @brief Reads a list of strings associated with the supplied key
     * from the currently open CSV database file.
     *
     * @pre !key.empty()
     *
     * @param[in] key The key associated with the list of strings to be
     *                read.
     * @param[out] data The allocated vector of strings to be read.
     * @param[in] append True if append to the list, otherwise overwite.
     */
    void
    getStringList(
        const std::string& key,
        std::vector<std::string> &data, 
        bool append);

private:
    /**
     * @brief Unimplemented copy constructor.
     */
    CSVDatabase(
        const CSVDatabase& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    CSVDatabase&
    operator = (
        const CSVDatabase& rhs);

};

}

#endif
