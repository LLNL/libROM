/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
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
#include <complex>

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
     * @brief Writes an array of integers associated with the supplied filename.
     *
     * @pre !file_name.empty()
     * @pre data != nullptr
     * @pre nelements > 0
     *
     * @param[in] file_name The filename associated with the array of values to be
     *                written.
     * @param[in] data The array of integer values to be written.
     * @param[in] nelements The number of integers in the array.
     */
    virtual
    void
    putIntegerArray(
        const std::string& file_name,
        const int* const data,
        int nelements);

    /**
     * @brief Writes an array of doubles associated with the supplied filename.
     *
     * @pre !file_name.empty()
     * @pre data != nullptr
     * @pre nelements > 0
     *
     * @param[in] file_name The filename associated with the array of values to be
     *                written.
     * @param[in] data The array of double values to be written.
     * @param[in] nelements The number of doubles in the array.
     */
    virtual
    void
    putDoubleArray(
        const std::string& file_name,
        const double* const data,
        int nelements);

    /**
     * @brief Writes a vector of doubles associated with the supplied filename to
     * the currently open CSV database file.
     *
     * @pre !file_name.empty()
     * @pre data != nullptr
     * @pre nelements > 0
     *
     * @param[in] file_name The filename associated with the vector of values to be
     *                written.
     * @param[in] data The vector of double values to be written.
     * @param[in] nelements The number of doubles in the vector.
     */
    virtual
    void
    putDoubleVector(
        const std::string& file_name,
        const std::vector<double>& data,
        int nelements);

    /**
     * @brief Writes a vector of complex doubles associated with the supplied filename.
     *
     * @pre !file_name.empty()
     * @pre data != nullptr
     * @pre nelements > 0
     *
     * @param[in] file_name The filename associated with the vector of values to be
     *                written.
     * @param[in] data The vector of complex double values to be written.
     * @param[in] nelements The number of complex doubles in the vector.
     */
    virtual
    void
    putComplexVector(
        const std::string& file_name,
        const std::vector<std::complex<double>>& data,
        int nelements);

    /**
     * @brief Writes a vector of strings associated with the supplied filename.
     *
     * @pre !file_name.empty()
     * @pre data != nullptr
     * @pre nelements > 0
     *
     * @param[in] file_name The filename associated with the vector of values to be
     *                written.
     * @param[in] data The vector of strings to be written.
     * @param[in] nelements The number of strings in the vector.
     */
    virtual
    void
    putStringVector(
        const std::string& file_name,
        const std::vector<std::string>& data,
        int nelements);

    /**
     * @brief Reads an array of integers associated with the supplied filename.
     *
     * @pre !file_name.empty()
     * @pre data != nullptr || nelements == 0
     *
     * @param[in] file_name The filename associated with the array of values to be
     *                read.
     * @param[out] data The allocated array of integer values to be read.
     * @param[in] nelements The number of integers in the array.
     */
    virtual
    void
    getIntegerArray(
        const std::string& file_name,
        int* data,
        int nelements);

    /**
     * @brief Reads a vector of integers associated with the supplied filename.
     *
     * @pre !file_name.empty()
     *
     * @param[in] file_name The filename associated with the vector of values to be
     *                read.
     * @param[out] data The allocated vector of integer values to be read.
     * @param[in] append If true, append to the vector, otherwise overwrite.
     */
    void
    getIntegerVector(
        const std::string& file_name,
        std::vector<int> &data,
        bool append = false);

    /**
     * @brief Count the number of elements in
     * an array of doubles associated with the supplied filename.
     *
     * @pre !file_name.empty()
     *
     * @param[in] file_name The filename associated with the array of values to be
     *                read.
     */
    virtual
    int
    getDoubleArraySize(const std::string& file_name)
    {
        std::vector<double> tmp;
        getDoubleVector(file_name, tmp);
        return tmp.size();
    }

    /**
     * @brief Reads an array of doubles associated with the supplied filename.
     *
     * @pre !file_name.empty()
     * @pre data != nullptr || nelements == 0
     *
     * @param[in] file_name The filename associated with the array of values to be
     *                read.
     * @param[out] data The allocated array of double values to be read.
     * @param[in] nelements The number of doubles in the array.
     */
    virtual
    void
    getDoubleArray(
        const std::string& file_name,
        double* data,
        int nelements);

    /**
     * @brief Reads a sub-array of doubles associated with the supplied filename.
     *
     * @pre !file_name.empty()
     * @pre data != nullptr || nelements == 0
     *
     * @param[in] file_name The filename associated with the array of values to be
     *                read.
     * @param[out] data The allocated sub-array of double values to be read.
     * @param[in] nelements The number of doubles in the full array.
     * @param[in] idx The set of indices in the sub-array.
     */
    virtual
    void
    getDoubleArray(
        const std::string& file_name,
        double* data,
        int nelements,
        const std::vector<int>& idx);

    /**
     * @brief Reads an array of doubles associated with the supplied filename.
     *
     * @pre !file_name.empty()
     * @pre data != nullptr || nelements == 0
     *
     * @param[in] file_name The filename associated with the array of values to be
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
        const std::string& file_name,
        double* data,
        int nelements,
        int offset,
        int block_size,
        int stride);

    /**
     * @brief Reads a vector of doubles associated with the supplied filename.
     *
     * @pre !file_name.empty()
     *
     * @param[in] file_name The filename associated with the vector of values to be
     *                read.
     * @param[out] data The allocated vector of double values to be read.
     * @param[in] append If true, append to the vector, otherwise overwrite.
     */
    void
    getDoubleVector(
        const std::string& file_name,
        std::vector<double> &data,
        bool append = false);

    /**
     * @brief Reads a vector of strings associated with the supplied filename.
     *
     * @pre !file_name.empty()
     *
     * @param[in] file_name The filename associated with the vector of strings to be
     *                read.
     * @param[out] data The allocated vector of strings to be read.
     * @param[in] append If true, append to the vector, otherwise overwrite.
     */
    void
    getStringVector(
        const std::string& file_name,
        std::vector<std::string> &data,
        bool append = false);

    /**
     * @brief Count the number of lines of CSV database file.
     *
     * @pre !file_name.empty()
     *
     * @param[in] file_name The filename associated with the vector of strings to be
     *                read.
     */
    int
    getLineCount(
        const std::string& file_name);

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
