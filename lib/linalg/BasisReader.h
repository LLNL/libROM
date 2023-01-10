/******************************************************************************
 *
 * Copyright (c) 2013-2023, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: A class that reads basis vectors from a file.

#ifndef included_BasisReader_h
#define included_BasisReader_h

#include "utils/Utilities.h"
#include "utils/Database.h"
#include <string>
#include <vector>

namespace CAROM {

class Matrix;
class Vector;
class Database;

/**
 * Class BasisReader reads the basis vectors from a file written by class
 * BasisWriter.
 *
 * @see BasisWriter
 */
class BasisReader {
public:
    /**
     * @brief The constructor for the BasisReader class takes the base part
     * of the name of the files holding the basis vectors and the file
     * format.
     *
     * @pre !base_file_name.empty()
     *
     * @param[in] base_file_name The base part of the name of the files
     *                           holding the basis vectors.
     * @param[in] db_format Format of the file to read.
     *                      One of the implemented file formats defined in
     *                      Database.
     */
    BasisReader(
        const std::string& base_file_name,
        Database::formats db_format = Database::HDF5);

    /**
     * @brief Destructor.
     */
    ~BasisReader();

    /**
     * @brief Returns true if the basis vectors at requested time are
     * different from the last requested basis vectors.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     *
     * @param[in] time Time at which we are interested in the basis vectors.
     *
     * @return True if the basis vectors at the requested time are different
     *         from the last requested basis vectors.
     */
    bool
    isNewBasis(
        double time)
    {
        CAROM_VERIFY(0 < numTimeIntervals());
        CAROM_VERIFY(0 <= time);
        bool result = false;
        if (d_last_basis_idx == -1) {
            result = true;
        }
        else {
            int num_time_intervals = numTimeIntervals();
            int i;
            for (i = 0; i < num_time_intervals-1; ++i) {
                if (d_time_interval_start_times[i] <= time &&
                        time < d_time_interval_start_times[i+1]) {
                    break;
                }
            }
            result = i != d_last_basis_idx;
        }
        return result;
    }

    /**
     *
     * @brief Returns the spatial basis vectors for the requested time as a
     *        Matrix.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     *
     * @param[in] time Time for which we want the basis vectors.
     *
     * @return The spatial basis vectors for the requested time.
     */
    Matrix*
    getSpatialBasis(
        double time);

    /**
     *
     * @brief Returns the first n spatial basis vectors for the requested time
     *        as a Matrix.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     * @pre 0 < n <= numColumns()
     *
     * @param[in] time Time for which we want the basis vectors.
     * @param[in] n    The number of spatial basis vectors desired.
     *
     * @return The spatial basis vectors for the requested time.
     */
    Matrix*
    getSpatialBasis(
        double time,
        int n);

    /**
     *
     * @brief Returns spatial basis vectors from start_col to end_col for the
     *        requested time as a Matrix.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     * @pre 0 < start_col <= numColumns()
     * @pre start_col <= end_col <= numColumns()
     *
     * @param[in] time         Time for which we want the basis vectors.
     * @param[in] start_col    The starting column desired.
     * @param[in] end_col      The starting column desired.
     *
     * @return The spatial basis vectors for the requested time.
     */
    Matrix*
    getSpatialBasis(
        double time,
        int start_col,
        int end_col);

    /**
     *
     * @brief Returns the first n spatial basis vectors for the requested time
     *        as a Matrix that capture the given energy fraction.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     * @pre 0 <= ef <= 1.0
     *
     * @param[in] time Time for which we want the basis vectors.
     * @param[in] ef   The desired energy fraction.
     *
     * @return The spatial basis vectors for the requested time.
     */
    Matrix*
    getSpatialBasis(
        double time,
        double ef);

    /**
     *
     * @brief Returns the temporal basis vectors for the requested time as
     *        a Matrix.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     *
     * @param[in] time Time for which we want the basis vectors.
     *
     * @return The temporal basis vectors for the requested time.
     */
    Matrix*
    getTemporalBasis(
        double time);

    /**
     *
     * @brief Returns the first n temporal basis vectors for the requested time
     *        as a Matrix.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     * @pre 0 < n <= numColumns()
     *
     * @param[in] time Time for which we want the basis vectors.
     * @param[in] n    The number of temporal basis vectors desired.
     *
     * @return The temporal basis vectors for the requested time.
     */
    Matrix*
    getTemporalBasis(
        double time,
        int n);

    /**
     *
     * @brief Returns temporal basis vectors from start_col to end_col for the
     *        requested time as a Matrix.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     * @pre 0 < start_col <= numColumns()
     * @pre start_col <= end_col <= numColumns()
     *
     * @param[in] time         Time for which we want the basis vectors.
     * @param[in] start_col    The starting column desired.
     * @param[in] end_col      The starting column desired.
     *
     * @return The temporal basis vectors for the requested time.
     */
    Matrix*
    getTemporalBasis(
        double time,
        int start_col,
        int end_col);

    /**
     *
     * @brief Returns the first n temporal basis vectors for the requested time
     *        as a Matrix that capture the given energy fraction.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     * @pre 0 <= ef <= 1.0
     *
     * @param[in] time Time for which we want the basis vectors.
     * @param[in] ef   The desired energy fraction.
     *
     * @return The temporal basis vectors for the requested time.
     */
    Matrix*
    getTemporalBasis(
        double time,
        double ef);

    /**
     *
     * @brief Returns the singular values for the requested time.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     *
     * @param[in] time Time for which we want the basis vectors.
     *
     * @return The temporal basis vectors for the requested time.
     */
    Vector*
    getSingularValues(
        double time);

    /**
     *
     * @brief Returns the largest singular values for the requested time
     *        that capture the given energy fraction.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     * @pre 0 <= ef <= 1.0
     *
     * @param[in] time Time for which we want the basis vectors.
     * @param[in] ef   The desired energy fraction.
     *
     * @return The temporal basis vectors for the requested time.
     */
    Vector*
    getSingularValues(
        double time,
        double ef);

    /**
     *
     * @brief Returns the dimension of the system on this processor.
     *
     * @return The dimension of the system on this processor.
     */
    int
    getDim(
        const std::string kind,
        double time);

    /**
     *
     * @brief Returns the number of samples (columns) in file.
     *
     * @return The number of samples in file.
     */
    int
    getNumSamples(
        const std::string kind,
        double time);

    /**
     *
     * @brief Returns the snapshot matrix for the requested time.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     *
     * @param[in] time Time for which we want the basis vectors.
     *
     * @return The snapshot matrix for the requested time.
     */
    Matrix*
    getSnapshotMatrix(
        double time);

    /**
     *
     * @brief Returns the first n columns of the snapshot matrix for the requested time.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     * @pre 0 < n <= numColumns()
     *
     * @param[in] time Time for which we want the basis vectors.
     * @param[in] n    The number of basis vectors desired.
     *
     * @return The snapshot matrix for the requested time.
     */
    Matrix*
    getSnapshotMatrix(
        double time,
        int n);

    /**
     *
     * @brief Returns the snapshot matrix from start_col to end_col for the requested time.
     *
     * @pre 0 < numTimeIntervals()
     * @pre 0 <= time
     * @pre 0 < start_col <= numColumns()
     * @pre start_col <= end_col <= numColumns()
     *
     * @param[in] time Time for which we want the basis vectors.
     * @param[in] start_col    The starting column desired.
     * @param[in] end_col      The starting column desired.
     *
     * @return The snapshot matrix for the requested time.
     */
    Matrix*
    getSnapshotMatrix(
        double time,
        int start_col,
        int end_col);

private:
    /**
     * @brief Unimplemented default constructor.
     */
    BasisReader();

    /**
     * @brief Unimplemented copy constructor.
     */
    BasisReader(
        const BasisReader& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    BasisReader&
    operator = (
        const BasisReader& rhs);

    /**
     * @brief Number of time intervals.
     *
     * @return The number of time intervals.
     */
    int
    numTimeIntervals()
    {
        return static_cast<int>(d_time_interval_start_times.size());
    }

    /**
     * @brief The start time of each time interval.
     */
    std::vector<double> d_time_interval_start_times;

    /**
     * @brief The database being read from.
     */
    Database* d_database;

    /**
     * @brief Base file name stored for consistency between reading and writing.
     */
    std::string base_file_name_;

    /**
     * @brief Full file name of database including rank.
     */
    std::string full_file_name;

    /**
     * @brief The last time at which basis vectors were requested.
     */
    int d_last_basis_idx;
};

}

#endif
