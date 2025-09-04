/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
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
#include <memory>

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
     * @param[in] dim Number of rows of basis that will be read from a file.
     *                If negative, will use the dimension from the rank-specific local file.
     */
    BasisReader(
        const std::string& base_file_name,
        Database::formats db_format = Database::formats::HDF5,
        const int dim = -1,
        MPI_Comm comm = MPI_COMM_WORLD);

    /**
     * @brief Destructor.
     */
    ~BasisReader();

    /**
     *
     * @brief Returns the spatial basis vectors as a Matrix.
     *
     * @return The spatial basis vectors.
     */
    std::unique_ptr<Matrix>
    getSpatialBasis();

    /**
     *
     * @brief Returns the first n spatial basis vectors for the requested time
     *        as a Matrix.
     *
     * @pre 0 < n <= numColumns()
     *
     * @param[in] n    The number of spatial basis vectors desired.
     *
     * @return The spatial basis vectors for the requested time.
     */
    std::unique_ptr<Matrix>
    getSpatialBasis(
        int n);

    /**
     *
     * @brief Returns spatial basis vectors from start_col to end_col for the
     *        requested time as a Matrix.
     *
     * @pre 0 < start_col <= numColumns()
     * @pre start_col <= end_col <= numColumns()
     *
     * @param[in] start_col    The starting column desired.
     * @param[in] end_col      The final column desired.
     *
     * @return The spatial basis vectors for the requested time.
     */
    std::unique_ptr<Matrix>
    getSpatialBasis(
        int start_col,
        int end_col);

    /**
     *
     * @brief Returns the first n spatial basis vectors for the requested time
     *        as a Matrix that capture the given energy fraction.
     *
     * @pre 0 <= ef <= 1.0
     *
     * @param[in] ef   The desired energy fraction.
     *
     * @return The spatial basis vectors for the requested time.
     */
    std::unique_ptr<Matrix>
    getSpatialBasis(
        double ef);

    /**
     *
     * @brief Returns the temporal basis vectors for the requested time as
     *        a Matrix.
     *
     * @return The temporal basis vectors for the requested time.
     */
    std::unique_ptr<Matrix>
    getTemporalBasis();

    /**
     *
     * @brief Returns the first n temporal basis vectors for the requested time
     *        as a Matrix.
     *
     * @pre 0 < n <= numColumns()
     *
     * @param[in] n    The number of temporal basis vectors desired.
     *
     * @return The temporal basis vectors for the requested time.
     */
    std::unique_ptr<Matrix>
    getTemporalBasis(
        int n);

    /**
     *
     * @brief Returns temporal basis vectors from start_col to end_col for the
     *        requested time as a Matrix.
     *
     * @pre 0 < start_col <= numColumns()
     * @pre start_col <= end_col <= numColumns()
     *
     * @param[in] start_col    The starting column desired.
     * @param[in] end_col      The final column desired.
     *
     * @return The temporal basis vectors for the requested time.
     */
    std::unique_ptr<Matrix>
    getTemporalBasis(
        int start_col,
        int end_col);

    /**
     *
     * @brief Returns the first n temporal basis vectors for the requested time
     *        as a Matrix that capture the given energy fraction.
     *
     * @pre 0 <= ef <= 1.0
     *
     * @param[in] ef   The desired energy fraction.
     *
     * @return The temporal basis vectors for the requested time.
     */
    std::unique_ptr<Matrix>
    getTemporalBasis(
        double ef);

    /**
     *
     * @brief Returns the singular values for the requested time.
     *
     * @return The temporal basis vectors for the requested time.
     */
    std::unique_ptr<Vector>
    getSingularValues();

    /**
     *
     * @brief Returns the largest singular values for the requested time
     *        that capture the given energy fraction.
     *
     * @pre 0 <= ef <= 1.0
     *
     * @param[in] ef   The desired energy fraction.
     *
     * @return The temporal basis vectors for the requested time.
     */
    std::unique_ptr<Vector>
    getSingularValues(
        double ef);

    /**
     *
     * @brief Returns the dimension of the system on this processor.
     *
     * @param[in] kind  Type of matrix whose row-dimension is returned.
     *                  "basis" - local dimension of spatial basis
     *                  "snapshot" - local dimension of snapshot matrix
     *                  "temporal_basis" - global dimension of temporal basis
     *
     * @return The dimension of the system on this processor.
     */
    int
    getDim(
        const std::string kind);

    /**
     *
     * @brief Returns the number of samples (columns) in file.
     *
     * @return The number of samples in file.
     */
    int
    getNumSamples(
        const std::string kind);

    /**
     *
     * @brief Returns the snapshot matrix for the requested time.
     *
     * @return The snapshot matrix for the requested time.
     */
    std::unique_ptr<Matrix>
    getSnapshotMatrix();

    /**
     *
     * @brief Returns the first n columns of the snapshot matrix for the requested time.
     *
     * @pre 0 < n <= numColumns()
     *
     * @param[in] n    The number of basis vectors desired.
     *
     * @return The snapshot matrix for the requested time.
     */
    std::unique_ptr<Matrix>
    getSnapshotMatrix(
        int n);

    /**
     *
     * @brief Returns the snapshot matrix from start_col to end_col for the requested time.
     *
     * @pre 0 < start_col <= numColumns()
     * @pre start_col <= end_col <= numColumns()
     *
     * @param[in] start_col    The starting column desired.
     * @param[in] end_col      The final column desired.
     *
     * @return The snapshot matrix for the requested time.
     */
    std::unique_ptr<Matrix>
    getSnapshotMatrix(
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
     * @brief The database being read from.
     */
    Database* d_database;

    /**
     * @brief The database being read from.
     */
    Database::formats d_format;

    /**
     * @brief Base file name stored for consistency between reading and writing.
     */
    std::string base_file_name_;

    /**
     * @brief Full file name of database including rank.
     */
    std::string full_file_name;

    /**
     * @brief Dimension of the basis on this processor.
     *
     * If negative, use the dimension from the rank-specific local file.
     */
    const int d_dim;

    /**
     * @brief Dimension of the basis on this processor.
     *
     * If negative, use the dimension from the rank-specific local file.
     */
    int d_global_dim;

    bool d_distributed;
};

}

#endif
