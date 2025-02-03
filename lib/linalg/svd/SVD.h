/******************************************************************************
 *
 * Copyright (c) 2013-2024, Lawrence Livermore National Security, LLC
 * and other libROM project developers. See the top-level COPYRIGHT
 * file for details.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 *****************************************************************************/

// Description: An abstract class defining the interface to the generic SVD
//              algorithm.

#ifndef included_SVD_h
#define included_SVD_h

#include "linalg/Matrix.h"
#include "linalg/Options.h"
#include <vector>

namespace CAROM {

/**
 * Class SVD defines the interface to the generic SVD algorithm.  The API is
 * intentionally small.  One may collect the samples, compute the SVD, and get
 * the dimension of the system.
 */
class SVD
{
public:
    /**
     * @brief Constructor.
     *
     * @param[in] options The struct containing the options for this abstract
     * SVD class.
     * @see Options
     */
    SVD(
        Options options);

    /**
     * @brief Collect the new sample, u_in at supplied time.
     *
     * @pre u_in != 0
     *
     * @param[in] u_in The new sample.
     * @param[in] add_without_increase If true, the addLinearlyDependent is invoked.
     *                                 This only applies to incremental SVD.
     *
     * @return True if the sampling was successful.
     */
    virtual
    bool
    takeSample(
        double* u_in,
        bool add_without_increase) = 0;

    /**
     * @brief Returns the dimension of the system on this processor.
     *
     * @return The dimension of the system on this processor.
     */
    int
    getDim() const
    {
        return d_dim;
    }

    /**
     * @brief Returns the basis vectors for the current time interval.
     *
     * @return The basis vectors for the current time interval.
     */
    virtual
    std::shared_ptr<const Matrix>
    getSpatialBasis() = 0;

    /**
     * @brief Returns the temporal basis vectors for the current time interval.
     *
     * @return The temporal basis vectors for the current time interval.
     */
    virtual
    std::shared_ptr<const Matrix>
    getTemporalBasis() = 0;

    /**
     * @brief Returns the singular values for the current time interval.
     *
     * @return The singular values for the current time interval.
     */
    virtual
    std::shared_ptr<const Vector>
    getSingularValues() = 0;

    /**
     * @brief Returns the singular values for the current time interval.
     *
     * @return The singular values for the current time interval.
     */
    virtual
    std::shared_ptr<const Matrix>
    getSnapshotMatrix() = 0;

    /**
     * @brief Get the number of samples taken.
     *
     */
    int getNumSamples() const
    {
        return d_num_samples;
    }

    /**
     * @brief Get the maximum number of samples that can be taken.
     *        SVD class will return an error if the number of samples exceeds the maximum.
     *
     */
    int getMaxNumSamples() const
    {
        return d_max_num_samples;
    }

protected:
    /**
     * @brief Returns true if the next sample will result in a new time
     * interval.
     *
     * @return True if the next sample results in the creation of a new time
     *         interval.
     */
    bool
    isFirstSample() const
    {
        return (d_num_samples == 0);
    }

    /**
     * @brief Dimension of the system.
     */
    const int d_dim;

    /**
     * @brief Number of samples stored for the current time interval.
     */
    int d_num_samples;

    /**
     * @brief Number of rows in right singular matrix.
     */
    int d_num_rows_of_W;

    /**
     * @brief The maximum number of samples.
     */
    const int d_max_num_samples;

    /**
     * @brief The globalized basis vectors for the current time interval.
     *
     * The basis vectors are large and each process owns all of the basis
     * vectors.
     */
    std::shared_ptr<Matrix> d_basis;

    /**
     * @brief The globalized right basis vectors for the current time interval.
     *
     * Depending on the SVD algorithm, it may be  distributed across all
     * processors or each processor may own all of U.
     */
    std::shared_ptr<Matrix> d_basis_right;

    /**
     * @brief The matrix U which is large.
     *
     * Depending on the SVD algorithm, d_U may be  distributed across all
     * processors or each processor may own all of U.
     */
    std::shared_ptr<Matrix> d_U;

    /**
     * @brief The matrix U which is large.
     *
     * Depending on the SVD algorithm, d_W may be  distributed across all
     * processors or each processor may own all of U.
     */
    std::shared_ptr<Matrix> d_W;

    /**
     * @brief The vector S which is small.
     *
     * For all SVD algorithms, S is not distributed and the entire vector
     * exists on each processor.
     */
    std::shared_ptr<Vector> d_S;

    /**
     * @brief The globalized snapshot vectors for the current time interval.
     *
     * The snapshot vectors are large and each process owns all of the snapshot
     * vectors.
     */
    std::shared_ptr<Matrix> d_snapshots;

    /**
     * @brief Flag to indicate if results of algorithm should be printed for
     * debugging purposes.
     */
    bool d_debug_algorithm;

private:
    /**
     * @brief Unimplemented default constructor.
     */
    SVD();

    /**
     * @brief Unimplemented copy constructor.
     */
    SVD(
        const SVD& other);

    /**
     * @brief Unimplemented assignment operator.
     */
    SVD&
    operator = (
        const SVD& rhs);
};

}

#endif
